#!/bin/bash -l
# Do C3D10 and U5+U6 converge to the SAME answer?
#
# Both elements are consistent -- each passes a patch test -- so on a sequence
# of refined meshes they must approach one limit. If they do, that limit is the
# physical answer and the C3D4 result is the locking error. If they do not, one
# of them is wrong and the agreement at any single mesh is a coincidence.
#
# The geometry is frozen with SPAX_SAVE_PACKING/SPAX_LOAD_PACKING so only
# L_mesh changes; the slabs and bridges are deterministic from the deck, but
# the spherical population is not.
set -eu
cd "$(dirname "$0")/../.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
ROOT=${ROOT:-out_u6conv}
L=${L:-0.10}
MESHES=${MESHES:-"0.010 0.008 0.006"}
NSLABS=${NSLABS:-2}
BRIDGEF=${BRIDGEF:-0.2929}
NBRIDGES=${NBRIDGES:-2}
SLABVOF=${SLABVOF:-0.1000}
CPUS=${CPUS:-4}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MKL_THREADING_LAYER=GNU
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_SEED=20260822
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT/pack"
DISP=$($PY -c "print(0.01 * $L)")
HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion,n_slabs,slab_vof,bridge_fraction,n_bridges,slab_axis,bridge_correlation'
# NO SPHERE POPULATION.  The u4 version of this test carried VoF_sphere=0.0325
# and relied on SPAX_SAVE_PACKING to freeze it, but the freeze does not hold
# across a re-mesh: at L_mesh 0.010/0.008/0.006 the order-2 (C3D10) run drew a
# DIFFERENT packing from the order-1 runs -- porosity 0.0178 against 0.0119,
# flipping which arm got which -- so C3D10 read 4.7e9 against 6.1e9 and the
# "convergence" was two different cells being compared.  The slabs and bridges
# ARE deterministic from the deck, so dropping the spheres makes every arm at
# every mesh the identical geometry by construction, which is the only way the
# refinement column means anything.  It also leaves exactly the physics under
# test: brine layers at K/G ~ 5000.
mkrow () {
    printf '%s\nLAY,%s,%s,Composite,9.43e9,0.33,0.0000,0.020,0.005,Uniaxial Tension X,%s,Uniaxial Tension Z,%s,0.0000,0.0000,2.2e9,0.48,0.80,0.1,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,440029.33528897085,%s,%s,%s,%s,x,0.0\n' \
        "$HDR" "$L" "$1" "$DISP" "$DISP" "$NSLABS" "$SLABVOF" \
        "$BRIDGEF" "$NBRIDGES" > "$2"
}

FIRST=$(echo $MESHES | awk '{print $1}')
mkrow "$FIRST" "$ROOT/seed.csv"
if [ ! -f "$ROOT/pack/LAY.npy" ]; then
    echo "==== freezing the packing ===="
    SPAX_SAVE_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
        "$ROOT/seed.csv" "$ROOT/seedgen" > "$ROOT/seed.log" 2>&1 || true
fi

for lm in $MESHES; do
    tag="m${lm/./}"
    mkrow "$lm" "$ROOT/$tag.csv"
    for order in 1 2; do
        gen="$ROOT/${tag}_o$order"
        if [ ! -f "$gen/Job-LAY-utx.inp" ]; then
            echo "==== $lm order $order: meshing ===="
            SPAX_MESH_ORDER=$order SPAX_LOAD_PACKING="$ROOT/pack" \
                "$PY" -u SpaX_Standalone.py "$ROOT/$tag.csv" "$gen" \
                > "$ROOT/$tag.o$order.gen.log" 2>&1 || { echo " gen FAILED"; continue; }
        fi
    done
    # C3D10 from the order-2 mesh, U4 from the order-1 mesh
    # c3d4   -- the locking baseline, order 1
    # c3d10  -- the trusted non-locking reference (README: "plain C3D10, and
    #           nothing else needs building"), order 2, SAME frozen packing
    # u6     -- U5+U6 with patches over the BRINE ONLY, as the campaign runs it
    # u6all  -- U5+U6 with patches over BOTH phases.  nodalbbar's own comment
    #           says brine-only leaves interface nodes with one-sided patches,
    #           and in a thin slab most brine nodes ARE interface nodes.  The
    #           README lists the interface as a design assertion, not a
    #           measurement.  This arm is that measurement.
    for cfg in "c3d4 1 plain" "c3d10 2 plain" "u6 1 u6" "u6all 1 u6all"; do
        set -- $cfg
        name="$1"; ord="$2"; kind="$3"
        w="$ROOT/${tag}_$name"
        [ -f "$ROOT/${tag}_$name.out.csv" ] && continue
        rm -rf "$w"; mkdir -p "$w"
        cp "$ROOT/${tag}_o$ord"/Job-LAY-ut*.inp "$w/" 2>/dev/null || continue
        echo "  -- $lm $name"
        SPAX_CCX=${CCX:-ccx_u6} python3 SpaX_CalculiX.py convert "$w" > /dev/null
        if [ "$kind" = "u6" ] || [ "$kind" = "u6all" ]; then
            if [ "$kind" = "u6all" ]; then
                ES="--elset Sphere_Only --elset Matrix_Only"
            else
                ES="--elset Sphere_Only"
            fi
            for f in "$w"/Job-LAY-ut*-ccx.inp; do
                SPAX_BBAR_SOLVER="${BBAR_SOLVER:-PARDISO}" "$PY" \
                    elements_ccx/nodalbbar.py "$f" "$f.u6" $ES \
                    >> "$ROOT/$tag.$name.u6gen.log" 2>&1 || { echo "  gen FAILED"; continue 2; }
                mv "$f.u6" "$f"
            done
        fi
        SPAX_CCX=${CCX:-ccx_u6} python3 SpaX_CalculiX.py solve "$w" \
            --cpus "$CPUS" --jobs 1 | sed 's/^/     /'
        python3 SpaX_PostProcess.py "$ROOT/$tag.csv" "$w" \
            "$ROOT/${tag}_$name.out.csv" > "$ROOT/${tag}_$name.post.log" 2>&1
    done
done

echo
$PY - "$ROOT" "$MESHES" <<'PYEOF'
import csv, os, sys
root, meshes = sys.argv[1], sys.argv[2].split()


def g(tag, name, col):
    p = os.path.join(root, '%s_%s.out.csv' % (tag, name))
    if not os.path.isfile(p):
        return float('nan')
    for r in csv.DictReader(open(p)):
        try:
            return float(r[col])
        except (KeyError, ValueError, TypeError):
            return float('nan')
    return float('nan')


print('Do C3D10 and U5+U6 converge to the same E_x?  (frozen packing, order-1')
print('meshes for C3D4/U6, the order-2 mesh of the SAME cell for C3D10)')
print()
arms = [('c3d4', 'C3D4'), ('c3d10', 'C3D10'), ('u6', 'U5+U6 brine'),
        ('u6all', 'U5+U6 both')]
print('%-9s %-13s %14s %11s %11s' % ('L_mesh', 'element', 'E_x', 'vs C3D10', 'gap'))
for lm in meshes:
    tag = 'm' + lm.replace('.', '')
    ref = g(tag, 'c3d10', 'E_x')
    for key, label in arms:
        v = g(tag, key, 'E_x')
        if v != v:
            continue
        rel = 100 * (v - ref) / ref if ref == ref and ref else float('nan')
        print('%-9s %-13s %14.6e %10.2f%% %11.1e'
              % (lm, label, v, rel, g(tag, key, 'equilibrium_gap')))
    print()
print('C3D10 is the trusted non-locking reference on this cell.')
print('An element that is RIGHT closes on C3D10 as the mesh refines.')
print('C3D4 should stay high (locking).  A too-soft element undershoots and')
print('keeps going -- that is how U4 CAPPED failed, and it only showed up')
print('under refinement, never at a single mesh.')
PYEOF
