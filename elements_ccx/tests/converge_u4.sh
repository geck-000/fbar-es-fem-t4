#!/bin/bash -l
# Do C3D10 and U4 converge to the SAME answer?
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
ROOT=${ROOT:-out_u4conv}
L=${L:-0.10}
MESHES=${MESHES:-"0.010 0.008 0.006"}
NSLABS=${NSLABS:-2}
SLABVOF=${SLABVOF:-0.1000}
CPUS=${CPUS:-4}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MKL_THREADING_LAYER=GNU
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_SEED=20260822
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT/pack"
DISP=$($PY -c "print(0.01 * $L)")
HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion,n_slabs,slab_vof,bridge_fraction,n_bridges,slab_axis,bridge_correlation'
mkrow () {
    printf '%s\nLAY,%s,%s,Composite,9.43e9,0.33,0.0325,0.020,0.005,Uniaxial Tension X,%s,Uniaxial Tension Z,%s,0.0100,0.0225,2.2e9,0.48,0.80,0.1,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,440029.33528897085,%s,%s,0.2929,2,x,0.0\n' \
        "$HDR" "$L" "$1" "$DISP" "$DISP" "$NSLABS" "$SLABVOF" > "$2"
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
    for cfg in "c3d10 2 plain" "u4 1 u4"; do
        set -- $cfg
        name="$1"; ord="$2"; kind="$3"
        w="$ROOT/${tag}_$name"
        [ -f "$ROOT/${tag}_$name.out.csv" ] && continue
        rm -rf "$w"; mkdir -p "$w"
        cp "$ROOT/${tag}_o$ord"/Job-LAY-ut*.inp "$w/" 2>/dev/null || continue
        echo "  -- $lm $name"
        SPAX_CCX=ccx_u4 python3 SpaX_CalculiX.py convert "$w" > /dev/null
        if [ "$kind" = "u4" ]; then
            for f in "$w"/Job-LAY-ut*-ccx.inp; do
                python3 elements_ccx/u4ify.py "$f" "$f.u4" > /dev/null
                mv "$f.u4" "$f"
            done
        fi
        SPAX_CCX=ccx_u4 python3 SpaX_CalculiX.py solve "$w" \
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


print('Do C3D10 and U4 converge to the same E_x?')
print('%-9s %14s %14s %10s %11s %11s'
      % ('L_mesh', 'C3D10', 'U4', 'U4/C3D10', 'gap C3D10', 'gap U4'))
for lm in meshes:
    tag = 'm' + lm.replace('.', '')
    a, b = g(tag, 'c3d10', 'E_x'), g(tag, 'u4', 'E_x')
    print('%-9s %14.6e %14.6e %9.2f %% %11.1e %11.1e'
          % (lm, a, b, 100 * (b - a) / a if a else float('nan'),
             g(tag, 'c3d10', 'equilibrium_gap'), g(tag, 'u4', 'equilibrium_gap')))
print()
print('Converging together = both right and C3D4 is the locking error.')
print('Diverging = one of them is wrong.')
PYEOF
