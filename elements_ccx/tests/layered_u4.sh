#!/bin/bash -l
# U4 on a real periodic layered RVE: C3D4 vs C3D10 vs U4 on ONE geometry.
#
# The cell is deliberately small. The mixed system is symmetric indefinite, so
# it must go through SPOOLES with pivoting, and a saddle point pivots into far
# more fill-in than an SPD system of the same size -- the campaign cell at 218k
# equations passed 5.9 GB without finishing. Small enough to solve is the point;
# there is no Abaqus reference at this size, so the comparison is internal.
#
# Expectation from the cantilever benchmark: C3D4 reads far too stiff across
# the layers, U4 much softer, C3D10 somewhere between.
set -eu
cd "$(dirname "$0")/../.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
ROOT=${ROOT:-out_u4small}
L=${L:-0.10}
LMESH=${LMESH:-0.008}
NSLABS=${NSLABS:-2}
SLABVOF=${SLABVOF:-0.1000}
CPUS=${CPUS:-6}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_SEED=20260822
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT"
DISP=$($PY -c "print(0.01 * $L)")
HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion,n_slabs,slab_vof,bridge_fraction,n_bridges,slab_axis,bridge_correlation'
printf '%s\nLAY,%s,%s,Composite,9.43e9,0.33,0.0325,0.020,0.005,Uniaxial Tension X,%s,Uniaxial Tension Z,%s,0.0100,0.0225,2.2e9,0.48,0.80,0.1,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,440029.33528897085,%s,%s,0.2929,2,x,0.0\n' \
    "$HDR" "$L" "$LMESH" "$DISP" "$DISP" "$NSLABS" "$SLABVOF" > "$ROOT/und.csv"

for order in 1 2; do
    gen="$ROOT/o$order"
    if [ ! -f "$gen/Job-LAY-utx.inp" ]; then
        echo "==== generating order $order ===="
        SPAX_MESH_ORDER=$order "$PY" -u SpaX_Standalone.py "$ROOT/und.csv" "$gen" \
            > "$ROOT/o$order.gen.log" 2>&1 || { echo "  gen FAILED"; exit 1; }
    fi
    grep -m1 "Element types" "$ROOT/o$order.gen.log" | sed 's/^ */  /'
done

run () {   # $1 tag, $2 gen dir, $3 binary, $4 = u4|plain
    w="$ROOT/$1"; rm -rf "$w"; mkdir -p "$w"
    cp "$2"/Job-LAY-ut*.inp "$w/"
    echo "  -- $1"
    SPAX_CCX="$3" python3 SpaX_CalculiX.py convert "$w" > /dev/null
    if [ "$4" = "u4" ]; then
        for f in "$w"/Job-LAY-ut*-ccx.inp; do
            python3 elements_ccx/u4ify.py "$f" "$f.u4" | sed 's/^/     /'
            mv "$f.u4" "$f"
        done
    fi
    SPAX_CCX="$3" python3 SpaX_CalculiX.py solve "$w" --cpus "$CPUS" --jobs 1 \
        | sed 's/^/     /'
    python3 SpaX_PostProcess.py "$ROOT/und.csv" "$w" "$ROOT/$1.out.csv" \
        > "$ROOT/$1.post.log" 2>&1
}

run c3d4  "$ROOT/o1" ccx_spax plain
run c3d10 "$ROOT/o2" ccx_spax plain
run u4    "$ROOT/o1" ccx_u4   u4

echo
$PY - "$ROOT" <<'PYEOF'
import csv, os, sys
root = sys.argv[1]


def g(tag, col):
    p = os.path.join(root, tag + '.out.csv')
    if not os.path.isfile(p):
        return float('nan')
    for r in csv.DictReader(open(p)):
        try:
            return float(r[col])
        except (KeyError, ValueError, TypeError):
            return float('nan')
    return float('nan')


print('%-28s %14s %14s %12s' % ('', 'E_x (across)', 'E_z (in-plane)', 'eq gap'))
for tag, lab in (('c3d4', 'C3D4  (displacement)'),
                 ('c3d10', 'C3D10 (displacement)'),
                 ('u4', 'U4    (mixed, MINI)')):
    print('%-28s %14.6e %14.6e %12.2e'
          % (lab, g(tag, 'E_x'), g(tag, 'E_z'), g(tag, 'equilibrium_gap')))
b, u = g('c3d4', 'E_x'), g('u4', 'E_x')
if b == b and u == u and b:
    print()
    print('U4 is %+.2f %% against C3D4 across the layers.' % (100 * (u - b) / b))
    print('equilibrium_gap is a real check on the U4 row: the stiffness and the')
    print('force recovery share u4mat, unlike the B-bar patch.')
PYEOF
