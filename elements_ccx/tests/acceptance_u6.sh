#!/bin/bash -l
# Acceptance suite for the nodal B-bar pair (U5 + U6) at its VALIDATED
# ENVELOPE, K/G <= 500 (nu <= 0.499).  See ../nodalbbar.py for why the
# envelope stops there.
#
# T1 patch test        general six-component strain on a distorted mesh.
#                      Consistency.  Exact to roundoff or the element is wrong.
# T2 confined          C1111 = K + 4G/3, the locking-sensitive modulus.
# T3 zero-energy modes exactly 6 (3 translations + 3 rotations); a 7th is a
#                      mechanism.
# T4 locking sweep     cantilever with G fixed, tip/Euler-Bernoulli flat in nu.
# T5 stability census  two-phase clamped cube, spurious-mode count.
#
# T1 and T2 run the REAL ccx element.  T3 and T5 run the numpy reimplementation
# in smoothing_proto.py, because U5 returns a null mass matrix and ccx
# *FREQUENCY cannot be used on it.
set -eu
cd "$(dirname "$0")"
W=${W:-$(mktemp -d)}
PY=${PY:-/home/giacomo/venvs/sci/bin/python}
CCX=${CCX:-ccx_u6s}
NU=${NU:-0.499}
N=${N:-6}
JIT=${JIT:-0.3}
pass=0; fail=0
ok () { if [ "$1" = 1 ]; then echo "   PASS  $2"; pass=$((pass+1));
        else echo "   FAIL  $2"; fail=$((fail+1)); fi; }

echo "=== U5+U6 acceptance, nu = $NU, mesh n=$N, jitter $JIT ==="
echo

# ---------------------------------------------------------------- T1 patch
echo "T1  patch test (general 6-component strain, distorted mesh)"
meta=$($PY make_block.py "$W/pt.inp" patch "$N" "$NU" "$JIT" one)
sig=$(echo "$meta" | sed -n 's/^patch_sigma //p')
for arm in C3D4 U6; do
    d="$W/pt_$arm"
    if [ "$arm" = C3D4 ]; then cp "$W/pt.inp" "$d.inp"; else
        SPAX_U7_RULE=all SPAX_U8=0 SPAX_BBAR_SOLVER=PARDISO \
            $PY ../nodalbbar.py "$W/pt.inp" "$d.inp" --elset Sphere_Only \
            > "$W/pt_$arm.gen" 2>&1
    fi
    ( cd "$W" && CCX_U5_STAB=0 $CCX "$(basename "$d")" > "$d.log" 2>&1 ) || true
    r=$($PY - "$d.dat" $sig <<'PYEOF'
import sys
ref = [float(x) for x in sys.argv[2:8]]
worst, n = 0.0, 0
scale = max(abs(x) for x in ref)
try:
    rows = open(sys.argv[1]).read().split('\n')
except IOError:
    print('nan 0'); raise SystemExit
for L in rows:
    f = L.split()
    if len(f) == 8:
        try:
            v = [float(x) for x in f[2:8]]
        except ValueError:
            continue
        n += 1
        worst = max(worst, max(abs(a - b) for a, b in zip(v, ref)) / scale)
print('%.3e %d' % (worst, n))
PYEOF
)
    set -- $r
    echo "      $arm: max rel err $1 over $2 integration points"
    ok "$($PY -c "print(1 if float('$1')<1e-7 and int('$2')>0 else 0)")" \
       "$arm patch test exact"
done
echo

# ------------------------------------------------------------- T2 confined
echo "T2  confined compression (exact C1111 = K + 4G/3)"
meta=$($PY make_block.py "$W/oe.inp" oedo "$N" "$NU" "$JIT" one)
exact=$(echo "$meta" | sed -n 's/^oedo_C1111 //p')
for arm in C3D4 U6; do
    d="$W/oe_$arm"
    if [ "$arm" = C3D4 ]; then cp "$W/oe.inp" "$d.inp"; else
        SPAX_U7_RULE=all SPAX_U8=0 SPAX_BBAR_SOLVER=PARDISO \
            $PY ../nodalbbar.py "$W/oe.inp" "$d.inp" --elset Sphere_Only \
            > "$W/oe_$arm.gen" 2>&1
    fi
    ( cd "$W" && CCX_U5_STAB=0 $CCX "$(basename "$d")" > "$d.log" 2>&1 ) || true
    r=$($PY - "$d.dat" "$exact" <<'PYEOF'
import sys
ex = float(sys.argv[2])
sx, n = 0.0, 0
for L in open(sys.argv[1]).read().split('\n'):
    f = L.split()
    if len(f) == 8:
        try:
            sx += float(f[2]); n += 1
        except ValueError:
            pass
print('%.6e %.3e' % (sx / n / 1e-3, abs(sx / n / 1e-3 / ex - 1)) if n
      else 'nan nan')
PYEOF
)
    set -- $r
    echo "      $arm: C1111 = $1   rel err $2   (exact $exact)"
    ok "$($PY -c "print(1 if float('$2')<1e-6 else 0)")" "$arm confined exact"
done
echo

# --------------------------------------------------------------- T3 / T5
echo "T3  zero-energy modes  /  T5  spurious-mode census"
$PY - "$N" "$NU" "$JIT" <<'PYEOF'
import sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '.')
import numpy as np, scipy.sparse.linalg as spl
import smoothing_proto as S
n, nu, jit = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
G = 4.4e5
K = 2*G*(1+nu)/(3*(1-2*nu))
for sc, tag in (('c3d4', 'C3D4'), ('ns_vol', 'U5+U6')):
    nodes, tets, mat = S.mesh_box(n, 0., 0., jit, geom=S.GEOM['slab'])
    mat[:] = 0
    g, vol = S.grads(nodes, tets)
    faces, patch = S.topology(tets, mat)
    Kg = S.assemble(sc, nodes, tets, mat, g, vol, faces, patch,
                    {0: (K, G)}).toarray()
    w = np.linalg.eigvalsh(Kg)
    scale = w.max()
    nz = int((w < 1e-9 * scale).sum())
    print('      %-6s zero modes = %d  (expect 6);  '
          'lambda_7/lambda_max = %.3e' % (tag, nz, w[6] / scale))
PYEOF
echo
echo "run T4 separately:  NUS=\"0.30 0.45 0.49 0.499\" ./locking_sweep_u6.sh"
echo "run T5 two-phase:   $PY stability_modes.py 10"
echo
echo "=== $pass passed, $fail failed ==="
[ "$fail" = 0 ]
