#!/bin/bash -l
# Volumetric-locking sweep for the NODAL B-BAR pair (U5 + U6 patches),
# alongside the stock C3D4 and MINI (U4) already covered by locking_sweep.sh.
#
# Same construction: a tip-loaded cantilever of linear tets swept in Poisson's
# ratio with the SHEAR modulus held fixed at G = 1 MPa, so the true compliance
# only falls 13% from nu=0.3 to the incompressible limit.  Reported as
# tip/Euler-Bernoulli, which divides out the coarse-mesh discretisation error
# and leaves the nu-dependence alone.
#
#   flat across nu      -> no volumetric locking
#   collapsing to zero  -> locking (C3D4 reaches 0.0177 at nu=0.4999)
#
# STAB is the graded stabilisation coefficient; the arm is run at each value
# given so the locking cost of stabilising is visible rather than assumed.
set -eu
cd "$(dirname "$0")"
WORK=${WORK:-$(mktemp -d)}
NX=${NX:-20}; NY=${NY:-4}; NZ=${NZ:-4}
NUS=${NUS:-"0.30 0.45 0.49 0.499 0.4999 0.49999"}
STABS=${STABS:-"0 0.07 1.0"}
CCX=${CCX:-ccx_u6s}
PY=${PY:-/home/giacomo/venvs/sci/bin/python}

echo "cantilever ${NX}x${NY}x${NZ}, G fixed at 1 MPa, tip load 1 kN"
echo "tip / Euler-Bernoulli;  flat = no locking,  -> 0 = locking"
echo
hdr=$(printf '%-9s %11s' 'nu' 'C3D4')
for s in $STABS; do hdr="$hdr$(printf '%14s' "U5+U6 s=$s")"; done
echo "$hdr"

tipval () {   # $1 = .dat
    $PY - "$1" <<'PYEOF'
import sys
v, n = 0.0, 0
try:
    rows = open(sys.argv[1]).read().split('\n')
except IOError:
    print('nan'); raise SystemExit
for r in rows:
    f = r.split()
    if len(f) in (4, 5):
        try:
            v += float(f[2]); n += 1
        except ValueError:
            pass
print('%.6e' % (v / n) if n else 'nan')
PYEOF
}

for nu in $NUS; do
    line=$(printf '%-9s' "$nu")
    base="$WORK/b_${nu}"
    eb=$($PY make_beam.py "$base.inp" C3D4 "$nu" "$NX" "$NY" "$NZ" \
         | sed 's/.*EB_tip=//')
    ( cd "$WORK" && "$CCX" "$(basename "$base")" > /dev/null 2>&1 ) || true
    d=$(tipval "$base.dat")
    line="$line$($PY -c "
try: print('%11.4f' % (float('$d')/float('$eb')))
except Exception: print('%11s' % 'FAILED')")"
    # the U5+U6 arm: same mesh, whole beam retyped, patches over EALL
    u="$WORK/u_${nu}"
    SPAX_U7_RULE=graded SPAX_U8=0 SPAX_BBAR_SOLVER=SPOOLES \
        $PY ../nodalbbar.py "$base.inp" "$u.inp" --elset EALL \
        > "$WORK/gen_${nu}.log" 2>&1 || true
    for s in $STABS; do
        ( cd "$WORK" && CCX_U5_STAB="$s" "$CCX" "$(basename "$u")" \
            > /dev/null 2>&1 ) || true
        du=$(tipval "$u.dat")
        line="$line$($PY -c "
try: print('%14.4f' % (float('$du')/float('$eb')))
except Exception: print('%14s' % 'FAILED')")"
    done
    echo "$line"
done
