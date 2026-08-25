#!/bin/bash -l
# Volumetric-locking test for U4 against the stock C3D4.
#
# A tip-loaded cantilever of linear tets, swept in Poisson's ratio with the
# SHEAR modulus held fixed at G=1 MPa. Holding G is what makes the sweep a
# clean test: the true compliance then goes as 1/E = 1/(2G(1+nu)), so the tip
# deflection falls by only 13% from nu=0.3 to the incompressible limit. An
# element that locks volumetrically stiffens without bound instead and its
# deflection collapses.
#
# Reported as delta_num / delta_EulerBernoulli, so the fixed discretisation
# error of a coarse linear-tet beam divides out and only the nu-dependence is
# left. Flat across the sweep = no locking. Collapsing = locking.
set -eu
cd "$(dirname "$0")"
WORK=${WORK:-$(mktemp -d)}
NX=${NX:-20}; NY=${NY:-4}; NZ=${NZ:-4}
NUS=${NUS:-"0.30 0.45 0.49 0.499 0.4999 0.49999"}
CCX=${CCX:-ccx_u4}

echo "cantilever ${NX}x${NY}x${NZ}, G fixed at 1 MPa, tip load 1 kN"
echo
printf '%-10s | %-24s | %-24s\n' 'nu' 'C3D4' 'U4  (MINI P1+bubble/P1)'
printf '%-10s | %11s %12s | %11s %12s\n' '' 'tip' 'vs EB' 'tip' 'vs EB'
for nu in $NUS; do
    line=$(printf '%-10s |' "$nu")
    for et in C3D4 U4; do
        tag="$WORK/b_${et}_${nu}"
        eb=$(python3 make_beam.py "$tag.inp" "$et" "$nu" "$NX" "$NY" "$NZ" \
             | sed 's/.*EB_tip=//')
        ( cd "$WORK" && "$CCX" "$(basename "$tag")" > /dev/null 2>&1 ) || true
        d=$(python3 - "$tag.dat" <<'PYEOF'
import sys
v, n = 0.0, 0
try:
    rows = open(sys.argv[1]).read().split('\n')
except IOError:
    print('nan'); raise SystemExit
for r in rows:
    f = r.split()
    # C3D4 prints node + 3 dof; U4 prints node + 4 (the fourth is pressure).
    if len(f) in (4, 5):
        try:
            v += float(f[2]); n += 1
        except ValueError:
            pass
print('%.6e' % (v / n) if n else 'nan')
PYEOF
)
        line="$line$(python3 -c "
d,eb='$d','$eb'
try: print('%12.5e %11.4f' % (float(d), float(d)/float(eb)))
except Exception: print('%12s %11s' % ('FAILED','-'))
") |"
    done
    echo "$line"
done
echo
echo "vs EB flat across nu = no volumetric locking."
echo "vs EB collapsing toward zero = locking."

# --- second test: refinement at fixed near-incompressible nu ---------------
#
# Volumetric locking is the discretisation pathology that refinement does NOT
# cure: the inf-sup constant does not improve with h, so a locking element
# stays stuck at the wrong answer however fine the mesh. A stable element
# converges. This is the same signature the campaign cells show -- ccx C3D4
# moves 1.06% between mesh levels on the undrained layered cell while Abaqus
# C3D4H moves 10.55% (see ../../calculix/README.md).
#
NUFIX=${NUFIX:-0.4999}
echo
echo "refinement at nu=$NUFIX (locking does not cure with h; stability does)"
printf '%-12s %14s %14s %10s\n' 'mesh' 'C3D4 tip' 'U4 tip' 'U4/C3D4'
for m in "8 2 2" "16 3 3" "24 4 4" "32 6 6"; do
    set -- $m
    r=""
    for et in C3D4 U4; do
        t="$WORK/r_${et}"
        python3 make_beam.py "$t.inp" "$et" "$NUFIX" "$1" "$2" "$3" > /dev/null
        ( cd "$WORK" && "$CCX" "$(basename "$t")" > /dev/null 2>&1 ) || true
        v=$(python3 - "$t.dat" <<'PYEOF'
import sys
tot, n = 0.0, 0
try:
    rows = open(sys.argv[1]).read().split('\n')
except IOError:
    print('nan'); raise SystemExit
for line in rows:
    f = line.split()
    if len(f) in (4, 5):
        try:
            tot += float(f[2]); n += 1
        except ValueError:
            pass
print('%.5e' % (tot / n) if n else 'nan')
PYEOF
)
        r="$r $v"
    done
    printf '%-12s %14s %14s %10s\n' "${1}x${2}x${3}" $r \
        "$(python3 -c "
a,b='$(echo $r | awk '{print $1}')','$(echo $r | awk '{print $2}')'
try: print('%.1fx' % (float(b)/float(a)))
except Exception: print('-')
")"
done
