#!/bin/bash -l
# Self-determining stabilisation: solve, detect, pin, repeat.
#
# The only criterion is physical and parameter-free: in a linear RVE with a
# prescribed macroscopic strain, no displacement FLUCTUATION can exceed the
# applied macroscopic displacement.  Plain C3D4 produces none; the nodal B-bar
# null mode produces thousands, at up to 19x.  So any node above the applied
# value is carrying the spurious mode, and the tets touching it must keep their
# own volumetric term instead of delegating to a patch that cannot see it.
#
# Iterating that to a fixed point lets each cell determine its own
# stabilisation from its own solution -- no coefficient, no threshold.
#
#   SRC=<c3d4 deck> DISP=<applied displacement> [ITERS=5] selfstab.sh <workdir>
set -eu
SRC=${SRC:?set SRC to the C3D4 deck}
DISP=${DISP:?set DISP to the applied macroscopic displacement}
ITERS=${ITERS:-5}
CCX=${CCX:-ccx_u6s}
PY=${PY:-/home/giacomo/venvs/sci/bin/python}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
W=${1:?workdir}
mkdir -p "$W"; : > "$W/pin.txt"
for it in $(seq 0 $((ITERS-1))); do
    SPAX_U7_RULE=all SPAX_U8=0 SPAX_PIN="$W/pin.txt" SPAX_BBAR_SOLVER=PARDISO \
        "$PY" "$ROOT/elements_ccx/nodalbbar.py" "$SRC" "$W/it$it.inp" \
        --elset Sphere_Only > "$W/gen$it.log" 2>&1
    $PY - "$W/it$it.inp" "$W/j.inp" <<'PYEOF'
import sys
t=open(sys.argv[1]).read(); i=t.upper().rindex('*END STEP')
open(sys.argv[2],'w').write(t[:i]+'*NODE FILE\nU\n'+t[i:])
PYEOF
    ( cd "$W" && CCX_U5_STAB=1.0 OMP_NUM_THREADS=4 CCX_NPROC_STIFFNESS=4 $CCX j > "run$it.log" 2>&1 )
    n=$($PY - "$W/j.frd" "$DISP" "$W/pin.txt" <<'PYEOF'
import sys
frd,disp,out=sys.argv[1],float(sys.argv[2]),sys.argv[3]
old=set()
try: old=set(int(x) for x in open(out).read().split() if x.strip())
except IOError: pass
bad=set(); inb=False; mx=0.0
for L in open(frd):
    if 'DISP' in L: inb=True; continue
    if not inb: continue
    if L.startswith(' -3'): break
    if L.startswith(' -1'):
        try:
            n=int(L[3:13])
            v=max(abs(float(L[13:25])),abs(float(L[25:37])),abs(float(L[37:49])))
        except ValueError: continue
        mx=max(mx,v)
        if v>disp: bad.add(n)
allb=old|bad
open(out,'w').write('\n'.join(str(x) for x in sorted(allb)))
print('%d %d %.4e'%(len(bad),len(allb),mx))
PYEOF
)
    set -- $n
    echo "  iter $it: $1 node(s) over the applied displacement (max |u| = $3); $2 pinned in total"
    cp "$W/j.dat" "$W/it$it.dat" 2>/dev/null || true
    if [ "$1" = 0 ]; then echo "  converged: no fluctuation exceeds the applied displacement"; break; fi
done
