#!/bin/bash -l
# Same-geometry mesh sweep for the F-barES-FEM-T4 arm.
#
# The Abaqus reference is NOT a convergence study: every one of its 16 rows is
# a separate packing (und and drn differ in porosity even within a seed), one
# seed is non-monotone, and neither seed's R increments shrink.  So it cannot
# say whether R converges or wanders.
#
# This can, for OUR element: the deck row is reused verbatim with only L_mesh
# overridden, SPAX_SEED fixed, and the drained twin is the SAME MESH with one
# elastic card rewritten.  R therefore carries NO packing noise at all -- the
# only thing that changes between points is h.
set -eu
cd /home/giacomo/SpaX/.claude/worktrees/CalculiX
LM=$1
TAG=$(echo "$LM" | tr -d '.')
ROOT=out_meshconv/$TAG
PY=/home/giacomo/venvs/sci/bin/python
export SPAX_CCX=ccx_fbar OMP_NUM_THREADS=8
export SPAX_MESH_TIMEOUT=14400 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821
# Without this the iterative solver converges loosely and the drained
# control comes back with a transverse reaction of 3.7e-03 instead of
# 2e-07 -- which looks like a broken deck.  Patch 0001 reads it.
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}
export SPAX_MESH_ORDER=1
mkdir -p "$ROOT"

$PY - params/rve_layermesh.csv "$LM" "$ROOT/LAY.und.csv" "$ROOT/LAY.drn.csv" <<'PYEOF'
import csv, sys
deck, lm, outu, outd = sys.argv[1:5]
rows = list(csv.DictReader(open(deck)))
r = dict(next(x for x in rows if x['run_id'] == 'LMESH_m0p0240_und_s1'))
r['run_id'] = 'LAY'; r['L_mesh'] = lm
for path, K in ((outu, None), (outd, '2.2e+06')):
    q = dict(r)
    if K: q['K_inclusion'] = K
    w = csv.DictWriter(open(path, 'w', newline=''), fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerow(q)
print('L_mesh=%s L=%s n_slabs=%s slab_vof=%s' % (lm, r['L'], r['n_slabs'], r['slab_vof']))
PYEOF

gen="$ROOT/gen"
[ -f "$gen/Job-LAY-utx.inp" ] || $PY -u SpaX_Standalone.py "$ROOT/LAY.und.csv" "$gen" > "$ROOT/gen.log" 2>&1
grep -m1 "Done: LAY" "$ROOT/gen.log" | sed 's/^ */  /'

for st in und drn; do
  w="$ROOT/$st"; rm -rf "$w"; mkdir -p "$w"; cp "$gen"/Job-LAY-utx.inp "$w/"
  if [ "$st" = drn ]; then
    $PY - "$w/Job-LAY-utx.inp" <<'PYEOF'
import sys
K, G = 2.2e6, 440029.33528897085
E = 9.0*K*G/(3.0*K+G); nu = (3.0*K-2.0*G)/(2.0*(3.0*K+G))
p = sys.argv[1]; L = open(p).readlines(); out = []; i = 0
while i < len(L):
    out.append(L[i])
    if L[i].strip().lower().startswith('*material, name=mat_inclusion'):
        out.append(L[i+1]); out.append('%r, %r\n' % (E, nu)); i += 3; continue
    i += 1
open(p, 'w').writelines(out)
PYEOF
  fi
  python3 SpaX_CalculiX.py convert "$w" > /dev/null
  # Direct solver for BOTH twins: R is a ratio of two moduli and an iterative
  # tolerance lands straight in it.
  sed -i 's/^\*STATIC.*/*STATIC, SOLVER=PARDISO/' "$w/Job-LAY-utx-ccx.inp"
done

# drained: plain C3D4, the denominator, exactly as report_abaqus_ratio expects
SPAX_CCX_REAL=ccx_fbar SPAX_CCX_MEMMAX=24G OMP_NUM_THREADS=8 \
  ccx_capped "$ROOT/drn/Job-LAY-utx-ccx" > "$ROOT/drn/solve.log" 2>&1
python3 SpaX_PostProcess.py "$ROOT/LAY.drn.csv" "$ROOT/drn" "$ROOT/drn.out.csv" \
  > "$ROOT/drn.post.log" 2>&1

# undrained: F-barES-FEM-T4 c=1
wu="$ROOT/und_fbar1"; rm -rf "$wu"; mkdir -p "$wu"
$PY elements_ccx/fbares.py "$ROOT/und/Job-LAY-utx-ccx.inp" \
    "$wu/Job-LAY-utx-ccx.inp" --elset Sphere_Only --cycles 1 2>&1 | tail -3
# IN CORE unless CCX_PARDISO_OOC is set deliberately.  And run from inside
# $wu: MKL_PARDISO_OOC_PATH is relative to the WORKING DIRECTORY, and pointing
# it at a directory that does not exist made PARDISO fail -- which stock ccx
# ignored, returning E_und 3.83e+11 against 5.52e+09.  ccx now stops on a
# PARDISO error instead, but the cwd still has to be right.
( cd "$wu" && CCX_FBAR_C=1 ${CCX_PARDISO_OOC:+CCX_PARDISO_OOC=$CCX_PARDISO_OOC} \
    ${CCX_PARDISO_OOC:+MKL_PARDISO_OOC_PATH=./ooc_temp} \
    SPAX_CCX_REAL=ccx_fbar SPAX_CCX_MEMMAX=${FBAR_MEM:-24G} OMP_NUM_THREADS=8 \
    ccx_capped Job-LAY-utx-ccx > solve.log 2>&1 ) || true
SPAX_CCX_SIGMA_FROM_RF=1 python3 SpaX_PostProcess.py "$ROOT/LAY.und.csv" "$wu" \
  "$ROOT/und_fbar1.out.csv" > "$ROOT/und_fbar1.post.log" 2>&1 || true

$PY - "$ROOT" "$LM" <<'PYEOF'
import csv, sys, os
root, lm = sys.argv[1:3]
def ex(p):
    if not os.path.isfile(p): return float('nan')
    for r in csv.DictReader(open(p)):
        try: return float(r['E_x'])
        except Exception: return float('nan')
    return float('nan')
u = ex(os.path.join(root, 'und_fbar1.out.csv')); d = ex(os.path.join(root, 'drn.out.csv'))
print('RESULT L_mesh=%s  E_und=%.6e  E_drn=%.6e  R=%.4f' % (lm, u, d, u/d if d else float('nan')))
PYEOF
