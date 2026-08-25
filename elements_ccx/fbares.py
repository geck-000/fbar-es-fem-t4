#!/usr/bin/env python3
"""fbares.py -- rewrite a C3D4 deck as F-barES-FEM-T4.

Onishi, Iida & Amaya, Int. J. Comput. Methods 15(7) 1845003 (2018).  The
formulation, its small-strain reduction and the verification ladder are in
elements_ccx/docs/fbar_es_fem_t4.md; this script is only the deck side.

    K = K_dev   per EDGE (U2): V_h Bt^T D_dev Bt          eq. (1), (4), (13)
      + K_vol   per EDGE (U3): (K V_h) tbar^T sbar        eq. (6)-(11), (17)

The base tets stay in the deck retyped to U4 -- the null base tetrahedron --
contributing nothing: U2 and U3 carry the whole stiffness and both read the
tets for geometry through the 'U4' node->element map.

WHY THE CONNECTIVITY IS WRITTEN OUT AT ALL.  The elements recompute their own
weights from geometry -- which subsets of a ring form tets is not recoverable
from a node list -- but ccx has to know the stencil to size the element and to
build the matrix structure, so the node SET must be in the deck and must match
what u3vol walks exactly.  u3vol stops and names the element if it meets a
node the connectivity does not carry, so a mismatch is loud.

RING x SUPPORT: THE U3 ELEMENT MATRIX IS NOT DENSE OVER ITS STENCIL.

K_vol = (K V_h) tbar^T sbar is rank one, and its two factors have DIFFERENT
supports.  tbar is the UNSMOOTHED edge divergence of eq. (1) -- u3vol builds
it only from the tets that contain the edge -- so its support is exactly the
U2 ring, ~6.3 nodes.  Only sbar carries the wide E A^c support, ~33.7 nodes at
c = 1.  So s(ii,jj) can be nonzero only for i in the RING and j in the
SUPPORT: the (support - ring) x (support - ring) block is identically zero.

That block is most of the element, and mastruct.c was allocating structure for
all of it.  The connectivity is therefore written RING FIRST, the ring size is
carried in the type label, and mastruct.c / mafillsmas.f skip the outer block.
Measured at c = 1: insertions fall 3.3x and LMESH_m0p0120 becomes buildable.
Nothing about the assembled matrix changes -- the entries dropped are zero.

STENCIL WIDTH IS THE BINDING CONSTRAINT.  Measured on LMESH_m0p0240's soft
phase (117437 tets, 36323 nodes, 184572 edges):

    U2 deviatoric      mean  6.3 nodes/edge   max  14   ->   42 DOF
    U3 volumetric c=1  mean 33.7              max 173   ->  519 DOF
    U3 volumetric c=2  mean 93.6              max 494   -> 1482 DOF

c=2 cannot be expressed at all: ccx carries a user element's node count in the
single byte lakon(8:8) and userelements.f rejects NODES > 255.  This script
refuses it rather than letting the solver truncate.  c=1 needs the element
matrix widened from 150 to 520 DOF along the whole e_c3d_u* path.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import scipy.sparse as sp

LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
# Two suffix characters, LETTERS ONLY.  ccx's built-in element dispatch keys on
# digits in the label (elements.f: label(4:4).eq.'4' -> nope=4, '10' -> 10,
# '20' -> 20), so a type named U214 is claimed by the nope=4 rule before the
# *USER ELEMENT lookup and ccx reads 4 nodes instead of 14.  This showed up
# on only 3 of ~36000 elements.
# THE REAL CEILING ON c IS mastruct, NOT MEMORY AND NOT THE 255-NODE LIMIT.
#
# mastruct.c pushes one entry per off-diagonal (dof,dof) pair of the upper
# triangle of EVERY element onto `mast1`, with NO deduplication -- insert.c
# is a linked-list append, and the compression happens afterwards.  It costs
# 8 bytes each (mast1 and next, one ITG apiece) and the index is a 32-bit ITG
# in a stock build, so past 2^31 the 1.1x growth in insert.c overflows and
# ccx dies in u_realloc with a NEGATIVE allocation size.
#
# With the ring x support reduction above the count is
#
#     sum over edges of  T(3 n_h) - T(3 (n_h - r_h)) + T(3 r_h),  T(d)=d(d+1)/2
#
# still QUADRATIC in the stencil, so this stays the real ceiling on c --
# reached before memory is, and before the 255-node limit.
INS_LIMIT = 2 ** 31

NAMES = [a + b for a in LETTERS for b in LETTERS]


def iskw(line):
    s = line.lstrip()
    return s.startswith('*') and not s.startswith('**')


def parse(path):
    """nodes, tets, elset membership, elset -> material, material -> (E, nu)."""
    nodes, tets = {}, {}
    elset_of, mat_of, elastic = defaultdict(list), {}, {}
    mode, cur, curmat = None, None, None
    for ln in open(path):
        s = ln.strip()
        if not s or s.startswith('**'):
            continue
        if iskw(ln):
            u = s.upper().replace(' ', '')
            mode = None
            if u.startswith('*NODE') and 'OUTPUT' not in u:
                mode = 'n'
            elif u.startswith('*ELEMENT'):
                if 'TYPE=C3D4' in u:
                    mode = 'e'
                    cur = next((p.split('=')[1] for p in s.split(',')
                                if p.strip().upper().startswith('ELSET=')),
                               'ALL').strip()
            elif u.startswith('*SOLIDSECTION'):
                es = next((p.split('=')[1] for p in s.split(',')
                           if p.strip().upper().startswith('ELSET=')), None)
                mt = next((p.split('=')[1] for p in s.split(',')
                           if p.strip().upper().startswith('MATERIAL=')), None)
                if es and mt:
                    mat_of[es.strip()] = mt.strip()
            elif u.startswith('*MATERIAL'):
                curmat = next((p.split('=')[1] for p in s.split(',')
                               if p.strip().upper().startswith('NAME=')),
                              None)
                if curmat:
                    curmat = curmat.strip()
            elif u.startswith('*ELASTIC'):
                mode = 'el'
            continue
        f = [x.strip() for x in s.split(',') if x.strip()]
        if mode == 'n' and len(f) >= 4:
            nodes[int(f[0])] = (float(f[1]), float(f[2]), float(f[3]))
        elif mode == 'e' and len(f) >= 5:
            tets[int(f[0])] = tuple(int(x) for x in f[1:5])
            elset_of[cur].append(int(f[0]))
        elif mode == 'el' and curmat and len(f) >= 2:
            elastic.setdefault(curmat, (float(f[0]), float(f[1])))
            mode = None
    return nodes, tets, elset_of, mat_of, elastic


def stencils(conn, ncyc, chunk=20000):
    """Edge list and, per edge, the U2 ring and the U3 = E A^c node support.

    Built with boolean sparse products rather than a Python walk: on the real
    cells this is ~500k edges and the walk is the whole runtime.  The PATTERN
    is all that is needed -- the elements recompute the weights themselves --
    and verify_u8_chain.py has already checked that the weighted walk agrees
    with the prototype operator to 1e-15.

    conn is (ne, 4) LOCAL node indices for one material's tets.
    """
    ne = len(conn)
    nn = int(conn.max()) + 1 if ne else 0
    inc = sp.coo_matrix((np.ones(4 * ne, dtype=np.int8),
                         (np.repeat(np.arange(ne), 4), conn.ravel())),
                        shape=(ne, nn)).tocsr()
    inc.data[:] = 1

    # edges of this material, and the elements at each
    ekey = {}
    for e in range(ne):
        t = conn[e]
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = (t[i], t[j]) if t[i] < t[j] else (t[j], t[i])
                ekey.setdefault((int(a), int(b)), []).append(e)
    keys = sorted(ekey)
    rows = np.repeat(np.arange(len(keys)),
                     [len(ekey[k]) for k in keys])
    cols = np.fromiter((e for k in keys for e in ekey[k]), dtype=np.int64)
    E = sp.coo_matrix((np.ones(len(cols), dtype=np.int8), (rows, cols)),
                      shape=(len(keys), ne)).tocsr()
    E.data[:] = 1

    # A = P Q as a PATTERN: elements sharing a node
    A = (inc @ inc.T).tocsr()
    A.data[:] = 1

    ring, supp = [], []
    for lo in range(0, len(keys), chunk):
        blk = E[lo:lo + chunk]
        r = (blk @ inc).tocsr()          # U2: nodes of the edge's own tets
        r.data[:] = 1
        for i in range(r.shape[0]):
            ring.append(r.indices[r.indptr[i]:r.indptr[i + 1]])
        s = blk
        for _ in range(ncyc):
            s = (s @ A).tocsr()
            s.data[:] = 1
        w = (s @ inc).tocsr()            # U3: nodes of the E A^c support
        w.data[:] = 1
        for i in range(w.shape[0]):
            supp.append(w.indices[w.indptr[i]:w.indptr[i + 1]])
    return keys, ring, supp


def card(eid, conn_nodes):
    """One *ELEMENT card, wrapped.

    At most 15 fields per line and NO trailing comma.  textpart in elements.f
    is dimensioned (16); a trailing comma adds a field and overflows it, after
    which ccx reads the continuation as a fresh element and reports
    'element N is already defined'.  ccx continues on node count instead.
    """
    row = [str(eid)] + [str(x) for x in conn_nodes]
    return [','.join(row[c:c + 15]) for c in range(0, len(row), 15)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('dst')
    ap.add_argument('--elset', action='append', default=None,
                    help='element set to treat; repeat for more than one. '
                         'Default: every C3D4 set in the deck.')
    ap.add_argument('--cycles', type=int,
                    default=int(os.environ.get('CCX_FBAR_C', 1)),
                    help='c, the number of cyclic smoothings of J, eq. (6)-(7)')
    ap.add_argument('--solver', default=os.environ.get('SPAX_FBAR_SOLVER',
                                                       'PARDISO'))
    a = ap.parse_args()

    nodes, tets, elset_of, mat_of, elastic = parse(a.src)
    sets = a.elset or sorted(elset_of)
    for es in sets:
        if es not in elset_of:
            raise SystemExit('fbares: no C3D4 elements in ELSET=%s '
                             '(deck has %s)' % (es, sorted(elset_of)))
        if es not in mat_of:
            raise SystemExit('fbares: ELSET=%s has no *SOLID SECTION' % es)

    if a.cycles < 0 or a.cycles > 3:
        raise SystemExit('fbares: --cycles must be 0..3')

    lines = open(a.src).read().splitlines()
    maxel = max(tets) if tets else 0

    # --- build the stencils, per material -------------------------------
    #
    # PER MATERIAL, always.  The smoothing of eqs. (1) and (6)-(8) must not
    # cross a phase boundary: averaging a divergence across a 1000x modulus
    # contrast is not the method, and u6patch.f records that a shared patch
    # took equilibrium_gap from 2.5e-3 to 2.6e-1 on the finer layered cell.
    # An interface edge legitimately gets one smoothing domain per phase.
    u2, u3 = {}, {}
    stats = {}
    for es in sets:
        ids = sorted(elset_of[es])
        loc = {}
        conn = np.empty((len(ids), 4), dtype=np.int64)
        for i, e in enumerate(ids):
            for j, n in enumerate(tets[e]):
                if n not in loc:
                    loc[n] = len(loc)
                conn[i, j] = loc[n]
        back = np.empty(len(loc), dtype=np.int64)
        for n, i in loc.items():
            back[i] = n
        keys, ring, supp = stencils(conn, a.cycles)
        u2[es] = [(int(back[k[0]]), int(back[k[1]]), back[r]) 
                  for k, r in zip(keys, ring)]
        u3[es] = [(int(back[k[0]]), int(back[k[1]]), back[s])
                  for k, s in zip(keys, supp)]
        rs = np.array([len(r) for r in ring])
        ss = np.array([len(s) for s in supp])
        stats[es] = (len(ids), len(loc), len(keys), rs, ss)

    # --- the mastruct insertion wall ------------------------------------
    #
    # mastruct.c builds the matrix structure by pushing ONE entry per
    # (dof, dof) pair of every element onto `mast1` and compressing
    # afterwards, and its counter is a 32-bit ITG in a stock build.
    #
    # A U3 element is NOT dense over its whole stencil, and exploiting that
    # is what makes c = 1 reachable on a production cell -- see RING x
    # SUPPORT at the top of this file.  Only the first nring rows of the
    # element matrix can be nonzero, so the deck costs
    #
    #     sum over edges of  T(3 n_h) - T(3 (n_h - r_h))    (U3)
    #                      + T(3 r_h)                       (U2)
    #
    # with T(d) = d(d+1)/2, instead of T(3 n_h) + T(3 r_h) for the dense
    # block.  Still quadratic in the stencil, and still the real ceiling on
    # c, but a factor ~3.3 lower at c = 1.
    def _t(d):
        return d * (d + 1) // 2

    ins = 0
    for st in stats.values():
        rr = 3 * st[3].astype('int64')          # U2 ring dofs
        ds = 3 * st[4].astype('int64')          # U3 support dofs
        do = ds - rr                            # U3 dofs OUTSIDE the ring
        ins += int(_t(rr).sum() + (_t(ds) - _t(do)).sum())
    if ins >= INS_LIMIT:
        raise SystemExit(
            'fbares: this deck would make mastruct.c push %.2e (dof,dof) '
            'entries onto mast1, past the %.2e a 32-bit ITG can index -- ccx '
            'dies in u_realloc with a NEGATIVE allocation size while growing '
            'the list. Use --cycles %d, '
            'or coarsen the mesh -- the count is quadratic in the stencil, so '
            'one fewer smoothing cycle is worth about an order of magnitude.'
            % (ins, float(INS_LIMIT), max(a.cycles - 1, 0)))

    # --- the 255-node wall ----------------------------------------------
    worst = max((s[4].max(), es) for es, s in stats.items())
    if worst[0] > 255:
        raise SystemExit(
            'fbares: the c=%d volumetric stencil reaches %d nodes in '
            'ELSET=%s. ccx stores a user element\'s node count in the single '
            'byte lakon(8:8) and userelements.f rejects NODES > 255, so this '
            'cannot be expressed as an element at all. Use --cycles 1, or '
            'assemble the volumetric term outside the element loop (see '
            'docs/fbar_es_fem_t4.md section 8).'
            % (a.cycles, worst[0], worst[1]))

    # --- connectivity, ORDERED: ring first ------------------------------
    #
    # nl ALREADY contains both edge nodes -- it is the node set of the tets at
    # the edge (U2) or of the E A^c support (U3), and both contain the edge
    # itself.  Declaring len(nl)+2 made ccx read two fields past the end of
    # every card and swallow the next element's id, which showed up as
    # duplicate ids and a connectivity carrying its own successor.
    #
    # For U3 the ORDER now matters as well as the set: the first nring nodes
    # must be exactly the edge ring, because that is what mastruct.c and
    # mafillsmas.f use to skip the identically-zero outer block.  Getting it
    # wrong is not silent -- e_c3d_u3 checks that tbar vanishes past nring and
    # stops with the element number if it does not.
    groups = {}
    for es in sets:
        for na, nb, nl in u2[es]:
            rest = sorted(int(x) for x in nl if x != na and x != nb)
            groups.setdefault(('U2', es, len(nl), len(nl)), []).append(
                [na, nb] + rest)
        for (ra, rb, rl), (sa, sb, sl) in zip(u2[es], u3[es]):
            assert (ra, rb) == (sa, sb), 'u2/u3 edge lists out of step'
            ring = set(int(x) for x in rl)
            supp = set(int(x) for x in sl)
            if not ring <= supp:
                raise SystemExit(
                    'fbares: the edge ring of %d-%d is not contained in its '
                    'E A^c support -- the ring-first ordering the element '
                    'relies on is not well defined' % (ra, rb))
            inner = [ra, rb] + sorted(ring - {ra, rb})
            outer = sorted(supp - ring)
            groups.setdefault(('U3', es, len(supp), len(inner)), []).append(
                inner + outer)

    # --- type names -----------------------------------------------------
    #
    # U2:  U2<xy>              xy just disambiguates (kind, elset, size).
    # U3:  U3<r><xy>           r ENCODES nring: LETTERS[nring-1], read back by
    #                          mastruct.c as lakon(3:3) and by e_c3d_u3 and
    #                          mafillsmas.f the same way.  elements.f keys the
    #                          *USER ELEMENT lookup on label(2:5), so 'U' plus
    #                          four characters is all ccx can carry and three
    #                          suffix letters is the most that fits.
    suffix = {}
    for k in sorted(x for x in groups if x[0] == 'U2'):
        i = len(suffix)
        if i >= len(NAMES):
            raise SystemExit('fbares: more U2 (elset, size) groups than the '
                             '%d type names available' % len(NAMES))
        suffix[k] = NAMES[i]
    per = defaultdict(int)
    for k in sorted(x for x in groups if x[0] == 'U3'):
        nr = k[3]
        if not 1 <= nr <= len(LETTERS):
            raise SystemExit(
                'fbares: an edge ring spans %d nodes, and the ring size is '
                'carried in ONE label character (lakon(3:3), A..Z = 1..%d) so '
                'that mastruct.c can skip the zero outer block. Beyond that '
                'the reduction cannot be expressed; remesh, or drop the '
                'reduction and accept the dense structure.'
                % (nr, len(LETTERS)))
        j = per[nr]
        per[nr] += 1
        if j >= len(NAMES):
            raise SystemExit('fbares: more than %d U3 (elset, size) groups at '
                             'ring size %d' % (len(NAMES), nr))
        suffix[k] = LETTERS[nr - 1] + NAMES[j]

    # --- emit -------------------------------------------------------------
    out, u5decl, done_step = [], False, False
    drop, sawnp = False, any(
        iskw(l) and l.upper().replace(' ', '').startswith(('*NODEPRINT',
                                                           '*NODEFILE'))
        for l in lines)
    eid = maxel
    for ln in lines:
        if not iskw(ln):
            # a data line of a dropped block goes with it
            if drop:
                continue
        else:
            u = ln.upper().replace(' ', '')
            # NO ELEMENT IN AN F-bar DECK CARRIES STRESS.  The base tets are
            # U4 and return a null matrix; U2 and U3 are smoothing domains
            # with no shape function of their own and write nothing to stx.
            # Left in, *EL PRINT would report a column of exact zeros as if
            # it were the answer.  Read the result from displacements and
            # reactions instead -- which is why a *NODE PRINT is added below
            # if the deck has none.
            drop = u.startswith('*ELPRINT') or u.startswith('*ELFILE')
            if drop:
                out.append('** SPAX fbares: dropped '
                           + ln.split(',')[0].strip()
                           + ' -- no element in an F-bar deck carries stress')
                continue
            if u.startswith('*ENDSTEP') and not sawnp:
                sawnp = True
                # *NODE FILE, not *NODE PRINT: printing needs an NSET and
                # ccx does not define NALL itself (cgx does), so a
                # *NODE PRINT,NSET=NALL is answered with 'node set NALL does
                # not exist' and an EMPTY .dat.  *NODE FILE needs no set and
                # writes the .frd the post-processing already reads.
                out.append('*NODE FILE')
                out.append('U,RF')
            if u.startswith('*ELEMENT') and any(
                    ('ELSET=' + e).upper() in u for e in sets):
                if not u5decl:
                    out.append('*USER ELEMENT,TYPE=U4,NODES=4,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3')
                    u5decl = True
                out.append(ln.replace('C3D4', 'U4').replace('c3d4', 'U4'))
                continue
            if u.startswith('*STATIC'):
                # The volumetric operator of eq. (17) is NOT symmetric, so the
                # assembled matrix is not either and incomplete-Cholesky PCG
                # does not apply.  PARDISO takes it through the asymmetric
                # path (nasym -> mafillsmas -> mtype=11).
                ln = '*STATIC, SOLVER=' + a.solver
            if u.startswith('*STEP') and not done_step:
                done_step = True
                # ALL *USER ELEMENT declarations first, then all *ELEMENT
                # blocks.  ccx sorts the deck into per-keyword chains
                # (keystart.f: *USER ELEMENT is position 3, *ELEMENT is 4), so
                # interleaving them puts a multi-line element's continuation
                # at a chain boundary, where ccx reads it as a fresh element
                # and reports 'element N is already defined'.
                out.append('** SPAX F-barES-FEM-T4(c=%d): %d edge domains'
                           % (a.cycles, sum(len(u2[e]) for e in sets)))
                for k in sorted(groups):
                    out.append('*USER ELEMENT,TYPE=%s%s,NODES=%d,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3'
                               % (k[0], suffix[k], k[2]))
                for k in sorted(groups):
                    kind, es, sz, nr = k
                    tag = 'FBD' if kind == 'U2' else 'FBV'
                    out.append('*ELEMENT,TYPE=%s%s,ELSET=%s_%s'
                               % (kind, suffix[k], tag, es))
                    for conn in groups[k]:
                        eid += 1
                        # konl(1), konl(2) ARE THE EDGE NODES, in both U2 and
                        # U3; u2edge/u3vol identify the edge from them and
                        # find the tets themselves.  For U3, konl(1..nring)
                        # is the edge ring -- see above.
                        out.extend(card(eid, conn))
                # Their OWN elset + *SOLID SECTION, never the phase's: an
                # element in the phase elset joins its *EL PRINT set, and
                # printoutelem.f would try to integrate a smoothing domain
                # that has no material volume of its own -- and it would
                # corrupt the volume-averaged stress the homogenisation reads.
                for k in sorted(groups):
                    kind, es, sz, nr = k
                    tag = 'FBD' if kind == 'U2' else 'FBV'
                    out.append('*SOLID SECTION,ELSET=%s_%s,MATERIAL=%s'
                               % (tag, es, mat_of[es]))
        out.append(ln)

    open(a.dst, 'w').write('\n'.join(out) + '\n')

    print('fbares: %s -> %s   (c = %d)' % (a.src, a.dst, a.cycles))
    for es in sets:
        nel, nnd, ned, rs, ss = stats[es]
        E, nu = elastic.get(mat_of[es], (float('nan'), float('nan')))
        kg = (1.0 / (3.0 * (1.0 - 2.0 * nu))) / (1.0 / (2.0 * (1.0 + nu)))
        print('  %-14s %7d tets  %6d nodes  %7d edges   material %-14s '
              'K/G = %.0f' % (es, nel, nnd, ned, mat_of[es], kg))
        print('      U2 ring   mean %5.1f  p99 %4d  max %4d nodes '
              '(%d DOF)' % (rs.mean(), np.percentile(rs, 99), rs.max(),
                            3 * (rs.max() + 0)))
        print('      U3 stencil mean %5.1f  p99 %4d  max %4d nodes '
              '(%d DOF)' % (ss.mean(), np.percentile(ss, 99), ss.max(),
                            3 * (ss.max() + 0)))
    ndof = 3 * max(s[4].max() for s in stats.values())
    print('  %d element types; widest element %d DOF; %.2e (dof,dof) '
          'insertions into mastruct = %.1f GB (limit %.2e)'
          % (len(groups), ndof, ins, ins * 8.0 / 2.0 ** 30, float(INS_LIMIT)))
    if ndof > 765:
        print('  NOTE: the e_c3d_u* family and mafillsm.f hold 765 DOF '
              '(255 nodes, the lakon(8:8) encoding limit).  %d DOF cannot '
              'be a ccx user element at all.' % ndof)
    print('  run with CCX_FBAR_C=%d -- the same value the deck was '
          'generated with' % a.cycles)


if __name__ == '__main__':
    main()
