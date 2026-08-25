"""Nodal-averaged B-bar for one phase of a CalculiX deck.

    python3 nodalbbar.py <in-ccx.inp> <out.inp> [--elset Sphere_Only]

Splits the phase's response in two so the bulk modulus never appears in an
element-local operator and cannot lock:

    K = K_dev   per tet, element type U5 (deviatoric only)
      + K_vol   per node, built here out of ccx primitives

The nodal volumetric strain is
    V_a     = sum over the phase's elements at node a of V_e/4
    theta_a = (1/V_a) sum over those elements of (V_e/4) div(u)|_e
            = sum over the 1-ring of  b_n^a . u_n
and its energy is (1/2) K V_a theta_a^2.

Rather than build a patch element, theta is carried as one extra dof on a
dummy node, tied by an *EQUATION and given energy by a grounded SPRING1 --
all existing ccx features, so no change to the matrix-structure code.

The scaling is folded into the equation so every spring is identical:
with  t_a = sqrt(K V_a) theta_a  the energy is exactly (1/2) t_a^2, i.e. unit
stiffness for all of them, and one *SPRING card covers the lot. Emitting one
card per node would otherwise mean tens of thousands of element sets.

    *EQUATION :  t_a  -  sum_n sqrt(K V_a) b_n^a . u_n  =  0

t_a is written first because ccx eliminates the leading term as the dependent
dof, and t_a is the one that must go.
"""
import os
import sys
import numpy as np


def main():
    argv = sys.argv[1:]
    # --elset may be given more than once. Covering only the soft phase gives
    # its interface nodes ONE-SIDED patches, and in a slab 2-3 elements thick
    # most of its nodes are interface nodes -- the averaging then has almost
    # nothing to average over, which is why the brine-only version tracked
    # plain C3D4 under refinement while the same element gained 50x on a
    # homogeneous cantilever. Pass every elset to average across the interface.
    elsets = []
    while '--elset' in argv:
        i = argv.index('--elset')
        elsets.append(argv[i + 1])
        del argv[i:i + 2]
    if not elsets:
        elsets = ['Sphere_Only']
    elset = elsets[0]
    # --min-body: leave small DISCONNECTED soft bodies as plain C3D4.
    #
    # A patch spans the centre node's 1-ring, roughly two element layers.  In a
    # body only a few elements across, that is a large fraction of the whole
    # body, so the volumetric constraint is averaged over most of it and the
    # body can breathe internally at no cost -- an incompressible inclusion
    # then behaves like a void.  Measured on BRKB_b280: the cell carries brine
    # inclusion spheres ~3.9 elements across, a patch covers ~51% of one, and
    # U6 read E_und 12.57% below C3D4 (R 8.97% below Abaqus C3D4H).  Deleting
    # the spheres from the same cell dropped that to 0.95%.
    #
    # The threshold has to act on ELEMENTS, not patches: U5 is deviatoric-only,
    # so a node left without a patch loses its volumetric stiffness entirely.
    # Small bodies therefore stay C3D4 -- complete elements, nothing removed.
    minbody = 0
    dropped_ids = set()
    while '--min-body' in argv:
        i = argv.index('--min-body')
        minbody = int(argv[i + 1])
        del argv[i:i + 2]
    src, dst = argv[0], argv[1]
    lines = open(src).read().split('\n')

    def iskw(l):
        return l.startswith('*') and not l.startswith('**')

    # --- parse nodes, the phase's elements, and the material -------------
    co, els, maxnode = {}, [], 0
    allel = []
    mode = None
    cur_elset = None
    for ln in lines:
        if iskw(ln):
            u = ln.upper().replace(' ', '')
            if u.startswith('*NODE') and 'OUTPUT' not in u and 'PRINT' not in u \
               and 'FILE' not in u:
                mode = 'n'
            elif u.startswith('*ELEMENT') and any(
                    ('ELSET=' + e).upper() in u for e in elsets):
                mode = 'e'
                cur_elset = next(e for e in elsets
                                 if ('ELSET=' + e).upper() in u)
            elif u.startswith('*ELEMENT') and 'TYPE=C3D4' in u:
                mode = 'a'          # the other phase: needed for anchoring
            else:
                mode = None
            continue
        f = [x.strip() for x in ln.split(',') if x.strip()]
        if mode == 'n' and len(f) >= 4:
            n = int(f[0])
            co[n] = np.array([float(f[1]), float(f[2]), float(f[3])])
            maxnode = max(maxnode, n)
        elif mode == 'e' and len(f) >= 5:
            els.append(([int(x) for x in f[1:5]], cur_elset, int(f[0])))
        elif mode == 'a' and len(f) >= 5:
            allel.append((int(f[0]), [int(x) for x in f[1:5]]))
    if not els:
        raise SystemExit('nodalbbar: no elements in ELSET=%s' % elset)

    if minbody > 0:
        # connected components of the soft phase, by shared nodes
        n2e = {}
        for i, (nl, _, _eid) in enumerate(els):
            for n in nl:
                n2e.setdefault(n, []).append(i)
        par = list(range(len(els)))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for n, idx in n2e.items():
            r = find(idx[0])
            for j in idx[1:]:
                rj = find(j)
                if rj != r:
                    par[rj] = r
        size = {}
        for i in range(len(els)):
            r = find(i)
            size[r] = size.get(r, 0) + 1
        keep = [i for i in range(len(els)) if size[find(i)] >= minbody]
        dropped = len(els) - len(keep)
        nb = sum(1 for r, s in size.items() if s < minbody)
        if dropped:
            print('  --min-body %d: %d of %d soft bodies (%d elements, %.1f%%)'
                  ' left as C3D4' % (minbody, nb, len(size), dropped,
                                     100.0 * dropped / len(els)))
        dropped_ids = set(els[i][2] for i in range(len(els))
                          if size[find(i)] < minbody)
        els = [els[i] for i in keep]
        if not els:
            raise SystemExit('nodalbbar: --min-body %d removed every element'
                             % minbody)

    # CCX_U5_STAB is applied entirely inside ccx -- e_c3d_u5.f gives the element
    # stab*K and u6patch.f scales the patch by (1-stab) -- because the patch K
    # is read from the material card at solve time, not written into the deck.
    # Nothing here has to change; it is validated and echoed so a deck carries a
    # record of the value it was generated alongside.
    stab = float(os.environ.get('CCX_U5_STAB', '0') or 0)
    # 1.0 is INCLUSIVE and is the parameter-free setting: a U7 tet then keeps
    # its full volumetric term and drops out of every patch, which is exactly
    # a plain C3D4.  Excluding it here rejected the default configuration.
    if not 0.0 <= stab <= 1.0:
        raise SystemExit('nodalbbar: CCX_U5_STAB must be in [0,1], got %r' % stab)

    # U7: a retyped tet all of whose nodes are INTERIOR to the retyped set.
    #
    # A node that also touches an element which is NOT retyped is anchored --
    # that element keeps its full volumetric stiffness and holds the node.  A
    # node with no such neighbour is restrained volumetrically by its patch
    # alone, and the patch operator has a null space: a pattern with theta_a=0
    # over a cluster of such nodes costs only G-level energy and, for brine
    # (K/G = 5000), grows until it dominates.  Measured on BRKB_b280: 2467
    # nodes at up to 19x the applied displacement, 100% of them interior and
    # 0% on an interface, forming sheets through the brine slab mid-surfaces
    # -- which is why it bled E_x (sheets normal to x) and left E_z alone.
    #
    # ccx stabilises exactly these (CCX_U5_STAB, e_c3d_u5.f).  Anchored tets
    # stay U5 and unstabilised, because stabilising them re-locks the element
    # where it was already correct.
    retyped = set(e for _, _, e in els)
    anchored = set()
    for eid, nl in allel:
        if eid not in retyped:
            anchored.update(nl)
    # SPAX_U7_RULE: all (default) | graded | any
    #
    # DEFAULT IS THE PARAMETER-FREE RULE.  With CCX_U5_STAB=1 the 'all' rule
    # says only: *a tet none of whose nodes is anchored keeps its own
    # volumetric term; a tet with at least one anchored node delegates to the
    # nodal patch.*  There is no scalar in that -- it is a statement about when
    # the patch operator can be trusted, and stab=1 on a U7 tet makes it
    # numerically identical to a plain C3D4 (deviatoric + full K, and it drops
    # out of every patch).  Measured on BRKB: 4 of 5 inside the Abaqus seed
    # floor, the miss being b150 at +1.63% on a 1.54% floor built from three
    # seeds, and b280 at -0.60% -- BETTER than the tuned graded rule's -1.35%.
    #
    # 'graded' (stab scaled by interior-node count) reaches 5 of 5 but only at
    # a fitted CCX_U5_STAB=0.07, with 0.09 points of margin at b150.  A fitted
    # constant that close to the floor is not worth the extra case; keep
    # 'graded' available for comparison, not as the default.
    #
    # SPAX_U8 (face-jump stabiliser) defaults OFF.  It works at gamma=0.01 but
    # is a weaker discriminator -- it penalises the jump to ONE neighbour,
    # while the deficit form penalises deviation from the volume-weighted
    # 1-ring average (measured discrimination 0.40 vs 0.98 between a clean and
    # a failing cell) -- and it adds ~254k elements for no gain.
    #
    # A tet is stabilised in proportion to how many of its 4 nodes are INTERIOR
    # -- i.e. touch no element that keeps its own volumetric term.  Those are
    # the only nodes the spurious mode lives on (measured on BRKB_b280: 100% of
    # its 2467 nodes interior, 0% on an interface), because an interface node
    # is already anchored by the neighbouring C3D4.
    #
    # The two all-or-nothing rules bracket the answer and neither works:
    #   all  -- only 4/4 tets, 21% of the phase.  Fixes b150 but a tet with
    #           three interior nodes still carries the mode, so b280 stalls at
    #           -2.29% even at stab = 1.
    #   any  -- >=1 interior node, 85-98% of the phase, i.e. the global case
    #           again, which re-locks b150 to +1.64%.
    # graded emits U7A..U7D for 1..4 interior nodes; ccx scales stab by n/4.
    # Count 0 stays plain U5, unstabilised.
    #
    # LETTERS, not digits: ccx's built-in element dispatch keys on digits in
    # the label (elements.f), which is what made the U6 types letters-only.
    _rule = os.environ.get('SPAX_U7_RULE', 'all').lower()
    retyped = set(e for _, _, e in els)
    anchored = set()
    for eid, nl in allel:
        if eid not in retyped:
            anchored.update(nl)
    nint = {e: sum(1 for n in nl if n not in anchored) for nl, _, e in els}
    if _rule == 'ring':
        # PARAMETER-FREE, and finer than 'all'.  'all' asks whether a tet's
        # nodes are anchored; the sharper question is whether each node's
        # PATCH reaches an anchor, because a patch that touches the interface
        # cannot host the null mode even when its centre node is interior.
        # A node is trusted if it is anchored or any node of its 1-ring is
        # (BFS depth <= 1); the tet is then weighted by the fraction of its
        # nodes that are NOT trusted.  Both the criterion and the weights come
        # from mesh topology -- there is no scalar to choose.
        ring_of = {}
        for nl, _, _e in els:
            for a in nl:
                ring_of.setdefault(a, set()).update(nl)
        trusted = set()
        for a, rg in ring_of.items():
            if a in anchored or (rg & anchored):
                trusted.add(a)
        nun = {e: sum(1 for n in nl if n not in trusted) for nl, _, e in els}
        u7grade = {e: k for e, k in nun.items() if k >= 1}
    elif _rule == 'uniform':
        # Domain-wide stabilisation, no topological targeting at all.  This is
        # the shape a local-projection / Dohrmann-Bochev stabiliser has: the
        # coefficient comes from the material (order 2G in the pressure-deficit
        # space, i.e. CCX_U5_STAB ~ 2G/K), not from asking which tets look
        # dangerous.  Carried so the coefficient can be probed WITHOUT the
        # targeting rule confounding it.
        u7grade = {e: 4 for e in nint}
    elif _rule == 'all':
        u7grade = {e: 4 for e, k in nint.items() if k == 4}
    elif _rule == 'any':
        u7grade = {e: 4 for e, k in nint.items() if k >= 1}
    else:
        u7grade = {e: k for e, k in nint.items() if k >= 1}
    # --pin FILE: node ids (one per line) whose patches have been SHOWN to be
    # untrustworthy by the parameter-free detector -- a displacement
    # fluctuation cannot exceed the applied macroscopic displacement, and C3D4
    # produces none, so any node above it is carrying the null mode.  Every tet
    # touching such a node keeps its own volumetric term.
    #
    # This closes the loop: solve, detect, pin, re-solve until the detector
    # reads zero.  The cell decides how much stabilisation it needs from its
    # own solution, so there is no constant anywhere -- neither a coefficient
    # nor a topological threshold.
    _pin = os.environ.get('SPAX_PIN', '')
    if _pin and os.path.isfile(_pin):
        pinned = set(int(x) for x in open(_pin).read().split() if x.strip())
        npin = 0
        for nl, _, e in els:
            if any(n in pinned for n in nl):
                if u7grade.get(e, 0) != 4:
                    npin += 1
                u7grade[e] = 4
        print('  --pin: %d node(s) flagged by the detector; %d further tet(s) '
              'keep their own K' % (len(pinned), npin))
    u7 = set(u7grade)
    if u7:
        import collections as _c
        h = _c.Counter(u7grade.values())
        print('  U7 rule=%s: %d of %d retyped tets stabilised %s; %d plain U5'
              % (_rule, len(u7), len(els),
                 dict(sorted(h.items())), len(els) - len(u7)))

    # SPAX_U8: face-jump stabilisers.  One rank-one element per INTERNAL face
    # of the soft phase, penalising gamma*K*V*(div|e1 - div|e2)^2 over the five
    # nodes of the two tets sharing that face.  Connectivity is ordered
    # (shared face, apex1, apex2) because e_c3d_u8 relies on that to rebuild
    # the two tets.  Faces on the phase boundary are skipped: there is no
    # second soft tet, and the neighbouring C3D4 already anchors those nodes.
    u8 = []
    if os.environ.get('SPAX_U8', '0') != '0':
        face = {}
        for nl, _, eid in els:
            for om in range(4):
                f = tuple(sorted(nl[:om] + nl[om + 1:]))
                face.setdefault(f, []).append((eid, nl[om]))
        for f, own in face.items():
            if len(own) == 2:
                u8.append((list(f), own[0][1], own[1][1]))
        print('  %d U8 face-jump stabilisers (%d internal faces of the soft '
              'phase)' % (len(u8), len(u8)))

    # material of each elset and its bulk modulus. Each patch spans one
    # material, so each needs its own section card and its own K.
    mat_of = {}
    for ln in lines:
        u = ln.upper().replace(' ', '')
        if u.startswith('*SOLIDSECTION'):
            es = mt = None
            for fld in ln.split(','):
                t = fld.strip()
                if t.upper().startswith('ELSET='):
                    es = t.split('=', 1)[1].strip()
                elif t.upper().startswith('MATERIAL='):
                    mt = t.split('=', 1)[1].strip()
            if es and mt:
                mat_of[es] = mt
    missing = [e for e in elsets if e not in mat_of]
    if missing:
        raise SystemExit('nodalbbar: no *SOLID SECTION for %s' % missing)

    def elastic_of(mt):
        for i2, ln2 in enumerate(lines):
            if ln2.strip().upper().replace(' ', '').startswith(
                    '*MATERIAL,NAME=' + mt.upper()):
                for j2 in range(i2 + 1, min(i2 + 6, len(lines))):
                    f2 = [x.strip() for x in lines[j2].split(',') if x.strip()]
                    if len(f2) == 2 and not lines[j2].startswith('*'):
                        return float(f2[0]), float(f2[1])
                break
        return None, None

    Kof = {}
    for e in elsets:
        E, nu = elastic_of(mat_of[e])
        if E is None:
            raise SystemExit('nodalbbar: no elastic card for material %s'
                             % mat_of[e])
        Kof[e] = E / (3.0 * (1.0 - 2.0 * nu))

        # VALIDATED ENVELOPE: K/G <= 500, i.e. nu <= 0.499.
        #
        # This element is NOT inf-sup stable.  It is locking-free -- that part
        # is solid, and it is what C3D4 is not:
        #
        #   cantilever tip/Euler-Bernoulli at nu=0.4999   0.6580 vs C3D4 0.0205
        #   two-phase census, lambda_1/G at K/G=5000      0.545  vs C3D4 47.78
        #
        # but the nodal patch operator has a null space, and the spurious mode
        # it admits grows in proportion to K/G, because it lets the soft phase
        # escape a K-level restraint at a G-level cost.  Measured on the
        # prototype sphere cell (elements_ccx/tests/smoothing_proto.py), max
        # fluctuation over the applied displacement:
        #
        #   K/G = 5000   5.69      K/G = 500   1.05      C3D4   0.20
        #
        # and on the real BRKB_b280 cell, 19x at K/G = 5000.  At K/G = 500 the
        # two-phase eigenvalue census finds NO mode with volumetric content
        # above 0.3, so this is the largest ratio at which the element passes
        # locking and stability together.
        #
        # Raising the ratio is a modelling choice, not a loss of physics: brine
        # is a liquid, G = 0, and 4.4e5 is already a regularisation.  The G
        # sweep in params/rve_brine.csv puts the cost of nu 0.4999 -> 0.499 at
        # +0.17% (b150) and +0.21% (b280) on R, against seed floors of 1.5-2.1%.
        # nu = 0.49 costs +1.73% / +1.90%, about a full floor, which is why the
        # envelope stops at 0.499 and not lower.  Use SPAX_BRINE_KG=500 when
        # generating the deck.
        #
        # SPAX_BBAR_KGMAX overrides the ceiling; 0 disables the check.
        _G = E / (2.0 * (1.0 + nu))
        _kg = Kof[e] / _G if _G > 0 else float('inf')
        _cap = float(os.environ.get('SPAX_BBAR_KGMAX', 500.0) or 0.0)
        if _cap and _kg > _cap * 1.001:
            raise SystemExit(
                'nodalbbar: %s has K/G = %.0f (nu = %.5f), outside this '
                'element\'s validated envelope of K/G <= %.0f (nu <= 0.499).\n'
                '  The nodal B-bar is locking-free but not inf-sup stable, and '
                'its spurious\n  mode grows as K/G: 1.05x the applied '
                'displacement at 500, 5.69x at 5000\n  on the prototype cell, '
                '19x on BRKB_b280.\n'
                '  Generate the deck with SPAX_BRINE_KG=500 (raises G at fixed '
                'K; costs\n  ~0.2%% on R), or set SPAX_BBAR_KGMAX to override '
                'deliberately.'
                % (mat_of[e], _kg, (3 * _kg - 2) / (2 * (3 * _kg + 1)), _cap))

    # --- nodal volumes and the averaged divergence operator --------------
    Va = {}
    b = {}                       # b[a][n] = 3-vector
    for e, es, _eid in els:
        p = np.array([co[n] for n in e])
        J = np.array([p[1] - p[0], p[2] - p[0], p[3] - p[0]]).T
        det = np.linalg.det(J)
        vol = abs(det) / 6.0
        if vol <= 0.0:
            continue
        Ji = np.linalg.inv(J)
        gr = np.zeros((4, 3))
        gr[1:, :] = Ji                      # rows of J^-1 are grad L1..L3
        gr[0, :] = -gr[1:, :].sum(axis=0)   # grads sum to zero
        w = vol / 4.0
        for a in e:
            key = (a, es)
            Va[key] = Va.get(key, 0.0) + w
            ba = b.setdefault(key, {})
            for k, n in enumerate(e):
                ba[n] = ba.get(n, np.zeros(3)) + w * gr[k]
    for key in b:
        for n in b[key]:
            b[key][n] /= Va[key]

    # --- emit ----------------------------------------------------------
    #
    # The volumetric term goes in as U6 PATCH ELEMENTS, not as *EQUATION +
    # SPRING1. That first design was verified correct in every part -- and in
    # ccx it still came out 1.55x too stiff on a confined block, because
    # adjacent patches overlap so each mesh dof appears as an independent term
    # in ~24 equations and ccx's MPC machinery does not survive that. One
    # patch alone reproduced an independent Python assembly exactly; the full
    # set did not. As elements the MPCs are gone.
    #
    # *USER ELEMENT fixes NODES= per TYPE, and ring sizes vary, so one type is
    # declared per distinct ring size.
    keys = sorted(Va, key=lambda k: (k[1], k[0]))
    ring = {k: [k[0]] + sorted(n for n in b[k] if n != k[0]) for k in keys}
    bysize = {}
    for k in keys:
        bysize.setdefault((k[1], len(ring[k])), []).append(k)
    # One *USER ELEMENT type per (material, ring size): NODES= is fixed per
    # type. Two suffix characters give 36^2 names -- the label only has to
    # start 'U6' for the dispatcher, and positions 7-8 carry ndof and nope.
    # LETTERS ONLY. ccx's built-in element dispatch keys on digits in the
    # label: elements.f:361 has "label(4:4).eq.'4' -> nope=4", :358 has
    # "label(4:5).eq.'10' -> nope=10", and :320 "'20' -> nope=20". A user type
    # named U614 is therefore claimed by the nope=4 rule before the
    # *USER ELEMENT lookup, so ccx reads 4 nodes instead of 18 and treats the
    # continuation line as a new element -- reported as "element N is already
    # defined" only when that line's first field happens to collide with a real
    # element id, which is why just 3 of ~36000 showed up. 26x26 names is ample.
    # A ring is capped at 255 nodes, and the cap is structural rather than
    # chosen: ccx carries a user element's node count in the single label byte
    # lakon(8:8) (userelements.f rejects NODES > 255 outright), and the element
    # matrix along the whole e_c3d_u* path is dimensioned 765 = 3*255.  Meshes
    # get closer to this than one would guess -- a sphere packing produced a
    # node carried by 248 tets, giving a 127-node ring -- so fail here, where
    # the offending node can be named, rather than inside the solver.
    toobig = [k for k in keys if len(ring[k]) > 255]
    if toobig:
        worst = max(toobig, key=lambda k: len(ring[k]))
        raise SystemExit(
            'nodalbbar: %d patch(es) exceed the 255-node ring limit; the '
            'worst is node %s with %d nodes. ccx stores a user element\'s '
            'node count in one byte, so this cannot be expressed at all. '
            'Remesh to remove the pathological valence.'
            % (len(toobig), worst[0], len(ring[worst])))

    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    names = [a + c for a in letters for c in letters]
    if len(bysize) > len(names):
        raise SystemExit('nodalbbar: %d (material, ring size) groups, more '
                         'than the U6 type names available' % len(bysize))
    suffix = {g: names[i] for i, g in enumerate(sorted(bysize))}

    # Patch element ids continue from the deck's real maximum. Starting at
    # 1e8 makes ccx's ne = max(ne, id) run to 100 million, and every
    # element-sized array (ipkon, lakon, ielmat, ...) is allocated to that.
    maxel = 0
    mode2 = None
    for ln in lines:
        if iskw(ln):
            mode2 = 'e' if ln.upper().replace(' ', '').startswith('*ELEMENT') \
                else None
            continue
        if mode2 == 'e':
            f2 = [x.strip() for x in ln.split(',') if x.strip()]
            if f2:
                try:
                    maxel = max(maxel, int(f2[0]))
                except ValueError:
                    pass
    out, done_step, u5decl = [], False, False
    # PER-BLOCK bookkeeping.  Retyped elements are held back and re-emitted
    # after *STEP, and they must go back under THEIR OWN *ELEMENT header.
    # Keying them by grade alone merged the elsets: with two --elset arguments
    # every U7 element from both was flushed under whichever header came last,
    # so on a two-phase deck the ice landed in ELSET=Sphere_Only and took the
    # brine's *SOLID SECTION.  Symptom was u6patch dying with 'node not in its
    # own connectivity', because the patches were built per material while the
    # elements were not.  Single-elset decks were never affected -- there is
    # only one header -- so brine-only results are unchanged.
    in_soft, blk, c3d4_back = False, -1, {}
    u7decl2 = False
    hdr_of = {}
    u7_back = {}                 # (block, grade) -> lines
    eid = maxel
    for ln in lines:
        if iskw(ln):
            u = ln.upper().replace(' ', '')
            if u.startswith('*ELEMENT') and any(
                    ('ELSET=' + e).upper() in u for e in elsets):
                if not u5decl:
                    out.append('*USER ELEMENT,TYPE=U5,NODES=4,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3')
                    u5decl = True
                if u7 and not u7decl2:
                    for _g in sorted(set(u7grade.values())):
                        out.append('*USER ELEMENT,TYPE=U7%s,NODES=4,'
                                   'INTEGRATIONPOINTS=1,MAXDOF=3'
                                   % 'ABCD'[_g - 1])
                    u7decl2 = True
                out.append(ln.replace('C3D4', 'U5').replace('c3d4', 'U5'))
                blk += 1
                hdr_of[blk] = ln
                # Elements in bodies below --min-body keep their volumetric
                # term and so must stay COMPLETE C3D4.  Retyping them to U5
                # without giving them patches would strip the volumetric
                # stiffness altogether -- far worse than the problem being
                # fixed -- so they are routed into a second block below.
                in_soft = True
                continue
            in_soft = False
            if u.startswith('*STEP') and not done_step:
                done_step = True
                if u8:
                    out.append('*USER ELEMENT,TYPE=U8,NODES=5,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3')
                for _b, _g in sorted(u7_back):
                    out.append('** SPAX: %d/4 nodes interior -- stabilised at '
                               '%d/4 of CCX_U5_STAB' % (_g, _g))
                    out.append(hdr_of[_b]
                               .replace('C3D4', 'U7' + 'ABCD'[_g - 1])
                               .replace('c3d4', 'U7' + 'ABCD'[_g - 1]))
                    out.extend(u7_back[(_b, _g)])
                for _b in sorted(c3d4_back):
                    out.append('** SPAX --min-body: %d element(s) in bodies too'
                               ' small to patch, kept as complete C3D4'
                               % len(c3d4_back[_b]))
                    out.append(hdr_of[_b])
                    out.extend(c3d4_back[_b])
                out.append('** SPAX nodal-averaged B-bar: %d patches, '
                           'one material each' % len(keys))
                # All *USER ELEMENT declarations first, then all *ELEMENT
                # blocks. ccx sorts the deck into per-keyword chains
                # (keystart.f: *USER ELEMENT is position 3, *ELEMENT is 4), and
                # interleaving them puts a multi-line element's continuation at
                # a chain boundary, where ccx reads it as a fresh element and
                # reports "element N is already defined".
                for grp in sorted(bysize):
                    es, sz = grp
                    out.append('*USER ELEMENT,TYPE=U6%s,NODES=%d,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3'
                               % (suffix[grp], sz))
                for grp in sorted(bysize):
                    es, sz = grp
                    t = 'U6' + suffix[grp]
                    out.append('*ELEMENT,TYPE=%s,ELSET=BBARP_%s' % (t, es))
                    for k in bysize[grp]:
                        eid += 1
                        row = [str(eid)] + [str(x) for x in ring[k]]
                        # ccx takes at most 16 fields per line, and a continued
                        # card MUST end in a comma. Without it ccx reads the
                        # next line as a fresh element and reports
                        # "element N is already defined" -- the first field of
                        # the continuation is a node number that collides with
                        # a real element id. Rings here reach 37 nodes, so
                        # every patch over 15 was silently malformed.
                        # textpart in elements.f is dimensioned (16), so a
                        # line must not exceed 16 fields; a trailing comma adds
                        # one more and overflows it, after which ccx reads the
                        # continuation as a fresh element and reports
                        # "element N is already defined". Keep to 15 fields and
                        # no trailing comma -- ccx continues on node count.
                        for c in range(0, len(row), 15):
                            out.append(','.join(row[c:c + 15]))
                for es in elsets:
                    if any(g[0] == es for g in bysize):
                        out.append('*SOLID SECTION,ELSET=BBARP_%s,MATERIAL=%s'
                                   % (es, mat_of[es]))
        if iskw(ln) and ln.strip().upper().startswith('*STATIC'):
            # PARDISO by default. The assembled system is SPD, so incomplete-
            # Cholesky PCG is valid -- but the patch rows are ~5x denser than
            # C3D4's and carry a K/mu = 5000 contrast spread across each node
            # patch, which the preconditioner does not like: on LMESH_m0p0240
            # PARDISO finished while ITERATIVE CHOLESKY was still running after
            # 24 minutes. Override with SPAX_BBAR_SOLVER if memory is tighter
            # than time.
            ln = '*STATIC, SOLVER=' + os.environ.get('SPAX_BBAR_SOLVER',
                                                     'PARDISO')
        if in_soft:
            f2 = [x.strip() for x in ln.split(',') if x.strip()]
            try:
                e2 = int(f2[0])
                if e2 in dropped_ids:
                    c3d4_back.setdefault(blk, []).append(ln)
                    continue
                if e2 in u7:
                    u7_back.setdefault((blk, u7grade[e2]), []).append(ln)
                    continue
            except (ValueError, IndexError):
                pass
        out.append(ln)

    if u8:
        # ids continue after the U6 patches.  Their OWN elset + *SOLID SECTION,
        # never the phase's: a stabiliser put in the phase elset joins its
        # *EL PRINT set, and printoutelem.f would then try to integrate a
        # 5-node element that has no material volume -- the same trap the U6
        # patches are kept out of.  It also corrupts the volume-averaged stress
        # the homogenisation reads.  A separate section card gives them K
        # without putting them in the printed set.
        out2 = []
        out2.append('*ELEMENT,TYPE=U8,ELSET=BBARJ_%s' % elset)
        for f, a1, a2 in u8:
            eid += 1
            out2.append('%d, %d, %d, %d, %d, %d'
                        % (eid, f[0], f[1], f[2], a1, a2))
        out2.append('*SOLID SECTION,ELSET=BBARJ_%s,MATERIAL=%s'
                    % (elset, mat_of[elset]))
        istep = next(i for i, l in enumerate(out)
                     if l.strip().upper().startswith('*STEP'))
        out[istep:istep] = out2

    open(dst, 'w').write('\n'.join(out))
    rsz = [len(ring[k]) for k in keys]
    print('nodalbbar: %s -> %s' % (src, dst))
    print('  %d elements retyped to U5 in ELSET=%s' % (len(els), elset))
    print('  %d U6 patches (one material each), ring %d..%d (mean %.1f), '
          '%d types' % (len(keys), min(rsz), max(rsz),
                        sum(rsz) / len(rsz), len(bysize)))
    for e in elsets:
        print('  %-14s material %-16s K = %.4e'
              % (e, mat_of[e], Kof[e]))
    print('  volumetric term delivered as elements -- no MPCs')


if __name__ == '__main__':
    main()
