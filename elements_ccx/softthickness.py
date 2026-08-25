"""Local thickness of the soft phase, in ELEMENTS, by BFS over the tet graph.

WHY THIS METRIC. A nodal B-bar patch spans the centre node's 1-ring -- roughly
two element layers. Where the soft phase is locally thinner than that, the
patch averages the volumetric constraint over the full thickness of the body
and the material can breathe internally at no cost: an incompressible phase
then behaves like a compliant one. On BRKB_b280 that cost E_und 12.57 % against
plain C3D4 and put R 8.97 % below Abaqus C3D4H, four times outside the seed
floor.

WHY NOT A SIMPLER MEASURE. Every single-variable geometric metric tried on that
cell came back FLAT across a series whose error grew 8x: patch ring size
(11.3 -> 11.4), interface node fraction (79.2 % -> 71.9 %, the wrong way),
connected-body size (-17 %), tet quality (no slivers). Deleting the sphere
population fixed the error, but deleting the voids alone (-1.34 %) or the brine
inclusions alone (-2.11 %) did not -- only both together (-12.57 %), which is
strongly superadditive. No single population controls the quantity that
matters. Local thickness is the one field that all of them act on at once, so
it is measured directly rather than inferred from any one of them.

METHOD. A soft element is on the boundary if one of its four faces is not
shared with another soft element -- that face abuts the other phase or a free
surface. Those elements are depth 0; BFS then propagates inward through shared
faces. Depth d means d element layers of soft material lie between the element
and the nearest boundary, so the local thickness is about 2d+1 elements.

    python3 softthickness.py <deck.inp> [--elset NAME] [--warn 1.5]

Exit status 1 if the warn threshold is breached, so it can gate a campaign.
"""
import collections
import re
import sys


def read(inp, elsets):
    conn, member, cur, gen = {}, collections.defaultdict(set), None, False
    for ln in open(inp):
        u = ln.strip().upper().replace(' ', '')
        if ln.strip().startswith('*'):
            cur, gen = None, False
            if u.startswith('*ELEMENT') and ('TYPE=C3D4' in u or 'TYPE=U5' in u):
                m = re.search(r'ELSET=([^,]+)', u)
                cur = ('el', m.group(1) if m else None)
            elif u.startswith('*ELSET'):
                m = re.search(r'ELSET=([^,]+)', u)
                cur = ('es', m.group(1) if m else None)
                gen = 'GENERATE' in u
            continue
        s = ln.strip()
        if not s or cur is None:
            continue
        f = [x for x in s.replace(',', ' ').split() if x]
        if cur[0] == 'el' and len(f) >= 5:
            conn[int(f[0])] = tuple(int(x) for x in f[1:5])
            member[cur[1]].add(int(f[0]))
        elif cur[0] == 'es':
            try:
                v = [int(x) for x in f]
            except ValueError:
                continue
            if gen and len(v) >= 2:
                sp = v[2] if len(v) > 2 else 1
                member[cur[1]].update(range(v[0], v[1] + 1, sp))
            else:
                member[cur[1]].update(v)
    soft = set()
    for e in elsets:
        soft |= member.get(e.upper(), set())
    return conn, soft & set(conn)


FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def depths(conn, soft):
    """BFS depth in element layers from the soft phase's own boundary."""
    face_owner = collections.defaultdict(list)
    for e in soft:
        n = conn[e]
        for f in FACES:
            face_owner[tuple(sorted((n[f[0]], n[f[1]], n[f[2]])))].append(e)
    nbr = collections.defaultdict(list)
    seed = set()
    for f, owners in face_owner.items():
        if len(owners) == 2:
            a, b = owners
            nbr[a].append(b)
            nbr[b].append(a)
        else:
            # unshared face: the other phase, or the outside of the cell
            seed.update(owners)
    d = {e: 0 for e in seed}
    q = collections.deque(seed)
    while q:
        e = q.popleft()
        for m in nbr[e]:
            if m not in d:
                d[m] = d[e] + 1
                q.append(m)
    for e in soft:                       # a body with no boundary at all
        d.setdefault(e, -1)
    return d


def main(argv):
    inp = argv[1]
    elsets = []
    warn = None
    a = argv[2:]
    while '--elset' in a:
        i = a.index('--elset'); elsets.append(a[i + 1]); del a[i:i + 2]
    if '--warn' in a:
        i = a.index('--warn'); warn = float(a[i + 1]); del a[i:i + 2]
    if not elsets:
        elsets = ['Sphere_Only']
    conn, soft = read(inp, elsets)
    if not soft:
        raise SystemExit('softthickness: no elements in %s' % ','.join(elsets))
    d = depths(conn, soft)
    # thickness in elements ~ 2*depth + 1
    t = sorted(2 * d[e] + 1 for e in soft)
    n = len(t)
    hist = collections.Counter(t)
    print('%s' % inp)
    print('  soft phase: %d elements in ELSET=%s' % (n, ','.join(elsets)))
    print('  local thickness (elements)   count      %%     cum%%')
    cum = 0
    for k in sorted(hist):
        cum += hist[k]
        print('      %-4s %22d %6.1f%% %7.1f%%'
              % (k, hist[k], 100.0 * hist[k] / n, 100.0 * cum / n))
    thin = sum(v for k, v in hist.items() if k <= 1)
    med = t[n // 2]
    print('  median %d   |  <=1 element thick: %d (%.1f%%)'
          % (med, thin, 100.0 * thin / n))
    if warn is not None:
        frac = 100.0 * thin / n
        if med < warn or frac > 25.0:
            print('  *** WARNING: soft phase is under-resolved for nodal B-bar.')
            print('      A patch spans ~2 element layers, so where the local')
            print('      thickness is 1 element the volumetric constraint is')
            print('      averaged over the whole thickness and leaks.')
            return 1
        print('  OK for nodal B-bar at threshold %.2f' % warn)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
