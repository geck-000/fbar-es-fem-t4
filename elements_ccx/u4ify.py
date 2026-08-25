"""Convert one phase of a CalculiX deck to the U4 mixed tetrahedron.

    python3 u4ify.py <in-ccx.inp> <out.inp> [--elset Sphere_Only]

Three edits:

1. The target element set changes from C3D4 to U4, and a *USER ELEMENT card is
   emitted for it. (ccx reorders keywords internally -- keystart.f puts
   *USER ELEMENT at position 3, ahead of *ELEMENT -- so the card only has to
   precede the first *STEP.)

2. PERIODIC PRESSURE EQUATIONS. This is the part the generator cannot know
   about. Its periodic constraints cover dofs 1-3 and carry the macroscopic
   strain jump through a reference node:

       *EQUATION
       3
       <master>, 1,  1
       <slave>,  1, -1
       <refnode>,1, -1

   Pressure is a scalar fluctuation with no macroscopic jump, so its periodic
   condition is simply p(master) = p(slave) -- two terms, no reference node.
   Without these the pressure field is unconstrained across the cell faces,
   which matters precisely in the layered decks, where the brine slab spans the
   cell and reaches the periodic boundary.

   Pairs are taken from the existing displacement equations: within one block,
   the terms whose node carries a pressure dof (i.e. belongs to a U4 element)
   are the master/slave pair; the reference node is in no element and so is
   excluded automatically. Each pair is emitted once, not once per dof.

3. SOLVER=SPOOLES. The pressure block enters negative, so the assembled system
   is symmetric indefinite and the incomplete-Cholesky PCG this repository
   defaults to cannot be used on it.
"""
import os
import sys


def is_kw(line):
    return line.startswith('*') and not line.startswith('**')


def kwname(line):
    return line.split(',')[0].strip().upper().replace(' ', '')


def main():
    argv = sys.argv[1:]
    elset = 'Sphere_Only'
    if '--elset' in argv:
        i = argv.index('--elset')
        elset = argv[i + 1]
        del argv[i:i + 2]
    src, dst = argv[0], argv[1]

    lines = open(src).read().split('\n')

    # --- pass 1: retype the target elset and collect its nodes -------------
    out = []
    u4nodes = set()
    n_u4 = 0
    i = 0
    while i < len(lines):
        ln = lines[i]
        if is_kw(ln) and kwname(ln) == '*ELEMENT':
            up = ln.upper()
            if ('ELSET=' + elset).upper() in up.replace(' ', ''):
                out.append('*USER ELEMENT,TYPE=U4,NODES=4,'
                           'INTEGRATIONPOINTS=1,MAXDOF=4')
                out.append(ln.replace('C3D4', 'U4').replace('c3d4', 'U4'))
                i += 1
                while i < len(lines) and not is_kw(lines[i]):
                    f = lines[i].split(',')
                    if len(f) >= 5:
                        for x in f[1:5]:
                            u4nodes.add(int(x))
                        n_u4 += 1
                    out.append(lines[i])
                    i += 1
                continue
        out.append(ln)
        i += 1

    if not n_u4:
        raise SystemExit('u4ify: no C3D4 elements found in ELSET=%s' % elset)

    # --- pass 2: harvest periodic pairs from the displacement equations ----
    pairs, seen = [], set()
    i = 0
    while i < len(out):
        if is_kw(out[i]) and kwname(out[i]) == '*EQUATION':
            i += 1
            if i >= len(out):
                break
            try:
                nterm = int(out[i].split(',')[0])
            except (ValueError, IndexError):
                continue
            i += 1
            vals = []
            while i < len(out) and len(vals) < 3 * nterm:
                if is_kw(out[i]):
                    break
                vals += [v.strip() for v in out[i].split(',') if v.strip()]
                i += 1
            terms = [(int(vals[k]), int(vals[k + 1]), float(vals[k + 2]))
                     for k in range(0, min(len(vals), 3 * nterm), 3)]
            hit = [(n, c) for n, d, c in terms if n in u4nodes]
            if len(hit) == 2:
                (a, ca), (b, cb) = hit
                if ca * cb < 0:
                    key = (min(a, b), max(a, b))
                    if key not in seen:
                        seen.add(key)
                        pairs.append((a, b, ca, cb))
            continue
        i += 1

    # --- pass 3: emit, inserting the pressure equations before *STEP -------
    final, done = [], False
    for ln in out:
        if not done and is_kw(ln) and kwname(ln) == '*STEP':
            final.append('** SPAX: periodic pressure continuity, p(m)=p(s).')
            final.append('** %d node pairs on the periodic faces carry a U4 '
                         'pressure dof.' % len(pairs))
            for a, b, ca, cb in pairs:
                final.append('*EQUATION')
                final.append('2')
                final.append('%d, 4, %s' % (a, repr(ca)))
                final.append('%d, 4, %s' % (b, repr(cb)))
            done = True
        if is_kw(ln) and kwname(ln) == '*STATIC':
            ln = '*STATIC, SOLVER=' + os.environ.get('SPAX_U4_SOLVER', 'PARDISO')
        final.append(ln)

    open(dst, 'w').write('\n'.join(final))
    print('u4ify: %s -> %s' % (src, dst))
    print('  %d elements retyped to U4 in ELSET=%s' % (n_u4, elset))
    print('  %d nodes carry a pressure dof' % len(u4nodes))
    print('  %d periodic pressure equations added' % len(pairs))
    print('  solver set to %s (the mixed system is indefinite)'
          % os.environ.get('SPAX_U4_SOLVER', 'PARDISO'))


if __name__ == '__main__':
    main()
