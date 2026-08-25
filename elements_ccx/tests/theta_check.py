"""Volumetric-strain smoothness for the U5/U6 nodal B-bar element.

The displacement-method analogue of pressure_check.py. With U6 there is no
theta degree of freedom to read -- the patch stiffness is condensed into the
element -- so theta is recovered from the solved displacement field using the
same operator the element uses:

    theta_a = sum over the 1-ring of b_n^a . u_n
    b_n^a   = (1/V_a) sum over U5 elements at a of (V_e/4) grad(L_n)|_e

Smoothness is measured as each node's deviation from the mean of its
neighbours, relative to the field's own variation. The displacement components
are the on-mesh controls: a nodal average should be no rougher than the
displacement fluctuations it is built from.

    python3 theta_check.py <deck-ccx.inp> <deck-ccx.dat>

Needs *NODE PRINT (U) covering the U5 nodes.
"""
import math
import sys
from collections import defaultdict

import numpy as np


def read_deck(inp):
    co, els, mode = {}, [], None
    for ln in open(inp):
        s = ln.strip()
        if s.startswith('*'):
            u = s.upper().replace(' ', '')
            if u.startswith('*NODE') and 'PRINT' not in u and 'SET' not in u \
               and 'FILE' not in u and 'OUTPUT' not in u:
                mode = 'n'
            elif u.startswith('*ELEMENT') and 'TYPE=U5' in u:
                mode = 'e'
            else:
                mode = None
            continue
        f = [x.strip() for x in s.split(',') if x.strip()]
        if mode == 'n' and len(f) >= 4:
            co[int(f[0])] = np.array([float(f[1]), float(f[2]), float(f[3])])
        elif mode == 'e' and len(f) >= 5:
            els.append([int(x) for x in f[1:5]])
    return co, els


def read_disp(dat):
    out, ind = {}, False
    for line in open(dat):
        low = line.lower()
        if 'displacements' in low and 'for set' in low:
            ind = True
            continue
        if line.strip().startswith(('stresses', 'forces', 'volume')):
            ind = False
            continue
        if not ind:
            continue
        f = line.split()
        if len(f) in (4, 5):
            try:
                out[int(f[0])] = np.array([float(f[1]), float(f[2]), float(f[3])])
            except ValueError:
                pass
    return out


def osc(field, nodes, nbr):
    have = [n for n in nodes if n in field and nbr[n]]
    if len(have) < 10:
        return float('nan')
    v = [field[n] for n in have]
    m = sum(v) / len(v)
    rf = math.sqrt(sum((x - m) ** 2 for x in v) / len(v))
    d = []
    for n in have:
        nb = [field[k] for k in nbr[n] if k in field]
        if nb:
            d.append(field[n] - sum(nb) / len(nb))
    ro = math.sqrt(sum(x * x for x in d) / len(d))
    return ro / rf if rf else float('nan')


def main():
    inp, dat = sys.argv[1], sys.argv[2]
    co, els = read_deck(inp)
    u = read_disp(dat)
    Va, b, nbr, nodes = {}, {}, defaultdict(set), set()
    for e in els:
        p = np.array([co[n] for n in e])
        J = np.array([p[1] - p[0], p[2] - p[0], p[3] - p[0]]).T
        vol = abs(np.linalg.det(J)) / 6.0
        if vol <= 0:
            continue
        Ji = np.linalg.inv(J)
        gr = np.zeros((4, 3))
        gr[1:, :] = Ji
        gr[0, :] = -gr[1:, :].sum(axis=0)
        w = vol / 4.0
        nodes.update(e)
        for i in range(4):
            for j in range(4):
                if i != j:
                    nbr[e[i]].add(e[j])
        for a in e:
            Va[a] = Va.get(a, 0.0) + w
            ba = b.setdefault(a, {})
            for k, n in enumerate(e):
                ba[n] = ba.get(n, np.zeros(3)) + w * gr[k]
    theta = {}
    for a in b:
        t = 0.0
        ok = True
        for n, v in b[a].items():
            if n not in u:
                ok = False
                break
            t += float(v @ u[n])
        if ok:
            theta[a] = t / Va[a]
    ctrl = [osc({n: u[n][c] for n in u}, nodes, nbr) for c in range(3)]
    t = osc(theta, nodes, nbr)
    print('U5+U6 volumetric strain smoothness')
    for c, name in enumerate(('u_x', 'u_y', 'u_z')):
        print('  %-22s %7.3f' % (name + ' (control)', ctrl[c]))
    hi = max(x for x in ctrl if x == x)
    print('  %-22s %7.3f   %s'
          % ('theta (vol strain)', t,
             'OK -- within controls' if t <= hi else 'ABOVE CONTROLS'))
    print('  patches: %d' % len(theta))


if __name__ == '__main__':
    main()
