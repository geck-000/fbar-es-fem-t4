"""Checkerboard diagnostic for the U4 pressure field.

An under-stabilised mixed element does not usually announce itself in the
homogenised modulus -- it shows up as a pressure field oscillating node to
node. This reads the nodal pressure (dof 4) out of a solved U4 deck and
measures that oscillation directly:

    jump ratio = RMS(p_i - p_j) over mesh edges  /  RMS(p - mean p)

A smooth field has a jump ratio well below 1: neighbours differ by much less
than the field varies overall. A checkerboard has neighbours differing by as
much as the whole field, so the ratio approaches (and can exceed) 1.

    python3 pressure_check.py <deck-ccx.inp> <deck-ccx.dat> [label]
"""
import sys
import math


def read_edges(inp):
    """Element edges of the U4 (and C3D4) tets in the deck."""
    edges, nodes = set(), set()
    inblock = False
    for line in open(inp):
        s = line.strip()
        if s.startswith('*'):
            u = s.upper().replace(' ', '')
            inblock = u.startswith('*ELEMENT') and ('TYPE=U4' in u)
            continue
        if not inblock or not s:
            continue
        f = [x.strip() for x in s.split(',') if x.strip()]
        if len(f) < 5:
            continue
        n = [int(x) for x in f[1:5]]
        nodes.update(n)
        for a in range(4):
            for b in range(a + 1, 4):
                edges.add((min(n[a], n[b]), max(n[a], n[b])))
    return edges, nodes


def read_pressure(dat):
    """{node: pressure} from the 4th component of a nodal displacement block."""
    p, indisp = {}, False
    for line in open(dat):
        low = line.lower()
        if 'displacements' in low and 'for set' in low:
            indisp = True
            continue
        if line.strip().startswith(('stresses', 'forces', 'volume')):
            indisp = False
            continue
        if not indisp:
            continue
        f = line.split()
        if len(f) == 5:
            try:
                p[int(f[0])] = float(f[4])
            except ValueError:
                pass
    return p


def main():
    inp, dat = sys.argv[1], sys.argv[2]
    label = sys.argv[3] if len(sys.argv) > 3 else inp
    edges, nodes = read_edges(inp)
    p = read_pressure(dat)
    have = [n for n in nodes if n in p]
    if len(have) < 10:
        raise SystemExit('%s: only %d brine nodes have a printed pressure -- '
                         'add *NODE PRINT for the brine node set'
                         % (label, len(have)))
    vals = [p[n] for n in have]
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    rms_field = math.sqrt(var)
    js, nj = 0.0, 0
    for a, b in edges:
        if a in p and b in p:
            js += (p[a] - p[b]) ** 2
            nj += 1
    rms_jump = math.sqrt(js / nj) if nj else float('nan')
    print('%-12s nodes=%-7d edges=%-8d mean p=%12.5e' % (label, len(have), nj, mean))
    print('%-12s RMS(p-mean)=%11.4e  RMS(edge jump)=%11.4e  ratio=%6.3f'
          % ('', rms_field, rms_jump, rms_jump / rms_field if rms_field else float('nan')))


if __name__ == '__main__':
    main()
