"""Build a tet-meshed cantilever for the volumetric-locking test.

The diagnostic is a Poisson sweep with the SHEAR modulus held fixed. Then the
true solution converges to a finite incompressible limit as nu -> 1/2 (Euler-
Bernoulli compliance goes as 1/E = 1/(2G(1+nu)), a 13% change from nu=0.3 to
nu=0.5), while a volumetrically locking element stiffens without bound and its
tip deflection collapses toward zero. No analytic reference is needed: the
shape of delta(nu) is the answer.

    python3 make_beam.py <out.inp> <element> <nu> [nx ny nz]

element is C3D4 or U4.
"""
import sys


def build(nx, ny, nz, lx, ly, lz):
    nid = {}
    nodes = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                nid[(i, j, k)] = len(nodes) + 1
                nodes.append((i * lx / nx, j * ly / ny, k * lz / nz))
    # 6-tet subdivision of each hex, oriented positive
    hexmap = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
              (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
    corner = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    els = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = [nid[(i + a, j + b, k + c2)] for a, b, c2 in corner]
                for t in hexmap:
                    n = [c[t[0]], c[t[1]], c[t[2]], c[t[3]]]
                    p = [nodes[x - 1] for x in n]
                    v = 0.0
                    a1 = [p[1][d] - p[0][d] for d in range(3)]
                    a2 = [p[2][d] - p[0][d] for d in range(3)]
                    a3 = [p[3][d] - p[0][d] for d in range(3)]
                    v = (a1[0] * (a2[1] * a3[2] - a2[2] * a3[1])
                         - a1[1] * (a2[0] * a3[2] - a2[2] * a3[0])
                         + a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]))
                    if v < 0:
                        n[2], n[3] = n[3], n[2]
                    els.append(n)
    return nodes, els, nid


def main():
    out, etype, nu = sys.argv[1], sys.argv[2], float(sys.argv[3])
    nx, ny, nz = (int(x) for x in (sys.argv[4:7] or (20, 4, 4)))
    lx, ly, lz = 10.0, 1.0, 1.0
    G = 1.0e6                      # held fixed across the sweep
    E = 2.0 * G * (1.0 + nu)
    nodes, els, nid = build(nx, ny, nz, lx, ly, lz)

    fixed = [nid[(0, j, k)] for k in range(nz + 1) for j in range(ny + 1)]
    tip = [nid[(nx, j, k)] for k in range(nz + 1) for j in range(ny + 1)]
    ftot = 1.0e3

    L = ["*NODE"]
    for i, (x, y, z) in enumerate(nodes, 1):
        L.append("%d, %.12e, %.12e, %.12e" % (i, x, y, z))
    if etype == "U4":
        L.append("*USER ELEMENT,TYPE=U4,NODES=4,"
                 "INTEGRATIONPOINTS=1,MAXDOF=4")
    L.append("*ELEMENT,TYPE=%s,ELSET=EALL" % etype)
    for i, n in enumerate(els, 1):
        L.append("%d, %d, %d, %d, %d" % (i, n[0], n[1], n[2], n[3]))
    L += ["*NSET,NSET=NTIP"] + [",".join(str(x) for x in tip[i:i + 8])
                                for i in range(0, len(tip), 8)]
    L += ["*SOLID SECTION,ELSET=EALL,MATERIAL=M",
          "*MATERIAL,NAME=M", "*ELASTIC", "%.12e, %.12e" % (E, nu),
          "*STEP", "*STATIC,SOLVER=SPOOLES", "*BOUNDARY"]
    for n in fixed:
        L.append("%d,1,3,0." % n)
    L += ["*CLOAD"]
    for n in tip:
        L.append("%d,2,%.12e" % (n, ftot / len(tip)))
    L += ["*NODE PRINT,NSET=NTIP", "U", "*END STEP"]
    open(out, "w").write("\n".join(L) + "\n")

    I = ly * lz ** 3 / 12.0
    print("%s nu=%.5f  E=%.6e  nodes=%d els=%d  EB_tip=%.6e"
          % (etype, nu, E, len(nodes), len(els),
             ftot * lx ** 3 / (3.0 * E * I)))


if __name__ == "__main__":
    main()
