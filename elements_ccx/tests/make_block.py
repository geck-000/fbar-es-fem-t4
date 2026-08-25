"""Emit a CalculiX deck for the U5+U6 acceptance suite.

    make_block.py OUT MODE n nu [jitter] [phase]

MODE
  patch    every boundary node prescribed u = eps.X for a general eps (all six
           components non-zero).  The exact answer is a uniform stress state,
           so any consistent element must reproduce it to roundoff at every
           integration point -- and a general eps, on a DISTORTED mesh, is
           what catches a shape-function-gradient bug that an axis-aligned
           unit cube cannot see.
  oedo     confined compression: u_x prescribed on x=0,1, roller elsewhere.
           Exact C1111 = K + 4G/3 for a homogeneous cell.
  modes    completely free block, *FREQUENCY.  Expect exactly 6 zero
           eigenvalues (3 translations + 3 rotations) and nothing else near
           zero; a 7th is a spurious mechanism.

phase: 'one' (homogeneous soft) or 'two' (soft sphere r=0.3 in ice).
"""
import sys

import os

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import smoothing_proto as S                                  # noqa: E402

EPS = np.array([[1.0e-3, 3.0e-4, 2.0e-4],
                [3.0e-4, -5.0e-4, 1.5e-4],
                [2.0e-4, 1.5e-4, 7.0e-4]])                   # symmetric


def main():
    out, mode, n, nu = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
    jit = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
    phase = sys.argv[6] if len(sys.argv) > 6 else 'one'

    G = 4.4e5
    K = 2.0 * G * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
    E = 9.0 * K * G / (3.0 * K + G)
    Ei, ni = 9.37e9, 0.33

    two = phase == 'two'
    # Jitter INTERIOR nodes only.
    #
    # smoothing_proto.mesh_box displaces boundary nodes too -- it has to, so
    # that periodic image pairs keep matching.  These decks use Dirichlet
    # boundary conditions instead, selected by x <= tol / x >= 1-tol, so a
    # displaced boundary node is no longer ON the face and never gets
    # constrained.  The first run of the suite failed everything including
    # C3D4's own patch test (max rel err 2.36, confined C1111 9.9e6 against an
    # exact 2.2e8) because the block was essentially unrestrained.  C3D4
    # passing T1 and T2 is this harness's control; if it does not, the harness
    # is what is broken.
    #
    # Distorting the interior is what the patch test actually needs -- it is
    # what exercises the global shape-function gradients -- so nothing is lost
    # by leaving the faces flat.
    nodes, tets, mat = S.mesh_box(n, 0.0, 0.30, 0.0,
                                  geom=S.GEOM['sphere' if two else 'slab'])
    if jit:
        # QUALITY-CONTROLLED JITTER.
        #
        # A free 0.3/h jitter of the Freudenthal split makes slivers: at n=6 it
        # produced a tet of volume 5.0e-05 against a mean of 7.7e-04, 15x
        # small.  In exact arithmetic the affine patch solution is still
        # exact, but B ~ 1/h for a sliver, so at K/G = 500 the solver's
        # round-off is amplified into a stress error of order the stress
        # itself -- C3D4 read -2.17e+05 where the answer is +2.64e+05.  That
        # is the MESH failing the test, not the element, and it would have
        # been reported as an element defect.
        #
        # So the amplitude is backed off until the worst tet is at least
        # QMIN of the mean volume.  Distortion is the point of the test; ill
        # conditioning is not.
        qmin = float(os.environ.get('SPAX_BLOCK_QMIN', 0.15))
        rng = np.random.default_rng(11)
        tol0 = 1e-12
        inter = ~((nodes <= tol0) | (nodes >= 1.0 - tol0)).any(axis=1)
        d = (rng.random(nodes.shape) - 0.5) * 2.0 / n
        base = nodes.copy()
        amp = jit
        for _ in range(40):
            nodes = base.copy()
            nodes[inter] += amp * d[inter]
            p4 = nodes[tets]
            vv = np.linalg.det(p4[:, 1:, :] - p4[:, :1, :]) / 6.0
            if np.abs(vv).min() >= qmin * np.abs(vv).mean() and vv.min() > 0:
                break
            amp *= 0.8
        sys.stderr.write('  jitter %.3f -> %.3f (min/mean tet volume %.3f)\n'
                         % (jit, amp, np.abs(vv).min() / np.abs(vv).mean()))
    if not two:
        mat[:] = 1
    S.grads(nodes, tets)                       # fixes any inverted tets

    L = ['*NODE']
    for i, p in enumerate(nodes):
        L.append('%d, %.12e, %.12e, %.12e' % (i + 1, p[0], p[1], p[2]))
    soft = [e for e in range(len(tets)) if mat[e] == 1]
    hard = [e for e in range(len(tets)) if mat[e] == 0]
    L.append('*ELEMENT,TYPE=C3D4,ELSET=Sphere_Only')
    for e in soft:
        L.append('%d, %s' % (e + 1, ', '.join(str(x + 1) for x in tets[e])))
    if hard:
        L.append('*ELEMENT,TYPE=C3D4,ELSET=Matrix_Only')
        for e in hard:
            L.append('%d, %s' % (e + 1, ', '.join(str(x + 1) for x in tets[e])))
    L += ['*SOLID SECTION,ELSET=Sphere_Only,MATERIAL=Mat_Inclusion',
          '*MATERIAL,NAME=Mat_Inclusion', '*ELASTIC', '%.12e, %.12e' % (E, nu)]
    if hard:
        L += ['*SOLID SECTION,ELSET=Matrix_Only,MATERIAL=Mat_Matrix',
              '*MATERIAL,NAME=Mat_Matrix', '*ELASTIC',
              '%.12e, %.12e' % (Ei, ni)]

    # ccx does not define NALL -- cgx does.  Without this, every
    # *NODE PRINT,NSET=NALL below is answered with 'node set NALL does not
    # exist' and the .dat comes back with no displacements and no reactions
    # at all, silently.
    L += ['*NSET,NSET=NALL,GENERATE', '1, %d, 1' % len(nodes)]

    tol = 1e-9
    onb = np.flatnonzero(((nodes <= tol) | (nodes >= 1.0 - tol)).any(axis=1))

    if mode == 'modes':
        L += ['*STEP', '*FREQUENCY,SOLVER=SPOOLES', '12', '*NODE PRINT,NSET=NALL',
              'U', '*END STEP']
        L.insert(len(L) - 5, '*DENSITY') if False else None
        # a density is required for an eigenvalue analysis
        L.insert(L.index('*ELASTIC'), '*DENSITY')
        L.insert(L.index('*DENSITY') + 1, '9.2e2')
    elif mode == 'oedo':
        L.append('*BOUNDARY')
        for i, p in enumerate(nodes):
            if p[0] <= tol:
                L.append('%d, 1, 1, 0.0' % (i + 1))
            if p[0] >= 1.0 - tol:
                L.append('%d, 1, 1, %.12e' % (i + 1, 1.0e-3))
            if p[1] <= tol or p[1] >= 1.0 - tol:
                L.append('%d, 2, 2, 0.0' % (i + 1))
            if p[2] <= tol or p[2] >= 1.0 - tol:
                L.append('%d, 3, 3, 0.0' % (i + 1))
        L += ['*STEP', '*STATIC,SOLVER=PARDISO',
              '*EL PRINT,ELSET=Sphere_Only', 'S', '*NODE PRINT,NSET=NALL',
              'U,RF', '*END STEP']
    else:                                       # patch
        L.append('*BOUNDARY')
        for i in onb:
            u = EPS @ nodes[i]
            for c in range(3):
                L.append('%d, %d, %d, %.12e' % (i + 1, c + 1, c + 1, u[c]))
        L += ['*STEP', '*STATIC,SOLVER=PARDISO',
              '*EL PRINT,ELSET=Sphere_Only', 'S', '*END STEP']

    open(out, 'w').write('\n'.join(L) + '\n')
    # the analytic uniform stress for the patch test
    tr = EPS.trace()
    sig = 2.0 * G * (EPS - tr / 3.0 * np.eye(3)) + K * tr * np.eye(3)
    print('K=%.6e G=%.6e nu=%.6f KG=%.1f' % (K, G, nu, K / G))
    print('patch_sigma %.9e %.9e %.9e %.9e %.9e %.9e'
          % (sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[0, 2], sig[1, 2]))
    print('oedo_C1111 %.9e' % (K + 4.0 * G / 3.0))
    print('nodes %d tets %d soft %d' % (len(nodes), len(tets), len(soft)))


if __name__ == '__main__':
    main()
