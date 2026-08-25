# U4 patch tests

Both prescribe a uniform strain state on **every** node of a single tet, so the
exact answer is a constant stress and any consistent element must reproduce it
to roundoff. Run with `ccx_u4 u4test` / `ccx_u4 u4brine` and read the `.dat`.

| Deck | Geometry | Material | Point |
|---|---|---|---|
| `u4test.inp` | unit tet on the axes | E = 9.43e9, ν = 0.33 | compressible baseline |
| `u4brine.inp` | distorted, off-axis | K = 2.2 GPa, G = 0.44 MPa, ν = 0.4999 | the real case, and it exercises the global-derivative path |

Results, against closed form:

| | σxx | σyy = σzz | shear |
|---|---|---|---|
| `u4test` analytic | 1.397192e+07 | 6.881690e+06 | 0 |
| `u4test` U4 | **1.397192E+07** | **6.881690E+06** | ≤ 2e-10 |
| `u4brine` analytic | 2.200587e+06 | 2.199707e+06 | 0 |
| `u4brine` U4 | **2.200587E+06** | **2.199707E+06** | ≤ 7e-12 |

Exact at all 15 integration points in both.

`u4brine` is the one that matters twice over. It is at the brine's own Poisson
ratio, where a displacement element locks; and it is deliberately not
axis-aligned, because on a unit axis-aligned tet `dN/dxi` equals `dN/dx` and the
test cannot see whether the shape-function derivatives were converted to global
coordinates. They were not, at first — `shape4tet` returns at `iflag=2` before
applying the inverse Jacobian, and only `iflag=3` gives global derivatives.
A unit-tet-only test would have shipped that bug.

**What these do not test:** the inf-sup behaviour. A patch test is satisfied by
plenty of unstable elements — it says the element is consistent, not that it is
stable. Stability is what MINI's bubble is for, and demonstrating it needs a
real constrained problem, not one tet.

# Volumetric locking: `locking_sweep.sh`

The patch tests establish consistency. This one establishes that the bubble is
doing its job — the thing a patch test cannot see, because unstable elements
pass patch tests too.

A tip-loaded cantilever of linear tets, swept in Poisson's ratio **with the
shear modulus held fixed** at G = 1 MPa. Holding G is what makes it clean: the
true compliance then goes as 1/E = 1/(2G(1+ν)), so the tip deflection falls by
only 13 % from ν = 0.3 to the incompressible limit. An element that locks
volumetrically stiffens without bound instead. Normalising by the
Euler–Bernoulli tip deflection divides out the fixed discretisation error of a
coarse linear-tet beam and leaves only the ν-dependence.

16×3×3:

| ν | C3D4 tip | vs EB | U4 tip | vs EB |
|---|---|---|---|---|
| 0.30 | 8.663e-01 | 0.5631 | 9.693e-01 | 0.6300 |
| 0.45 | 5.607e-01 | 0.4065 | 8.490e-01 | 0.6155 |
| 0.49 | 2.411e-01 | 0.1796 | 8.017e-01 | 0.5973 |
| 0.499 | 5.142e-02 | 0.0385 | 7.878e-01 | 0.5905 |
| 0.4999 | 2.361e-02 | 0.0177 | 7.863e-01 | 0.5897 |
| 0.49999 | 2.066e-02 | 0.0155 | 7.862e-01 | **0.5896** |

**C3D4 collapses by 36×.** Textbook volumetric locking. **U4 drifts 6 % and
settles** — 0.5905, 0.5897, 0.5896 — converging to a finite incompressible
limit, which is the defining behaviour of a locking-free element. At ν =
0.49999 U4 is 38× softer than C3D4, i.e. C3D4 is 38× too stiff.

## Refinement at fixed ν = 0.4999

Locking is the discretisation pathology refinement does *not* cure — the
inf-sup constant does not improve with h — so this separates a locking element
from a merely coarse one.

| mesh | C3D4 tip | U4 tip | U4/C3D4 |
|---|---|---|---|
| 8×2×2 | 2.813e-02 | 3.999e-01 | 14.2× |
| 16×3×3 | 2.361e-02 | 7.863e-01 | 33.3× |
| 24×4×4 | 2.428e-02 | 9.936e-01 | 40.9× |
| 32×6×6 | 3.339e-02 | 1.119e+00 | 33.5× |

**C3D4 does not move** across a 4× refinement — 0.028 → 0.033, stuck at roughly
1/40 of the answer. **U4 climbs monotonically**, 0.400 → 1.119, converging.
(U4 is still rising at the finest mesh: linear tets also *shear*-lock in
bending, and that part does cure with h. The volumetric pathology is the one
being tested here and it is gone.)

This is the same signature the campaign cells show. On the undrained layered
cell, refining 0.0240 → 0.0120 moves Abaqus C3D4H by −10.55 % and CalculiX
C3D4 by −1.06 %: the displacement element is stuck for exactly the reason
demonstrated here.
