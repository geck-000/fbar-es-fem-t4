# New CalculiX elements

Fortran sources added to `ccx` to give it the mixed displacement/pressure
element it does not ship — the one Abaqus calls `C3D4H`. Copy into
`<ccx_2.23>/src`, add them to `SCCXF` in `Makefile.inc`, apply the dispatcher
wiring in `../patches_ccx/`, and rebuild.

| File | Role |
|---|---|
| `u4mat.f` | the element operator: builds and bubble-condenses the mixed blocks, in closed form |
| `e_c3d_u4.f` | turns that operator into the element stiffness (`mafillsm.f` calls it) |
| `resultsmech_u4.f` | turns the *same* operator into stresses and internal forces (`resultsmech.f` calls it) |

Both consumers go through `u4mat`, deliberately. Splitting them is what leaves
`../patches_ccx/0002-bbar-mean-dilatation.patch` reporting an `equilibrium_gap`
of 1.0: patch the stiffness, recover the forces with the unpatched operator,
and the reaction cross-check — the only reference-free convergence evidence
this repository has — stops meaning anything.

## Using it

```
*USER ELEMENT,TYPE=U4,NODES=4,INTEGRATIONPOINTS=1,MAXDOF=4
*ELEMENT,TYPE=U4,ELSET=Brine
...
*STATIC,SOLVER=SPOOLES
```

`MAXDOF=4` is what makes it work at all: `allocation.f:1110` raises `mi(2)` to
it, which gives `mastruct.c` room for a fourth nodal DOF. DOFs are node-major —
1,2,3 displacement, 4 pressure.

**`SOLVER=SPOOLES` is not optional in this build.** The pressure block enters
negative, so the assembled system is symmetric *indefinite*, and
incomplete-Cholesky PCG — this repository's production solver, and the only one
whose memory fits a large cell — requires positive definiteness.

And SPOOLES is worse on a saddle point than its SPD numbers suggest, because
pivoting drives fill-in: the campaign cell at 218k equations passed **5.9 GB
without finishing**, against 1.28 GB for a 167k-equation SPD system in
`../calculix/README.md`. This is the binding constraint on using U4 at
production scale, not the element.

## Getting an efficient solver

Three routes, in increasing order of effort.

**1. Compile in PARDISO. ccx already has the interface and it is already set up
for exactly this matrix.** `src/pardiso.c` sets `mtype=-2` — MKL PARDISO's code
for *real symmetric indefinite* — and `src/pastix.c` exists beside it. Both are
guarded by `-DPARDISO` / `-DPASTIX`, and this build sets neither:

```
$ grep -oE '\-D[A-Z]+' Makefile | sort -u
-DARCH -DARPACK -DMATRIXSTORAGE -DNETWORKOUT -DSPOOLES
```

So no source needs writing — install MKL (it is not on this machine;
`ldconfig -p | grep mkl` is empty), add `-DPARDISO` and the MKL libraries to
the Makefile, rebuild, and `*STATIC,SOLVER=PARDISO` handles U4 with proper
indefinite ordering and threading. **This is the cheapest path by a wide
margin** and should be tried before anything below.

**2. Eliminate the pressure and recover a positive-definite system.** Because
`A_lb = 0` (see `u4mat.f`), the assembled system is exactly

```
[ A    Bᵀ ]        C' = C + B_b A_bb⁻¹ B_bᵀ
[ B   −C' ]
```

`C` is SPD for any finite `K` and the added term is PSD, so `C'` is invertible
and the Schur complement

```
A_eff = A + Bᵀ C'⁻¹ B
```

is symmetric **positive definite** — incomplete-Cholesky PCG applies, memory
returns to the ~370 MB regime, and there are no extra global DOFs. Lumping `C'`
to its row sums makes `C'⁻¹` diagonal and the assembly direct.

The catch is sparsity. `B` couples a pressure node to the displacement nodes of
its own element, so `Bᵀ C'⁻¹ B` couples displacement nodes that share a
*pressure node* — the node patch, wider than the element stencil. `mastruct.c`
has to build that larger pattern. That is the same structural change
`../calculix/README.md` identifies for nodal-averaged B-bar, and the two turn
out to be the same piece of work.

**3. Keep the mixed form and precondition it properly** (block-diagonal or
Uzawa). Most code, least reuse of what ccx has. Not recommended while route 1
is uncompiled.

## Closed-form integration, and why A_lb vanishes

`u4mat.f` uses no quadrature. Every integral is a monomial in the barycentric
coordinates over a straight tet, where
`∫ L₁^a L₂^b L₃^c L₄^d dV = 6V·a!b!c!d!/(a+b+c+d+3)!`, so with `gᵢ = ∇Lᵢ`
constant and `Σᵢ gᵢ = 0`:

| integral | value |
|---|---|
| `∫ Lᵢ dV` | `V/4` |
| `∫ Lᵢ Lⱼ dV` | `V/10` (i=j), `V/20` (i≠j) |
| `∫ Lᵢ ∂b/∂x_d dV` | `−(256V/840)·gᵢ[d]` |
| `∫ ∂b/∂x_c ∂b/∂x_d dV` | `(256²V/15120)·Σᵢ gᵢ[c]gᵢ[d]` |

**The linear–bubble block `A_lb` is exactly zero.** The bubble vanishes on every
face, so `∫∇b dV = ∮ b n dS = 0`, and every linear–bubble deviatoric term
carries a factor of it. Checked against the earlier 15-point implementation:
`max|A_lb| / max|A_ll| = 7e-16`, and the two agree on `A_ll` to 13 digits.

That collapses the condensation to one line — the bubble's *only* effect is a
stabilisation added to the pressure block:

```
A_ll' = A_ll ,   B_l' = B_l ,   C' = C + B_b A_bb⁻¹ B_bᵀ
```

which is the standard characterisation of MINI as P1/P1 plus a parameter-free
stabilisation, and it is what makes route 2 above possible.

Two practical consequences. The element declares **`INTEGRATIONPOINTS=1`**,
which matters more than it sounds: ccx sizes `sti`, `eme`, `xstiff` and `stx`
at `mi(1)` for *every* element in the model, so a 15-point brine element taxed
the whole mesh — **1.74 GB against 0.12 GB** on a 345k-element cell. And the
single output point carries the *exact* element volume average, because the
bubble contributes nothing to the mean strain; that removes a real trap, since
ccx's `.dat` reader collapses integration points by arithmetic mean, which
equals the volume average only for equal quadrature weights — true for C3D4 and
C3D10, false for the 15-point tet rule.

## Why MINI (P1⊕bubble/P1) and not the P1/P0 that Abaqus documents

`C3D4H` is documented as a linear tet with *constant* pressure, and the obvious
implementation is to make that pressure element-internal and statically
condense it, which needs no new global DOF and no solver change. That is the
recipe in the note this work was checked against, and **it produces an element
bit-for-bit identical to plain `C3D4`.**

The reason is short. A linear tet has one integration point, so the strain is
constant over the element and an element-constant pressure represents the
divergence *exactly*. Condensing it is then algebraically the mean-dilatation
B-bar operator, and B-bar on a one-point tet is the identity. This is not a
prediction: `../patches_ccx/0002-bbar-mean-dilatation.patch` implements exactly
that operator and measures `E_x` on an undrained layered cell as bit-identical
with it on and off.

So a faithful P1/P0 clone cannot close any gap. Making a linear tet
locking-free requires enriching it:

* **P1/P1 unenriched** violates the inf-sup condition — checkerboard pressure.
* **MINI (P1⊕bubble/P1)** restores inf-sup parameter-free. The bubble is
  element-internal and condenses locally; the pressure stays nodal, continuous
  and global. Implemented here.
* **P1/P1 with Brezzi–Pitkäranta stabilisation** would also work and is less
  code, but introduces a mesh-dependent tuning constant — unwanted when the
  whole purpose is validating against Abaqus.

The cost of the choice is exactly the thing the condensation route was avoiding:
a global pressure DOF, hence an indefinite system, hence SPOOLES.

## What is verified, and what is not

**Verified: the patch test passes exactly, at brine incompressibility.** See
`tests/`. A uniform strain state prescribed on every node of one tet has a
constant-stress exact solution, and U4 reproduces it to roundoff at all 15
integration points:

| | σxx | σyy = σzz | shear |
|---|---|---|---|
| E = 9.43e9, ν = 0.33, unit tet — analytic | 1.397192e+07 | 6.881690e+06 | 0 |
| U4 | **1.397192E+07** | **6.881690E+06** | ≤ 2e-10 |
| K = 2.2 GPa, G = 0.44 MPa, ν = 0.4999, distorted tet — analytic | 2.200587e+06 | 2.199707e+06 | 0 |
| U4 | **2.200587E+06** | **2.199707E+06** | ≤ 7e-12 |

The formulation is the standard mixed one, `σ = 2μ dev(ε) + p I` with
`div u − p/K = 0`, assembled as `[[A, Bᵀ],[B, −C]]`. λ never appears, which is
the point — it is the term that blows up as ν → ½.

**Not verified:**

* **Inf-sup stability.** A patch test is passed by plenty of unstable elements;
  it establishes consistency, not stability. The bubble is there to supply
  stability, and showing it needs a genuinely constrained problem.
* **Anything at RVE scale.** No layered cell has been run with U4, so nothing
  here yet says whether it closes the 8–12 % gap against Abaqus `C3D4H` that
  `../calculix/README.md` measures.
* **Periodic pressure coupling** — see below. On the layered decks this has to
  be solved before a U4 run means anything.

One bug worth recording because the obvious test misses it: `shape4tet` returns
at `iflag=2` *before* applying the inverse Jacobian, so `shp(1:3,i)` are
derivatives with respect to the local coordinates. Only `iflag=3` gives global
ones. On a unit axis-aligned tet the two coincide, so the first patch test
passed with the bug in place; the distorted tet in `tests/u4brine.inp` is what
exposes it.

## Known gap for the layered decks

The periodic `*EQUATION` constraints the generator writes cover DOFs 1–3 only.
In the spherical decks the brine is interior and that is fine. In the **layered**
decks the brine slab spans the cell and reaches the periodic faces, so its
pressure field needs periodic coupling too. Running U4 on those cells without
adding pressure equations across the face pairs would leave the pressure
unconstrained where it matters most.

# STATUS: U4 fails the Abaqus validation. Do not use it yet.

Verified: consistency (patch test exact at ν = 0.4999 on a distorted tet) and
absence of *shear*-regime locking (cantilever, C3D4 collapses 36×, U4 holds).

**Not correct on a layered RVE.** Against the campaign cell `LMESH_m0p0240`,
using the drained twin that is C3D4 in both codes:

| | R = E_x(und)/E_x(drn) | vs Abaqus |
|---|---|---|
| Abaqus C3D4H / C3D4 | 2.4701 (seed spread 0.61 %) | — |
| CalculiX C3D4 / C3D4 | 2.5094 | +1.59 % |
| **CalculiX U4 / C3D4** | **1.2751** | **−48.4 %** |

`equilibrium_gap` was 1.9e-07, so this is a clean solve of a wrong model.

**The diagnosis, and one wrong turn on the way to it.** The stabilisation is
439× the compressibility term in a single matrix entry, which looks like the
answer and is not: it annihilates uniform pressure *exactly* (row-sum ratio
0.99999999999999534), which is precisely what MINI should do. The element is
structurally right.

What is wrong is magnitude. The bubble penalises pressure *gradients* with a
coefficient ~h²/μ, derived for Stokes where μ is O(1). Relative to the physical
compressibility ~h³/K the ratio is

```
S/C  ~  0.05 K/mu        -- independent of h
```

and the brine has K/μ = 5000. The consequence, measured on the small cell:

| | E_x (across layers) |
|---|---|
| C3D4 undrained, K = 2.2 GPa | 6.349e9 |
| C3D10 undrained, K = 2.2 GPa | 4.738e9 |
| **U4 undrained, K = 2.2 GPa** | **3.600e9** |
| C3D4 **drained**, K = 2.2 MPa | 3.446e9 |

**U4's undrained answer is within 4.5 % of a genuinely drained cell.** A factor
of a thousand in bulk modulus contributes essentially nothing.

**Why the earlier tests missed it, which is the lesson worth keeping.** Both
passing tests are structurally blind to this:

* the **patch test** prescribes a uniform field, so the pressure is uniform,
  `B_bᵀp = 0`, the bubble amplitude is identically zero, and the stabilisation
  is not exercised at all;
* the **cantilever** is deviatoric-dominated, so the bulk modulus barely enters.

Consistency and shear behaviour were verified; the bulk response under a
*varying* pressure field was not, and that is the only regime the confined
brine occupies. A test that passes is not evidence for behaviour it cannot see.

**The remedy under test.** `CCX_U4_STAB=CAPPED` scales the stabilisation by
`θ = tr(C)/(tr(C)+tr(S))`, so it can never exceed the physical compressibility
it perturbs: `θ → 1` when the material is compressible (recovering textbook
MINI exactly), and `S_eff → C` when it is not. Parameter-free, and it leans on
the pressure mass term being nonzero — true here, because the brine is
ν = 0.49993 rather than exactly ½. Default remains unscaled MINI.

This is a deliberate departure from the textbook element and is **not yet
validated**. It must be measured against Abaqus C3D4H before use, exactly as
the unscaled version was.

# CAPPED also fails, for the opposite reason. Use C3D10.

`CCX_U4_STAB=CAPPED` was meant to stop the stabilisation swamping the physics.
It does — and it breaks stability instead.

| element | L_mesh | R | vs converged (1.9897) |
|---|---|---|---|
| Abaqus C3D4H | 0.0240 | 2.4701 | +24.1 % |
| Abaqus C3D4H | 0.0120 | 2.3786 | +19.5 % |
| CalculiX C3D4 | 0.0240 | 2.5094 | +26.1 % |
| CalculiX C3D4 | 0.0120 | 2.5852 | +29.9 % |
| **CalculiX C3D10** | **0.0120** | **2.1249** | **+6.8 %** |
| CalculiX U4 MINI | 0.0240 | 1.2751 | −35.9 % |
| CalculiX U4 CAPPED | 0.0240 | 2.2627 | +13.7 % |
| CalculiX U4 CAPPED | 0.0120 | 1.7020 | −14.5 % |

CAPPED passes *through* the right answer between the two meshes and keeps
softening. The coarse-mesh +13.7 % was a crossing, not agreement.

## The pressure field says why — and the first metric for it was wrong

The obvious checkerboard measure, RMS edge-to-edge jump over RMS field
variation, gives 0.991 for MINI and 1.281 for CAPPED, which looks conclusive
and is not. Controls on the same mesh: `u_y`, a physical displacement
fluctuation, scores **1.316** — higher than either pressure. On an unstructured
tet mesh that ratio is dominated by edge length, not oscillation.

The measure that does discriminate is each node's deviation from the mean of
its own neighbours, relative to the field variation:

| | u_x | u_y | u_z | **pressure** |
|---|---|---|---|---|
| MINI | 0.367 | 0.794 | 0.780 | **0.754** |
| CAPPED | 0.614 | 0.825 | 0.761 | **1.034** |

MINI's pressure (0.754) sits *inside* the range of the physical fluctuation
fields — no checkerboard, exactly as inf-sup stability promises. CAPPED's
(1.034) exceeds every displacement control on the same mesh and is 37 % above
MINI's. **CAPPED is under-stabilised**, which is why it is too soft and why it
gets worse under refinement.

So the two failures are complementary and neither is fixable by tuning:

* **MINI** — pressure stable, but the Stokes-scaled stabilisation destroys the
  brine's bulk stiffness (R = 1.28, within 4.5 % of a drained cell).
* **CAPPED** — bulk stiffness restored, stability lost.

## What to use instead

**Plain C3D10, and nothing else needs building.** At `L_mesh` = 0.0120 it is
+6.8 % from the converged answer, against Abaqus's own C3D4H at +19.5 % and
CalculiX C3D4 at +29.9 % on the same mesh. No new element, no MKL, no mixed
formulation — order 2 simply does not lock enough to matter here.

The cost is 8× the equations (601k → 4.9M on this cell) and it is not free, but
it is correct, and it is the recommendation until someone builds a properly
inf-sup stable element scaled for K/μ ~ 5000. `pressure_check.py` and the
neighbour-deviation measure are the acceptance test any such element must pass,
alongside the refinement sequence.

# Next: nodal-averaged B-bar, and why it escapes the trap

## Why U4 cannot be rescued

`S = B_b A_bb⁻¹ B_bᵀ` is *simultaneously* the stabilisation and the spurious
compliance. `CAPPED` scaled S down; `STIFFB` stiffened `A_bb`, which also
shrinks S. Two different-looking fixes, one identical lever:

| variant | u_x | u_y | u_z | pressure | R @0.0240 |
|---|---|---|---|---|---|
| MINI | 0.367 | 0.794 | 0.780 | **0.754** OK | 1.2751 |
| CAPPED | 0.614 | 0.825 | 0.761 | **1.034** above controls | 2.2627 |
| STIFFB | 0.653 | 0.834 | 0.781 | **1.061** above controls | 2.3241 |

Any pressure stabilisation for a P1/P1-type space scales as h²/μ on dimensional
grounds, while the physical compressibility scales as h²/K. At K/μ = 5000 the
stabilisation exceeds the physics by ~K/μ **regardless of mesh**. Shrink it to
protect the bulk stiffness and you fall below inf-sup. No window exists at this
contrast — a property of the linear tet, not of this implementation, and why
Abaqus's own C3D4H is only +19.5 % at `L_mesh` = 0.0120.

## Why nodal-averaged B-bar is different

Stability comes from **patch averaging**, a geometric construction, not a
1/μ-scaled penalty, so the scaling argument does not apply. It carries **no
pressure DOF**: a pure displacement method, so the matrix stays symmetric
**positive definite** and incomplete-Cholesky PCG works. PARDISO becomes an
option rather than a requirement.

```
V_a = sum over elements at node a of V_e/4            (nodal volume)
θ_a = (1/V_a) sum over those elements of (V_e/4) div(u)|_e
K   = K_dev (element-wise, as now) + K_vol (built from θ)
```

## Cost, measured on `LMESH_m0p0120` (1 211 410 elements, 214 539 nodes)

| | equations | nnz |
|---|---|---|
| C3D4 | 601k | 26.5 M |
| **nodal B-bar** | **601k** | **126.5 M** (4.8× C3D4) |
| C3D10 | 4.9M | 343 M |

2.7× fewer nonzeros than C3D10 and 8× fewer equations. Elements per node 22.6;
1-ring 14.7 mean / **27 max**; 2-ring 70.1 mean / 132 max.

## Implementation route — no `mastruct.c` change needed

The objection to nodal B-bar is that `K_vol` couples a whole patch, widening
the stencil past the element graph. Not needed here: `*USER ELEMENT` accepts up
to **255 nodes** (`userelements.f:83`) and `mastruct.c:137` reads the node count
from the element label, so a patch *is* expressible as an element. Our worst
patch is 27 nodes.

Two pure-displacement user elements, both `MAXDOF=3`:

| type | nodes | role |
|---|---|---|
| `U5` | 4 | linear tet, **deviatoric only** — C3D4 minus its volumetric term |
| `U6` | the 1-ring of one node (≤27 here) | that node's volumetric patch stiffness |

One `U6` per mesh node, one `U5` per tet. Every piece of plumbing already
exists and is exercised: label-driven `nope`, `mafillsm` dispatch,
`resultsmech_u` recovery, and the `printoutelem.f` volume fix. A generator
computes the 1-rings and writes the `U6` connectivity, as `u4ify.py` writes the
periodic pressure equations.

Acceptance tests unchanged — they caught every failure so far: patch test
exact, refinement tracking Abaqus toward R = 1.9897, and the *volumetric strain*
field checked for oscillation against the displacement controls on the same
mesh, the way `pressure_check.py` checks pressure.

# U5 + nodal-averaged B-bar

The route that replaces U4. Two pieces, and only one of them is an element:

| file | role |
|---|---|
| `e_c3d_u5.f`, `resultsmech_u5.f` | `U5` — linear tet carrying the **deviatoric** stiffness only |
| `nodalbbar.py` | builds the volumetric half out of existing ccx features |

```
K = K_dev   per tet (U5)
  + K_vol   per node: theta_a tied to the surrounding displacements by an
            *EQUATION, given energy by a grounded SPRING1
```

The bulk modulus never enters an element-local operator, so it cannot lock.
Stability comes from averaging over the node patch — geometry, with no 1/μ
scaling anywhere — which is exactly what U4 could not have: there the
stabilisation and the spurious compliance were the same term.

And there is **no pressure unknown**, so the matrix stays symmetric positive
definite: `SOLVER=ITERATIVE CHOLESKY` works, PARDISO optional rather than
mandatory.

## Why no patch element was needed

`U6` was going to be a patch element — one per node, spanning its 1-ring.
`*USER SECTION` can attach constants to user elements, but only the *same*
constants to a whole set, so per-patch data would have meant ~36 000
single-element sections. Carrying θ as one extra DOF on a dummy node instead
needs no new element at all, and ccx's MPC cascade handles the widened stencil
that `mastruct.c` would otherwise have to build.

Two details that make it practical:

* **Every spring is identical.** With `t_a = √(K·V_a)·θ_a` the energy is exactly
  `½ t_a²`, so one `*SPRING` of unit stiffness covers every patch and `V_a`
  lives in the `*EQUATION` coefficients. Otherwise each node needs its own
  element set.
* **The θ node's unused DOFs must be constrained.** It carries three DOFs and
  only DOF 1 is used; left free the matrix is singular and SPOOLES *stops
  silently* after "Factoring the system of equations" — no error, no
  `Job finished`, an empty `.dat`.

## Stress output uses a different operator from the stiffness, deliberately

`resultsmech_u5` reports `σ = 2μ dev(ε) + K div(u)` from the **element's** own
divergence, while the stiffness contains no volumetric term at all. That is
exact for the homogenisation because θ_a is the V-weighted mean of the
surrounding element divergences and each element feeds four nodes:

```
sum_a V_a theta_a  =  sum_e V_e div(u)|_e
```

so the volume-*integrated* volumetric stress is identical computed either way.
Internal forces stay deviatoric, matching the stiffness; the springs supply the
rest of the reaction, so `equilibrium_gap` stays a real check rather than the
~1.0 it reads under the B-bar patch.

## Verification

| test | U5 + nodal B-bar | for contrast |
|---|---|---|
| patch test, ν = 0.33, unit tet | **exact** 1.397192E+07 | — |
| patch test, ν = 0.4999, distorted | **exact** 2.200587E+06 | — |
| locking sweep, vs EB at ν = 0.49999 | **0.5632**, settling | C3D4 0.0155 |
| volumetric strain smoothness | **0.456**, controls 0.412/0.620/0.608 | U4 pressure: MINI 0.754, CAPPED 1.034, STIFFB 1.061 |

θ is smoother than two of the three displacement controls, which is what a
patch average should be. This is the measure that failed both U4 repairs.

**Still open:** the layered RVE against Abaqus at both mesh levels. CAPPED
cleared the coarse point and then shot past the converged answer under
refinement, so a single good R proves nothing — the 0.0120 point is the test.
Also untested: the phase interface. Only the brine is split into U5 + patches;
the ice stays ordinary C3D4, and patches at interface nodes cover brine
elements only, so no 1000× modulus contrast is smeared. That is a design
assertion, not yet a measurement, and a bad interface treatment would show up
as an offset in R rather than in any of the four tests above.

## STATUS: the formulation is right, the ccx delivery is not

The RVE match looked perfect and is not to be believed. On `LMESH_m0p0240`,
`R = 2.4677` against Abaqus's 2.4701 — 0.10 %, inside the 0.61 % seed spread —
but `equilibrium_gap = 6.93e-01` and `E_z` came out 3.6 % *stiffer* than C3D4,
which replacing volumetric stiffness with an average can never do.

A confined-compression block, where the answer is closed form
(`M = K + 4G/3`), isolates it:

| mesh | free displacement DOFs | reaction / exact |
|---|---|---|
| 1×1×1 (6 tets, 8 nodes) | none | **1.0000** |
| 2×2×2 (48 tets, 27 nodes) | few | 1.5485 |
| 4×4×4 (384 tets, 125 nodes) | many | 3.0515 |

The error appears **only when displacement DOFs are free** and grows with the
interior-node count.

**Every ingredient is exact in isolation:**

| checked | result |
|---|---|
| U5 deviatoric stiffness, confined block | reaction ratio **1.0000** |
| nodal volumes | Σ V_a = mesh volume exactly |
| divergence operator vs an analytic field | θ_a = 1e-4 = div(u) at every node |
| `SPRING1` convention | k = 1, u = 1 → RF = 1.0 |
| **the whole operator, assembled in Python** | **reaction ratio 1.0000** |

That last row is the one that settles it. The same `K_dev + Σ_a K V_a b⊗b`
assembled directly reproduces the closed form exactly; routed through ccx's
`*EQUATION` + `SPRING1` it does not. **Nodal B-bar is sound here — the delivery
mechanism is at fault.**

It is also provably a bug rather than a limitation: nodal averaging satisfies
`Σ_a V_a θ_a² ≤ Σ_e V_e θ_e²` by Cauchy–Schwarz, so this method can only ever be
*softer* than element-wise volumetric stiffness. Coming out 3× stiffer is
impossible for a correct assembly.

**Prime suspect**: the θ DOF carries the spring *and* is the dependent DOF of
its own `*EQUATION`, so its stiffness entry has both indices MPC-dependent.
`mafillsm.f` has a distinct branch for that case and it is the one path none of
the passing tests exercises.

**Routes out**, in order of preference:

1. Give the volumetric term to a real element so no MPC is involved — the
   original `U6` patch element. `*USER SECTION` blocks the obvious per-patch
   data route, but ccx's substructure path (`matrix2userelem.f`,
   `writesubmatrix.f`) reads an externally assembled stiffness matrix as a user
   element, and the Python assembly above already produces exactly that matrix.
2. Keep the MPC but move the spring off the dependent DOF.
3. Fix the both-DOFs-dependent branch in `mafillsm.f`, if that is genuinely
   where it is.

## What the verification ladder did and did not catch

Patch test, locking sweep and volumetric-strain smoothness (0.456, inside
controls) **all passed** on the broken assembly. The patch test cannot see it
because prescribing every displacement means the springs never carry load; the
other two are insensitive to an overall volumetric scale error.

Only **reaction against a closed-form modulus on a confined block** found it,
and that test should run before any RVE, not after. It is now
`tests/oedometer_m2.inp`.

### Root cause: overlapping MPCs, not the method

Every candidate was eliminated by test, not by argument:

| hypothesis | test | verdict |
|---|---|---|
| U5 stiffness wrong | confined block, deviatoric only | exact, 1.0000 |
| nodal volumes wrong | Σ V_a vs mesh volume | exact |
| divergence operator wrong | θ_a vs analytic field | exact, 1e-4 everywhere |
| emitted equations wrong | θ reconstructed from the deck | exact, 1e-4 everywhere |
| `SPRING1` convention wrong | k=1, u=1 | RF = 1.0 |
| MPC + spring condensation broken | 2 nodes, θ = 3u, k = 1, F = 1 | u = 1/9 **exact** |
| equations over-constrain u | equations with springs removed | deviatoric exactly, 1.0000 |
| **the operator itself** | **same operator assembled in Python** | **exact, 1.0000** |

The discriminator:

| patches | Python | ccx |
|---|---|---|
| 1 | 7.069748e+01 | **7.069748e+01** |
| 27 (all) | 2.200587e+05 | 3.407531e+05 |

**One patch is exact; the full set is 1.55× too stiff.** ccx handles a single
`*EQUATION` + `SPRING1` pair perfectly and degrades once many equations share
the same displacement DOFs — which is intrinsic here, because adjacent node
patches overlap by construction. Each mesh DOF appears as an independent term
in up to ~24 equations.

This cannot be fixed in the generator. The volumetric term has to reach the
matrix without going through ccx's MPC machinery.

### The fix: make the patch a real element after all

Give `K_vol^a = K V_a b⊗b` to a `U6` element spanning the patch's 1-ring —
the original design, which the `*USER SECTION` per-set-constants limit pushed
me away from. Two ways to supply the per-patch data:

1. **Compute it inside the element.** Make node 1 of the `U6` connectivity the
   patch centre and have `e_c3d_u6.f` build a node→element map once (cached,
   built in a serial phase to stay safe under the OpenMP assembly), then form
   `b` from the surrounding tets exactly as `nodalbbar.py` does. No deck
   plumbing at all.
2. **Feed the matrices in.** ccx's substructure path (`matrix2userelem.f`,
   `writesubmatrix.f`) reads an externally assembled stiffness matrix as a user
   element, and the Python assembly above already produces exactly these
   matrices — but one file per patch makes this impractical at 36 000 patches.

Route 1 is the one to build. `*USER ELEMENT` already permits 255 nodes
(worst patch here is 27), `mastruct.c` reads the count from the label, and U5
plus the dispatcher and volume fixes are all in place — the remaining work is
one element routine.

## Per-material patches: the design, and an unresolved deck-format failure

Applying U5+U6 to the brine alone leaves interface brine nodes with one-sided
patches, and in a slab 2-3 elements thick most brine nodes *are* interface
nodes -- so the averaging has little to average over. Measured: R goes
2.3843 -> 2.5602 under refinement, i.e. back to plain C3D4, while the *same
element* gains 8.9x -> 50.9x over C3D4 on a homogeneous cantilever as the mesh
refines. The element is fine; treating one phase is not.

**The ice needs patches too.** Locking follows isochoric *deformation*, not the
material's own ν: ice beside a near-incompressible slab is forced to deform
almost isochorically, so C3D4 ice locks despite ν = 0.33.

**But patches must not span the interface.** The stress recovery relies on
`Σ_a V_a θ_a = Σ_e V_e div(u)|_e`, which holds for ONE K. Across a 1000×
modulus jump it fails, and the reported stress stops matching the transmitted
force. So: one patch per (node, material). `u6patch` filters contributing tets
by the patch element's own material; the generator emits one `*USER ELEMENT`
type per (material, ring size).

### The deck-format failure: found and fixed

**Cause: ccx's built-in element dispatch claims type names containing
digits.** `elements.f:361` has `elseif(label(4:4).eq.'4') then nope=4`, `:358`
has `'10' -> nope=10`, `:320` has `'20' -> nope=20`. A user element named
`U614` therefore gets `nope=4` from the digit rule *before* the `*USER ELEMENT`
lookup runs, so ccx reads 4 nodes instead of 18 and treats the continuation
line as a fresh element card.

Only 3 of ~36 000 continuation-needing patches reported an error, because
`"element N is already defined"` fires only when the misread line's first field
collides with a real element id. The rest were corrupted silently.

The fix is to name types from **letters only** (`U6AA`…`U6ZZ`, 676 available
against 58 needed), which avoids every digit-keyed rule.

Ruled out on the way, all by test rather than argument: a missing or extra
trailing comma on continued cards; more than 16 fields overflowing
`textpart(16)`; interleaving `*USER ELEMENT` with `*ELEMENT` across a
keyword-chain boundary; a cap on the number of user element types (`nuel_` is
dynamic); the deck being internally malformed (a validator found 0 mismatches
over all 93 904 elements); and continuation being broken in general (a minimal
18-node deck whose continuation line deliberately starts with a colliding
element id parses correctly).

A second real bug found alongside it: patch element ids started at 1e8, so
ccx's `ne = max(ne, id)` sized every element array for 100 million elements.
Ids now continue from the deck's real maximum.

### The equilibrium gap: NOT the constraints -- the internal forces vanish

Corrects an earlier conclusion in this file. With **no constraints of any kind**
-- the real RVE mesh, uniaxial displacement on the two x-faces, lateral faces
traction-free, so the average-stress theorem gives `<σxx>·A = F` exactly:

| | `<σxx>·A` | reaction | ratio |
|---|---|---|---|
| C3D4 | 5.8196e3 | 5.8179e3 | **1.000291** |
| U5+U6 | 5.1345e3 | **−1.2985e2** | **−39.5** |

The U5+U6 reaction is essentially zero while its volume-averaged stress is
sane. Internal forces are not reaching the constrained nodes. This reproduces
without periodic constraints, so everything ruled out earlier about MPCs,
three-term equations and 3D chaining was ruled out correctly but was never the
issue -- the failure needs only the real mesh and free DOFs.

What is now known:

| condition | result |
|---|---|
| structured block, free DOFs, any constraint form tested | **exact** |
| real RVE mesh, every node prescribed | stress **exact**, but no free DOFs to exercise the stiffness |
| real RVE mesh, free DOFs, no constraints | **reaction ≈ 0** |

So the discriminator is the *mesh*, not the boundary conditions. The structured
block has patches built from a regular 6-tets-per-hex subdivision; the RVE mesh
is unstructured, with rings from 4 to 32 nodes and a wide spread of element
sizes. A reaction near zero with a plausible stress means `fn` is being
accumulated for only a small part of the model.

### ROOT CAUSE: patches larger than 20 nodes overflow ccx's element matrix

`mafillsm.f:69` and `e_c3d_u.f:68` declare **`s(60,60)`** — 60 DOFs, i.e. 20
nodes at 3 DOF/node. U6 patch rings reach **32 nodes on the small cell and 37
on the campaign cell**, so `e_c3d_u6` writes `s(ii,jj)` with indices up to 96
into a 60×60 array. Fortran does not bounds-check: the writes land in other
columns or past the array, and `mafillsm` reads the same out-of-bounds entries
when assembling.

It accounts for every observation:

| observation | explanation |
|---|---|
| structured block exact under every BC | its rings are ≤15 nodes = 45 DOFs, inside the limit |
| real RVE mesh fails with free DOFs | rings to 32 nodes = 96 DOFs, past the limit |
| stress output exact under prescribed strain | stress comes from `u6patch`, never from `s` |
| stiffness and recovery return **bit-identical** `va`, `kva`, `b` | both call `u6patch`, which is correct — the corruption is downstream, in `s` |
| reaction ≈ 0 with a plausible `<σxx>` | the assembled stiffness is corrupted, so the solution is too soft, while U5 still reports `K_e·div(u)` at full stiffness |

Two things this clears, both of which looked like suspects and were not:

* **U5-only giving RF ≈ 1.1e-11 is physically correct**, not a missing
  assembly: deviatoric-only means `K = 0`, so `E = 9KG/(3K+G) = 0` and a
  laterally-free bar has genuinely zero axial stiffness.
* **`va = 2.8e-10` on a 4-node ring is correct too.** Its centre node sits in 7
  tets but only one of that material, so per-material patching is working as
  designed.

**The fix** is to enlarge the element matrix along the whole user-element path —
`s`, `sm` and `ff` in `mafillsm.f`, `e_c3d_u.f` and the U5/U6 routines — to at
least `3 × max_ring`. Capping the ring instead is not an option: dropping nodes
from a patch breaks `Σ_a V_a θ_a = Σ_e V_e div(u)|_e`, and the rank-one patch
stiffness cannot be split across two elements.

A cheap interim check, before touching ccx: regenerate with only patches of ≤20
nodes given to U6 and the rest left as plain C3D4. That is not the right
method, but if the equilibrium gap collapses it confirms the diagnosis end to
end.

### Two earlier claims to correct

* The shared-patch (both-phase) run gave `equilibrium_gap` 2.5e-3 at 0.0240 and
  2.6e-1 at 0.0120, and I attributed that to the K-identity breaking across the
  interface. **Those decks were also malformed** by the same continuation issue,
  so the attribution is unproven -- the algebra stands, the demonstration does
  not.
* `R = 2.1568` (+8.4 % from converged) from that run was reported as a large
  improvement before its equilibrium gap was checked. It is not trustworthy.

# STATUS 2026-08-23: the element is sound. Cross-code breadth is what is thin.

Two more delivery bugs found and fixed (`../patches_ccx/0008-*.patch`), then
the element re-tested from the bottom up. Every test with a definite right
answer now passes, including the two that killed the previous elements.

## The two bugs

1. **The 150-DOF capacity was still a guess, and still too small.** Rings of
   4..40 had been seen, so 50 nodes looked ample. A sphere-packing cell
   contains a node carried by **248 tets** -> a **127-node ring**, 381 DOF. The
   guard from 0006 caught it cleanly. The limit is now 765 = 3x255, which is
   not another guess: ccx stores a user element's node count in the single
   `lakon(8:8)` byte and `userelements.f` rejects `NODES > 255`, so this is the
   largest patch ccx can express at all.

2. **`e_c3d_u6` zeroed a fixed 60x60 block of a larger array.** The stiffness
   loop assigns the whole `3*nope` square so `s` survived, but the element
   never writes `ff`, and `mafillsm` reads `ffval = ffu(jj)` UNCONDITIONALLY,
   before it tests `rhsi`. For any patch over 20 nodes, `ff(61..3*nope)` was
   the previous element's leftovers. Layered rings reach 40 and sphere rings
   127, so every U6 number taken before this was re-run.

## Acceptance ladder, re-run against the fixed binary

| test | result | what it rules out |
|---|---|---|
| patch test, distorted mesh, general strain **with shear**, nu = 0.4999 | **exact**, all 6 components, all 162 ip | inconsistency |
| confined compression vs `M = K + 4G/3`, N = 1,2,4,8 | **1.0000** every level | the previous delivery gave 1.0000 / 1.5485 / **3.0515** here |
| **zero-energy modes**, structured and 0.5/0.8-jittered | **exactly 6**, every mesh | spurious mechanisms |
| U5 alone, same meshes | **exactly 7** | confirms the deviatoric split is exact: rigid body + uniform dilatation, and U6 restores that one mode and no more |
| softest genuine mode vs C3D4 | 0.976-0.978 structured, 0.54 jittered | over-softness (see below) |
| refinement vs C3D10, 4 meshes, frozen geometry | **monotone from above, no crossing** | the CAPPED failure |

**On over-softness, which was the live worry.** U6 agrees with C3D4 on the
softest genuine mode to **2.4 %** on structured meshes. Jittering the mesh
stiffens *C3D4's* by **82 %** while U6 moves 5 %. A physical mode does not
stiffen because nodes moved, so U6 departs from C3D4 only where C3D4 is
provably wrong. This is the measurement that separates "converges faster" from
"too soft", and the eigenvalue count is what makes it a measurement rather than
an inference.

## Refinement: the test that killed CAPPED

Frozen deterministic geometry (slabs and bridges only -- the sphere population
is NOT reproducible across a re-mesh, see the trap below), slab thickness
0.0125 matched to the campaign cell, C3D10 as the trusted reference.

| L_mesh | slab/mesh | C3D4 | U5+U6 brine | U5+U6 both |
|---|---|---|---|---|
| 0.024 | 0.52 | +2.97 % | +2.76 % | **+2.04 %** |
| 0.012 | 1.04 | +1.66 % | +1.10 % | **+0.70 %** |
| 0.008 | 1.56 | +1.13 % | +0.59 % | **+0.34 %** |
| 0.006 | 2.08 | +0.85 % | +0.38 % | **+0.22 %** |

Four refinements spanning the whole range of the Abaqus reference sequence.
Every arm approaches C3D10 from ABOVE and none crosses; U6's error halves each
refinement. CAPPED was already 14.5 % *below* by its second mesh. C3D10 itself
is settled -- last three meshes within 0.04 %.

## Two traps in the comparison methodology, both of which bit

**1. The drained twin is the mesh-equivalence CONTROL, not just a denominator.**
Both codes use plain C3D4 for the drained cell, so any disagreement there is
mesh, not element. It passes at 0.0240 (-0.07 %) and FAILS at 0.0120
(+3.58 %, against a 1.5 % seed spread): at the same nominal `L_mesh` the
CalculiX mesh is effectively coarser. Abaqus's drained value is flat from
0.0120 on (~2.33) while CalculiX's is still descending (2.517 -> 2.424). So the
0.0120 R comparison is not like-for-like in either direction, and the U6 vs
C3D4H disagreement there is unresolved, not resolved. **Gate every R row on the
drained check.**

**2. R's denominator must be plain C3D4 for every arm.** Abaqus's R is
`E_und(C3D4H) / E_drn(C3D4)` -- it substitutes the hybrid element in the
UNDRAINED cell only, so the element under test must be substituted only there
too. This is not a formality: patching both phases softens the ice ~4 %, and in
the drained cell the load runs through thin ice bridges, so `E_drn` moves 4.3 %
straight into R. Dividing by each arm's own drained value made the SAME element
read +3.49 % at b020 and +0.62 % at b040 purely from which denominator existed
yet -- a reporting artefact that looked exactly like a physical flip.

**Corollary: do not use the both-phase arm on drained cells.** It moves E_drn
to -5.30 % against Abaqus at b040, outside the 2.24 % floor, where plain C3D4
sits at -1.03 % inside it. Whether that softening is an improvement (linear
tets genuinely are stiff in thin ligaments) or an overshoot needs C3D10 on a
drained cell. Untested.

## The frozen-packing trap in `converge_u4.sh`

`SPAX_SAVE_PACKING` does not survive a re-mesh. At `L_mesh` 0.010/0.008/0.006
the order-2 (C3D10) run drew a DIFFERENT sphere population from the order-1
runs -- porosity 0.0178 against 0.0119, flipping which arm got which -- so
C3D10 read 4.7e9 against 6.1e9 and the "convergence" was two different cells.
The slabs and bridges ARE deterministic, so `converge_u6.sh` drops the sphere
population entirely. **Check `porosity`/`phi_soft_total` agree across arms
before reading any modulus.**

## What is still thin

* **Cross-code breadth.** Only `L_mesh` 0.0240/0.0330 has verified mesh
  equivalence, so that is where the Abaqus comparison lives. 42 bracket stems
  with und/drn pairs and 3 seeds each sit there unused -- that is the
  validation, not more refinement of one cell.
* **Which patch coverage.** `u6all` wins at all four refinement meshes; on the
  bracket cells the two are within 0.5 % of each other and both well inside the
  seed floor. Not decided.
* **`m0p0080` and finer with U6 is not reachable in 32 GB.** The patch couples
  every node of a 1-ring to every other, so nz/eq goes 21.5 -> 49.2 on the same
  mesh; `m0p0120` peaked at 13.5 GB and `m0p0080` projects to ~68 GB.

# THE SPURIOUS MODE: mechanism found, fix partial. Do not run the cluster yet.

BRKB (5 layered und/drn conditions, Abaqus R from 2.9 to 19.4, all at L_mesh
0.0240 where the drained mesh-equivalence gate PASSES) put U5+U6 inside the
Abaqus seed floor at b020/b040/b080/b150 where plain C3D4 was outside it. At
**b280 it failed: R 8.97 % below Abaqus, four times the 2.07 % floor.**
Reproduced on a second packing (-12.35 % vs -12.57 %), so it is systematic.

## What it is

The U5 tet is deviatoric-only, and for brine G = K/5000, so ALL volumetric
restraint comes from the nodal patch. A brine node ON THE ICE INTERFACE is also
held by the neighbouring C3D4, which keeps its own volumetric term. A brine
node in the INTERIOR has no such anchor, and the patch operator b^a has a null
space: a displacement pattern with theta_a = 0 at every node of an interior
cluster costs only G-level energy and grows by ~K/G.

Measured on b280 (`tmp/leak/`, from the solved displacement field):

* volumetric energy deficit `1 - E_nodal/E_elem` = **97.9 %** (39.6 % with the
  spheres deleted); 90 % of it in 5.7 % of the elements
* **max |u| = 9.4e-02 against an applied 5.0e-03** -- 19x. The C3D4 control on
  the SAME mesh has max |u| = 4.8e-03 and NOT ONE node over the applied value;
  U6 has 2467
* the mode's support is **100 % brine-interior nodes, 0 % interface, 0 % ice,
  0 % triple points**
* those 2467 nodes form 22 clusters, and every large one is a SHEET (singular
  values 0.045/0.032/0.006, bbox 0.022 x 0.248 x 0.200) -- the mid-surfaces of
  the brine slabs, thin in x
* a hot node's own patch is 68 % hot (12 % baseline); 849 patches are >80 % hot

Sheets normal to x explains the direction: E_x fell 12.6 %, **E_z only 0.2 %**.
And equilibrium_gap stayed at 3.9e-07 throughout, because the mode is very
nearly free -- **the gap cannot see this failure.**

## Five hypotheses that were wrong, and why that matters

Every geometric metric was FLAT across a series whose error grew 8x, so none of
them is the discriminator. Recorded so they are not re-tried:

| hypothesis | measured | verdict |
|---|---|---|
| patch ring too large | 11.3 -> 11.4 | flat |
| interface / one-sided patches | 79.2 % -> 71.9 % | wrong direction |
| brine cut into small pockets | -17 % | too small |
| slivers / tet quality | no slivers, b280 the best | no |
| local thickness (BFS, `softthickness.py`) | best cell is the THINNEST | anti-correlates |

The last one is the instructive failure: thin brine is FINE, because a thin band
has every node on an interface and therefore anchored. Thick brine is the
hazard. `--min-body` in `nodalbbar.py` was written on the "small inclusions get
polluted patches" story and moved the answer by **0.01 %** -- the story was
clean and wrong. Deleting the spheres from the cell fixed it (-0.95 %), but
deleting only the voids (-1.34 %) or only the inclusions (-2.11 %) did not: the
populations are superadditive because together they fatten the brine enough for
interior sheets to exist.

## Refinement makes it WORSE

| L_mesh | elements across the brine slab | dE_und | R vs Abaqus |
|---|---|---|---|
| 0.024 | 0.66 | -12.57 % | -8.97 % |
| 0.016 | 1.00 | -10.13 % | -4.24 % |
| 0.012 | 1.33 | **-21.60 %** | **-15.79 %** |

So there is no resolution criterion to hide behind: refining adds interior
nodes and gives the mode more room. This also explains the LMESH 0.0120 result
(R 15 % low) that had been put down to mesh inequivalence -- same effect.

## The fix, and where it stands

`CCX_U5_STAB` returns a fraction of K to the ELEMENT (`e_c3d_u5.f`) and takes
the same fraction out of the patch (`u6patch.f`), so the two halves still sum
to K, the mean-dilatation identity is untouched and the patch test stays exact
for any value. The element charges `stab*K*div(u)^2`, which is exactly what the
patch average cannot see, so it penalises the null modes directly.
`resultsmech_u5.f` carries the same term -- adding it to the stiffness alone
took equilibrium_gap from 2e-07 to 1.6e-01.

GLOBAL stabilisation trades one case against another -- the CAPPED failure
again:

| case | floor | C3D4 | stab 0 | stab 0.05 | stab 0.10 |
|---|---|---|---|---|---|
| b020 | 7.99 % | +0.72 | -0.75 | -0.18 | +0.01 |
| b150 | 1.54 % | +3.62 | **-0.43** | +1.67 LOCK | +2.17 LOCK |
| b280 | 2.07 % | +4.11 | **-8.97** | **-0.88** | +0.71 |

LOCALISED stabilisation (element type **U7** = interior tets, stabilised; U5 =
anchored, not) removes the trade at b150 but does not fully close b280:

| U7 rule | tets marked | b150 | b280 |
|---|---|---|---|
| `all` (4/4 nodes interior) | 21 % | inside at every stab | only -2.29 % at stab 1.0 |
| `any` (>=1 node interior) | 85-98 % | re-locks +1.64 % | -0.88 % |

`all` is too narrow -- a tet with three interior nodes carries the mode and is
left unstabilised. `any` is so broad it is the global case again. **The rule
must be graded: stabilise each tet in proportion to how many of its nodes are
interior.** That is the open work. Note also that equilibrium_gap grows with
stab on MIXED patches (1e-03 at stab 1.0), so the U5/U7 bookkeeping inside a
shared patch needs checking as part of it.

Plumbing for U7 is done: `e_c3d_u.f` and `resultsmech_u.f` dispatch '7',
`e_c3d_u5.f`/`resultsmech_u5.f` stabilise only '7', `u6patch.f` scales only U7
contributions and includes U7 in its map, `mastruct.c`/`mastructcs.c` accept
'7' in their user-element whitelists (4 sites -- omitting them makes ccx treat
the element as a substructure and abort in `mastructread`), and `nodalbbar.py`
emits U7 with `SPAX_U7_RULE=all|any`. `CCX_U5_STAB=0` reproduces the
unstabilised element bit-for-bit.

## GRADED stabilisation: this is the fix

`SPAX_U7_RULE=graded` (now the default) emits **U7A..U7D** for tets with 1..4
nodes interior to the soft phase; ccx scales `CCX_U5_STAB` by n/4. Tets with no
interior node stay plain U5 and unstabilised. The grade rides in `lakon(3:3)`
as a LETTER, because ccx's built-in dispatch keys on digits in the label.

At `CCX_U5_STAB=0.07`, the whole BRKB family is inside the Abaqus seed floor:

| case | R_abq | floor | C3D4 | U6 raw | **U6 graded** |
|---|---|---|---|---|---|
| b020 | 19.3558 | 7.99 % | +0.72 % | -0.75 % | **-0.28 %** |
| b040 | 12.9543 | 2.95 % | +2.67 % | +0.94 % | **+1.59 %** |
| b080 | 7.5979 | 2.32 % | +1.81 % | +0.08 % | **+0.68 %** |
| b150 | 6.2903 | 1.54 % | +3.62 % | -0.43 % | **+1.43 %** |
| b280 | 2.9150 | 2.07 % | +4.11 % | **-8.97 %** | **-1.41 %** |

C3D4 outside on 2 of 5, raw U6 outside on 1 of 5, graded U6 outside on none.

And it survives refinement, which is the test raw U6 failed outright:

| L_mesh | el/slab | R raw | vs abq | R graded | vs abq |
|---|---|---|---|---|---|
| 0.024 | 0.66 | 2.6534 | -8.97 % | 2.8738 | -1.41 % |
| 0.016 | 1.00 | 2.7913 | -4.24 % | 2.9850 | +2.40 % |
| 0.012 | 1.33 | 2.4546 | **-15.79 %** | 2.9172 | **+0.08 %** |

### What is still not settled

* **0.07 is a tuned constant, and b150 has 0.11 points of margin** (+1.43 % on a
  1.54 % floor). Unlike CAPPED this is aimed at a measured mechanism and the
  response is monotone rather than a crossing -- but it is still a constant
  fitted to five cells. It needs the other tight-floor families (nbridges,
  density) before it can be called general.
* **equilibrium_gap rises from 2e-07 to ~1e-04** once stab is on. It falls with
  refinement (2e-4 / 1e-4 / 5e-5) and does not move E_x at the percent level,
  but it means the U5/U7 bookkeeping inside a SHARED patch is not exact: a
  patch spanning tets of different grades gives up a blended fraction of K
  while each element charges its own. That should be made exact.
* The 0.016 point is +2.40 %, just outside; the sequence -1.41 / +2.40 / +0.08
  is not monotone. Each refinement level re-meshes with a different sphere
  packing, so a few percent of that is packing noise, but it has not been
  separated.
* The acceptance ladder has NOT been re-run against the graded element: patch
  test, zero-energy modes, and the C3D10 convergence sweep were all verified on
  the unstabilised element. `CCX_U5_STAB=0` reproduces it bit-for-bit, so the
  old results still stand for stab=0, but they say nothing about stab>0.

## 2026-08-24 — the linear tet is exhausted: MINI falsified, and why no scalar can work

Three measurements taken together close off every conforming-linear-tet route
that has been tried, including the one this file has been recommending.

### 1. `CCX_U5_STAB` is the wrong operator, not a mis-scaled one

Local-projection / Dohrmann–Bochev theory does not leave the stabilisation
coefficient free: the pressure-deficit space should carry stiffness of order
`2G`, i.e. `CCX_U5_STAB ~ 2G/K`. For brine that is `4e-4`, three orders of
magnitude below the fitted `0.07`. A new rule `SPAX_U7_RULE=uniform` (every
retyped tet graded D, no topological targeting) was added so the coefficient
could be probed without the targeting rule confounding it.

On BRKB_b280 undrained, applied displacement 5.0e-03:

| `CCX_U5_STAB` | max abs u | vs applied | nodes over |
|---|---|---|---|
| 0 | 9.44e-02 | 18.9x | 2467 |
| **4e-4** (theory) | 9.09e-02 | **18.2x** | 2420 |
| 4e-3 | 7.28e-02 | 14.6x | 2142 |
| 4e-2 | 2.97e-02 | 5.9x | 1235 |
| 0.07 | 2.06e-02 | 4.1x | 930 |
| 1.0 | 5.24e-03 | 1.0x | 2 |

The theory-scaled value dents a 19x mode by 4%. The mode only dies at
`stab = 1.0`, which *is* C3D4 (elements carry all of K, patches carry none).
So the trade between locking and the spurious mode is monotone in the scalar
and there is no interior point that removes both. That kills the calibrated,
graded, and self-determining variants at once -- they differ only in *where*
the scalar is applied, never in the fact that it is one.

### 2. MINI (`U4`) is catastrophic, not marginal

`U4` has been in the tree since the beginning and was set aside on a layered
cell. Re-measured on the bracket cells, undrained arm only, denominator plain
C3D4 as always:

| | R | vs Abaqus C3D4H |
|---|---|---|
| **BRKB_b150** (R_abq 6.2903, floor 1.54%) | | |
| C3D4 | 6.5181 | +3.62% |
| U5+U6 raw | 6.2631 | -0.43% |
| **U4 MINI** | **1.8906** | **-69.9%** |
| **BRKB_b280** (R_abq 2.9150, floor 2.07%) | | |
| C3D4 | 3.0348 | +4.11% |
| U5+U6 raw | 2.6534 | -8.97% |
| **U4 MINI** | **1.2743** | **-56.3%** |

`R -> 1` means the undrained cell reads the same as the drained one: the
brine's 2.2 GPa bulk modulus is gone entirely. The cause is the one `u4mat.f`
already documents -- the P1 pressure cannot see the bubble's dilatation, so the
bubble opens a compliance channel worth K/G = 5000.

### 3. `STIFFB` shows the bubble cannot be both stabiliser and honest

`CCX_U4_STAB=STIFFB` gives the bubble its own volumetric energy, closing that
channel. It also drops the condensed stabilisation `S = B_b A_bb^-1 B_b^T` from
order `h^2/G` to order `V/K` -- the same scaling as the physical compressibility
term it sits next to, which is no stabilisation at all. Measured on b280:

| U4 variant | max abs u | vs applied | nodes over |
|---|---|---|---|
| MINI | 2.05e-02 | 4.10x | 2160 |
| CAPPED | 5.73e-02 | 11.46x | 4870 |
| STIFFB | 5.57e-02 | 11.14x | 4651 |

STIFFB lands back with the unstabilised element. The bubble is simultaneously
the stabiliser and the leak, exactly as the U5 header says killed U4 the first
time; giving it volumetric energy trades one for the other with no gain.

### What this means

Every formulation tried is P1 displacement with a pressure the mesh cannot
support: element-wise P0 (C3D4, locking), nodal P1 (U5+U6, spurious mode), or
P1 with a stabilisation (`CCX_U5_STAB`, MINI) that interpolates between them.
The inf-sup requirement is that the stabilisation scale as `h^2/G`; the physics
requires the constraint compliance to be `V/K`. At `K/G = 5000` on a mesh with
0.66-1.33 elements per brine slab those two demands are a factor of thousands
apart, and nothing chosen inside the P1 family can satisfy both.

**The fix must leave P1.** Quadratic displacement is inf-sup stable with *no*
stabilisation and therefore no parameter and no leak -- either Taylor-Hood
P2/P1 or P2 with reduced volumetric integration (P2/P0). Order 2 costs about
7.4x the DOFs on the BRKB cells (211k -> 1.68M with the periodic equations).
On the frozen-packing convergence cell stock C3D10 is already mesh-converged to
0.33% over a 2.5x refinement while C3D4 sits 1.5-3% high, which is the first
evidence that the stock quadratic tet needs no custom element at all.

## 2026-08-24, later — that conclusion was wrong: the element works at nu = 0.499

The section above ends "the fix must leave P1".  **That is superseded.**  It was
written before the element had ever been put through a locking test or an
eigenvalue census, and both of those changed the answer.  What was actually
wrong was not the element but the *stabilisation* and the *operating point*.

### The stabilisation was the locking, not the cure

`elements_ccx/tests/locking_sweep_u6.sh` — cantilever, shear modulus held
fixed, tip deflection over Euler-Bernoulli.  Flat in nu = no locking.

| nu | C3D4 | U5+U6 s=0 | U5+U6 s=0.07 | U5+U6 s=1.0 |
|---|---|---|---|---|
| 0.30 | 0.6754 | 0.7419 | 0.7366 | 0.6754 |
| 0.45 | 0.5261 | 0.7238 | 0.7033 | 0.5261 |
| 0.49 | 0.2630 | 0.7000 | 0.6098 | 0.2630 |
| 0.499 | 0.0560 | **0.6799** | 0.3002 | 0.0560 |
| 0.4999 | 0.0205 | **0.6580** | **0.0685** | 0.0205 |
| 0.49999 | 0.0166 | **0.5947** | 0.0221 | 0.0166 |

Unstabilised U5+U6 is genuinely locking-free.  `CCX_U5_STAB=0.07` -- carried for
weeks as "the fix" because it was 5/5 inside the Abaqus floor on BRKB -- reads
0.0685 at the brine's own Poisson ratio, against C3D4's 0.0205.  It re-creates
almost exactly the pathology the element exists to remove.  Obvious in
hindsight: 0.07*K at K/G = 5000 is 350*G, not a small perturbation.  It had
simply never been measured, because the only tests ever run on the graded
element were RVE ratios, and R cannot separate locking from over-softening.

### The eigenvalue census says the same thing

`elements_ccx/tests/stability_modes.py`, two-phase clamped cube (a homogeneous
one cannot show the instability -- every node patch reaches the clamped
boundary, so none is unanchored).  lambda_1/G is the softest mode over the
shear scale; it grows with K under locking and stays flat when there is none.

| K/G | C3D4 | U5+U6 s=0 | U5+U6 s=0.07 |
|---|---|---|---|
| 50 | 0.943 | 0.296 | 0.380 |
| 500 | 6.10 | 0.377 | 1.027 |
| 5000 | **47.78** | **0.545** | 6.489 |

C3D4's softest mode costs 47.8x the shear scale at the brine ratio -- locking,
quantified.  U5+U6 s=0 is flat over a hundredfold change in K/G.

### The operating point: K/G = 500, nu = 0.499

The element is locking-free but NOT inf-sup stable, and its spurious mode grows
as K/G because the soft phase escapes a K-level restraint at a G-level cost.
Prototype sphere cell, max fluctuation over the applied displacement:
5.69x at K/G = 5000, **1.05x at 500**, against C3D4's 0.20x.  At K/G = 500 the
census finds NO mode with volumetric content above 0.3.

That ratio is a modelling choice, not data -- brine is a liquid, G = 0, and
4.4e5 is already a regularisation.  Cost of nu 0.4999 -> 0.499, measured with
C3D4 on both (same element, so the shift is pure physics): **+0.17%** on b150,
**+0.21%** on b280, against seed floors of 1.5-2.1%.  nu = 0.49 costs +1.73% /
+1.90%, about a full floor, which is why the envelope stops at 0.499.

`SPAX_BRINE_KG=500` (SpaX_Standalone.py) raises G at fixed K; `nodalbbar.py`
refuses a deck outside the envelope unless `SPAX_BBAR_KGMAX` overrides it.

### Acceptance suite

`elements_ccx/tests/acceptance_u6.sh`, at nu = 0.499:

| test | U5+U6 s=0 | C3D4 |
|---|---|---|
| T1 patch, 6 components, 1296 ip, distorted | 0.0e+00 | 0.0e+00 |
| T2 confined, exact C1111 = K + 4G/3 | 0.0e+00 | 0.0e+00 |
| T3 zero-energy modes | 6 | 6 |
| T4 cantilever tip/EB | 0.6799 | 0.0560 |
| T5 census lambda_1/G, spurious count | 0.377, 0 | 6.10, 0 |

**C3D4 is the harness's control, and it earned its place.**  The suite's first
two runs reported four failures including C3D4 failing its own patch test --
both times the harness, not the element.  First, `mesh_box` jitters boundary
nodes (it must, to keep periodic image pairs matching), so the `x <= tol`
boundary detection found nothing and the block was unrestrained.  Second, a
free 0.3/h jitter of the Freudenthal split makes slivers -- at n=6, a tet of
volume 5.0e-05 against a mean of 7.7e-04 -- and since B ~ 1/h, at K/G = 500 the
solver's round-off is amplified into a stress error of order the stress itself
(C3D4 read -2.17e+05 where the answer is +2.64e+05).  Both would have been
written up as element defects without the control.  `make_block.py` now jitters
interior nodes only and backs the amplitude off until the worst tet is at least
`SPAX_BLOCK_QMIN` (0.15) of the mean volume.
