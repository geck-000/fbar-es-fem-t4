# F-barES-FEM-T4: formulation, small-strain reduction, and what it commits us to

Source: Y. Onishi, R. Iida, K. Amaya, *F-barES-FEM-T4 for viscoelastic problems*,
Int. J. Comput. Methods **15**(7) 1845003 (2018).  Equation numbers below are the
paper's.  `~/Downloads/onishi2019.pdf`; text dump kept alongside the run logs.

The method it builds on is Onishi, Iida & Amaya (2017), which is where the
recommended `c` per Poisson ratio actually lives.  We do not have that paper.
What the 2018 paper states about the envelope is quoted verbatim in section 5.

---

## 1. Notation

| symbol | meaning |
|---|---|
| `ᵉ(·)` | element quantity, constant over a T4 |
| `ⁿ(·)` | node quantity |
| `ᵸ(·)` | edge quantity, constant over an edge smoothing domain |
| `ᵉV`   | element volume (reference) |
| `ᵸV`   | edge smoothing volume, `Σ_{e∈ᵸE} ᵉV/6` |
| `ⁿV`   | node smoothing volume, `Σ_{e∈ⁿE} ᵉV/4` |
| `ᵸE`   | elements attached to edge `h`; `ⁿE` elements attached to node `n`; `ᵉN` the 4 nodes of `e` |

The `/6` and `/4` are the edge and node counts of a T4, so each element
distributes its whole volume across its 6 edges and its 4 nodes.  Every
smoothing operator below therefore has **row sums equal to 1**, which is what
makes the small-strain reduction in section 3 exact rather than approximate.

## 2. The finite-strain chain

**(a) Edge-smoothed deformation gradient — ES-FEM, applied once.** Eq. (1):

```
ᵸF̃ = Σ_{e∈ᵸE} ᵸw_e ᵉF ,        ᵸw_e = (ᵉV/6) / ᵸV
```

equivalently through smoothed shape-function gradients, eq. (3):
`ᵸF̃_ij = ᵸÑ_{P,j} x_{P:i}`.

**(b) Isovolumetric split.** Eqs. (4)-(5):

```
ᵸJ̃ = det(ᵸF̃) ,      ᵸF̃^iso = (ᵸJ̃)^(-1/3) ᵸF̃
```

**(c) c-time cyclic smoothing of J — NS-FEM, applied c times.** Eqs. (6)-(7):

```
ⁿJ̄ = (1/ⁿV) Σ_{e∈ⁿE} ᵉJ̃ (ᵉV/4)          (6)   element -> node
ᵉJ̄ = (1/4) Σ_{n∈ᵉN} ⁿJ̄                   (7)   node -> element
```

repeated `c` times, with `ᵉJ̄` fed back in as `ᵉJ̃` on the second and later
passes.  Then once, to the edges, eq. (8):

```
ᵸJ̄ = (1/ᵸV) Σ_{e∈ᵸE} ᵉJ̄ (ᵉV/6)           (8)   element -> edge
```

**(d) Volumetric gradient and recombination.** Eqs. (9), (11):

```
ᵸF^vol = (ᵸJ̄)^(1/3) I
ᵸF̄     = ᵸF^vol · ᵸF̃^iso = (ᵸJ̄ / ᵸJ̃)^(1/3) ᵸF̃
```

`det(ᵸF̄) = ᵸJ̄` by construction: the edge keeps ES-FEM's shape change and takes
its volume change from the cyclically smoothed field.

**(e) Stress.** Eqs. (12)-(13), Hencky strain `H̄` of `ᵸF̄`:

```
T^hyd = K tr(H̄) I ,     T^dev = 2G₀ (H̄^dev - Σᵢ gᵢ H^vᵢ)
```

For us the Prony terms are absent, so `T^dev = 2G H̄^dev`.  Two identities
matter: `tr(H̄) = ln det(ᵸF̄) = ln ᵸJ̄`, and `H̄^dev = H̃^dev` because `F̄` and `F̃`
differ only by a spherical factor.  **So the deviatoric stress never sees the
cyclic smoothing, and the hydrostatic stress never sees anything else.**

**(f) Internal force — and the one line that defines F-bar.** Eq. (17), with
the paper's own note:

> "Note that the stretching tensor in this equation, `D̃`, is **not** the
> deformation rate of `ᵸF̄` in Eq. (11) but that of `ᵸF̃` in Eq. (1) due to the
> adoption of the F-bar method."

Stress from the **modified** gradient, virtual work paired with the
**unmodified** one.  This is Petrov-Galerkin and the tangent, eq. (19), is
**non-symmetric**.  That is not an incidental detail — see section 4.

## 3. Small-strain reduction (what we actually implement)

Let `ε_e = sym ∇u|_e` and `θ_e = tr ε_e = div u|_e`, both constant per T4.
Collect the three smoothing operators as sparse matrices:

```
Q : elem -> node ,  Q[n,e] = (ᵉV/4)/ⁿV
P : node -> elem ,  P[e,n] = 1/4
E : elem -> edge ,  E[h,e] = (ᵉV/6)/ᵸV
A = P Q             one cycle, elem -> elem
```

All three have unit row sums, so on `ᵉJ ≈ 1 + θ_e` they act affinely and the
constant passes through untouched: the whole chain is **linear in θ**.

```
ᵸθ̃ = (E θ)_h                 edge strain trace        (ES-FEM, once)
ᵸθ̄ = (E A^c θ)_h             cyclically smoothed      (eqs. 6-8)
ᵸε̃ = (E ε)_h
ᵸε̄ = dev(ᵸε̃) + (1/3) ᵸθ̄ I    the F-bar strain, eq. (11)
```

Stress at the edge: `σ_h = 2G dev(ᵸε̃) + K ᵸθ̄ I`.  Virtual work paired with
`δᵸε̃` (section 2f).  With `B̃_h = Σ_e ᵸw_e B_e` the edge strain-displacement
matrix and `D_div` the element divergence operator:

```
K = Σ_h ᵸV B̃_hᵀ D_dev B̃_h            deviatoric: symmetric ES-FEM
  + K_bulk (E D_div)ᵀ W (E A^c D_div)   volumetric: NON-SYMMETRIC
```
with `W = diag(ᵸV)`.  Test space `E D_div`, trial space `E A^c D_div`.

At `c = 0` the two spaces coincide and the scheme collapses to plain selective
ES-FEM-T4, symmetric.  That is the only `c` for which symmetry is legitimate.

## 4. The Petrov-Galerkin form, and what it does not change

The obvious-looking Galerkin assembly `B̄ᵀ D B̄` is not what section 2f
specifies, and the difference is measurable.  What it is *not* is a change in
the constraint count: `K_vol = G_testᵀ W G_trial` has rank `min(rank G_test,
rank G_trial)` and the same right null space `ker(G_trial)` in either form, so
both admit exactly the same set of volumetric-free displacement modes.  The
counts, measured directly by `verify_fbar.py` V4 on a brine sphere with 1104
elements, 311 nodes and 1630 edges:

| c | rank `E` (test) | rank `E A^c` (trial) | constraints vs 3·n_node = 933 |
|---|---|---|---|
| 0 | 1011 | 1011 | ratio 0.92 — over-constrained, still locking territory |
| 1 | 1011 |  308 | ratio 3.03 — the Q1P0 / MINI count |
| 2 | 1011 |  308 | ratio 3.03 |

So the cyclic smoothing does exactly what it is designed to do: it relieves the
edge-volumetric constraint from over-constrained at `c = 0` to the textbook
count at `c ≥ 1`.  An ideal constraint *count* is not inf-sup stability,
though, which is the open question in section 6.

What the Petrov-Galerkin form changes is the distribution of internal force and
energy across those modes, not which modes exist.  It is required by eq. (17)
regardless of how large the difference turns out to be.

**Consequence for a CalculiX implementation.** ccx calls PARDISO with
`mtype = -2`, symmetric indefinite, and stores one triangle.  A correct
F-barES-FEM-T4 cannot be assembled into that path.  It needs `mtype = 11` and
a full-matrix storage scheme — a solver-level change, not a `*USER ELEMENT`.
This is a second, independent blocker on top of the `lakon(8:8)` 255-node
connectivity cap, which the `E A^c` stencil (~2c+1 element rings) breaches at
`c ≥ 2` and is already tight at `c = 1`.

## 5. The envelope the paper claims for itself

Verbatim, section 3.1:

> "According to our previous research for hyperelastic materials [Onishi et al.
> (2017)], the recommended c of F-barES-FEM-T4 for this material (**Poisson's
> ratio ν = 0.49 at most**) is 1 or 2."

and the validation material is `E = 1 MPa, ν = 0.3` instantaneous relaxing to
`ν ≈ 0.49` long-term.  **ν = 0.49 is K/G ≈ 50.**  Our operating point is
K/G = 500 (ν = 0.499); brine as specified is K/G = 5000 (ν = 0.4999).  So we
are asking the method for one to two orders of magnitude beyond anything the
authors published, with `c` stated to be ν-dependent and no formula given for
it.  That is a reason to measure our own operating point rather than to assume
either outcome; section 7 does.  It is not, on the evidence there, a wall.

## 6. Verification

`elements_ccx/tests/verify_fbar.py`, all passing:

* **V1** unit row sums on `Q`, `P`, `E`, `R` and on the composed chain `S`, to
  1e-15, for `c = 0..3` and both readings of eq. (6).  A constant `J` field
  survives the chain untouched, which is what makes section 3 exact.
* **V2** patch test: homogeneous block, uniform strain, `C1111 = K + 4G/3`
  recovered to 4e-16 relative for every `c` and both readings.
* **V3** `S = E` identically at `c = 0`, so the scheme collapses to selective
  ES-FEM-T4 exactly, as it must.
* **V4** the constraint-count table in section 4.

`elements_ccx/tests/verify_fbar_nl.py` checks the finite-strain element of
section 2 against closed-form answers, all passing:

* **N1** `f(0) = 0`; **N2** `f(rigid translation) = 0`.
* **N3a** unjittered box, free surfaces, uniform stretch: interior nodal
  residual vanishes to 1e-13 relative, and the resultant on the `x = 1` face
  matches the analytic Hencky Cauchy stress `(K + 4G/3) ln λ` to 1e-13, for
  `c = 0..3` at both `λ = 1.001` and `λ = 1.2`.  A 20% stretch, so this
  exercises eqs. (1)-(18) well away from the linear regime.
* **N3b** jittered periodic cell, homogeneous: the affine field is recovered
  with fluctuation 2e-12.
* **N4** the finite-strain element converges to its own small-strain reduction
  at O(eps): 6e-6 relative in C1111 at `eps = 1e-5`, for `c = 0, 1, 2`.

The reading ambiguity in eq. (6) — the input to the first cycle is written
`ᵉJ̃`, with the tilde the paper otherwise reserves for smoothed quantities, so
it is either the raw element Jacobian `det(ᵉF)` or an element restriction of
the already edge-smoothed `ᵸJ̃` — is implemented both ways and switchable with
`SPAX_FBAR_JIN=elem|edge`.  It changes the numbers by a few per cent and
changes no conclusion below.  `elem` is the default.

## 7. What it does, measured

**A retraction first.**  An earlier version of this section reported that the
method comes apart above K/G = 50 — fluctuation growing monotonically with `c`
and 29% of C1111 lost by `c = 3` — and concluded that it could not reach our
operating point.  That was wrong, and the cause was not the element.  At
`jitter = 0.3` the prototype's structured mesh **tangles**: one inverted tet at
`n = 6`, fourteen at `n = 16`.  `grads()` silently flipped their node order,
which restored positive volumes and hid it, but the tets still overlapped their
neighbours, so the mesh was not a partition of the cell.  On that mesh nothing
passes its own patch test — C3D4 included, at `fluc` 1.4e-3 with a residual of
0.21 against the exact affine field, at interior nodes.  `mesh_box` now shrinks
the jitter amplitude until no tet is inverted and the worst is at least
`SPAX_MESH_QMIN` (0.15) of the mean, the same guard `make_block.py` has carried
since the acceptance suite hit this.  Every number below is remeasured.

Brine sphere `r = 0.30` in ice, periodic cell, requested jitter 0.3 (shrunk to
0.64 of that at `n = 8` by the quality guard), `n = 8`.  `fluc` is max
displacement fluctuation per unit applied strain; `p-sch` is the scheme's own
pressure field, jump across brine-brine faces over mean `|p|`.

| K/G | scheme | C1111 | fluc | p-sch |
|---|---|---|---|---|
| 500  | c3d4   | 2.6926e+08 | 0.44 | 0.193 |
| 500  | ns_vol | 2.5067e+08 | 0.45 | 0.009 |
| 500  | fbar_0 | 2.5788e+08 | 0.54 | 0.049 |
| 500  | fbar_1 | 2.4946e+08 | 0.51 | 0.005 |
| 500  | fbar_2 | 2.4841e+08 | 0.44 | 0.002 |
| 500  | fbar_3 | 2.4814e+08 | 0.41 | 0.001 |
| 5000 | c3d4   | 2.5110e+09 | 0.38 | 0.139 |
| 5000 | ns_vol | 2.4141e+09 | 0.92 | 0.004 |
| 5000 | fbar_0 | 2.4527e+09 | 0.67 | 0.035 |
| 5000 | fbar_1 | 2.4102e+09 | 1.07 | 0.002 |
| 5000 | fbar_2 | 2.4046e+09 | 0.61 | 0.001 |
| 5000 | fbar_3 | 2.4032e+09 | 0.47 | 0.000 |

Both of the paper's claims now hold at our operating point and beyond it:

1. **Pressure.**  `p-sch` falls monotonically with `c` at every K/G, by two
   orders of magnitude from C3D4 by `c = 3`.  This is Fig. 6 of the paper.
2. **Displacement.**  `fluc` is not degraded by `c`; at K/G = 5000 it is
   non-monotone and `c = 3` (0.47) is *better* than `c = 0` (0.67) and better
   than NS-FEM's 0.92.  C1111 moves only −2.0% from `c = 0` to `c = 3` at
   K/G = 5000, against −29% on the tangled mesh.  This is the paper's "regardless
   of the number of cyclic smoothings".

At K/G = 5000, `fbar_3` is the best arm in the prototype: it is 4.3% softer
than C3D4 (which locks, so softer is the right direction), it carries the
smallest fluctuation of any smoothed scheme, and its pressure field is
effectively oscillation-free.  It beats the incumbent `ns_vol` on all three.

The paper's stated envelope (section 5, `ν = 0.49 at most`) is a statement
about where the authors validated `c`, not a wall we have found.

## 8. The CalculiX implementation

Delivered. Patches `0009` (elements, dispatchers, build) and `0010` (the
asymmetric path) are in `../../patches_ccx/`; the binary is `ccx_fbar`.

**Correcting section 4.** It said a non-symmetric tangent needs a solver-path
change ccx cannot make. That is wrong: ccx already carries a full asymmetric
path — element routines raise `nasym`, `mafillsmasmain.c`/`mafillsmas.f`
assemble it, and `pardiso.c` selects `mtype = 11`. The Petrov-Galerkin form of
eq. (17) is delivered as specified.

**The names are U2 and U3, not U7 and U8.** `nodalbbar.py` already emits
`U7<letter>` for its graded U5 variants and `U8` for its 5-node theta carrier,
and `e_c3d_u8.f` exists in the ccx tree as the face-jump stabiliser. `U2` and
`U3` were the free family digits.

### What decides the architecture: stencil width, measured

On the campaign's `LMESH_m0p0240` soft phase (117 437 tets, 36 323 nodes,
184 572 edges; element valence per node mean 12.9, max 74):

| operator | nodes/edge mean | p99 | max | DOF at max |
|---|---|---|---|---|
| U2 deviatoric, eq. (1) | 6.3 | 10 | 14 | 42 |
| U3 volumetric, `c = 1` | 33.7 | 78 | 173 | **519** |
| U3 volumetric, `c = 2` | 93.6 | 259 | 494 | **1482** |

Consequences, and they are hard limits rather than judgement calls:

* **`c = 1` is deliverable as an element.** 173 nodes is inside
  `*USER ELEMENT`'s 255-node limit, and patch `0008` had already widened the
  user-element matrix to 765 DOF, so nothing had to be widened for this.
* **`c = 2` cannot be an element on that cell.** `userelements.f:83` rejects
  `NODES > 255` and `mastruct.c` reads the node count from a single byte of
  the label. It would need a direct global assembly pass.
* **`*EQUATION` + `SPRING1` is not an option**, and this is settled rather
  than untried: the U6 header records that the campaign built exactly that,
  verified every piece in isolation, and still got a confined-compression
  reaction **1.55x** the closed-form answer once patches overlapped.

Since `c = 1` measures within noise of `c = 2` and `c = 3` on our cells
(section 7, and n = 16: `fluc` 0.34 vs 0.36 vs 0.34 at K/G 500), `c = 1` is the
production point and the element route is enough.

### Files

| file | role |
|---|---|
| `u2edge.f` | edge ring, `V_h` and the smoothed gradients of eq. (1) |
| `e_c3d_u2.f` | `U2`, deviatoric `V_h Bt^T D_dev Bt` |
| `u3vol.f` | the eq. (6)-(8) chain as a row walk; returns `sbar`, `tbar`, `K V_h` |
| `e_c3d_u3.f` | `U3`, volumetric `K V_h tbar^T sbar`, raises `nasym` |
| `resultsmech_u2.f` / `resultsmech_u3.f` | internal-force recovery, through the *same* `u2edge` / `u3vol` calls |
| `fbares.py` | generator: `U2` rings and `U3` stencils per material |
| `tests/verify_fbares_deck.py` | independent walk of the generated deck |

`verify_u8_chain.py` mirrors `u3vol.f`'s walk in Python and checks it against
the prototype operator `S = E A^c`: agreement to 1e-15 for `c = 1, 2`, on
jittered and unjittered meshes, with unit row sums preserved. (Its name still
carries the old U8 label.) The element-side algorithm was therefore verified
before it compiled, not after it produced a suspicious number.

Two details the generator must respect, both of which are silent failures
otherwise:

* `konl(1)` and `konl(2)` are the two edge nodes, in both `U2` and `U3`.
* The `U3` connectivity must be the *same* `E A^c` support the element walks.
  `u3vol` stops with a message naming the element if it meets a node the
  connectivity does not carry.

### The base tets are `U5Z`

The tets stay in the model because `u2edge` and `u3vol` read their geometry
through the `'U5'` node->element map, and they must contribute nothing because
`U2` and `U3` supply the whole stiffness. `fbares.py` retypes them to `U5Z`,
and `e_c3d_u5` / `resultsmech_u5` return on suffix `Z`. It is a deck-level
switch rather than an environment variable on purpose: an environment variable
would also null the tets of a nodal-B-bar (`U5` + `U6`) deck, where `U5`
carries the deviatoric itself.

### What an F-bar deck cannot report

No element carries stress: `U5Z` is null, and `U2`/`U3` are smoothing domains
with no shape function of their own. `fbares.py` therefore drops `*EL PRINT`
and `*EL FILE` — left in, they would print a column of exact zeros as if it
were the answer — and adds a `*NODE FILE` if the deck has none. Read the result
from displacements and reactions.

For the same reason there is no mass matrix, so ccx `*FREQUENCY` is unavailable,
exactly as for `U5` + `U6`. The zero-energy-mode census runs in
`smoothing_proto.py` instead, and for `c >= 1` it must use **singular values**
of `K`, not eigenvalues: the symmetric part of the Petrov-Galerkin operator is
indefinite, and an eigenvalue census silently reports the wrong thing.

### Two defects the wiring exposed

Both are in stock ccx and both were silent.

1. **The asymmetric mechanical assembly loop covered only contact elements.**
   `mafillsmasmain.c` partitions that range over `ne0+1..ne`, and with no
   contact it sets `nea = ne0+1 > neb = ne0` — empty. U3 lives among the
   ordinary elements, so the whole volumetric half of the scheme vanished from
   the matrix while `neq` and the sparsity stayed correct. The patch test came
   back with interior displacements **200x** the applied field.

2. **`nope = lakon[8*i+7]` read the node count as a signed char.** Over 127
   nodes it went negative and `mastruct` built no structure for that element;
   `add_sm_st_as` then stopped with *coefficient should be 0*. Patch `0008` had
   argued the ceiling was 255 from the Fortran side, which was right about the
   Fortran and wrong about the C. U6's widest measured ring is 127 exactly, so
   it never hit it; a `c = 2` stencil of 138 nodes did.

### Acceptance, real ccx (n = 6, jitter 0.3, nu = 0.499, K/G = 500)

C3D4 is the control. It is what caught the tangled prototype mesh (section 7)
and it is what says the harness is sound.

| test | C3D4 | c = 0 | c = 1 | c = 2 |
|---|---|---|---|---|
| patch, `max abs(u - eps.X) / max abs(u)` | 3.165e-07 | 3.165e-07 | 3.165e-07 | 3.165e-07 |
| confined `C1111` (exact 2.204400000e+08) | 2.204400e+08 | 2.204400e+08 | 2.204400e+08 | 2.204400e+08 |
| `abs(sum RF) / max abs(RF)` | 1.6e-15 | 7.1e-15 | 3.5e-15 | 5.8e-15 |

The three F-bar arms agree with the control to the last digit printed; 3.2e-07
is this deck's PARDISO round-off floor at K/G = 500, not the elements'. The
force balance is the cross-check that matters most: the reactions come from
`resultsmech_u2`/`u3` and the matrix from `mafillsmas`, through completely
separate code paths.

Zero-energy modes, `smoothing_proto.py`, clamped-free block, singular values:
null dimension 6 for `c = 0, 1, 2` at both K/G = 500 and 5000, with
`sigma_7/sigma_max` at the G scale — no spurious mechanism and no locking.
Asymmetry `max abs(K - K^T) / max abs(K)` is 3.3e-02 at `c = 1` and 2.4e-02 at
`c = 2` against 2.7e-16 at `c = 0`, so the Petrov-Galerkin structure of eq. (17)
is live in the assembled matrix and not just in the element.

### `V_n` is per material

Eq. (6)'s nodal volume must be summed over the smoothing material only. Caching
it once over all elements at a node and then skipping foreign elements in the
walk divides by an inflated `V_n`, breaks the unit row sum of `Q`, and stops
the chain preserving a constant `J` — a patch-test failure that appears only at
the brine/ice interface. `u3vol` recomputes it per material where it is used.

## 9. Against Abaqus C3D4H on the layered cells

The measurement the whole exercise is for. `calculix/layered_abaqus_ratio.sh`
compares

    R = E_x(undrained) / E_x(drained)

against the stored Abaqus campaign, where the undrained cell was meshed with
C3D4H and the drained one with plain C3D4. The denominator is always plain
C3D4 in both codes, so only the undrained element differs, and Abaqus's own
seed-to-seed spread in R sets the noise floor the difference has to beat.
Locking makes the undrained cell too stiff and so inflates R.

| cell | element | R_ccx | R_abq | floor | excess | verdict |
|---|---|---|---|---|---|---|
| LMESH_m0p0240 | C3D4 | 2.5074 | 2.4701 | 0.61% | **+1.51%** | locking |
| | U5+U6 brine | 2.4572 | | | −0.52% | inside |
| | F-barES `c=0` | 2.4893 | | | +0.78% | locking |
| | F-barES `c=1` | 2.4448 | | | −1.03% | below ref |
| LMESH_m0p0120 | C3D4 | 2.5736 | 2.3786 | 0.75% | **+8.20%** | locking |
| | U5+U6 brine | 2.1718 | | | −8.70% | below ref |
| | F-barES `c=0` | 2.4933 | | | +4.82% | locking |
| | F-barES `c=1` | 2.0551 | | | −13.60% | below ref |

`c=1` on the finer cell needs patch `0011` (ring x support) to build at all and
patch `0012` (out-of-core PARDISO) to factor: 600 662 equations, 68.6M
nonzeros, 2143 s.

### The same-mesh reference is not converged, and that changes the reading

The stored campaign carries two FINER Abaqus meshes, and R falls hard with
refinement:

| L_mesh | Abaqus C3D4H R, per seed | mean |
|---|---|---|
| 0.0240 | 2.4626, 2.4776 | 2.4701 |
| 0.0120 | 2.3875, 2.3697 | 2.3786 |
| 0.0080 | 2.1113, 1.8651 | 1.9882 |
| 0.0060 | 1.9650, 2.0145 | 1.9897 |

So **the C3D4H value each arm is being scored against is itself 24% (0.0240)
and 20% (0.0120) above the level Abaqus's own finest meshes reach.** The
table below scores against that level -- but see the two sections after it:
~1.989 is NOT a converged value, and nothing here should be quoted as one.

| L_mesh | arm | R | vs converged 1.9890 |
|---|---|---|---|
| 0.0240 | CalculiX C3D4 | 2.5074 | +26.1% |
| | F-barES `c=0` | 2.4893 | +25.2% |
| | Abaqus C3D4H | 2.4701 | +24.2% |
| | U5+U6 brine | 2.4572 | +23.5% |
| | F-barES `c=1` | 2.4448 | +22.9% |
| 0.0120 | CalculiX C3D4 | 2.5736 | +29.4% |
| | F-barES `c=0` | 2.4933 | +25.4% |
| | Abaqus C3D4H | 2.3786 | +19.6% |
| | U5+U6 brine | 2.1718 | +9.2% |
| | **F-barES `c=1`** | **2.0551** | **+3.3%, INSIDE the fine-mesh scatter** |

At `0.0240` mesh error swamps everything -- every arm sits within 23-26% and
the ordering means little. At `0.0120` they separate, and F-barES-FEM-T4 at
`c = 1` is the only arm that lands in the same region as Abaqus's two finest
meshes, on a mesh 1.5-2x coarser.

### 1.989 IS NOT A CONVERGED VALUE, and this table must not be read as if it were

An earlier version of this section called 1.989 "converged" and scored the
arms against it. That was wrong. **The stored Abaqus sequence is not a
convergence study and cannot be made into one.**

*Every one of its 16 rows is a separate packing.* Not just seed to seed: the
undrained and drained cells of the SAME seed at the SAME mesh have different
porosity (`m0p0060_und_s1` 0.011083 against `m0p0060_drn_s1` 0.009905). So
each R is a ratio of two different random cells, and refining the mesh
re-packs as well as re-meshes.

*Neither seed shows asymptotic behaviour.* Per seed, R over
`0.0240, 0.0120, 0.0080, 0.0060`:

| seed | R | increments | verdict |
|---|---|---|---|
| s1 | 2.4626, 2.3875, 2.1113, 1.9650 | −0.075, −0.276, −0.146 | monotone, but increments **grow then shrink** |
| s2 | 2.4776, 2.3697, 1.8651, 2.0145 | −0.108, −0.505, **+0.150** | **not monotone** |

A converging sequence has increments that shrink as h shrinks. s1's largest
change is in the MIDDLE of the refinement; s2 reverses direction. Tracing it
to `E_x`: s2's undrained cell reads 4.313e9 at `0.0080` against 5.586e9 at
`0.0120` and 4.620e9 at `0.0060` -- one point out of line with both
neighbours, in a quantity whose drained twin is stable to 2% across the whole
sweep.

So what 1.989 actually is: the mean of the two finest mesh levels, whose two
seeds disagree by 12.4% at `0.0080`. What the data supports is only that **R
falls by about 20% between `0.0240` and `0.0060` and the two finest levels
both land near 2.0 rather than continuing to fall.** That is suggestive of an
asymptote somewhere near 2; it does not establish one, and no arm should be
scored to a tenth of a percent against it.

### What CAN be checked, and is

Our own points do not have the packing problem. `params/rve_layermesh.csv`
rows `LMESH_m0p0240_und_s1` and `LMESH_m0p0120_und_s1` differ **only** in
`run_id` and `L_mesh`; with `SPAX_SEED` fixed the packing is identical and
only h changes, and the drained twin is the same mesh with one elastic card
rewritten. So a CalculiX mesh sweep carries no packing noise at all, and it
is the one thing here that can answer whether R varies smoothly.

Over the two points available, the arms already separate by DIRECTION:

| arm | R(0.0240) | R(0.0120) | direction |
|---|---|---|---|
| CalculiX C3D4 | 2.5074 | 2.5736 | **up** +2.6% |
| F-barES `c=0` | 2.4893 | 2.4933 | flat +0.2% |
| U5+U6 brine | 2.4572 | 2.1718 | down −11.6% |
| F-barES `c=1` | 2.4448 | 2.0551 | down −15.9% |

The two arms that unlock are the two that move toward the region Abaqus's fine
meshes occupy; C3D4 moves away from it and `c = 0` does not move at all. Two
points cannot show convergence either, which is why
`elements_ccx/tests/meshconv.sh` fills in intermediate `L_mesh` values on the
same geometry.

**The finer end is out of reach for `c = 1`.** `0.0080` is 3.4x the elements of
`0.0120`, which puts the mastruct insertion count back over 2^31 even with the
ring x support reduction of patch `0011`. So our element cannot be taken to the
mesh where Abaqus's sequence flattens; the sweep can only establish whether the
trend is smooth over the range that is reachable.

### It does NOT converge over the meshes we can reach, and here is why

**Retracting the paragraph this section used to contain.** With three points
(`0.0240, 0.0180, 0.0120`) F-barES `c = 1` looked monotone with shrinking
increments and fitted `R = R_inf + C h^p` at `p = 1.54`, `R_inf = 1.851`. That
fit was EXACT -- three points, three unknowns, residual 4e-16 -- so it had no
redundancy and could not test the power law it asserted. A fourth point at
`h = 0.0360` destroys it:

| h = L_mesh | F-barES `c=1` | CalculiX C3D4 |
|---|---|---|
| 0.0360 | 2.4100 | 2.4712 |
| 0.0240 | 2.4448 | 2.5074 |
| 0.0180 | 2.2322 | 2.5407 |
| 0.0120 | 2.0551 | 2.5736 |
| steps | **+0.035**, −0.213, −0.177 | +0.036, +0.033, +0.033 |
| fit over all four | p = 0.00, `R_inf` = −1234, residual 1.1e−01 | p = 0.12, `R_inf` = 3.32 |

F-barES `c = 1` **turns over** between `0.0360` and `0.0240`, and no power law
fits the four points at all. C3D4 stays monotone but its four-point fit
degenerates too (`p = 0.12`, no meaningful limit) -- it is essentially linear
in h, which is not asymptotic behaviour either.

**The reason is that the mesh never resolves the microstructure.** The cell is
`L = 0.5` with `n_slabs = 4` and `slab_vof = 0.1`, so a brine slab is about
`0.5 x 0.1 / 4 = 0.0125` thick. Elements through one slab:

| h | 0.0360 | 0.0240 | 0.0180 | 0.0120 | 0.0060 (Abaqus's finest) |
|---|---|---|---|---|---|
| elements per slab | 0.35 | 0.52 | 0.69 | **1.04** | 2.08 |

Every mesh in this study is at or below ONE element through the soft layer that
carries the entire undrained effect. You cannot be in an asymptotic regime
there, and that -- not the elements -- is what makes both our sequence and
Abaqus's look the way they do. A converged R would want maybe 4-8 elements
across the slab, `h ~ 0.002-0.003`, which is 64-216x the elements of `0.0120`:
out of reach for `c = 1` by a wide margin, and probably for anything.

**What survives, and it is not nothing.** Over the three finest meshes the
DIRECTIONS are unambiguous and opposite: F-barES `c = 1` falls (2.4448 ->
2.0551), C3D4 rises (2.5074 -> 2.5736), and Abaqus C3D4H falls too (2.4701 ->
2.3786 -> ~1.99). The F-bar element moves the way the hybrid element moves and
plain C3D4 moves away from both, getting worse with refinement -- which is what
volumetric locking looks like on this cell and is the thing the element is
there to remove. That is a statement about direction, on a controlled
comparison at fixed mesh. It is not a converged modulus and must not be quoted
as one.

**Which K/G these numbers are at, and it is not the target one.** The
campaign's undrained decks carry `*ELASTIC 1320000, 0.4999` for the brine, so
every row of the table above is at **K/G = 5000**, with the drained twin at
`nu = 0.406` (K/G = 5). That is the campaign's own definition of undrained --
brine's real bulk modulus, 2.2 GPa against G = 4.4e5 -- and it cannot be moved
without losing the comparison, because the stored Abaqus reference was run at
it. So this table sits one to two orders of magnitude above the element's
validated envelope (K/G = 500, `nu = 0.499`) and above the paper's own
(`nu = 0.49 at most`, section 5). The section 8 acceptance ladder IS at
K/G = 500 and is exact there. Read this table as a controlled
element-substitution comparison at an extreme point, not as evidence about
the operating point.

**What it says.** The cyclic smoothing does exactly what the paper says it
does, and on `m0p0240` it steps straight over the target: `c = 0` — plain
selective ES-FEM-T4 — still locks at +0.78%, one cycle moves R by −1.81
points to −1.03%, and the ±0.61% Abaqus floor falls in the gap between them.
`c` is an integer count of smoothings, so there is no setting in between. On
this cell U5+U6 lands inside the floor and F-barES-FEM-T4 does not.

On the finer `m0p0120` cell, where C3D4 locks by +8.20% and U5+U6 overshoots
by −8.70%, `c = 0` halves the locking to +4.82% and is the closest any
CalculiX arm has come — but it still locks, and `c = 1`, the arm that would
close it, cannot be built at all (below).

**`c = 1` at `m0p0120` is out of reach, and the reason is `mastruct`, not
memory and not the 255-node limit.** `insert.c` appends one entry per
off-diagonal upper-triangle `(dof, dof)` pair of every element to `mast1`
with no deduplication, at 8 bytes each (`mast1` and `next`), and grows the
list by 1.1x with a 32-bit ITG index. A U3 element is dense over its whole
stencil, so the count is `sum_h d_h(d_h+1)/2` with `d_h = 3 n_h` — quadratic
in the stencil:

| cell | `c` | insertions | structure memory | outcome |
|---|---|---|---|---|
| LMESH_m0p0240 | 0 | 7.23e7 | 0.5 GB | 6.76M nnz, 40 s |
| LMESH_m0p0240 | 1 | 1.12e9 | 8.3 GB | 30.05M nnz, 154 s |
| LMESH_m0p0120 | 0 | 2.74e8 | 2.0 GB | 24.94M nnz, 471 s |
| LMESH_m0p0120 | 1 | 5.16e9 | 38.4 GB | `*ERROR in u_realloc: size(bytes)=-8589934592` |

Past `2^31` the 1.1x growth overflows and ccx asks for a negative allocation.
`fbares.py` now computes the count and refuses the deck with that number in
the message, rather than letting a two-hour generation end in a `u_realloc`
crash. Raising it needs a `long long` ITG build *and* more RAM than this
machine has, or the direct global assembly route of section 8.

**How E_x is read.** An F-bar deck carries no element stress, so `sigma_bar`
comes from the reference-point reaction alone
(`SPAX_CCX_SIGMA_FROM_RF=1` in `SpaX_CalculiX.py`) — the other of the two
independent measurements the campaign normally cross-checks against each
other, and the one Abaqus's macroscopic modulus is defined by. The
`equilibrium_gap` column is therefore empty for these arms by construction.
What stands in for it: on `m0p0240` the transverse reference points came back
with reactions of 2.7e-07 against 1.5e+07 at the driven one — traction-free
to 1.7e-14 — and the fixed centre carried 1.7e-08, so the periodic constraint
set and the element are consistent to round-off.


## 10. Mesh convergence on a cell that actually resolves the slab

Section 9 could not answer whether R converges, because no mesh in that study
put more than one element through the brine layer. `tests/make_slabconv.py`
makes the slab thick instead of the mesh fine: a brine slab at
`x in [0.4, 0.6]` pierced by square ice bridges, so `0.2 n` elements span it,
and every phase boundary is a multiple of 0.1 so a mesh with `n % 10 == 0`
resolves the geometry EXACTLY at every n. The soft fraction is 16.8% at every
n; nothing but h changes between points.

Loading is uniaxial STRESS across the slab -- rollers on the low faces only.
That choice is not cosmetic and cost three attempts to find:

* confined compression ACROSS the slab puts the layer in uniaxial STRAIN,
  where `div u` is uniform inside it and incompressibility never binds. C3D4
  and F-barES agreed to 0.1-0.2% at every bridge pattern and every K/G, while
  R itself ranged over 1.17 to 7.16.
* loading ALONG the slab is ice-dominated: R only 1.03 to 1.26.
* uniaxial stress lets the cell contract laterally, the ice (`nu = 0.33`) and
  the brine (`nu -> 0.5`) disagree about by how much, and the brine is forced
  to deform at nearly constant volume. That is the constraint that locks.

### K/G = 500, the operating point

| n | el/slab | C1111 drn | und C3D4 | und F-bar `c=1` | R C3D4 | R F-bar | C3D4 excess |
|---|---|---|---|---|---|---|---|
| 10 | 2 | 1.1577e+09 | 1.6158e+09 | 1.2560e+09 | 1.3957 | 1.0850 | **+28.64%** |
| 20 | 4 | 9.6249e+08 | 1.2224e+09 | 1.0338e+09 | 1.2700 | 1.0741 | +18.24% |
| 30 | 6 | 9.2762e+08 | 1.0971e+09 | 9.9705e+08 | 1.1827 | 1.0749 | +10.03% |
| 40 | 8 | 9.1347e+08 | 1.0686e+09 | 9.8273e+08 | 1.1698 | 1.0758 | +8.74% |
| 50 | 10 | 8.9982e+08 | 1.0205e+09 | 9.6840e+08 | 1.1341 | 1.0762 | +5.38% |
| 60 | 12 | 8.9317e+08 | 1.0039e+09 | 9.6170e+08 | 1.1240 | 1.0767 | +4.40% |

**F-barES-FEM-T4 is converged from two elements through the slab.** R reads
1.0850, 1.0741, 1.0749, 1.0758, 1.0762, 1.0767 -- flat to within 1% across a
six-fold refinement, with a 0.2% wiggle at the coarse end and monotone from
`n = 20` on.

**Plain C3D4 is not converged at twelve.** R falls 1.3957 -> 1.1240 and is
still falling, heading down toward the F-bar value, as it must: both elements
are consistent and share a limit.

The last column is the locking, isolated -- how much stiffer C3D4 reads than
the unlocked element on the SAME mesh. It falls 28.6, 18.2, 10.0, 8.7, 5.4 per
cent: a discretisation error on its way out, which is exactly what volumetric
locking is. F-barES-FEM-T4 reaches at `n = 10` an answer C3D4 has not reached
at `n = 60`, a 216-fold difference in element count.

This is the believability evidence section 9 could not supply, and it is at
K/G = 500 rather than 5000.

### Abaqus C3D4H is a fourth arm on the SAME mesh

`make_slabconv.py` emits `m_abq.inp` (C3D4H in the inclusion) beside
`m_ccx.inp` (C3D4) from one geometry -- same nodes, same elements, same
boundary conditions, only the element keyword differs. So C3D4H can be run on
exactly these meshes rather than being a reference from another campaign with
another packing, which is the weakness of every comparison in section 9.
`hpc/submit_slabconv_abaqus.sh` runs the sweep on CSC Roihu (`abaqus/2026`),
and `hpc/postprocess_slabconv.sh` (chained after the solve array) reads the
reactions back out of the `.odb` with `elements_ccx/tests/slabconv_extract.py`
-- the ODB route the other campaigns use, not the `*NODE PRINT` block.

The drained twins agree to 0.003% at every `n` -- both codes, both elements --
because nothing locks at `nu = 0.41`.  The undrained arm is the discriminator,
and `c = 0` (plain selective ES-FEM-T4) is included because `c` is an integer
count of cycles, so the question "does `c = 1` overshoot?" has to be answered
with the neighbouring setting:

| n | el/slab | C1111 und C3D4H | und F-bar `c=1` | R C3D4H | R F-bar `c=1` | R F-bar `c=0` | R C3D4 |
|---|---|---|---|---|---|---|---|
| 10 | 2 | 1.4802e+09 | 1.2560e+09 | 1.2787 | 1.0850 | 1.2335 | 1.3957 |
| 20 | 4 | 1.0466e+09 | 1.0338e+09 | 1.0875 | 1.0741 | 1.1235 | 1.2700 |
| 30 | 6 | 1.0031e+09 | 9.9705e+08 | 1.0815 | 1.0749 | 1.0985 | 1.1827 |
| 40 | 8 | 9.8715e+08 | 9.8273e+08 | 1.0807 | 1.0758 | 1.0978 | 1.1698 |
| 50 | 10 | 9.7128e+08 | 9.6840e+08 | 1.0795 | 1.0762 | 1.0897 | 1.1341 |
| 60 | 12 | 9.6396e+08 | 9.6170e+08 | 1.0793 | 1.0767 | 1.0882 | 1.1240 |

**F-barES-FEM-T4 `c=1` and Abaqus C3D4H share a limit.**  C1111 und agrees to
0.23% at `n = 60`, and a Richardson extrapolation (`R = R_inf + C h^2`,
`h = 1/n`, least squares over `n = 40, 50, 60`) puts the two limits at
**1.0773 and 1.0780 -- 0.06% apart**.  C3D4 is still 4% stiff at `n = 60`, and
its fit residual (3.9e-3) is 30x the others' because locking is not a clean
`h^2` error -- it is not in the asymptotic regime any of these meshes reaches.

**One smoothing cycle removes the last locking, and no more.**  `c = 0` still
reads 0.83% stiff at `n = 60` (R 1.0882, limit 1.0795); one cycle takes it to
0.23% (R 1.0767, limit 1.0773), and the reference C3D4H sits at 1.0793 /
1.0780 between them.  So `c = 1` lands inside the reference's own convergence
scatter, `c = 0` is still short, and there is no half-step available.

`analysis/slabconv_report.py` reproduces this table, the Richardson fits, and
the `C1111`-vs-`h` figure (`out_slabconv/<case>/slabconv_convergence.png`).
