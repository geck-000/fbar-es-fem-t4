# The F-barES-FEM-T4 element (U2 + U3 + U4)

Fortran sources added to `ccx` to give it the locking-free 4-node tetrahedron
that Abaqus calls `C3D4H`, delivered as three `*USER ELEMENT` types over the
underlying T4 connectivity:

| Type | Nodes | Role |
|---|---|---|
| `U2` | edge ring | deviatoric edge-smoothed stiffness `V_h B~^T D_dev B~` (ES-FEM) |
| `U3` | `E A^c` support | volumetric chain `(K V_h) tbar^T sbar`, the F-bar term |
| `U4` | 4 | null base tetrahedron -- carries geometry only, no stiffness |

It is a displacement formulation -- no pressure degree of freedom, hence no
saddle-point solve -- that suppresses volumetric locking by combining
edge-based strain smoothing (ES-FEM) with the F-bar dilatation treatment and a
cyclic smoothing of the volumetric strain whose count `c` is the free
parameter. See `docs/fbar_es_fem_t4.md` for the formulation, its small-strain
reduction, the verification ladder and the Abaqus `C3D4H` comparison.

| File | Role |
|---|---|
| `u2edge.f` | edge-ring geometry: `V_h` and the smoothed gradients of eq. (1) |
| `u3vol.f` | the eq. (6)-(8) chain as a row walk; returns `V_h`, `sbar`, `tbar` |
| `e_c3d_u2.f` | `U2` stiffness |
| `e_c3d_u3.f` | `U3` stiffness, raises `nasym` (the operator is not symmetric) |
| `e_c3d_u4.f` | `U4` null base tet (zeroes the element matrices) |
| `resultsmech_u2.f` / `resultsmech_u3.f` | internal-force recovery through the *same* `u2edge`/`u3vol` calls |
| `resultsmech_u4.f` | `U4` null base (no internal force) |
| `fbar_lock.c` | pthread mutex for the one-time node->element map build |
| `fbares.py` | generator: retypes the base tets to `U4`, emits the `U2`/`U3` elements |

## Building CalculiX with the element

1. Copy the `*.f` and `*.c` sources into `<ccx_2.23>/src`, add them to
   `SCCXF`/`SCCXC` in `Makefile.inc`.
2. Apply the patches in `../patches_ccx/` (see its `README.md`).
3. Rebuild.

Two hard limits follow from the stencil widths: `*USER ELEMENT` caps
connectivity at 255 nodes, so `c = 1` is deliverable but `c >= 2` is not
without a direct global-assembly pass; and the non-symmetric tangent needs
PARDISO's general-matrix path (`mtype = 11`). Both are in the patch set.

## Generating a deck and running

```
*USER ELEMENT,TYPE=U2,NODES=<ring size>,INTEGRATIONPOINTS=1,MAXDOF=3
*USER ELEMENT,TYPE=U3,NODES=<support size>,INTEGRATIONPOINTS=1,MAXDOF=3
*USER ELEMENT,TYPE=U4,NODES=4,INTEGRATIONPOINTS=1,MAXDOF=3
```

`fbares.py` writes the whole deck from a `C3D4` deck; see its docstring and
`tests/slabconv.sh` for the end-to-end driver. No element in an F-bar deck
carries stress (the smoothing domains have no shape function of their own and
`U4` is null), so read the result from displacements and reactions.

## Verification

`tests/verify_fbar.py` (operator checks V1--V4), `tests/verify_fbar_nl.py`
(finite-strain checks N1--N4), `tests/verify_u3_chain.py` (the Fortran walk
against the Python operator) and `tests/verify_fbares_deck.py` (the generator)
verify the element against closed-form answers. `tests/smoothing_proto.py`
carries the Python prototype the Fortran was checked against.
