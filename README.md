# F-barES-FEM-T4 for CalculiX

A locking-free 4-node tetrahedral element for [CalculiX](http://www.calculix.de/),
implementing the F-barES-FEM-T4 formulation of Onishi, Iida and Amaya
(*Accurate viscoelastic large deformation analysis using F-bar aided edge-based
smoothed finite element method for 4-node tetrahedral meshes (F-BarES-FEM-T4)*,
Int. J. Comput. Methods **15**(7) 1845003, 2018, DOI 10.1142/S0219876218450032).

The element is delivered as two `*USER ELEMENT` types over the underlying T4
connectivity:

| type | role |
|---|---|
| `U2` | deviatoric edge-ring stiffness (symmetric ES-FEM block) |
| `U3` | volumetric chain `E A^c` (the cyclically smoothed F-bar term) |

It is a **displacement formulation** --- no pressure degree of freedom, hence no
saddle-point solve --- that suppresses volumetric locking by combining
edge-based strain smoothing (ES-FEM) with the F-bar dilatation treatment and a
cyclic smoothing of the volumetric strain whose count `c` is the method's free
parameter. On a purpose-built cell that resolves the near-incompressible
inclusion layer, F-barES-FEM-T4 (`c=1`) and Abaqus `C3D4H` share a
Richardson-extrapolated limit to **0.06%**, while plain `C3D4` remains 4% stiff
at the finest mesh.

## Layout

| path | contents |
|---|---|
| `elements_ccx/` | Fortran element sources, the `fbares.py` deck generator, the verification/validation tests, and `docs/fbar_es_fem_t4.md` |
| `patches_ccx/` | patches to a stock `ccx_2.23` source tree |

## Building CalculiX with the element

1. Copy the `*.f` sources from `elements_ccx/` into the CalculiX `src/`
   directory and add them to `SCCXF` in `Makefile.inc`.
2. Apply the patches in `patches_ccx/` (see its `README.md`).
3. Rebuild.

The F-barES-FEM-T4 element itself is patches `0009` (elements, dispatchers,
build), `0010` (asymmetric Petrov--Galerkin path), `0011` (ring-x support
reduction) and `0012` (out-of-core PARDISO). The earlier patches are the
dependency chain and the precursor MINI (`U4`) and nodal-B-bar (`U5`/`U6`)
elements documented in `elements_ccx/README.md`.

## Generating a deck and running

`elements_ccx/fbares.py` rewrites a plain `C3D4` deck as F-barES-FEM-T4; see its
docstring and `elements_ccx/tests/slabconv.sh` for the end-to-end driver. Two
hard limits follow from the stencil widths: `*USER ELEMENT` caps connectivity at
255 nodes, so `c=1` is deliverable but `c>=2` is not without a direct
global-assembly pass; and the non-symmetric tangent needs PARDISO's
general-matrix path (`mtype = 11`).

## Verification and validation

* `elements_ccx/tests/verify_fbar.py` --- operator checks V1--V4 (unit row sums,
  patch test recovering `C1111 = K + 4G/3`, the `c=0` collapse to selective
  ES-FEM-T4, the volumetric-constraint rank).
* `elements_ccx/tests/verify_fbar_nl.py` --- finite-strain checks N1--N4.
* `elements_ccx/tests/verify_u8_chain.py` --- the Fortran volumetric walk
  against the reference Python operator `S = E A^c`.
* `elements_ccx/tests/locking_sweep.sh`, `stability_modes.py` --- locking and
  spurious-mode behaviour.
* `elements_ccx/tests/make_slabconv.py`, `slabconv.sh`, `report_slabconv.py`,
  `slabconv_extract.py` --- the mesh-convergence cell and the Abaqus `C3D4H`
  comparison.

The full formulation, implementation notes and validation results are in
`elements_ccx/docs/fbar_es_fem_t4.md`.
