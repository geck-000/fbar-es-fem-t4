# F-barES-FEM-T4 tests

| Test | What it verifies |
|---|---|
| `verify_fbar.py` | operator checks V1--V4: unit row sums on `Q,P,E,R` and the chain `S = E A^c`; the patch test recovering `C1111 = K + 4G/3`; the `c = 0` collapse to selective ES-FEM-T4; the volumetric-constraint rank |
| `verify_fbar_nl.py` | finite-strain checks N1--N4: `f(0) = 0`, rigid translation, a 20% stretch patch, and convergence to the small-strain reduction |
| `verify_u3_chain.py` | the Fortran `u3vol.f` walk against the reference Python operator `S = E A^c` |
| `verify_fbares_deck.py` | the `fbares.py` deck generator, walked independently |
| `smoothing_proto.py` | the Python prototype of the smoothing operators the Fortran was checked against |
| `stability_modes.py` | spurious-mode census: distinguishes locking, stable, and unstable elements |
| `make_slabconv.py` / `slabconv.sh` / `report_slabconv.py` / `slabconv_extract.py` | the mesh-convergence cell and the Abaqus `C3D4H` comparison |
| `meshconv.sh` | same-geometry mesh sweep for the F-bar arm |

The two-element patch tests establish *consistency* (a uniform strain state is
reproduced to roundoff). `locking_sweep.sh`-style behaviour and the inf--sup
question are covered by the stability census and the mesh-convergence cell,
because an unstable element passes a patch test.
