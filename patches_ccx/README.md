# Patches to CalculiX

Applied against a stock `ccx_2.23` source tree, in order:

```bash
cd <ccx_2.23>
git apply -p0 <patch>
cd src && make -j8
```

| Patch | What it does |
|---|---|
| `0001-wide-user-element-matrix.patch` | widens the user-element assembly path from 60 to 765 DOF so the U3 volumetric stencil (up to 173 nodes at `c = 1`) fits |
| `0002-fbar-es-fem-t4-elements.patch` | wires `U2`/`U3`/`U4` into the element and results dispatchers and the build; adds the user-element digit list and the EVOL branch |
| `0003-asymmetric-user-element-path.patch` | opens ccx's asymmetric assembly/solve path to user elements (`U3` raises `nasym`), and stops truncating the node count at 127 |
| `0004-ring-support-reduction.patch` | stops allocating the identically-zero `(support − ring) × (support − ring)` block of `U3` |
| `0005-pardiso-out-of-core.patch` | opt-in out-of-core PARDISO for factors that do not fit in RAM |

The element sources themselves (`u2edge.f`, `u3vol.f`, `e_c3d_u2.f`,
`e_c3d_u3.f`, `e_c3d_u4.f`, `resultsmech_u2.f`, `resultsmech_u3.f`,
`resultsmech_u4.f`, `fbar_lock.c`) live in `../elements_ccx/` and are copied
into `src/` before building; add them to `SCCXF`/`SCCXC` in `Makefile.inc`.
See `../elements_ccx/README.md`.
