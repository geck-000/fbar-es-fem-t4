# Patches to CalculiX

Applied against a stock `ccx_2.23` source tree:

```bash
cd <ccx_2.23>
git apply -p0 <patch>
cd src && make -j8
```

| Patch | What it does | Default behaviour |
|---|---|---|
| `0001-iterative-tolerance.patch` | exposes the iterative solvers' convergence constant as `CCX_ITER_TOL`, and the iteration cap as `CCX_ITER_MAXIT` | unset → byte-for-byte stock |
| `0002-bbar-mean-dilatation.patch` | element-level B-bar (mean dilatation) above a Poisson-ratio threshold, `CCX_BBAR_NU` | unset → bit-for-bit stock |
| `0003-user-element-u4.patch` | wires the `U4` mixed tetrahedron in `../elements_ccx/` into the element and results dispatchers and the build | inert unless a deck declares `TYPE=U4` |

Read `0001` before trusting any iterative result: out of the box the solvers
stop at 0.5 % of the mean load and return a quietly wrong answer.

Read the header of `0002` before trusting it at all — it is a verified no-op on
`C3D4`, it only moves `C3D10` where the mesh is too coarse to resolve the
feature under test, and while it is active `equilibrium_gap` is meaningless
because only the stiffness is patched.
