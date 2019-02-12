## Create samples with displacements and velocities from force constants

Input files: `force_constants.dat` is a 3Nx3N matrix written with `np.savetxt()`,
holding the force constants for the geometry file, `geometry.in.supercell` in this case.

From the input files, create new samples with, e.g.,

```
create_samples geometry.in.supercell 300 -fc force_constants.dat  --deterministic -n 2
```
