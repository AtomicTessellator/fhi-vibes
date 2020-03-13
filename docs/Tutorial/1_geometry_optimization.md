!!! danger
	In this Tutorial, we assume you are already familiar with running *FHI-aims* calculations, and that you have [installed](../README.md#installation) and [configured](../README.md#configuration) `FHI-vibes` successfully.
	
!!! info
	For vibrational studies, it is crucial to use structures that are accurately  relaxed. Before starting with actual phonon calculations, we thus learn how to perform a geometry optimization with *FHI-vibes*

## <a name="1_GeometryOptimization"></a> Optimize Your Geometry

### Define Inputs

As input, we use an fcc-diamond Silicon crystal that we create with *ase*:
```python
"""run this in a script or in a python shell"""

from ase.build import bulk

si = bulk("Si")

si.write("geometry.in", scaled=True)
```

Next, we generate an input file for running a relaxation via the command line interface (CLI) of `FHI-vibes`:

```
vibes template relaxation
```

This should write a file called `relaxation.in` to your working directory. Inspect the file with an editor of your choice. The file should contain

```
[geometry]
file =                   geometry.in

[control]
xc =                     pw-lda
k_grid =                 [4, 4, 4]
use_symmetric_forces =   True

[basissets]
default =                light

[relaxation]
driver =                 BFGS
fmax =                   0.001
unit_cell =              True
```

You can start the calculation with `vibes run relaxation`. We suggest pipe the output, e.g., like this:

```
vibes run relaxation > vibes.relaxation.log &
```

`vibes` will create a working directory with the default name `relaxation` and will handle running the `aims` calculations and using a [straightforward BFGS algorithm implemented in ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#bfgs). You will find the converged structure in `relaxation/geometry.in.next_step`, and a summary of the relaxtion path in `relaxation/relaxation.log`.