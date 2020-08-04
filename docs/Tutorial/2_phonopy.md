<a name="2_Phonopy"></a>

!!! info "Prerequisites"

	- For vibrational studies, it is crucial to use structures that are accurately  relaxed. Before starting with actual phonon calculations, make sure you are familiar with [geometry optimization](1_geometry_optimization.md).
	- Create a new working directory and copy over the `geometry.in.next_step` file you obtained from the previous geometry optimization as your new `geometry.in` file.



## Perform a phonon calculation

Setting up a `phonopy` calculation is similar to settings up a `relaxation` (or any other workflow supported by `FHI-vibes`). To ensure that our physical settings don't change, we will copy the `relaxation.in` obtained in the [previous part of the tutorial](1_geometry_optimization.md) to the new working directory and rename it to `phonopy.in`. Please delete the `relaxation` specific sections `[relaxation]` and `[relaxation.kwargs]` and add settings for a phonopy calculation by running

```
vibes template phonopy >> phonopy.in
```

??? info "`phonopy.in`"
	```
	[calculator]
    name:                          aims
    
    [calculator.parameters]
    xc:                            pw-lda
    
    [calculator.kpoints]
    density:                       2
    
    [calculator.basissets]
    default:                       light
    
    [calculator.socketio]
    port:                          12345
    
    [phonopy]
    supercell_matrix:              [1, 1, 1]
    displacement:                  0.01
    is_diagonal:                   False
    is_plusminus:                  auto
    symprec:                       1e-05
    q_mesh:                        [45, 45, 45]
    workdir:                       phonopy
    ```

Obviously the most important section in the `phonopy.in` input file is `[phonopy]` which containts information about how the supercells with displacements should be set up to compute the force constants from the [finite-differences method](0_intro.md#Phonons). An explanation for the full list of keywords is found in the [documentation](../Documentation/phonopy.md). The most important two are explaned in the following:

### Supercell Matrix (`supercell_matrix`)

The supercell matrix $M_{\rm S}$ given as `supercell_matrix` will be used to [generate the lattice of the supercell from the lattice of the primitive unitcell by matrix multiplication:](https://phonopy.github.io/phonopy/phonopy-module.html#supercell-matrix)

$$
\begin{align}
	\require{mediawiki-texvc}
	\def\t#1{\text{#1}}
	\begin{pmatrix}
		\mathbf a_\t{S}^\t{t} \\ \mathbf b_\t{S}^\t{t} \\ \mathbf c_\t{S}^\t{t}
	\end{pmatrix}
	=
	M_\t{S} \cdot
	\begin{pmatrix}
	\mathbf a_\t{u}^\t{t} \\ \mathbf b_\t{u}^\t{t} \\ \mathbf c_\t{u}^\t{t}
	\end{pmatrix}
	 ~.
	 \label{eq:smatrix}
\end{align}
$$

Here, $\mathbf a_\t{u}^\t{t}$, $\mathbf b_\t{u}^\t{t}$, $\mathbf c_\t{u}^\t{t}$ are the transposed lattice vectors (row-vectors) of the (primitive) unit cell and $\mathbf a_\t{S}$, $\mathbf b_\t{S}$, $\mathbf c_\t{S}$ label the lattice vectors of the supercell respectively. `supercell_matrix` can be given in any shape that lets itself transform trivially to a $3 \times 3$-matrix. For example, `[1, 1, 1]` gets transformed to the $3 \times 3$ unit matrix.

### Displacement (`displacement`)

The `displacement` tag will set the amplitude of the finite displacement in $\AA$. The same parameter is called [`DISPLACEMENT_DISTANCE` in `phonopy`](https://phonopy.github.io/phonopy/setting-tags.html#displacement-distance). In principle, this is a numerical parameter that needs to be optimized. A smaller `discplacement` results in a better approximation to the true second derivative of the potential. However, a too small displacement generates too small forces that can be severely affected by other sources of computational noise, e.g., finite grids etc. For production purposes, the default value of $d = 0.01\,\AA$ usually works quite well and with a properly set up force calculator, there is no need to increase the displacement further.

### Run the calculation

Let's stick to the default settings in `phonopy.in` for the moment and run the calculation with

```
vibes run phonopy | tee log.phonopy 
```

The calculation should take only a few seconds (depending on you computer).


CC: This Restart and submission settings are too hidden and also not really applicable here, given that there is only 1 displacement: 
I suggest to make separate section "Running/Restarting on Clusters" after the Anharmonicity/before the High-throughput Section and
just point to it from here and in other sections 

#### Restart a calculation

If you need to restart a calculation, e.g., because you're working on a cluster and your job does not fit into a walltime or stops for another reason before the total number of simulation steps is reached, you can simply re-run or re-submit `vibes run phonopy`. It will restart the calculation from the last completed step. 

#### Automatic restarts

Optionally, `vibes` can restart the job by itself using a `[restart]` section in `phonopy.in`. To this end, add

```
[restart]
command = sbatch submit.sh
```

to your `phonopy.in`, where `sbatch submit.sh` is the command you use to submit the phonon calculation to the queue.

### Postprocessing

The `vibes run` command takes care that all _ab initio_ calculations are performed, but some additional, rapid postprocessing
is needed to obtain the phonon-related quantities. The postprocessing  itself can be performed interactivley with

```
vibes output phonopy phonopy/trajectory.son --full
```


CC:UPDATE OUTPUT BELOW:

??? note "Terminal output"
	```
    [phonopy.postprocess] Start phonopy postprocess:
    [trajectory]   Parse `phonopy/trajectory.son`
    [son] read file:  phonopy/trajectory.son
    [son] process:    |||||||||||||||||||||||||||||||||||||  2/2
    [trajectory]   .. create atoms
    [progress]        |||||||||||||||||||||||||||||||||||||  1/1
    [trajectory]   .. done in 0.001s
    [phonopy.postprocess] .. done in 0.034s
    [phonopy.postprocess] 
    Extract phonopy results:
    [phonopy.postprocess] .. q_mesh:   [45, 45, 45]
    [phonopy.postprocess] .. write force constants
    [phonopy.postprocess] Extract basic results:
    [phonopy.postprocess] .. write primitive cell
    [phonopy.postprocess] .. write supercell
    [phonopy.postprocess] .. write force constants to FORCE_CONSTANTS
    [phonopy.postprocess] Extract bandstructure
    [phonopy.postprocess] .. write yaml
    [phonopy.postprocess] .. plot
    [phonopy.postprocess] .. all files written to phonopy/output in 1.113s
    * Message from file vibes/phonopy/postprocess.py, line 123, function check_negative_frequencies:
        --> Negative frequencies found at G = [0 0 0]:

    # Mode   Frequency
        1 -1.62419e-07 THz
        2 -1.15231e-07 THz
    [phonopy.postprocess] 
    Frequencies at Gamma point:
    q = [0. 0. 0.] (weight= 1)
    # Mode   Frequency
        1   -0.0000002 THz
        2   -0.0000001 THz
        3    0.0000002 THz
        4   15.7649077 THz
        5   15.7649078 THz
        6   15.7649079 THz
    ```
This will:

CC: I removed the internal work not interesting for beginners.

- Compute the phonon bandstructure along CC:WHICH PATH and save it in `phonopy/output/bandstructure.pdf`.
- Compute the density of states using a $45 \times 45 \times 45$ $\bf q$ point grid and the Tetrahedron method. 
  The density of states will be plotted alongside the bandstructure to a file `output/bandstructure_dos.pdf`, and written to a data file [`total_dos.dat`](https://phonopy.github.io/phonopy/output-files.html#total-dos-dat-and-projected-dos-dat).
  The q-grid can be adjusted by specifying it with an additional flag `--q_mesh`.
CC: My phonpy.in already contains   q_mesh:  [45, 45, 45] What is the preferred method that can be used to change qmesh?
- Compute the harmonic free energy $F^{\rm ha}$ and the harmonic heat capacity at constant volume, $C_V$, i.e., the thermal properties accessible in the harmonic approximation using the DOS and it q-point settings.  
  An overview plot is saved to `output/thermal_properties.pdf` and the detailed output is written to [`output/thermal_properties.yaml`](https://phonopy.github.io/phonopy/output-files.html#thermal-properties-yaml).
- CC: What about the animations?

CC: Add a sentence/paragraph that explains when/why the AI calculations need to be run again and when not (e.g. for the q-mesh)

CC: Add all discussed plots below:
??? info "Bandstructure"
	![image](bandstructure.png)
	
**Congratulations!** You have just performed a full (but not yet converged!) _ab initio_ phonon bandstructure calculation.

Note that the CLI also allows to only run a subset of the postprocessing, e.g.,
```
vibes output phonopy phonopy/trajectory.son -v -bs
```
only outputs the bandstructure.

## Choosing a supercell size

!!! info
	The ideal supercell size and shape depends on your problem at hand and it is difficult to give definite advice. In practice, the supercell size needs to be converged until the target property of interest is not changing anymore.  
        To facilitate this, there is a CLI tool that can help you creating supercells of different sizes.

There is a [CLI utility](../Documentation/cli.md#vibes-utils)  in`FHI-vibes` that can help you to find supercells of different sizes:

```
vibes utils make-supercell
```

For example

```
vibes utils make-supercell geometry.in -n 8
```

will find the conventional, cubic cell of silicon with 8 atoms:

```
...
Settings:
  Target number of atoms: 8

Supercell matrix:
 python:  [-1,  1,  1,  1, -1,  1,  1,  1, -1]
 cmdline: -1 1 1 1 -1 1 1 1 -1
 2d:
[[-1, 1, 1],
 [1, -1, 1],
 [1, 1, -1]]

Superlattice:
[[5.42906529 0.         0.        ]
 [0.         5.42906529 0.        ]
 [0.         0.         5.42906529]]

Number of atoms:  8
  Cubicness:         1.000 (1.000)
  Largest Cutoff:    2.715 AA
  Number of displacements: 1 (1)

Supercell written to geometry.in.supercell_8
```

It will tell you the supercell matrix that you can use in `phonopy.in` (`python:  [-1,  1,  1,  1, -1,  1,  1,  1, -1]`), the generated superlattice, a "cubicness" score based on the filling ratio of the largest sphere fitting into the cell, the largest cutoff in which any neighbor is not a periodic image of a closer neighbor to estimate boundary effects, and the number of supercells with displacements that  `phonopy` will create. It will also write the structure to `geometry.in.supercell_8` which you can inspect, e.g., with `jmol`.

To run a calculation for such a supercell, one just needs to...
CC: Describe and show results as above.

CC: Add Warning: Please save this results since they are needed for the latter tutorial Harmonic Sampling.

### Practical guideline

In practice, the convergence with supercell size needs always to be checked carefully, since it depends on the range of the interactions present in your system, e.g., long ranged unscreened van-der-Waals interactions require larger supercells than short-ranged covalent ones as here in Si. Along the same lines, the acceptable supercell size depends also on the properties you are interested in. Free energies and specific heats converge faster than individual frquenecies.  Using a cubic-as-possible supercell shape and playing around with `vibes utils make-supercell` and a little bit of experience will do the job. For the example at hand, it might for instance be instructive to check the convergence of different properties for the supercell sizes A, B, and C. CC: Add A, B, C
As a reference,a supercell of CxCxC gives very nicely converged band structures for the system inspected here:

CC: Add output/image for that!

