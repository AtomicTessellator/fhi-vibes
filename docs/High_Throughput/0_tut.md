# High-Throughput Calculations with FHI-vibes

During the course of the tutorial we primarily used individual files and the command line interface to set up and run all calculations.
While this is useful for in-depth studies it becomes impractical for high-throughput applications where one wants to run hundreds to thousands of materials using the same settings.
To accomplish this we developed a set of [FireWorks](https://materialsproject.github.io/fireworks/) workflows that automatically generates the correct calculations from a given set of input geometries and a workflow.in file.
After setting up FireWorks with the given guide, you will be able to use the high-throughput functionality within FHI-vibes.

Because of the cost of running Phonopy, Phono3py, and molecular dynamics calculations we find it easier to run FireWorks in the reservation mode.
Because of this we made some modifications to the basic FireWorks framework:

1. **Combined Launcher (`vibes claunch`): Run only some of the functions on queues**
    - A new launcher type that can run postporcessing steps locally on your machine, while using clusters to do all the force evaluations
    - Only functions launch_params.task2queue will be sent to clusters
    - Created a new function to launch a Firework remotely
2. **Allow GSS-API authorization with Kerberos**
    - Allows you to use GSS-API for passwordless login

## Setting up a workflow file

The input file for running a high-throughput workflow is similar to their command line counterparts, but with multiple task blocks defined within a single input file.
Because of this new sections need to be added in order to ensure a smooth operation of all the workflows, which are defined here

### The `fireworks` Block
Outside of the sections defined in `.fireworksrc` a few new terms need to be defined for the individual workflows:

* `fireworks.name`: name to prepend to all workflows defined in the workflow
* `fireworks.workdir.local`: Directory to run tasks on the local machine
* `fireworks.workdir.remote`: Directory to run tasks on remote hosts

### The `qadapter` Blocks
In addition to having multiple tasks to run, each task can optionally have its own `qadapter` section that dictates the computational resources that should be requested for it run properly.
These blocks include the following keywords:

* `nodes`: Number of nodes to request
* `ntasks_per_node`: Number of tasks per node
* `walltime`: Wallclock time for the calculation to run
* `queue`: The queue that you wish to run the jobs on
* `account`: The account to charge the hours too.

If any of these keywords is not included in the `qadapter` file then the defaults in the remote host's `my_qadapter.yaml` file will be used

### k-grid Optimization
Because these workflows are designed to run an many materials, optimizing the kgird for FHI-aims calculations for all materials maybe necessary. To accommodate this include the following section
```
[optimize_kgrid]
dfunc_min = 1e-9
```
Here `optimize_kgrid.dfunc_min` is the minimum change in energy necssary to say the k-grid is converged

### Relaxation
Relaxing the structures is largely the same as what was done with command line tools outside of two things:

 1. The option of using ASE or the electronic structure package relaxation algorithms
 2. The ability to do a multi step relaxation.

`relaxation.use_ase_relax` is the switch that controls which relaxation techniques are used.

Because the relaxation can be multi-step the relaxation terms are defined by numeric sub-features `relaxation.1`, `relaxation.2`, etc. each defining a particular step in the relaxation process in numerical order.
To do a single step relation, just use `relaxation.1`.
The terms defined in each block should either be code specific if `relaxation.use_ase_relax` is False, and the standard choices if `relaxation.use_ase_relax` is True.

An example of a single step calculation using ASE's optimizer is defined below, see `vibes/examples/high_throughput/workflow.in` for a multi-step example using FHI-aims.
```
[relaxation]
use_ase_relax = True

[relaxation.1]
driver = BFGS
fmax = 0.001
unit_cell = True
scalar_pressure = 0.01
hydrostatic_strain = True
```

### Phonopy Calculations
The inputs for phonopy are the same as those used in the command line utilities, with an additional term to control if supercell convergence: `phonopy.convergence`.
If `phonopy.convergence` is True then the supercell size will iteratively increase by until the Tanimoto similarity of the phonon density of states larger than 0.80.

To control either the convergence criteria or the base supercell matrix to iterate, make `phonopy.convergence` a block with the following terms
```
[phonopy.convergence]
minimum_similiarty_score= 0.85
sc_matrix_base = [-1,1,1,1,-1,1,1,1,-1]
```
where `minimum_similiarty_score` is the new criteria, and `sc_matrix_base` is the supercell matrix to determine the magnitude of the supercell size increase

For example of how the supercells increase if the following sections of an input file was given:
```
[phonopy]
supercell_matrix = [2, 2, 2]
convergence = True
```
Then the next supercell that will be tested would be `[4, 4, 4]`.
However if this was given instead
```
[phonopy]
supercell_matrix = [2, 2, 2]

[phonopy.convergence]
sc_matrix_base = [1, 1, 1]
```
Then the next supercell size calculated would be [3, 3, 3].
For a `sc_matrix_base` to be valid `supercell_matrix` must be an integer multiple of `sc_matrix_base`.

### Gruneisen Parameter Calculations
Once a phonopy calculation are complete, the Gruneisen parameter can be calculated from finite difference by including a gruneisen section
```
[gruneisen]
volume_factors = [0.99, 1.01]
```
where `volume_factors` describes the percentage difference in cell volume to calculate the Gruneisen parameters.

### Anharmonicity Quantification
The framework provides a means to quantify the anharmonicity from from statistical sampling using the `statistical_sampling` keyword
```
[statistical_sampling]
phonon_file = path/to/phonopy/trajectory.son
supercell_matrix = [-1,1,1,1,-1,1,1,1,-1]
temperatures = [300, 600]
debye_temp_fact = [1.0]
serial = True
n_samples = 1
plus_minus = True
mc_rattle = False
quantum = True
deterministic = True
zacharias = True
gauge_eigenvectors = True
ignore_negative = False
failfast = True
sobol = Faslse
random_seed = 13
propagate = False
```
Each of these key words corresponds to a term used to generate the samples from a phonopy calculation described in phonon_file.
If phonopy is also in the workflow then phonon_file is not needed and the model calculated would be used.
Each of the keywords correspond to the following terms

* `supercell_matrix`: The supercell matrix for the calculation, if not given use the one from the phonopy calculation. If the supercell matrix is different from the one in the phonon_file the phonopy force constants will be remapped onto the new supercell
* `temperatures`: list of temperatures to calculate the anharmonicity at
* `debye_temp_fact`: list of multipliers to add temperatures that are factors the materials Debye temperature
* `serial`: If True then do this in serial
* `n_samples`: number of samples to calculate for each temperature
* `plus_minus`: Use the deterministic sampling regime from Zacharias, et al
* `deterministic`: Populate all phonon modes to k_B T energy
* `gauge_eigenvectors`: Use a plus minus gauge for the eigenmodes
* `mc_rattle`: Rattle the structures using a Monte Carlo rattling method
* `quantum`: Populate phonon modes with a Bose-Einstein distribution
* `ignore_negative`: Ignore all imaginary modes
* `failfast`: If True fail if any imaginary modes are present or acustic modes are not near at Gamma
* `sobol`: If True use sobol quasi-random numbers
* `random_seed`: integer to seed the random number generator
* `propagate`: propagate the structure forward in time somewhat with ASE

### Molecular Dynamics
Finally a molecular dynamics calculation can be done using the same settings as the command line utilities.
There are only three minor changes to automatically generate the samples. All of these new parmeters are the same as the ones defined for the anharmonicity calculation.
```
[md]
phonon_file = path/to/phonopy/trajectory.son
temperatures = [300, 600]
supercell_matrix = [1, 1, 1]
```

## Adding and running workflows

To add a workflow to a particular LaunchPad use the following command
```
vibes fireworks add_wf -l LAUNCHPAD_FILE -w WORKFLOW_FILE
```
The default LAUNCHPAD_FILE is the same in your fw_config.py file (See FireWorks setup guide) and the default for WORKFLOW_FILE is workflow.in.

For running a workflow four options exist:

- `vibes fireworks claunch`: Run electronic structure calculations on clusters and everything else locally
- `vibes fireworks qlaunch`: Run all jobs on clusters using the queuing system
- `vibes fireworks rlaunch`: Run all jobs on locally
- The FireWorks utilities:

A more detailed description for each running option can be found in their respective documentation.
