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

The input file for running a high-throughput workflow is similar to their command line counterparts, but with an extra section titled general.
The general blocks represent a set of options that alters the workflow on a global level and includes the following keywords:

- **workdir_cluster:** The root directory to store files on a cluster
- **workdir_local:** The root directory to store files on your local machine
- **basisset:** The default basis set to use for all materials (Not stored in the same location as the individual step locations so that can change for pre-relaxation)
- **relax_structure:** If True relax all structures
- **use_tight_relax:** If True relax the structure with tight basis settings and use that structure for all future calculations.
- **opt_kgrid:** If True converge the k-grid with respect to the total energy of a material
- **kgrid_dfunc_min:** Threshold to converge the total energy change to

All of the working directories for the calculations represent parent directories with the actual files stored in workdir/empirical_chemical_formula/atoms_hash/task_name/ to prevent accidental overwriting of files.

Once the general section is added to the file, every step that you want to perform in the workflow must be directly added to workflows.
Additionally if you do not want to use the default job settings for running the jobs on a cluster corresponding task_qadapter settings are needed
For example if you want to compare the values of the anharmonicity quantification at 300 and 600 K for a set of FCC structures using the one-shot approach and molecular dynamics you'd have the following skeleton for your workflow.in file
```
[general]
workdir_cluster = /path/to/base/working/directory/on/cluster/
workdir_local = /path/to/base/working/directory/on/local/machine/
basisset = light
relax_structure = True

[geometry]
files = regular/expression/to/get/all/geometries.in

[control]
# Electronic structure settings
xc =                 pw-lda

[control_kpt]
# Standard k-point density to use for all calculations
density = 1.00

[phonopy]
supercell_matrix = [-2, 2, 2, 2, -2, 2, 2, 2 -2]
converge_phonons = True
conv_crit = 0.80
sc_matrix_original = [-1, 1, 1, 1, -1, 1, 1, 1 -1]
get_gruniesen = False

[statistical_sampling]
# phonon_file = /path/to/phonopy/trajectory.son
temperatures = [300, 600]
n_samples = 1
plus_minus = True

[md]
# phonon_file = /path/to/phonopy/trajectory.son
supercell_matrix = [-2, 2, 2, 2, -2, 2, 2, 2, -2]
driver = Langevin
timestep = 5
temperatures = [300, 600]
friction = 0.03
maxsteps = 2000
compute_stresses = 10
logfile = md.log

[light]

[phonopy_qadapter]
walltime = 4:00:00
nodes = 2

# [statistical_sampling_qadapter]
# walltime = 04:00:00
# nodes = 2

[md_qadapter]
walltime = 24:00:00
nodes = 2
```
This workflow has many of the same keywords as their command line counterparts, but phonopy has four additional keywords:

- **converge_phonons:** If True calculate the harmonic with a converged supercell as determined by the [Tanimoto similarity](https://en.wikipedia.org/wiki/Jaccard_index#Other_definitions_of_Tanimoto_distance) between successive supercell sizes
- **sc_matrix_original:** Base supercell matrix to scale at each step (supercell_matrix must be an integer scale multiple of sc_matrix_original)
- **conv_crit:** Minimum similarity score to be considered converged
- **get_gruniesen:** Automatically scale the volume and run the corresponding phonopy calculations to calculate the gruneisen parameters from finite differences

All of these keywords will automatically add new phonopy calculations with an updated qadatper for the increased computational cost.
If a phonopy is being calculated in the workflow, then the md/statistical sampling steps will automatically use the converged phonon model (if no supercell convergence just use the single calculation) otherwise you must pass the phonopy trajectories via the phonon_file keyword.
If no statistical_sampling_qadapter is given it will use the last phonopy_qadapter for all calculations.
All other keywords are the same as the command line utilities.

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
