<a name="Running Multiple Phonopy Calculations"></a>

??? info "Prerequisite"
    FireWorks dependencies installed and set up for vibes. See the [installation guide](../Installation/0_setup.md) for more information.


## Setup workflow.in file

Setting up a high-throughput workflow to perform multiple `phonopy` calculation is similar to settings up a single calculation, but with a few additional steps to ensure the calculations are done in a similar manner.
Because the high throughput workflows are designed to be flexible, there is no `vibes template` command to automatically generate them, but modifying the workflows you'll work with here and in [the multi-step tutorial](../1_multistep) should be good guide on how to get started.

In this case only two materials (Si-diamond and MgO-rock salt) will be calculated, but the workflow can be used to generate a harmonic model an arbitrary number of materials.
Also in order to ensure these tutorials can be run locally, we will not fully converge the results, but will explain how to do that in the relevant sections.

As the first step in the tutorial, get the starting structures for the materials from the NOMAD repository  [Si](https://encyclopedia.nomad-coe.eu/gui/#/material/17241/structure) and [MgO](https://encyclopedia.nomad-coe.eu/gui/#/material/17916) and store them in `Si/geometry.in` and `MgO/geometry.in`, respectively.
*Note: Ideally these structures would already be relaxed, we are purposefully not doing that here. Information on how to do that will be explained in [the multistep tutorial](../1_multistep). Please do that after this tutorial is finished.*

??? info "`Si/geometry.in`"
    ```
    lattice_vector      0.000      2.703      2.703
    lattice_vector      2.703      0.000      2.703
    lattice_vector      2.703      2.703      0.000
    atom_frac      0.00     0.00     0.00 Si
    atom_frac      0.25     0.25     0.25 Si

    ```

??? info "`MgO/geometry.in`"
    ```
    lattice_vector    0.000   2.104   2.104
    lattice_vector    2.104   0.000   2.104
    lattice_vector    2.104   2.104   0.000
    atom_frac    0.00    0.00    0.00  Mg
    atom_frac    0.50    0.50    0.50  O
    ```
Once the geometry files are added, create a phonopy workflow in `workflow.in` with the following contents.

??? info "`workflow.in`"
    ```
    [files]
    geometries:                    */geometry.in

    [fireworks]
    name:                          example_phonon_calculations

    [fireworks.workdir]
    local:                         analysis/
    remote:                        run/

    [calculator]
    name:                          aims

    [calculator.parameters]
    xc:                            pw-lda

    [calculator.kpoints]
    density:                       1

    [calculator.basissets]
    default:                       light

    [calculator.socketio]
    port:                          12345

    [phonopy]
    supercell_matrix:              [-2, 2, 2, 2, -2, 2, 2, 2, -2]
    displacement:                  0.01
    is_diagonal:                   False
    is_trigonal:                   False
    is_plusminus:                  auto
    symprec:                       1e-05
    q_mesh:                        [45, 45, 45]
    serial:                        True

    [phonopy.convergence]
    minimum_similiarty_score:      0.05
    sc_matrix_base:                [-1, 1, 1, 1, -1, 1, 1, 1, -1]

    [phonopy.qadapter]
    nodes:                         1
    walltime:                      00-04:00:00
    ```

This `workflow.in` is very similar to the `phonopy.in` file from the [phonopy tutorial](../../../Tutorial/2_phonopy), but with a few extra sections and keywords.

### FireWorks sections

The largest difference between this workflow and the `phonopy.in` from the command line interface is the additional sections defining FireWorks specific parameters.
The purpose of these sections is to organize both the LaunchPad and file structure of the machines running the workflows.
`[fireworks.workdir]` specifies the base working directory for where the workflows will run (`fireworks.workdir.remote`) and where all the post-processing will happen (`fireworks.workdir.local`).
The naming convention was chosen to reflect that in most cases the electronic structure calculations will be done on a remote cluster, while the postprocessing will be done locally on a lap/desktop computer.
If `fireworks.workdir.remote` is not set, both post-processing and running the jobs will be done in the base working directory defined in `fireworks.workdir.local`.
Furthermore to ensure the MgO and Si calculations do not overwrite each other these base working directories are appended with `{material_chemical_formula}/{atoms_hash}/` giving each material a unique working directory to store all of its data in.

### qadapter

Another new section in this workflow is `[phonopy.qadapter]`.
This section is used to define job/queue specific information (number of nodes and wallclock time limits in this case), for a more complete definition of what to include here see the [documentation](../../Documentation/2_qadapter).
Defaults for all of these parameters will be stored in the `my_qadapter.yaml` file.

### `[phonopy]`
In the `[phonopy]` section there are only two changes: no `workdir` keyword and the `serial` keyword.
For these workflows all working directories have been standardized, so each task does not require its own `workdir` keyword.
The `serial` keyword controls how FireWorks sets up the force evaluation calculations for phonopy.
If it is True then all calculations will be run serially on a single job, while if it is False then each force evaluation will be done as a separate calculations and tasks.
We recommend this is always True to make use of ASE's [socektio calculators](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html), but it is not necessary to run the calculations.

### Phonopy Convergence

This section is used to optionally converge the harmonic model with respect to the supercell size.
Because `phonopy` uses finite differences to calculate the harmonic approximation, if the supercell is too small an atomic displacement can interact with a periodic image of itself leading to errors.
The magnitude of these errors is material dependent, and is normally either ignored by using what is considered a sufficiently large supercell or checked via a manual inspection of the phonon density of states and band structures.
For high throughput applications manual inspection is impractical so we developed an automatic metric for determining supercell convergence using the [Tanimoto similarity score](https://en.wikipedia.org/wiki/Jaccard_index#Other_definitions_of_Tanimoto_distance) of two successively larger supercells to check convergence.
The score is calculated from the DOS on a 45x45x45 q-grid of the smaller, $\mathbf{d_{small}}$, and larger, $\mathbf{d_{large}}$, supercell
\begin{equation}
    score = \frac{\mathbf{d_{small}} \cdot \mathbf{d_{large}}}{\left|\mathbf{d_{small}}\right|^2 + \left|\mathbf{d_{large}}\right|^2 - \mathbf{d_{small}} \cdot \mathbf{d_{large}} }.
    \label{eq:tanimoto}
\end{equation}

If the similarity score is larger than `phonopy.convergence.minimum_similarty_score` then the harmonic model is considered converged and that part of the workflow ends.
If it is smaller than the threshold then the supercell is increased to
\begin{equation}
    M_\text{S, new} = \left(n\right) M_\text{S, base} + M_\text{0}
    \label{eq:update_phonon}
\end{equation}
where $M_\text{S, base}$ is defined in [`phonopy.convergence.sc_matrix_base`](../../Documentation/5_phonons/#sc_matrix_base), $n$ the current phonopy iteration, and $M_\text{0}$ is the original supercell defined by `phonopy.supercell_matrix`.
$M_\text{0}$ must be an integer scalar value of $M_\text{S, base}$, or the workflow will not be added to the `LaunchPad`.

A score of 0.80 is considered to be a good balance between getting fully converged results and not going to very large supercell sizes; however, you may want to increase it if very accurate results are needed or lower it if the unitcell of a material is already very large.
Here a significantly lower minimum of 0.05 is used to allow the workflows to be easily run on a laptop and so that none of the supercells get above 64 atoms.

In this case we set the initial supercell to be a 2x2x2 supercell of the conventional cell in order to reduce the over number of calculations needed to calculate the converged phonon model.
Instead of explicitly calculating the smaller supercells the force constants from the larger one are mapped onto them and used to check for the convergence.
In most cases if a 200 atom supercell is used in this scheme then a converged phonon model can be calculated from a single phonopy calculation.

## Running the Calculations

Unlike what is done with the command line interface running these calculations is done in two steps: first adding the workflow to the LaunchPad and then second running all jobs in the LaunchPad.
To add a calculation to the LaunchPad run:
```
vibes fireworks add_wf
```

and you should get the following output
??? info "Adding Workflow Output"
    ```
    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``

    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``

    [calculator]   Update aims k_grid with kpt density of 1 to [4, 4, 4]
    [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
    [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
    [calculator]   .. add `output_level: MD_light` to parameters (default)
    [calculator]   Calculator: aims
    [calculator]   settings:
    [calculator]     xc: pw-lda
    [calculator]     k_grid: [4, 4, 4]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     output_level: MD_light
    [calculator]     compute_forces: True
    [calculator]     use_pimd_wrapper: ('localhost', 12345)
    [calculator]     aims_command: mpiexec -n 1 /home/purcell/git/fhi_aims/bin/ipi.aims.190906.scalapack.mpi.x
    [calculator]     species_dir: /home/purcell/git/fhi_aims/species_defaults/light
    [fireworks]    Generating workflow for MgO
    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``

    * Message from file vibes/context.py, line 57, function workdir:
    --> workdir not set, return `workdir``

    [calculator]   Update aims k_grid with kpt density of 1 to [4, 4, 4]
    [calculator]   .. add `sc_accuracy_rho: 1e-06` to parameters (default)
    [calculator]   .. add `relativistic: atomic_zora scalar` to parameters (default)
    [calculator]   .. add `output_level: MD_light` to parameters (default)
    [calculator]   Calculator: aims
    [calculator]   settings:
    [calculator]     xc: pw-lda
    [calculator]     k_grid: [4, 4, 4]
    [calculator]     sc_accuracy_rho: 1e-06
    [calculator]     relativistic: atomic_zora scalar
    [calculator]     output_level: MD_light
    [calculator]     compute_forces: True
    [calculator]     use_pimd_wrapper: ('localhost', 12345)
    [calculator]     aims_command: mpiexec -n 1 /home/purcell/git/fhi_aims/bin/ipi.aims.190906.scalapack.mpi.x
    [calculator]     species_dir: /home/purcell/git/fhi_aims/species_defaults/light
    [fireworks]    Generating workflow for Si
    ```

If you are using a non-default launchpad file (normally `my_launchpad.yaml`) or a workflow not defined in `workflow.in` use the `-l` and `-w` flags respectively.

To confirm that both the workflows are in the LaunchPad run
```
lpad get_wflows
```
and you should see:
```
[
    {
        "state": "READY",
        "name": "example_phonon_calculations_Si_92b4ed53a4691c7621606216cab915fa8d5ac311--3",
        "created_on": "YYYY-MM-DDTHH:MM:SS",
        "states_list": "W-REA"
    },
    {
        "state": "READY",
        "name": "example_phonon_calculations_MgO_294b61fcf3f0b3b0c067796c380cade4e2153b30--1",
        "created_on": "YYYY-MM-DDTHH:MM:SS",
        "states_list": "W-REA"
    }
]
```
These workflows each have two jobs: one for setting up the calculation and one for performing the analysis.
The additional analysis step is done for creating a trajectory file if the calculations are running in parallel, and creating a trajectory locally if running with a combined launcher.

Once the workflows have been added to the LaunchPad you have multiple options to run it:

- `vibes fireworks claunch`: Run electronic structure calculations on clusters and everything else locally
- `vibes fireworks qlaunch`: Run all jobs on clusters using the queuing system
- `vibes fireworks rlaunch`: Run all jobs on locally
- The FireWorks utilities:

A more detailed description for each running option can be found in their respective documentation.

For this example we'll use `rlaunch`. After running
```
vibes fireworks rlaunch rapidfire
```
the output from
```
lpad get_wflows
```
should now be
```
[
    {
        "state": "COMPLETED",
        "name": "example_phonon_calculations_MgO_294b61fcf3f0b3b0c067796c380cade4e2153b30--1",
        "created_on": "2020-04-10T11:08:43.213000",
        "states_list": "C-C-C-C-C-C"
    },
    {
        "state": "COMPLETED",
        "name": "example_phonon_calculations_Si_92b4ed53a4691c7621606216cab915fa8d5ac311--3",
        "created_on": "2020-04-10T11:08:43.292000",
        "states_list": "C-C-C-C-C-C"
    }

]
```
Four additional tasks have now been added to each workflow corresponding to the force evaluation for the original supercell and three more tasks for the larger supercell for convergence tests.
Now all the workflows have been completed let's analyze the results.

## Analyzing the Results
The first step in analyzing the results is understanding the file structure of the directories.
First looking at the run directory there are the following folders:
```
run/MgO/294b61fcf3f0b3b0c067796c380cade4e2153b30:
sc_natoms_64

run/Si/92b4ed53a4691c7621606216cab915fa8d5ac311:
sc_natoms_64
```
The  `sc_natoms_64` folders each contain a `phonopy` directory where all the files generated from respectively running $n=2$ `phonopy` iterations.

Looking in the analysis folder there are the following folders
```
analysis/MgO/294b61fcf3f0b3b0c067796c380cade4e2153b30:
converged  sc_natoms_64

analysis/Si/92b4ed53a4691c7621606216cab915fa8d5ac311:
converged  sc_natoms_64
```
The `sc_natoms_64` folders each contain a `phonopy_analysis` directory that only contain the final phonopy `trajectory.son` file for those iterations.
Additionally the last `phonopy` iteration are stored in `converged` for easy access to the converged phonon calculations.
From here you can perform any analysis that is possible within `phonopy` on all the materials, and know that the results standardized with respect to all of the calculation parameters.
For example you can see the bandstructure and DOS of both materials by running
```
vibes output phonopy -bs -dos
```
in each of the `converged` folders.
Because of the large variety of possible analysis steps, there is no automated phonopy output scripts in the workflow, but the file structure can be used to easily make bash or python scripts to do all post-processing.

While we can now get converged phonon results, these workflows are incomplete because they require pre-relaxed structures to get physically relevant results.
It would be possible to separately relax all the structures and then use those in this workflow, but the easier solution would be to use [multi-step workflows](../1_multistep) in the next tutorial.
