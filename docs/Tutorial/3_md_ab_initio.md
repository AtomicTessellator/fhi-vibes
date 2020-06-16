# _ab initio_ Molecular Dynamics

!!! info
	We will now introduce _ab initio_ Molecular Dynamics simulations where the atomic forces come from a first-principles code. We will use `FHI-aims` as this calculator in the following. We assume:

	- You have familiarized yourself with running MD simulations with FHI-vibes for the LJ-Argon toy model. 
	- You are familiar with running `FHI-aims` calculations.
	- Optionally: You are familiar with running `FHI-aims` calculations on a workstation or supercomputer.

In principle, _ab initio_ molecular dynamics simulations are identical to MD simulations using an empirical force field like Lennard-Jones, only that the forces are computed from first principles, e.g., using density functional theory (DFT) with LDA or GGA xc-functional in the Born-Oppenheimer approximation. Thus,

$$
\begin{align}
\mathcal V ({\bf R}) = E_{\rm tot}^{\rm DFT} ({\bf R})~,
\label{eq:PES}
\end{align}
$$

where the potential energy $\mathcal V({\bf R})$ of a given atomic configuration $\bf R$ is given by the total energy  $E_{\rm tot}^{\rm DFT} ({\bf R})$ of the electronic and nuclear system computed for this structure.

There are however some technical differences, mostly because a single force evaluation is typically several orders of magnitude more expensive when performed by a DFT code compared to an empirical force field.[^footnote1]

## Setting up and ab initio MD

Setting up an _ab initio_ MD with FHI-vibes is in principle similar to setting up a force field MD. We will use an 8 atoms supercell of silicon with LDA xc-functional as previously in the tutorial on geometry optimization and phonon calculations. You can re-use the structure from [the tutorial on geometry optimization](1_geometry_optimization.md), as well as the calculator setup.

```
cp ../path/to/relaxation/geometry.in.next_step geometry.in.primitive
```

Create a supercell with 8 atoms:

```
vibes utils make-supercell geometry.in.primitive -n 8
mv geometry.in.primitive.supercell_8 geometry.in.supercell
```

Pre-thermalize supercell and create input file
```
vibes utils create-samples geometry.in.supercell -T 300
mv geometry.in.supercell.0300K geometry.in
```

Adjust your MD input file:

??? info "`md.in`"
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
    
    [md]
    driver =           Langevin
    timestep =         4
    temperature =      300
    friction =         0.02
    maxsteps =         2500
    
    [files]
    geometry:                      geometry.in
    primitive:                     geometry.in.primitive
    supercell:                     geometry.in.supercell
    ```

The flag `compute_stresses` in the section `[md]` will make FHI-aims compute the _ab initio_ stress every 10 steps during the MD simulation. This will provide access to pressure as earlier.[^footnote2]

## Run a calculation

This step is similar to [before](3_md_canonical_sampling.md#run-the-calculation), i.e., you run

```
vibes run md >> log.md &
```

### Submit calculation on a cluster

To efficiently perform _ab initio_ molecular dynamics simulations for systems larger than a few atoms, you will need a workstation or access to a supercomputer. To submit a `vibes` simulation to your supercomputer, follow these steps:

1. [Install `FHI-vibes` on your supercomputer](../../#installation),
2. set up a calculation as you have done earlier on your laptop,
3. submit the `vibes run` command to the queue.

??? info "Example `submit.sh` for `slurm` queue manager"
    ```
    #!/bin/bash -l

    #SBATCH -J md|vibes
    #SBATCH -o log/md.%j
    #SBATCH -e log/md.%j
    #SBATCH --mail-type=all
    #SBATCH --mail-user=your@mail.com
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=32
    #SBATCH --ntasks-per-core=1
    #SBATCH -t 24:0:00

    vibes run md md.in
    ```


## Postprocess
You can perform postprocessing of the pressure as explained [earlier](3_md_postprocess.md). Be aware that the simulation time is shorter when discarding the thermalization period.

??? info "reference pressure"
	For 8 atoms LDA-Silicon, you should get a potential pressure of $-0.61 \pm 0.06 {}$

## References
Running the calculation will take some time depending on the computer your working with. You find references [in our reference repository](https://gitlab.com/vibes-developers/vibes-tutorial-files/-/tree/master/3_molecular_dynamics/ab_initio). There you also find reference calculations for 64 and 216 atoms.

[^footnote1]: There are other differences like the inherent incompleteness of the SCF cycle in Kohn-Sham DFT schemes. Incomplete SCF convergence can introduce a systematic error that introduces energy drifts and other unphysical effects. Choosing the correct convergence settings is an important aspect of performing _ab initio_ MD simulations and is highly materials specific. Devising a strategy on how to choose these settings goes beyond the scope of this tutorial.

[^footnote2]: We compute the stress only every 10th step because computing the stress is numerically more expensive in _ab initio_ than computing the atomic forces only. Since we know that consecutive samples generated during MD are highly correlated, we don't loose valuable information by computing this quantity not in every single step.