# Tutorial

!!! warning "Warnings"

	- All tutorials assume you have a background in (_ab initio_) vibrational modeling. The background theory is written down to establish a common notation and cannot replace professional training by any means.
	- The settings used throughout the tutorials are chosen in order to allow for smooth calculations. They are _not_ sufficient for producing publication-ready scientific results.
	- We assume you have [installed](../Installation.md) and [configured](../Installation.md#configuration) `FHI-vibes` successfully.


In this tutorial, we introduce the functionality of `FHI-vibes` with hands-on examples.

### Outline

The following tutorials are available:

- [Geometry optimization](1_geometry_optimization.md)
- [Phonon calculations](2_phonopy_intro.md)
- [Molecular dynamics](3_md_intro.md)
- [Harmonic sampling](4_statistical_sampling.md)
- [Anharmonicity quantification](5_anharmonicity_quantification.md)
- [High-Throughput workflows](../High_Throughput/Tutorial/0_phonopy.md)

## Test systems for the tutorials

!!! info
	We assume that you are familiar with running *FHI-aims* for performing _ab initio_ calculations.

There are two test system available for running the tutorials:

1. fcc-Silicon with LDA exchange-correlation functional, and 
2. Lennard-Jones Argon.

Running the tutorial with LDA-Silicon will show you how to perform the calculations at full _ab initio_ quality. Running the tutorial with LJ-Argon is great to get a quick hands-on overview over the features provided by `FHI-vibes`.

### LDA-Silicon

??? info "Geometry input file `geometry.in`"
    ```
    lattice_vector 0.0000000000000000 2.7149999999999999 2.7149999999999999 
    lattice_vector 2.7149999999999999 0.0000000000000000 2.7149999999999999 
    lattice_vector 2.7149999999999999 2.7149999999999999 0.0000000000000000 
    atom_frac 0.0000000000000000 0.0000000000000000 0.0000000000000000 Si
    atom_frac 0.2500000000000000 0.2500000000000000 0.2500000000000000 Si
    ```

??? info "`calculator` section for task in put file"
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
    ```

### LJ-Argon

??? info "Geometry in put file `geometry.in`"
    ```
    lattice_vector 0.0000000000000000 2.6299999999999999 2.6299999999999999 
    lattice_vector 2.6299999999999999 0.0000000000000000 2.6299999999999999 
    lattice_vector 2.6299999999999999 2.6299999999999999 0.0000000000000000 
    atom 0.0000000000000000 0.0000000000000000 0.0000000000000000 Ar
    ```
??? info "`calculator` section for task input file"
    ```
        [calculator]
        name:                          lj

        [calculator.parameters]
        # parameters for LJ Argon
        sigma:    3.405
        epsilon:  0.010325 
        rc:       8.0
    ```
