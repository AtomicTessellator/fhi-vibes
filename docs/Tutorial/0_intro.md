# Tutorial

In this tutorial, we introduce some of the functionality of `FHI-vibes` with hands-on examples.

!!! warning
	The settings used throughout the tutorials are chosen in order to allow for smooth calculations. They are _not_ sufficient for producing publication-ready scientific results.

!!! info
	We assume that you are familiar with running *FHI-aims* calculations, and that you have [installed](../README.md#installation) and [configured](../README.md#configuration) `FHI-vibes` successfully.

There are two test system used in the tutorial:

1. Lennard-Jones Argon, and
2. fcc-Silicon with LDA exchange-correlation functional.

Running the tutorial with LJ-Argon is great to get a quick hands-on overview over the features provided by `FHI-vibes`. Running the tutorial with LDA-Silicon will show you how to perform the calculations at full _ab initio_ quality.

## Test systems

### LJ-Argon

??? info "`geometry.in`"
    ```
    lattice_vector 0.0000000000000000 2.6299999999999999 2.6299999999999999 
    lattice_vector 2.6299999999999999 0.0000000000000000 2.6299999999999999 
    lattice_vector 2.6299999999999999 2.6299999999999999 0.0000000000000000 
    atom 0.0000000000000000 0.0000000000000000 0.0000000000000000 Ar
    ```
??? info "`calculator` section"
    ```
        [calculator]
        name:                          lj

        [calculator.parameters]
        # parameters for LJ Argon
        sigma:    3.405
        epsilon:  0.010325 
        rc:       8.0
    ```
    
### LDA-Silicon

??? info "`geometry.in"
    ```
    lattice_vector 0.0000000000000000 2.7149999999999999 2.7149999999999999 
    lattice_vector 2.7149999999999999 0.0000000000000000 2.7149999999999999 
    lattice_vector 2.7149999999999999 2.7149999999999999 0.0000000000000000 
    atom_frac 0.0000000000000000 0.0000000000000000 -0.0000000000000000 Si
    atom_frac 0.2500000000000000 0.2500000000000000 0.2500000000000000 Si
    ```

??? info "`calculator` section"
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

### Outline

- [Geometry optimization](1_geometry_optimization.md)
- [Phonon calculations](2_phonopy_intro.md)
- [Molecular dynamics](3_md_intro.md)
- [Harmonic sampling](4_statistical_sampling.md)
- [Anharmonicity quantification](5_anharmonicity_quantification.md)
- [High-Throughput workflows](../High_Throughput/Tutorial/0_phonopy.md)
