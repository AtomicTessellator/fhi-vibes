---
title: 'FHI-vibes: _Ab initio_ Vibrational Simulations'
tags:
  - Python
  - Physics
  - phonons
  - transport
authors:
  - name: Florian Knoop
    orcid: 0000-0002-7132-039X
    affiliation: 1
  - name: Thomas A. R. Purcell
    orcid: 0000-0003-4564-7206
    affiliation: 1
  - name: Matthias Scheffler
    affiliation: 1
  - name: Christian Carbogno
    orcid: 0000-0003-0635-8364
    affiliation: 1
affiliations:
 - name: Fritz Haber Institute of the Max Planck Society, Berlin, Germany
   index: 1
date: XX March 2020
bibliography: paper.bib
---

# Summary

The vibrational motion of nuclei determines many important properties of materials, including their dynamical stability and transport [@Fultz2010]. Modeling the nuclear dynamics is therefore an important task for computational materials scientists in a broad range of sub-fields.

In _ab initio_ materials modeling, there are two main routes towards computing vibrational properties: _Lattice dynamics_ techniques start from the _harmonic approximation_ in which the interaction of atoms is approximated by a pairwise parabolic potential, allowing for an analytic solution of the equations of motion [@Dove1993]. On the other hand there are _molecular dynamics_ (MD) simulations in which the equations of motion are solved numerically by propagating the atoms at each time step with atomic forces determined from quantum chemistry methods _without_ making approximations to the interatomic potential. Physical properties can then be extracted as time averages of properly chosen observables [@Tuckerman2010].

The tasks of i) computing interatomic forces and ii) solving the nuclear equations of motion to obtain physical properties are often completely independent of each other. Moreover, task i) is typically several orders of magnitude more demanding computationally then task ii). It is therefore possible—and often advisable—to use different software packages for each task and ‘’glue‘’ them together with a flexible language like _python_. This philosophy has been adopted in the past by software packages like *ase* [@Larsen2017], *phonopy* [@Togo2015], and *i-pi* [@Kapil2019]. These packages implement methods for modeling vibrational properties of solids and molecules, but leave the computationally heavy calculation of interatomic forces to specialized force-field packages like *LAMMPS*[@Plimpton1993], or efficient quantum chemistry codes like *VASP*[@Kresse1996], *QuantumEspresso*[@Giannozzi2009], *FHI-aims*[@Blum2009], and many other codes. However, there are no frameworks that allow for an integrated production workflow for a range of vibrational simulation techniques with a unified syntax for inputs and output, independet of the force calculator. 

`FHI-vibes` is a `python` package that fills this gap by using `ase` as a backend in order to be able to represent materials and connect to force calculators, implement methods like geometry optimization and MD, and to connect to external codes like `phonopy` and `phono3py`[cite] that implement lattice dynamics techniques, while providing defined input and output files and a command line interface. It can help to set up and run calculations on local machines and clusters using the `slurm` submission system, provides tools for performing standard postprocessing and analyzing raw data like MD trajectories, and ships utilities for manipulating and sharing results. For advanced analysis, it provides a `python` API fully compatible either with ``ase`` [@Larsen2017], or with `numpy` [@Walt2011], `pandas` [@McKinney2011], and `xarray` [@Hoyer2017]. Furthermore, `FHI-vibes` provides a connection to *FireWorks*[@Jain2015], a workflow management system for running simulation workflows on extensive sets of materials in _high-throughput_ fashion. `FHI-vibes` is tightly integrated with *FHI-aims* to perform energy and force calculations, but extending the functionality to any calculator available via *ase* is straightforward.

`FHI-vibes` was used to produce the results in [@Knoop2020].

An extensive user guide including a comprehensive list of the available features, examples and a documentation of the API is available at [`flokno.gitlab.io/hilde/`](https://flokno.gitlab.io/hilde/).

# Acknowledgements

This work was supported by ERC ...

# References