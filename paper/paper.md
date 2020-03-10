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

The vibrational motion of nuclei determines many important properties of materials, including dynamical stability and transport [@Fultz2010]. Modeling the nuclear dynamics is therefore an important task for computational materials scientists in a broad range of sub-fields.

In _ab initio_ materials modeling, there are two main routes towards computing vibrational properties: _Lattice dynamics_ techniques start from the _harmonic approximation_ in which the interaction of atoms is approximated by a pairwise parabolic potential, allowing for an analytic solution of the equations of motion [@Dove1993]. On the other hand there are _molecular dynamics_ (MD) simulations in which the equations of motion are solved numerically by propagating the atoms at each time step with atomic forces determined from quantum chemistry methods, so that physical properties can be extracted as time averages of properly chosen observables [@Tuckerman2010].

There are codes that implement methods:

- phonopy [@Togo2015]
- i-pi [@Kapil2019]
- ase [@Larsen2017]
  
and quantum chemistry codes that provide energies and forces, but usually only a subset of implemented methods:

- vasp
- quantum espresso
- gpaw
- fhi-aims [@Blum2009]
- many others.

However, there are no frameworks that allow for an integrated production workflow with a unified syntax for inputs independet of the force calculator, and output. `FHI-vibes` is a `python` package that fills this gap by using `ase` as a backend in order to be able to represent materials and connect to force calculators, implement methods like relaxation and MD, and to connect to external codes like `phonopy`, using well defined input and output files. It can help to set up and run calculations on local machines and clusters using the `slurm` submission system, provides tools for performing standard postprocessing and analyzing raw data like MD trajectories, and ships utilities for manipulating and sharing results. For advanced analysis, it provides a `python` API fully compatible either with ``ase`` [@Larsen2017], or with `numpy` [@Walt2011], `pandas` [@McKinney2011], and `xarray` [@Hoyer2017].

`FHI-vibes` was used to produce the results in [@Knoop2020].

An extensive user guide including examples and a documentation of the API is available at
[`http://flokno.gitlab.io/hilde/`](https://flokno.gitlab.io/hilde/).

# 

# Acknowledgements

This work was supported by ERC ...

# References
