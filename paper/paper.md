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
date: July 2020
bibliography: paper.bib
---

# Summary

The vibrational motion of nuclei determines many important properties of materials, including their dynamical stability and transport [@Fultz2010]. Modeling the nuclear dynamics is therefore an important task for computational materials scientists in a broad range of sub-fields.

In _ab initio_ materials modeling, there are two main routes towards computing vibrational properties: _Lattice dynamics_ techniques that start from the _harmonic approximation_ and those that use _molecular dynamics_ (MD) simulations. In the _harmonic approximation_ the interaction of atoms is approximated by a pairwise parabolic potential, allowing for an analytic solution of the equations of motion [@Dove1993]. On the other hand the (MD) approach numerically solves the equations of motion by propagating the atoms at each time step with forces determined from quantum chemistry methods _without_ making approximations to the interatomic potential. Physical properties can then be extracted as time averages of properly chosen observables [@Tuckerman2010].

The tasks of i) computing interatomic forces and ii) solving the nuclear equations of motion to obtain physical properties are often completely independent of each other. Moreover, the first task is typically several orders of magnitude more computationally demanding than the second. It is therefore possible—and often advisable—to use different software packages for each task and ‘’glue‘’ them together with a flexible language like _python_. This philosophy has been adopted in the past by software packages like *ase* [@Larsen2017], *phonopy* [@Togo2015], and *i-pi* [@Kapil2019]. These packages implement methods for modeling vibrational properties of solids and molecules, but leave the computationally heavy calculation of interatomic forces to specialized force-field packages like *LAMMPS* [@Plimpton1993], or efficient quantum chemistry codes like *VASP* [@Kresse1996], *QuantumEspresso* [@Giannozzi2009], *FHI-aims* [@Blum2009], and many other codes. However, there are no frameworks that allow for an integrated production workflow for a range of vibrational simulation techniques with a unified syntax for inputs and output, independent of the force calculator.

`FHI-vibes` is a `python` package that fills this gap by using `ase` as a backend in order to be able to represent materials and connect to force calculators, implement methods like geometry optimization and MD, and connect to external codes like `phonopy` [@Togo2015] and `phono3py` [@Togo2015b] that implement lattice dynamics techniques. For all these tasks `FHI-vibes` provides defined input and output files and a command line interface. It can help to set up and run calculations on local machines and clusters using the `slurm` submission system, provides tools for performing standard postprocessing and analyzing raw data like MD trajectories, and ships utilities for manipulating and sharing results. For advanced analysis, it provides a `python` API fully compatible either with ``ase`` [@Larsen2017] or a combination of `numpy` [@Walt2011], `pandas` [@McKinney2011], and `xarray` [@Hoyer2017]. Furthermore, `FHI-vibes` provides a connection to *FireWorks* [@Jain2015], a workflow management system for running simulation workflows on extensive sets of materials in _high-throughput_ fashion. `FHI-vibes` is tightly integrated with *FHI-aims* [@Blum2009] to perform energy and force calculations, but extending the functionality to any calculator available via *ase* is straightforward.

An extensive user guide including a comprehensive list of the available features, examples and a documentation of the API is available at [`vibes.fhi-berlin.mpg.de`](http://vibes.fhi-berlin.mpg.de/).

`FHI-vibes` was used to produce the results in [@Knoop2020].

# Acknowledgements
T.P. would like to thank the Alexander von Humboldt Foundation for their support through the Alexander von Humboldt Postdoctoral Fellowship Program. This project was supported by TEC1p (the European Research Council (ERC) Horizon 2020 research and innovation programme, grant agreement No. 740233), BigMax (the Max Planck Society’s Research Network on Big-Data-Driven Materials-Science), and the NOMAD pillar of the FAIR-DI e.V. association.

# References
