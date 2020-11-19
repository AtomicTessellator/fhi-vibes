---
title: 'FHI-vibes: _Ab Initio_ Vibrational Simulations'
tags:
  - Python
  - Physics
  - Phonons
  - Transport
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

# Introduction

The vibrational motion of nuclei determines many important properties of materials, including their thermodynamic equilibrium and non-equilibrium properties. Accurately assessing the nuclear dynamics and the associated material properties is therefore an important task for computational materials scientists in a broad range of sub-fields. Of particular importance are simulation techniques that build on first-principles electronic-structure simulations and thereby allow to systematically investigate the virtually infinite space of materials, including those systems for which little or no experimental data is hitherto available [@Curtarolo2013]. This allows one
to design novel and improved materials with optimal properties for many applications, e.g., high-performance thermal insulators for gas and airplane turbines [@Evans2008], organic semicondcutors with long-term phase stabilities [@Salzillo2016], thermoelectric generators [@Snyder2008], and improved thermal management systems [@Tian2019].

Essentially, there are two distinct routes towards assessing vibrational properties:
In perturbative _lattice dynamics_ techniques, the potential-energy surface on which the nuclei move is approximated with a Taylor expansion around the equilibrium structure. 
Typically, one starts from a second-order expansion, i.e., the _harmonic approximation_, which allows for an analytic solution of the equations of motion [@Dove1993] 
and thus for a straightforward evaluation of observables (thermodynamic expectation values). Higher-order terms in the Taylor expansion can be accounted for perturbatively.
Conversely, _molecular dynamics_ (MD) based approaches account for the full, non-perturbative potential-energy surface _without_ approximating the actual interactions. This requires
one to solve the equations of motion numerically by propagating the atoms in time; physical properties can then be extracted as time  averages of properly chosen observables [@Tuckerman2010].
Although both _lattice dynamics_ and _molecular dynamics_ techniques aim at computing the same physical observables, the involved methodologies, formalisms, and challenges are quite different.
Accordingly, both methodologies also have different strengths and weaknesses: For instance, performing and analyzing MD simulations is typically computationally and conceptually more challenging,
whereas perturbative lattice dynamics calculations inherently rely on approximations that are hard to validate.

To date, a variety of different software packages exists at different degrees of sophistication in both fields. Prominent examples are the _phonopy_ code [@Togo2015] for performing _lattice dynamics_
calculations using Parlinski's finite-difference formalism [@Parlinski1997] and the _i-PI_ code [@Kapil2019] for performing classical MD and quantum-mechanical path-integral MD simulations.
Both packages interface with a variety of first-principles codes like *VASP* [@Kresse1996], *QuantumEspresso* [@Giannozzi2009], *Abinit* [@Gonze2020], *FHI-aims* [@Blum2009], and several others.


# Statement of need
To date, there is no software solution that allows for the seamless bridging and interlinking of _lattice dynamics_ and  _molecular dynamics_ based approaches, despite the fact that actual material science studies can profit in accuracy and efficiency by exploiting both approaches. For instance,
potential use cases include 
accelerating _MD_ calculations by starting from harmonic equilibrium configurations [@West2006], 
analyzing _MD_ simulations in terms of harmonic phonons [@Turney2009], 
investigating the range of validity of the perturbative expansion used in _lattice dynamics_ [@Knoop2020], 
and overcoming finite size- and time-effects in _ab initio_ Green Kubo simulations of the thermal conductivity [@Carbogno2016]. 
Given the wide opportunities for application, the aspect of _integration_, i.e., the ability to utilize different methodologies from distinct codes in a flexible fashion using a consistent user interface, is paramount. In particular, this is a prerequisite for automatizing these workflows to enable hierarchical high-throughput screening of whole material classes in a systematic fashion. For example, such a workflow would start from geometry optimizations followed by a study of harmonic properties for many materials, so to single out candidate materials for more involved, fully anharmonic aiMD simulation techniques. Along these lines, let us mention that providing descriptive input and output files is a prerequisite for sharing raw data and results in a transparent and interpretable way in the spirit of open science and the FAIR Principles [@Draxl2018]. On top of that, tracking the provenance [@AiiDA] across different codes facilitates the repurposing and analysis of the obtained data.

To facilitate the described scientific studies, _FHI-vibes 1.0_ offers the following main features:
*Add list from readme, maybe even longer, just make sure that we can tick all points/problems (NOT THE APPLICATIONS!) described above.*
More advanced features and workflows will be added in future releases.

# Summary

_FHI-vibes_ is a _python_ package that allows for such an integrated workflow. 
It uses the _Atomistic Simulation Environment (ASE)_ [@Larsen2017] as a backend in order to represent materials and to connect to various first-principles codes. Via _ASE_, _FHI-vibes_ provides a flexible framework for geometry optimization and MD, and connects to external codes like _spglib_ [@Togo2018], _phonopy_ [@Togo2015], _phono3py_ [@Togo2015b], and _hiphive_ [@Eriksson2019] that implement lattice dynamics techniques based on the harmonic approximation. For all these tasks, _FHI-vibes_ provides defined input files and a command line interface to set up and run calculations on local machines and clusters using the _slurm_ submission system. The output is organized in self-contained and descriptive output files that enable a straightforward exchange of the data obtained with different methodologies.
For advanced analysis, it provides an API fully compatible  with _ASE_ as well as _numpy_ [@Walt2011], _pandas_ [@McKinney2011], and _xarray_ [@Hoyer2017]; several user-friendly utilities allow to perform the most common postprocessing tasks within the command-line interface, such as providing comprehensive summaries of MD simulations or phonon calculations.

_FHI-vibes_ provides a connection to *FireWorks* [@Jain2015], a workflow management system for running simulation workflows on extensive sets of materials in _high-throughput_ fashion. _FHI-vibes_ is tightly integrated with *FHI-aims* [@Blum2009] to perform energy and force calculations, but extending the functionality to any calculator available via *ASE* is straightforward.

An extensive user guide including a comprehensive list of the available features, tutorials, and a reference documentation is available at [`vibes.fhi-berlin.mpg.de`](http://vibes.fhi-berlin.mpg.de/).

_FHI-vibes_ was used to produce the results in [@Knoop2020].

# Acknowledgements
The authors would like to thank Roman Kempt and Marcel Hülsberg for testing and providing valuable feedback. F.K. would like to thank Marcel Langer and Zhenkun Yuan for feedback and Ask Hjorth Larsen for valuable discussions. T.P. would like to thank the Alexander von Humboldt Foundation for their support through the Alexander von Humboldt Postdoctoral Fellowship Program. This project was supported by TEC1p (the European Research Council (ERC) Horizon 2020 research and innovation programme, grant agreement No. 740233), BigMax (the Max Planck Society’s Research Network on Big-Data-Driven Materials-Science), and the NOMAD pillar of the FAIR-DI e.V. association.

# References
