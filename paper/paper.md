---
title: 'FHI-vibes: _Ab initio_ Vibrational Simulations'
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
In perturbative _lattice dynamics_ techniques, the potential-energy surface on which the nuclei move is approximated with a Taylor expansion around an equilibrium structure.
Typically, one starts from a second-order expansion, i.e., the _harmonic approximation_, which allows for an analytic solution of the equations of motion [@Dove1993]
and thus for a straightforward evaluation of observables (thermodynamic expectation values). Higher-order terms in the Taylor expansion can be accounted for perturbatively.
Conversely, _molecular dynamics_ (MD) based approaches account for the full, non-perturbative potential-energy surface _without_ approximating the actual interactions. This requires
one to solve the equations of motion numerically by propagating the atoms in time; physical properties can then be extracted as time  averages of properly chosen observables [@Tuckerman2010].
Although both _lattice dynamics_ and _molecular dynamics_ techniques aim at computing the same physical observables, the involved methodologies, formalisms, and challenges are quite different.
Accordingly, both methodologies also have different strengths and weaknesses: For instance, performing and analyzing MD simulations is typically computationally and conceptually more challenging,
whereas perturbative lattice dynamics calculations inherently rely on approximations that are hard to validate.

To date, a variety of different software packages exists at different degrees of sophistication in both fields. Prominent examples are the _phonopy_ code [@Togo2015] for performing _lattice dynamics_
calculations using Parlinski's finite-difference formalism [@Parlinski1997] and the _i-PI_ code [@Kapil2019] for performing classical MD and quantum-mechanical path-integral MD simulations.
Note that both packages interface with a wide variety of first-principles codes like *VASP* [@Kresse1996], *QuantumEspresso* [@Giannozzi2009], *Abinit* [@Gonze2020], *FHI-aims* [@Blum2009], and many other. There is, however, no software solution that allows for the seamless bridging of _lattice dynamics_ and  _molecular dynamics_ based approaches, despite the fact that actual material science studies can profit in accuracy
and efficiency by exploiting both approaches. For instance, _MD_ simulations can be used to proof the validity of the approximations on which the _lattice dynamics_ approach relies [@Knoop2020]; similarly,
_MD_ simulations can be sped up and interpreted more efficiently using techniques and data from _lattice dynamics_ calculations [@West2006].


# Statement of need

Although bridging and interlinking between calculations performed with _lattice dynamics_ and _MD_ techniques is useful to accelerate, validate, and analyze material science studies, there are to date no frameworks that allow to tackle this task with integrated, customizable, and easily extendable workflows using a unified syntax for inputs and outputs. For instance, this allows to accelerate _MD_ calculations by starting from harmonic equilibrium configurations [@West2006], 
to analyze _MD_ simulations in terms of harmonic phonons [@Turney2009], to investigate the range of validity of the perturbative expansion used in _lattice dynamics_ [@Knoop2020], and to efficiently overcome finite size- and time-effects in _ab initio_ Green Kubo simulations of the thermal conductivity [@Carbogno2016].
The aspect of _integration_ becomes all the more important when one aims at (semi-)automatizing these workflows to enable hierarchical high-throughput screening of material classes in a systematic fashion. For example, such a workflow would start from a simple study of harmonic properties for many materials to single out candidate materials for more involved, (fully) anharmonic simulation techniques, thereby ensuring _accuracy_ and _efficiency_ at the same time.
Last but not least, providing descriptive input and output files is a prerequisite for sharing raw data and results in a transparent and interpretable way in the spirit of open science and the FAIR Principles [@Draxl2018]. Especially _ab initio_ MD (aiMD) simulations are a valuable resource beyond the scope of the studies they are performed for, as they provide a wealth of information about the potential-energy landscape. This information can be _repurposed_ by parameterizing models of the potential function by means of, e.g.,  (deep) Neural Networks [@Behler2007], or compressed sensing techniques [@Zhou2014], from which new properties _not_ included in the initial study can be obtained at a drastically reduced computational cost.

# Summary

_FHI-vibes_ is a _python_ package that allows for such an integrated workflow. It uses the _Atomistic Simulation Environment (ASE)_[@Larsen2017] as a backend in order to represent materials and to connect to various first-principles codes, to implement methods like geometry optimization and MD, and to connect to external codes like _phonopy_ [@Togo2015], _phono3py_ [@Togo2015b], and _hiphive_ [@Eriksson2019] that implement lattice dynamics techniques. For all these tasks, _FHI-vibes_ provides defined input files and a command line interface to set up and run calculations on local machines and clusters using the _slurm_ submission system. The output is organized in self-contained and descriptive output files that enable a straightforward exchange of calculated data, especially of aiMD trajectories in the spirit of the above discussion.
For advanced analysis, it provides an API fully compatible  with _ASE_ as well as _numpy_ [@Walt2011], _pandas_ [@McKinney2011], and _xarray_ [@Hoyer2017], while several utilities allow for standard postprocessing such as providing comprehensive summaries of MD simulations or phonon calculations.

_FHI-vibes_ provides a connection to *FireWorks* [@Jain2015], a workflow management system for running simulation workflows on extensive sets of materials in _high-throughput_ fashion. _FHI-vibes_ is tightly integrated with *FHI-aims* [@Blum2009] to perform energy and force calculations, but extending the functionality to any calculator available via *ASE* is straightforward.

An extensive user guide including a comprehensive list of the available features, tutorials, and a reference documentation is available at [`vibes.fhi-berlin.mpg.de`](http://vibes.fhi-berlin.mpg.de/).

_FHI-vibes_ was used to produce the results in [@Knoop2020].

# Acknowledgements
The authors would like to thank Roman Kempt and Marcel Hülsberg for testing and providing valuable feedback. F.K. would like to thank Marcel Langer and Zhenkun Yuan for feedback and Ask Hjorth Larsen for valuable discussions. T.P. would like to thank the Alexander von Humboldt Foundation for their support through the Alexander von Humboldt Postdoctoral Fellowship Program. This project was supported by TEC1p (the European Research Council (ERC) Horizon 2020 research and innovation programme, grant agreement No. 740233), BigMax (the Max Planck Society’s Research Network on Big-Data-Driven Materials-Science), and the NOMAD pillar of the FAIR-DI e.V. association.

# References
