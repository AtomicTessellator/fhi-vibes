FHI-vibes
===

Welcome to `FHI-vibes`, a `python` package for calculating, analyzing, and understanding the vibrational properties of solids from first principles. `FHI-vibes` is intended to seamlessly bridge between the harmonic approximation and fully anharmonic molecular dynamics simulations. To this end, `FHI-vibes` builds on several [existing packages](https://vibes-developers.gitlab.io/vibes/Credits/) and interfaces them in a consistent and user-friendly fashion.

In the documentation and tutorials, knowledge of first-principles electronic-structure theory as well as proficiency with _ab initio_ codes such as [FHI-aims](https://aimsclub.fhi-berlin.mpg.de/) and high-performance computing are assumed. Additional experience with python, the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/), or [Phonopy](https://atztogo.github.io/phonopy/) is helpful, but not needed.

`FHI-vibes` provides the following features:

- Geometry optimization via [ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-ase.optimize),
- harmonic phonon calculations via [Phonopy](https://atztogo.github.io/phonopy/),
- molecular dynamics simulations in [NVE](https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-nve-simulations-the-microcanonical-ensemble), [NVT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin), and [NPT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen) ensembles,
- [harmonic sampling](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.115504),
- [anharmonicity quantification](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.083809), and
- thermal transport with [ab initio Green-Kubo](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.224304) and [Phono3py](Documentation/phono3py.md)

Most of the functionality is high-throughput ready via [fireworks](https://materialsproject.github.io/fireworks/#).


## Overview

- [Installation](https://vibes-developers.gitlab.io/vibes/Installation)
- [Tutorial](https://vibes-developers.gitlab.io/vibes/Tutorial/0_intro)
- [Documentation](https://vibes-developers.gitlab.io/vibes/Documentation/0_intro)
- [Credits](https://vibes-developers.gitlab.io/vibes/Credits)
- [References](https://vibes-developers.gitlab.io/vibes/References)

## Changelog

#### v1.1.0

- Update ab initio Green Kubo method and [tutorials](https://vibes-developers.gitlab.io/vibes/Tutorial/7_green_kubo/). [(!100)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/100)
- Update Phono3py interface and [documentations](https://vibes-developers.gitlab.io/vibes/Documentation/phono3py/).
- Heatflux unit updated from `eV/AA^2/ps` to `eV/AA^2/fs`.

#### v1.0.5

- bugfix

#### v1.0.4

- looser dependencies
- [various bug fixes](https://gitlab.com/vibes-developers/vibes/-/merge_requests?scope=all&state=merged)

#### v1.0.3

- update dependencies to allow `phonopy` versions up to 2.8

#### v1.0.2

- First official release after passing the [JOSS review](https://github.com/openjournals/joss-reviews/issues/2671).
- Several additions to the documentation.

#### v1.0.0a10

- Enable conversion of trajectories to `ase.io.Trajectory` files for viewing with ASE [(!37)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/37)
- Important fix for running NPT dynamics [(!36)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/36)
- We have a changelog now!
