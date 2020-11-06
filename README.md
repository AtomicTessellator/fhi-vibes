FHI-vibes
===

Welcome to `FHI-vibes`, a `python` package for calculating, analyzing, and understanding the vibrational properties of anharmonic solids from first principles. `FHI-vibes` is intended to bridge between different methodologies, so to allow for a seamless assessment of vibrational properties with different approaches, ranging from the harmonic approximation to anharmonic MD. `FHI-vibes` builds on several [existing packages](https://vibes-developers.gitlab.io/vibes/Credits/) and interfaces them in a consistent and user-friendly fashion. 

Its main features are:

- Geometry optimization via [ASE](https://wiki.fysik.dtu.dk/ase/ase/optimize.html#module-ase.optimize),
- harmonic phonon calculations via [Phonopy](https://atztogo.github.io/phonopy/),
- molecular dynamics simulations in [NVE](https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-nve-simulations-the-microcanonical-ensemble), [NVT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin), and [NPT](https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.nptberendsen) ensembles,
- [harmonic sampling](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.115504), and
- [anharmonicity quantification](https://arxiv.org/abs/2006.14672).

Most of the functionality is high-throughput ready via [fireworks](https://materialsproject.github.io/fireworks/#).

## Overview

- [Installation](https://vibes-developers.gitlab.io/vibes/Installation/)
- [Tutorial](https://vibes-developers.gitlab.io/vibes/Tutorial/0_intro/)
- [Documentation](https://vibes-developers.gitlab.io/vibes/Documentation/0_intro/)
- [Credits](https://vibes-developers.gitlab.io/vibes/Credits/)
- [References](https://vibes-developers.gitlab.io/vibes/References/)


## News

- [Our anharmonicity measure got published.](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.083809)
- [the best is yet to come](https://www.youtube.com/watch?v=B-Jq26BCwDs)

## Changelog

#### v1.0.0a10

- Enable conversion of trajectories to `ase.io.Trajectory` files for viewing with ASE [(!37)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/37)
- Important fix for running NPT dynamics [(!36)](https://gitlab.com/vibes-developers/vibes/-/merge_requests/36)
- We have a changelog now!