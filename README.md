**Installation**

`python setup.py install --prefix .`

**Basic Setup**

`cp hilde.conf.template hilde.conf`

and edit according to system.

**New Features**
* Wrapper for `phono3py`
  * Preprocess and re-creation of Phono3py objects from precomputed force 
  constants, see examples
* Wrapper for `phonopy`
  * Preprocess and (some) postprocess, see examples
* Templates
  * `from hilde.templates.lammps import setup_lammps_si` to provide lammps calculator
* Brillouin zone helpers
  * `hilde.helpers.brillouinzone` features `get_paths`, `get_bands`, and
  `get_labels` to provide paths in the BZ that can be fed to `phonopy` via
  `phonon.set_bandstructure(bands)`, and
  `phonon.plot_band_structure(labels=labels)`.
  * These functions are used by `hilde.phonopy.plot_dos_and_bandstructure` to
  plot DOS and bandstructure in the working directory.
* Scripts:
  * `make_supercell`: create supercell from supercell matrix or
  target target
  * `geometry_info`: print geometry information for given input
  structure
* Symmetry Block Generation Functions
  * `AtomsInput`: A storage class that stores relevant information about a structure
  * `write_sym_constraints_geo`: Read any geometry.in file and use the list of `AtomInputs`
  to create a new supercell with a user defined symmetry block added to it
* FireWorks integration
  * Functions that can be used with PyTask to use FireWorks as a job manager

**Setup of Docker Image**
```
mkdir -p docker
cp Dockerfile requirements.txt docker
cd docker
docker login registry.gitlab.com
docker build -t registry.gitlab.com/floyd4k/hilde .
docker push registry.gitlab.com/floyd4k/hilde:latest   
```
