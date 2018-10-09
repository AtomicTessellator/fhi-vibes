**Installation**

`python setup.py install --prefix .`

**Basic Setup**

`cp hilde.conf.template hilde.conf`

and edit according to system.

**New Features**
* Templates
  * `from hilde.templates.lammps import setup_lammps_si` to provide lammps calculator
* Brillouin zone helpers
  * `hilde.helpers.brillouinzone` features `get_paths`, `get_bands`, and
  `get_labels` to provide paths in the BZ that can be fed to `phonopy` via
  `phonon.set_bandstructure(bands)`, and
  `phonon.plot_band_structure(labels=labels)`.
  * These functions are used by `hilde.phonopy.plot_dos_and_bandstructure` to
  plot DOS and bandstructure in the working directory.
