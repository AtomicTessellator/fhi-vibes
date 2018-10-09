""" Provide a readily set up lammps calculator """

import os
from pathlib import Path
from ase.calculators.lammpsrun import LAMMPS

def setup_lammps_si(workdir):
    """Set up an ASE lammps calculator for silicon with Tersoff potential """
    # LAMMPS context information
    lmp_path = Path(os.getenv("LAMMPS_PATH"))
    potential = str(lmp_path / "potentials" / "Si.tersoff")
    files = [potential]
    parameters = {"mass": ["* 1.0"],
                  "pair_style": "tersoff",
                  "pair_coeff": ['* * ' + potential + ' Si']}

    # Logging
    lammps = LAMMPS(parameters=parameters,
                    files=files,
                    tmp_dir=workdir / 'lammps')

    return lammps
