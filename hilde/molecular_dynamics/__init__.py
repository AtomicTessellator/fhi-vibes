""" prepare molecular dynamics simulations using the ASE classes """

from ase import units as u
from hilde.helpers.converters import input2dict, results2dict


def metadata2dict(atoms, calc, md):
    """Convert metadata information to plain dict

    Parameters
    ----------
    atoms: ASE Atoms Object
        Reference structure for MD calculations
    calc: ASE Calculator
        Calculator for the MD Run
    md: ASE Molecular Dynamics Obect
        MD propagator
    """

    md_dict = md.todict()
    # save time and mass unit
    md_dict.update({"fs": u.fs, "kB": u.kB, "dt": md.dt, "kg": u.kg})

    return {"MD": md_dict, **input2dict(atoms, calc)}


from .workflow import run_md
from .initialization import initialize_md
