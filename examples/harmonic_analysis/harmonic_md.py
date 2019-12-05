""" run harmonic md """

import os

from ase import units
from ase.calculators.lammpsrun import LAMMPS
from ase.md.verlet import VelocityVerlet

from vibes.io import read
from vibes.tdep.wrapper import parse_tdep_forceconstant
from vibes.helpers import progressbar

from vibes.molecular_dynamics.utils import FCCalculator, MDLogger


def lammps_si_tersoff_calculator(tmp_dir="./lammps"):
    """ create a lammps calculator for Si """
    lmp_path = os.getenv("LAMMPS_PATH")
    potential = os.path.join(lmp_path, "potentials", "Si.tersoff")
    files = [potential]
    parameters = {
        "mass": ["* 1.0"],
        "pair_style": "tersoff",
        "pair_coeff": ["* * " + potential + " Si"],
    }

    # New syntax introduces with https://gitlab.com/ase/ase/merge_requests/1000
    lammps = LAMMPS(parameters=parameters, files=files, tmp_dir=tmp_dir)

    return lammps


def run(
    maxsteps=1001,
    dt=5,
    harmonic=True,
    sample="geometry.in.supercell.300K",
    primitive="geometry.in.primitive",
    supercell="geometry.in.supercell",
    fc_file="infile.forceconstant",
    trajectory="trajectory.son",
):
    """ run Verlet MD, harmonic or force field """
    trajectory = trajectory
    atoms = read(sample)

    force_constants = parse_tdep_forceconstant(
        primitive=primitive,
        supercell=supercell,
        fc_filename=fc_file,
        two_dim=True,
        format="aims",
    )
    # force_constants.resize(2 * (3 * len(supercell),))

    supercell = read(supercell)
    if harmonic is True:
        calc = FCCalculator(supercell, force_constants)
    else:
        calc = lammps_si_tersoff_calculator()

    # generic md settings
    settings = {"atoms": atoms, "timestep": dt * units.fs}
    metadata = {"MD": {"fs": units.fs, "dt": dt}}

    md = VelocityVerlet(**settings)

    logger = MDLogger(atoms, trajectory, metadata=metadata, overwrite=True)

    atoms.calc = calc
    for _ in progressbar(range(maxsteps)):
        logger(atoms, info={"nsteps": md.nsteps, "dt": md.dt})
        md.run(1)


if __name__ == "__main__":
    run()
