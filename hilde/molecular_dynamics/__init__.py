""" prepare molecular dynamics simulations using the ASE classes """

from pathlib import Path
from ase import units as u
from hilde.trajectory.md import last_from_yaml
from .workflow import run_md


def setup_md(
    atoms,
    algorithm="Verlet",
    temperature=None,
    timestep=None,
    friction=None,
    logfile=None,
    workdir=".",
    trajectory=None,
    **kwargs,
):
    """ create and ase.md object with respective settings """

    if trajectory is None:
        trajectory = (Path(workdir) / "trajectory.yaml").absolute()
    else:
        trajectory = Path(trajectory).absolute()

    if not Path(workdir).is_dir():
        Path(workdir).mkdir()

    temp = temperature * u.fs
    dt = timestep * u.fs

    if "verlet" in algorithm.lower():
        from ase.md.verlet import VelocityVerlet

        md = VelocityVerlet(atoms, timestep=dt, logfile=logfile)

    elif "langevin" in algorithm.lower():
        from ase.md.langevin import Langevin

        md = Langevin(
            atoms, temperature=temp, timestep=dt, friction=friction, logfile=logfile
        )

    else:
        raise RuntimeError(f"Molecular dynamics mode {algorithm} is not suppported.")

    prepared = prepare_from_trajectory(atoms, md, trajectory)

    return atoms, md, prepared


def prepare_from_trajectory(atoms, md, trajectory="trajectory.yaml", **kwargs):
    """ Take the last step from trajectory and initialize atoms + md accordingly """

    trajectory = Path(trajectory).absolute()
    if trajectory.exists():
        last_atoms = last_from_yaml(trajectory)
        md.nsteps = last_atoms["MD"]["nsteps"]

        atoms.set_positions(last_atoms["atoms"]["positions"])
        atoms.set_velocities(last_atoms["atoms"]["velocities"])
        return True
    print(f"** {trajectory} does  not exist, nothing to prepare")
