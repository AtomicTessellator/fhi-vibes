""" prepare molecular dynamics simulations using the ASE classes """

from pathlib import Path
from ase import units as u
from hilde.trajectory.md import last_from_yaml
from .workflow import run_md
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, PhononHarmonics


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

    temp = temperature * u.kB
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

    if not prepared:
        atoms = initialize_md(atoms, temperature, **kwargs)

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


def initialize_md(atoms, temperature, force_constants=None, quantum=False, **kwargs):
    """ Either use Maxwell Boltzmann or PhononHarmonics to prepare the MD run """

    print(f"Prepare MD run at temperature {temperature}")

    if "force_constants" is not None:
        print("Initialize positions and velocities using force constants.")
        force_constants = np.loadtxt(force_constants)
        PhononHarmonics(
            atoms, force_constants, quantum=quantum, temp=temperature * u.kB
        )
    else:
        print("Initialize velocities according to Maxwell-Boltzmann distribution.")
        MaxwellBoltzmannDistribution(atoms, temp=temperature * u.kB)

    return atoms

