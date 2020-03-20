""" prepare molecular dynamics simulations using the ASE classes """

from pathlib import Path

import numpy as np
from ase import units as u
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, PhononHarmonics

from vibes.filenames import filenames
from vibes.helpers.fileformats import last_from_yaml
from vibes.helpers.warnings import warn


def setup_md(
    atoms,
    driver="Verlet",
    temperature=None,
    timestep=None,
    friction=None,
    logfile=None,
    workdir=".",
    trajectory_file=None,
    **kwargs,
):
    """Create and ase.md object with respective settings

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Reference structure for molecular dynamics
    driver: str
        Algorithm used to propagate the MD
    temperature: float
        Temperature used for Langevin propagator
    timestep: float
        MD timestep
    friction: float
        friction used for Langevin propagator
    logfile: str
        file to log in
    workdir: str or Path
        The working directory
    trajectory_file: str or Path
        The output trajectory file

    Returns
    -------
    atoms: ase.atoms.Atoms
        The refrence structure
    md: ase.md.MolecularDynamics
        The MD propagator
    prepared: bool
        True if prepared from trajectory

    Raises
    ------
    RuntimeError
        If driver is not supported
    """

    if trajectory_file is None:
        trajectory_file = (Path(workdir) / filenames.trajectory).absolute()
    else:
        trajectory_file = Path(trajectory_file).absolute()

    if not Path(workdir).is_dir():
        Path(workdir).mkdir()

    dt = timestep * u.fs
    md = None

    if "verlet" in driver.lower():
        from ase.md.verlet import VelocityVerlet

        md = VelocityVerlet(atoms, timestep=dt, logfile=logfile)

    elif "langevin" in driver.lower():
        from ase.md.langevin import Langevin

        if temperature is None:
            warn("temperature not set", level=2)

        if friction is None:
            warn("Friction not defined, set to 0.01", level=2)

        md = Langevin(
            atoms,
            temperature=temperature * u.kB,
            timestep=dt,
            friction=friction,
            logfile=logfile,
        )

    else:
        raise RuntimeError(f"Molecular dynamics mode {driver} is not suppported.")

    prepared = prepare_from_trajectory(atoms, md, trajectory_file)

    if md is None:
        raise RuntimeError("ASE MD algorithm has to be given")

    return atoms, md, prepared


def prepare_from_trajectory(atoms, md, trajectory_file=filenames.trajectory, **kwargs):
    """ Take the last step from trajectory and initialize atoms + md accordingly

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Reference structure for molecular dynamics
    md: ase.md.MolecularDynamics
        The MD propagator
    trajectory_file: str or Path
        The output trajectory file

    Returns
    -------
    bool
        True if prepared from the last step in a trajectory
    """

    trajectory_file = Path(trajectory_file).absolute()
    if trajectory_file.exists():
        last_atoms = last_from_yaml(trajectory_file)
        if "info" in last_atoms["atoms"]:
            md.nsteps = last_atoms["atoms"]["info"]["nsteps"]

            atoms.set_positions(last_atoms["atoms"]["positions"])
            atoms.set_velocities(last_atoms["atoms"]["velocities"])
            print(f"Resume MD from last step in\n  {trajectory_file}\n")
            return True

    print(f"** {trajectory_file} does not exist, nothing to prepare")
    return False


def initialize_md(
    atoms,
    temperature=None,
    force_constants=None,
    quantum=False,
    force_temp=True,
    deterministic=True,
    **kwargs,
):
    """Either use Maxwell Boltzmann or PhononHarmonics to prepare the MD run

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Reference structure for molecular dynamics
    temperature: float
        The temperature in Kelvin
    force_constants: str
        The filename of the file holding force constants for phonon rattle
    quantum: bool
        If True use Bose-Einstein distribution instead of Maxwell-Boltzmann
    force_temp: bool
        If True strictly enforce the initialization temperature
    deterministic: bool
        If True create sample deterministically

    Returns
    -------
    atoms: ase.atoms.Atoms
        Updated atoms with positions and velocities set by the initialization scheme
    """

    if temperature is None:
        return atoms

    print(f"Prepare MD run at temperature {temperature}")

    if force_constants is not None:
        print("Initialize positions and velocities using force constants.")
        force_constants = np.loadtxt(force_constants)
        PhononHarmonics(
            atoms,
            force_constants,
            quantum=quantum,
            temp=temperature * u.kB,
            plus_minus=deterministic,
            **kwargs,
        )
        if force_temp:
            temp0 = atoms.get_kinetic_energy() / len(atoms) / 1.5
            gamma = temperature * u.kB / temp0
            atoms.set_momenta(atoms.get_momenta() * np.sqrt(gamma))

    else:
        print("Initialize velocities according to Maxwell-Boltzmann distribution.")
        MaxwellBoltzmannDistribution(
            atoms, temp=temperature * u.kB, force_temp=force_temp
        )

    return atoms
