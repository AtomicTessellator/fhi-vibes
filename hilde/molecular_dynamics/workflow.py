""" run molecular dynamics simulations using the ASE classes """

from pathlib import Path

import numpy as np

from ase.calculators.socketio import SocketIOCalculator
from ase import md as ase_md
from ase import units as u

from hilde import son
from hilde.aims.context import AimsContext
from hilde.trajectory import step2file, metadata2file
from hilde.helpers.aims import get_aims_uuid_dict
from hilde.helpers.watchdogs import SlurmWatchdog as Watchdog
from hilde.helpers.paths import cwd
from hilde.helpers.socketio import get_port, get_stresses
from hilde.helpers.socketio import socket_stress_on, socket_stress_off
from hilde.helpers.compression import backup_folder as backup
from hilde.helpers.restarts import restart
from hilde.helpers import talk
from . import metadata2dict


_calc_dirname = "calculations"


def run_md(ctx):
    """ high level function to run MD """

    args = bootstrap(ctx)

    converged = run(**args)

    if not converged:
        restart(ctx.settings)
    else:
        talk("done.")


def bootstrap(ctx):
    """Load settings, prepare atoms, calculator and MD algorithm

    Parameters
    ----------
    ctx: MDContext
        Context for the workflow

    Returns
    -------
    dict
        The relevant information to run the MD with the following items

        atoms: ase.atoms.Atoms
            The reference structure
        calc: ase.calculators.calulator.Calculator
            The Calculator for the MD
        maxsteps: int
            Maximum number of steps for the MD
        compute_stresses: bool
            If True compute the stresses
        workdir: str
            working directory for the run
    """

    # read structure
    atoms = ctx.settings.get_atoms()

    # create aims from context
    aims_ctx = AimsContext(settings=ctx.settings)

    # make sure `compute_forces .true.` is set
    aims_ctx.settings.obj["compute_forces"] = True

    calc = aims_ctx.get_calculator()

    # create md from context
    obj = ctx.settings.obj
    # workdir has to exist
    Path(ctx.workdir).mkdir(exist_ok=True)
    md_settings = {
        "atoms": atoms,
        "timestep": obj.timestep * u.fs,
        "logfile": Path(ctx.workdir) / obj.logfile,
    }

    if "verlet" in obj.driver.lower():
        md = ase_md.VelocityVerlet(**md_settings)
    if "langevin" in obj.driver.lower():
        md = ase_md.Langevin(
            temperature=obj.temperature * u.kB, friction=obj.friction, **md_settings
        )

    compute_stresses = 0
    if "compute_stresses" in obj:
        compute_stresses = obj.compute_stresses

    # resume?
    prepare_from_trajectory(atoms, md, ctx.trajectory)

    return {
        "atoms": atoms,
        "calc": calc,
        "md": md,
        "maxsteps": ctx.maxsteps,
        "compute_stresses": compute_stresses,
        "workdir": ctx.workdir,
    }


def run(
    atoms,
    calc,
    md,
    maxsteps,
    compute_stresses=0,
    trajectory="trajectory.son",
    metadata_file="md_metadata.yaml",
    workdir=".",
    backup_folder="backups",
    **kwargs,
):
    """run and MD for a specific time

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Initial step geometry
    calc: ase.calculators.calulator.Calculator
        The calculator for the MD run
    md: ase.md.MolecularDynamics
        The MD propagator
    maxsteps: int
        Maximum number of steps
    compute_stresses: int or bool
        Number of steps in between each stress computation
            if False 0, if True 1, else int(compute_stresses)
    trajectory: Path or str
        trajectory file path
    metadata_file: str or Path
        File to store the metadata in
    workdir: Path or str
        Path to working directory
    backup_folder: str or Path
        Path to the back up folders

    Returns
    -------
    bool
        True if hit max steps or completed

    """

    # create watchdog
    watchdog = Watchdog()

    # create working directories
    workdir = Path(workdir)
    trajectory = (workdir / trajectory).absolute()
    calc_dir = workdir / _calc_dirname
    backup_folder = workdir / backup_folder

    # make sure compute_stresses describes a step length
    if compute_stresses is True:
        compute_stresses = 1
    elif compute_stresses is False:
        compute_stresses = 0
    else:
        compute_stresses = int(compute_stresses)

    # atomic stresses
    if calc.name == "aims":
        if compute_stresses:
            talk(f"Compute atomic stresses every {compute_stresses} steps")
            calc.parameters["compute_heat_flux"] = True

    # prepare the socketio stuff
    socketio_port = get_port(calc)
    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc
    atoms.calc = calc

    # does it make sense to start everything?
    if md.nsteps >= maxsteps:
        talk(f"MD run already finished, please inspect {workdir.absolute()}")
        return True

    # is the calculation similar enough?
    metadata = metadata2dict(atoms, calc, md)
    if trajectory.exists():
        old_metadata, _ = son.load(trajectory)
        check_metadata(metadata, old_metadata)

    # backup previously computed data
    backup(calc_dir, target_folder=backup_folder)

    talk("Start MD.")

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        calc_dir, mkdir=True
    ):

        if socketio_port is not None:
            atoms.calc = iocalc

        # log very initial step and metadata
        if md.nsteps == 0:
            metadata2file(metadata, file=trajectory)
            atoms.info.update({"nsteps": md.nsteps, "dt": md.dt})
            _ = atoms.get_forces()
            meta = get_aims_uuid_dict()
            step2file(atoms, atoms.calc, trajectory, metadata=meta)

        # store MD metadata locally
        metadata2file(metadata, file=metadata_file)

        while not watchdog() and md.nsteps < maxsteps:

            md.run(1)

            talk(f"Step {md.nsteps} finished, log.")

            if compute_stresses_now(compute_stresses, md.nsteps):
                stresses = get_stresses(atoms)
                atoms.calc.results["stresses"] = stresses

            # peek into aims file and grep for uuid
            atoms.info.update({"nsteps": md.nsteps, "dt": md.dt})
            meta = get_aims_uuid_dict()
            step2file(atoms, atoms.calc, trajectory, metadata=meta)

            if compute_stresses:
                if compute_stresses_next(compute_stresses, md.nsteps):
                    talk("switch stresses computation on")
                    socket_stress_on(iocalc)
                else:
                    talk("switch stresses computation off")
                    socket_stress_off(iocalc)

                continue

        talk("Stop MD.\n")

    # restart
    if md.nsteps < maxsteps:
        return False
    return True


def compute_stresses_now(compute_stresses, nsteps):
    """Return if stress should be computed in this step

    Parameters
    ----------
    compute_stresses: int
        Number of steps between each stress calculation
    nsteps: int
        Current step number

    Returns
    -------
    bool
        True if the stress should be computed at this step
    """
    return compute_stresses and (nsteps % compute_stresses == 0)


def compute_stresses_next(compute_stresses, nsteps):
    """Return if stress should be computed in the NEXT step

    Parameters
    ----------
    compute_stresses: int
        Number of steps between each stress calculation
    nsteps: int
        Current step number

    Returns
    -------
    bool
        True if the stress should be computed at the next step
    """
    return compute_stresses_now(compute_stresses, nsteps + 1)


def prepare_from_trajectory(atoms, md, trajectory):
    """Take the last step from trajectory and initialize atoms + md accordingly

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The initial geometry
    md: ase.md.MolecularDynamics
        The MD propagator
    trajectory: str or Path
        The trajectory file
    """

    trajectory = Path(trajectory)
    if trajectory.exists():
        last_atoms = son.last_from(trajectory)
        assert "info" in last_atoms["atoms"]
        md.nsteps = last_atoms["atoms"]["info"]["nsteps"]

        atoms.set_positions(last_atoms["atoms"]["positions"])
        atoms.set_velocities(last_atoms["atoms"]["velocities"])
        talk(f"Resume MD from step {md.nsteps} in\n  {trajectory}\n")
        return True

    talk(f"** {trajectory} does not exist, nothing to prepare")
    return False


def check_metadata(new_metadata, old_metadata):
    """Sanity check if metadata sets coincide

    Parameters
    ----------
    new_metadata: dict
        The metadata for this run
    old_metadata: dict
        The metadata for the run stored in the trajectory file

    Raises
    ------
    AssertionError
        If the metadata do not agree
    """
    om, nm = old_metadata["MD"], new_metadata["MD"]

    # check if keys coincide:
    # sanity check values:
    check_keys = ("md-type", "timestep", "temperature", "friction", "fs")
    keys = [k for k in check_keys if k in om.keys()]
    for key in keys:
        ov, nv = om[key], nm[key]
        if isinstance(ov, float):
            assert np.allclose(ov, nv, rtol=1e-10), f"{key} changed from {ov} to {nv}"
        else:
            assert ov == nv, f"{key} changed from {ov} to {nv}"

    # calculator
    om = old_metadata["calculator"]["calculator_parameters"]
    nm = new_metadata["calculator"]["calculator_parameters"]

    # sanity check values:
    for key in ("xc", "k_grid", "relativistic"):
        ov, nv = om[key], nm[key]
        if isinstance(ov, float):
            assert np.allclose(ov, nv, rtol=1e-10), f"{key} changed from {ov} to {nv}"
        else:
            assert ov == nv, f"{key} changed from {ov} to {nv}"
