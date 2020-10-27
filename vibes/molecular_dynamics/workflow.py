""" run molecular dynamics simulations using the ASE classes """

import numpy as np
from ase.calculators.socketio import SocketIOCalculator

from vibes import keys, son
from vibes.helpers import talk, warn
from vibes.helpers.aims import get_aims_uuid_dict
from vibes.helpers.backup import backup_folder as backup
from vibes.helpers.backup import default_backup_folder
from vibes.helpers.paths import cwd
from vibes.helpers.restarts import restart
from vibes.helpers.socketio import (
    get_socket_info,
    get_stresses,
    socket_stress_off,
    socket_stress_on,
)
from vibes.helpers.utils import Timeout
from vibes.helpers.watchdogs import SlurmWatchdog as Watchdog
from vibes.trajectory import metadata2file, step2file

from ._defaults import calculation_timeout, name


_calc_dirname = "calculations"
# _socket_timeout = 60
_prefix = name


def run_md(ctx, timeout=None):
    """ high level function to run MD """

    converged = run(ctx)

    if not converged:
        talk("restart", prefix=_prefix)
        restart(ctx.settings, trajectory_file=ctx.trajectory_file)
    else:
        talk("done.", prefix=_prefix)


def run(ctx, backup_folder=default_backup_folder):
    """run and MD for a specific time

    Args:
        ctx (MDContext): context of the MD
        backup_folder (str or Path): Path to the back up folders
    Returns:
        bool: True if hit max steps or completed
    """

    # extract things from context
    atoms = ctx.atoms
    calculator = ctx.calculator
    md = ctx.md
    maxsteps = ctx.maxsteps
    compute_stresses = ctx.compute_stresses
    settings = ctx.settings

    # create watchdog
    buffer = 3
    if compute_stresses > 0:
        buffer = 5
    watchdog = Watchdog(buffer=buffer)

    # create working directories
    workdir = ctx.workdir
    trajectory_file = ctx.trajectory_file
    calc_dir = workdir / _calc_dirname
    backup_folder = workdir / backup_folder

    # prepare the socketio stuff
    socketio_port, socketio_unixsocket = get_socket_info(calculator)
    if socketio_port is None:
        socket_calc = None
        # choose some 5 digit number
        socketio_port = np.random.randint(0, 65000)
    else:
        socket_calc = calculator
    atoms.calc = calculator

    # does it make sense to start everything?
    if md.nsteps >= maxsteps:
        msg = f"run already finished, please inspect {workdir.absolute()}"
        talk(msg, prefix=_prefix)
        return True

    # is the calculation similar enough?
    metadata = ctx.metadata
    if trajectory_file.exists():
        old_metadata, _ = son.load(trajectory_file)
        check_metadata(metadata, old_metadata)

    # backup previously computed data
    backup(calc_dir, target_folder=backup_folder)

    # back up settings
    if settings:
        with cwd(workdir, mkdir=True):
            settings.write()

    timeout = Timeout(calculation_timeout)
    with SocketIOCalculator(
        socket_calc, port=socketio_port, unixsocket=socketio_unixsocket
    ) as iocalc, cwd(calc_dir, mkdir=True):
        # make sure the socket is entered
        if socket_calc is not None:
            atoms.calc = iocalc

        # log very initial step and metadata
        if md.nsteps == 0:
            # carefully compute forces
            if not get_forces(atoms):
                return False
            # log metadata
            metadata2file(metadata, file=trajectory_file)
            # log initial structure computation
            atoms.info.update({keys.nsteps: md.nsteps, keys.dt: md.dt})
            meta = get_aims_uuid_dict()
            if compute_stresses:
                stresses = get_stresses(atoms)
                atoms.calc.results["stresses"] = stresses

            step2file(atoms, file=trajectory_file, metadata=meta)

        while not watchdog() and md.nsteps < maxsteps:

            # reset timeout
            timeout()

            if not md_step(md):
                break

            talk(f"Step {md.nsteps} finished, log.", prefix=_prefix)

            if compute_stresses_now(compute_stresses, md.nsteps):
                stresses = get_stresses(atoms)
                atoms.calc.results["stresses"] = stresses
            else:  # make sure `stresses` are not logged
                if "stresses" in atoms.calc.results:
                    del atoms.calc.results["stresses"]

            # peek into aims file and grep for uuid
            atoms.info.update({keys.nsteps: md.nsteps, keys.dt: md.dt})
            meta = get_aims_uuid_dict()
            step2file(atoms, atoms.calc, trajectory_file, metadata=meta)

            if compute_stresses:
                if compute_stresses_next(compute_stresses, md.nsteps):
                    talk("switch stresses computation on", prefix=_prefix)
                    socket_stress_on(atoms.calc)
                else:
                    talk("switch stresses computation off", prefix=_prefix)
                    socket_stress_off(atoms.calc)

        talk("Stop.\n", prefix=_prefix)

    # restart
    if md.nsteps < maxsteps:
        return False
    return True


def compute_stresses_now(compute_stresses, nsteps):
    """Return if stress should be computed in this step"""
    return compute_stresses and (nsteps % compute_stresses == 0)


def compute_stresses_next(compute_stresses, nsteps):
    """Return if stress should be computed in the NEXT step"""
    return compute_stresses_now(compute_stresses, nsteps + 1)


def check_metadata(new_metadata, old_metadata):
    """Sanity check if metadata sets coincide"""
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
        if key not in om and key not in nm:
            continue
        ov, nv = om[key], nm[key]
        if isinstance(ov, float):
            assert np.allclose(ov, nv, rtol=1e-10), f"{key} changed from {ov} to {nv}"
        else:
            assert ov == nv, f"{key} changed from {ov} to {nv}"


def get_forces(atoms):
    try:
        _ = atoms.get_forces()
        return True
    except OSError as error:
        warn(f"Error during force computation:")
        print(error, flush=True)
        return False


def md_step(md):
    try:
        md.run(1)
        return True
    except OSError as error:
        warn(f"Error during MD step:")
        print(error, flush=True)
        return False
