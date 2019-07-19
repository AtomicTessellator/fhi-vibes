""" run molecular dynamics simulations using the ASE classes """

import numpy as np
from ase.calculators.socketio import SocketIOCalculator
from hilde import son
from hilde.trajectory import step2file, metadata2file
from hilde.helpers.aims import get_aims_uuid_dict
from hilde.helpers.watchdogs import SlurmWatchdog as Watchdog
from hilde.helpers.paths import cwd
from hilde.helpers.socketio import get_port, get_stresses
from hilde.helpers.socketio import socket_stress_on, socket_stress_off
from hilde.helpers.compression import backup_folder as backup
from hilde.helpers.restarts import restart
from hilde.helpers import talk, warn
from . import metadata2dict


_calc_dirname = "calculations"
# _socket_timeout = 60
# _calculation_timeout = 30
_prefix = "md"


def run_md(ctx, timeout=None):
    """ high level function to run MD """

    converged = run(ctx)

    if not converged:
        talk("restart", prefix=_prefix)
        restart(ctx.settings, trajectory=ctx.trajectory)
    else:
        talk("done.", prefix=_prefix)


def run(ctx, backup_folder="backups"):
    """run and MD for a specific time

    Args:
        ctx (MDContext): context of the MD
        backup_folder (str or Path): Path to the back up folders
    Returns:
        bool: True if hit max steps or completed
    """

    # extract things from context
    atoms = ctx.atoms
    calc = ctx.calc
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
    trajectory = ctx.trajectory
    calc_dir = workdir / _calc_dirname
    backup_folder = workdir / backup_folder

    # prepare the socketio stuff
    socketio_port = get_port(calc)
    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc
    atoms.calc = calc

    # does it make sense to start everything?
    if md.nsteps >= maxsteps:
        msg = f"run already finished, please inspect {workdir.absolute()}"
        talk(msg, prefix=_prefix)
        return True

    # is the calculation similar enough?
    metadata = metadata2dict(atoms, calc, md)
    if trajectory.exists():
        old_metadata, _ = son.load(trajectory)
        check_metadata(metadata, old_metadata)

    # backup previously computed data
    backup(calc_dir, target_folder=backup_folder)

    # back up settings
    if settings:
        with cwd(workdir, mkdir=True):
            settings.obj["workdir"] = "."
            settings.write()

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        calc_dir, mkdir=True
    ):

        # make sure the socket is entered
        if socket_calc is not None:
            atoms.calc = iocalc

        # log very initial step and metadata
        if md.nsteps == 0:
            # carefully compute forces
            if not get_forces(atoms):
                return False
            # log metadata
            metadata2file(metadata, file=trajectory)
            # log initial structure computation
            atoms.info.update({"nsteps": md.nsteps, "dt": md.dt})
            meta = get_aims_uuid_dict()
            if compute_stresses:
                stresses = get_stresses(atoms)
                atoms.calc.results["stresses"] = stresses

            step2file(atoms, file=trajectory, metadata=meta, append_cell=False)

        while not watchdog() and md.nsteps < maxsteps:

            if not md_step(md):
                break

            talk(f"Step {md.nsteps} finished, log.", prefix=_prefix)

            if compute_stresses_now(compute_stresses, md.nsteps):
                stresses = get_stresses(atoms)
                atoms.calc.results["stresses"] = stresses

            # peek into aims file and grep for uuid
            atoms.info.update({"nsteps": md.nsteps, "dt": md.dt})
            meta = get_aims_uuid_dict()
            step2file(atoms, atoms.calc, trajectory, metadata=meta, append_cell=False)

            if compute_stresses:
                if compute_stresses_next(compute_stresses, md.nsteps):
                    talk("switch stresses computation on", prefix=_prefix)
                    socket_stress_on(iocalc)
                else:
                    talk("switch stresses computation off", prefix=_prefix)
                    socket_stress_off(iocalc)

                continue

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


# devel stuff from trying out pebble
# def launch_server(atoms, timeout=_socket_timeout):
#     """wraps `iocalc.launch_server` with pebble to kill the job if it takes too long
#
#     Args:
#         atoms: Atoms object with `iocalc` attached
#         timeout: the timeout in seconds
#     """
#     talk(f"Enter socket with timeout: {timeout}", prefix=_prefix)
#     iocalc = atoms.calc
#
#     @concurrent.process(timeout=timeout)
#     def launch(iocalc):
#         cmd = iocalc.calc.command.replace("PREFIX", iocalc.calc.prefix)
#         iocalc.calc.write_input(atoms)
#         iocalc.launch_server(cmd)
#         time.sleep(2)
#
#     return launch(iocalc).result()
#
#
# def get_forces(atoms, timeout=_calculation_timeout):
#     @concurrent.process(timeout=timeout)
#     def forces(atoms):
#         return atoms.get_forces()
#
#     try:
#         forces(atoms).result()
#         return True
#     except TimeoutError:
#         warn("TimeoutError raised while computing force", level=1)
#         warn("return False", level=1)
#         return False
#
#
# def md_step(md, timeout=_calculation_timeout):
#     @concurrent.process(timeout=timeout)
#     def step(md):
#         if md.nsteps > 2:
#             print("h")
#             time.sleep(20)
#         return md.run(1)
#
#     try:
#         step(md).result()
#         return True
#     except TimeoutError:
#         warn("TimeoutError raised while computing md step", level=1)
#         return False
#
