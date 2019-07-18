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
from hilde.helpers import talk, Timer
from . import metadata2dict


_calc_dirname = "calculations"


def run_md(ctx, timeout=None):
    """ high level function to run MD """

    converged = run(ctx)

    if not converged:
        restart(ctx.settings)
    else:
        talk("done.")


def run(ctx, backup_folder="backups", socket_timeout=60):
    """run and MD for a specific time

    Args:
        ctx (MDContext): context of the MD
        backup_folder (str or Path): Path to the back up folders
        socket_timeout (int): timeout for the socket communication
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
        talk(f"[md] run already finished, please inspect {workdir.absolute()}")
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

    msg = f"Enter socket with timeout: {socket_timeout}"
    socket_timer = Timer(msg, timeout=socket_timeout)

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        calc_dir, mkdir=True
    ):

        socket_timer("Socket entered")

        if socketio_port is not None:
            atoms.calc = iocalc

        # log very initial step and metadata
        if md.nsteps == 0:
            metadata2file(metadata, file=trajectory)
            atoms.info.update({"nsteps": md.nsteps, "dt": md.dt})
            _ = atoms.get_forces()
            meta = get_aims_uuid_dict()
            step2file(atoms, file=trajectory, metadata=meta, append_cell=False)

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
