""" run molecular dynamics simulations using the ASE classe """

from pathlib import Path
from warnings import warn
from ase.calculators.socketio import SocketIOCalculator
from hilde.trajectory import step2file, metadata2file, from_yaml
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.paths import cwd


def prepare_from_trajectory(atoms, md, trajectory="trajectory.yaml"):
    """ Take the last step from trajectory and initialize atoms + md accordingly """

    trajectory = Path(trajectory).absolute()
    if trajectory.exists():
        traj = from_yaml(trajectory)
        last_atoms = traj[-1]
        md.nsteps = last_atoms["MD"]["nsteps"]

        atoms.set_positions(last_atoms["atoms"]["positions"])
        atoms.set_velocities(last_atoms["atoms"]["velocities"])
    else:
        warn(f"{trajectory} does  not exist, nothing to prepare")


def run_md(
    atoms,
    calc,
    md,
    maxsteps=25000,
    trajectory="trajectory.yaml",
    metadata="md_metadata.yaml",
    socketio_port=None,
    walltime=1800,
    watchdog=None,
    workdir=".",
):
    """ run and MD for a specific time

    Args:
        atoms (Atoms): initial structure
        calc (Calculator): calculator to calculate forces (and stresses)
        md (MolecularDynamics): MD algorithm
        trajectory (str/Path): use to save trajectory as yaml
        socketio_port (int): use to run through socketio
        walltime (int, optional): Defaults to 1800.
        workdir (str, optional): Defaults to '.'.
    """

    if watchdog is None:
        watchdog = Watchdog(walltime=walltime)

    if metadata is not None:
        metadata2file(atoms, calc, md, file=metadata)

    trajectory = Path(trajectory).absolute()

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc
    atoms.calc = calc

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        workdir, mkdir=True
    ):

        if socketio_port is not None:
            atoms.calc = iocalc

        # log very initial step
        if md.nsteps == 0:
            step2file(atoms, atoms.calc, md, trajectory)

        # store MD metadata locally
        metadata2file(atoms, calc, md, file=metadata)

        while watchdog() == False and md.nsteps < maxsteps:
            md.run(1)
            step2file(atoms, atoms.calc, md, trajectory)

    return md.nsteps >= maxsteps
