""" run molecular dynamics simulations using the ASE classe """

from pathlib import Path
from ase.calculators.socketio import SocketIOCalculator
from hilde.trajectory import step2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog
from hilde.helpers.paths import cwd

def run_md(
    atoms,
    calc,
    md,
    trajectory="trajectory.yaml",
    socketio_port=None,
    walltime=1800,
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

    if Path(trajectory).exists():
        print('* restart not yet supported')


    watchdog = Watchdog(walltime=walltime)

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    with SocketIOCalculator(socket_calc, port=socketio_port) as calc, cwd(workdir):
        atoms.calc = calc

        while watchdog() == False:
            step2file(atoms, calc, trajectory)
            md.run(1)
