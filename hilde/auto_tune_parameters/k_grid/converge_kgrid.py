"""Function that performs a k-grid optimization for a structure"""
from pathlib import Path
import numpy as np

from ase.constraints import UnitCellFilter
from ase.calculators.socketio import SocketIOCalculator

from hilde.auto_tune_parameters.k_grid.kpointoptimizer import KPointOptimizer
from hilde.helpers.paths import cwd
from hilde.trajectory.relaxation import metadata2file, step2file
from hilde.watchdogs import WallTimeWatchdog as Watchdog


def converge_kgrid(
    atoms,
    calc,
    func=lambda x: x.calc.get_property("energy", x) / len(x),
    loss_func=lambda x: x,
    dfunc_min=1e-4,
    even=True,
    unit_cell=True,
    maxsteps=100,
    trajectory="kpt_trajectory.yaml",
    logfile="kpoint_conv.log",
    socketio_port=None,
    walltime=1800,
    workdir=".",
    kpts_density_init=1.0,
):
    """
    Converges the k-grid relative to some loss function
    Args:
        atoms: (ASE Atoms object) geometry of the system you are converging the k-grid on
        calc: (ASE Calculator object) calculator for the k-grid convergence
        func: (function) Function used to get the property the routine is trying to converge relative to the k-grid density
        loss_func: (function) Function used to transform the property obtained in func into a score to compare agsint
        dfunc_min: (float) Convergence criteria for the loss function
        even: (bool) If True kgrid must be even valued
        unit_cell: (bool) if True system is periodic
        maxsteps: (int) maximum steps to run the optimization over
        trajecotry: (str) file name to store the trajectory
        logfile: (str) file name for the log file
        socketio_port: (int) port number for interactions with the socket
        walltime: (int) length of the wall time for the job in seconds
        workdir: (str) working directory for the calculation
        kpts_density_init: (float) initial k-point density
    Returns: (bool) True if the convergence criteria is met
    """
    watchdog = Watchdog(walltime=walltime)

    workdir = Path(workdir).absolute()
    trajectory = Path(trajectory).absolute()

    kpt_settings = {
        "func": func,
        "loss_func": loss_func,
        "dfunc_min": dfunc_min,
        "even": even,
        "logfile": str(workdir / logfile),
    }

    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calc

    atoms.calc = calc
    opt_atoms = atoms

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        workdir / "calculation", mkdir=True
    ):
        if socketio_port is not None:
            atoms.calc = iocalc

        opt = KPointOptimizer(opt_atoms, **kpt_settings)

        # log very initial step and metadata
        if opt.nsteps == 0:
            metadata2file(atoms, calc, opt, trajectory)

        for ii, converged in enumerate(opt.irun(steps=maxsteps)):
            step2file(atoms, atoms.calc, opt, trajectory, unit_cell=unit_cell)
            if watchdog():
                break

    return converged, opt.kpts_density, calc
