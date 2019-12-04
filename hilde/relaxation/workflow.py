from pathlib import Path

from ase.calculators.socketio import SocketIOCalculator

from hilde.settings import Settings
from hilde.helpers.watchdogs import SlurmWatchdog as Watchdog
from hilde.helpers.paths import cwd
from hilde.helpers.socketio import get_port
from hilde.trajectory import step2file, metadata2file
from hilde.helpers.structure import clean_atoms
from hilde.helpers import talk, warn
from hilde.helpers.restarts import restart
from . import metadata2dict
from ._defaults import name


_prefix = name
_calc_dirname = "calculation"
_temp_geometry_filename = "geometry.in.next_step"


def run_relaxation(ctx):
    """ high level function to run relaxation"""

    converged = run(ctx)

    if not converged:
        talk("restart", prefix=_prefix)
        restart(ctx.settings, trajectory=ctx.trajectory)
    else:
        talk("done.", prefix=_prefix)


def run(ctx, backup_folder="backups"):
    """ run a relaxation with ASE"""

    watchdog = Watchdog()

    # extract things from context
    atoms = ctx.atoms
    calculator = ctx.calc
    opt = ctx.opt
    fmax = ctx.fmax

    workdir = ctx.workdir
    trajectory = ctx.trajectory
    calc_dir = workdir / _calc_dirname

    socketio_port = get_port(calculator)
    if socketio_port is None:
        socket_calc = None
    else:
        socket_calc = calculator

    atoms.calc = calculator

    opt_atoms = ctx.opt_atoms
    opt.atoms = opt_atoms
    opt.initialize()

    # is a filter used?
    filter = len(atoms) < len(opt_atoms)

    with SocketIOCalculator(socket_calc, port=socketio_port) as iocalc, cwd(
        calc_dir, mkdir=True
    ):
        if socketio_port is not None:
            atoms.calc = iocalc

        # log very initial step and metadata
        if opt.nsteps == 0:
            metadata2file(ctx.metadata, trajectory)

        talk(f"Start step {opt.nsteps}", prefix=_prefix)
        for ii, converged in enumerate(opt.irun(fmax=fmax)):
            if converged:
                talk("Relaxation converged.", prefix=_prefix)
                break

            forces = opt_atoms.get_forces()

            # residual forces (and stress)
            na = len(atoms)
            res_forces = (forces[:na] ** 2).sum(axis=1).max() ** 0.5 * 1000
            if filter:
                res_stress = (forces[na:] ** 2).sum(axis=1).max() ** 0.5 * 1000

            # log if it's not the first step from a resumed relaxation
            if not (ii == 0 and opt.nsteps > 0):
                talk(f"Step {opt.nsteps} finished.", prefix=_prefix)
                talk(f".. residual force:  {res_forces:.3f} meV/AA", prefix=_prefix)
                if filter:
                    talk(f".. residual stress: {res_stress:.3f} meV/AA", prefix=_prefix)

                talk("clean atoms before logging", prefix=_prefix)
                log_atoms = clean_atoms(atoms, decimals=ctx.decimals)
                log_atoms.info.update({"nsteps": opt.nsteps})

                talk(f".. log", prefix=_prefix)
                step2file(log_atoms, atoms.calc, trajectory)

                info_str = [
                    f"Relaxed with BFGS, fmax={fmax*1000:.3f} meV/AA",
                    f"nsteps = {opt.nsteps}",
                    f"residual force  = {res_forces:.6f} meV/AA",
                ]
                if filter:
                    info_str.append(f"residual stress = {res_stress:.6f} meV/AA")

                log_atoms.write(
                    workdir / _temp_geometry_filename,
                    format="aims",
                    scaled=False,
                    info_str=info_str,
                )

            if watchdog():
                break

    return converged
