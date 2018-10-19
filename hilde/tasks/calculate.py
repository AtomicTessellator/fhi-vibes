"""
    Functions to run several related calculations with using trajectories as
    cache
"""

import os
from contextlib import ExitStack
from pathlib import Path
from warnings import warn
from itertools import zip_longest
from ase.calculators.socketio import SocketIOCalculator
from ase.io import Trajectory

from hilde.helpers.paths import cwd
from hilde.structure import pAtoms

def calculate(atoms, calculator, workdir):
    """Short summary.
    Perform a dft calculation with ASE
    Args:
        atoms (Atoms or pAtoms): Structure.
        calculator (calculator): Calculator.
        workdir (folder): Folder to perform calculation in.

    Returns:
        int: 1 when error happened, otherwise None
    """
    calc_atoms = atoms.copy()
    calc_atoms.calc = calculator
    with cwd(workdir, mkdir=True):
        calc_atoms.write('geometry.in')
        calc_atoms.calc.calculate(calc_atoms)
        return calc_atoms


def calculators_coincide(calc1, calc2):
    """ Check if calculators coincide by checking the parameters """
    calcs_coincide = all(p1 == p2 for p1, p2 in
                         zip_longest(calc1.parameters,
                                     calc2.parameters))
    return calcs_coincide


def return_from_trajectory(cells, calculator, trajectory,
                           ordered=True):
    """ return pre computed atoms objects from trajectory """
    cells_computed = []
    if ordered:
        itrajectory = iter(trajectory)
        atoms = next(itrajectory)
        for cell in cells:
            if ((atoms == cell) and calculators_coincide(calculator,
                                                         atoms.calc)):
                cells_computed.append(atoms)
                atoms = next(itrajectory)
            else:
                cells_computed.append(cell)
    else:
        raise Exception('non ordered read from trajectory not yet implemented.')
    return cells_computed


def calculate_multiple(cells, calculator, workdir, trajectory_to=None,
                       trajectory_from=None):
    """Calculate several atoms object sharing the same calculator.

    Args:
        cells (Atoms): List of atoms objects.
        calculator (calculator): Calculator to run calculation.
        workdir (str/Path): working directory
        trajectory_to (Trajectory): store information to trajectory
        trajectory_from (Trajectory): read information from trajectory
        read (bool): Read from trajectory if True.

    Returns:
        list: Atoms objects with calculation performed.

    """

    # If trajectory_from is given, try to read pre-computed atoms from it
    cells_pre_calculated = cells
    if trajectory_from is not None:
        traj_file = os.path.join(workdir, trajectory_from)
        if Path(traj_file).exists():
            with Trajectory(traj_file, mode='r') as trajectory:
                cells_pre_calculated = return_from_trajectory(cells,
                                                              calculator,
                                                              trajectory)

    # Check if backup should be performed
    if trajectory_to is not None:
        traj_file = os.path.join(workdir, trajectory_to)
        if Path(traj_file).exists():
            Path(traj_file).rename(str(traj_file) + '.bak')

    workdirs = [Path(workdir) / f'{ii:05d}' for ii, _ in enumerate(cells)]

    cells_calculated = []
    for cell, wdir in zip(cells_pre_calculated, workdirs):
        if cell is None:
            cells_calculated.append(cell)
            continue
        if cell.calc is None:
            cell = calculate(cell, calculator, wdir)
        cells_calculated.append(cell)

        if trajectory_to is not None:
            with Trajectory(traj_file, mode='a') as trajectory:
                trajectory.write(cell)

    return cells_calculated


def calculate_multiple_socketio(cells, calculator, workdir, port,
                                trajectory='socketio.traj',
                                force=False,
                                log_file='socketio.log'):
    """
    Compute given list of atoms objects utilizing the SocketIOCalculator
    Args:
        cells: list of atoms objects
        calculator: ase.calculator for calculating the forces, needs to support socketio
        port: the port to be used for socket communication
        workdir: the workingdirectory for storing input and output files
        traj_file: file to store ase.Trajectory() object in. Basically a list of atom objects
        log_file: file to pipe status messages to

    Returns:
        list of forces for each atoms object in cells

    """

    if trajectory:
        raise Exception('Future Flo, implement me!')
        traj_file = Path(workdir) / trajectory
        traj, is_calculated = return_trajectory(cells, calculator,
                                                traj_file, force)
        if is_calculated and not force:
            return traj

    cells_calculated = []
    workdir.mkdir(exist_ok=True)
    with ExitStack() as stack, cwd(workdir, debug=True):
        calc = stack.enter_context(
            SocketIOCalculator(calculator, log=log_file.open('w'), port=port))
        atoms = cells[0].copy()
        atoms.calc = calc
        for _, cell in enumerate(cells):
            atoms.positions = cell.positions
            atoms.calc.calculate(atoms, system_changes=['positions'])
            cells_calculated.append(atoms)

            # prepare atoms to copy to Trajectory
            cell.calc = calculator
            cell.calc.results = atoms.calc.results
            traj.write(cell)
    return cells_calculated


def setup_multiple(cells, calculator, workdir):
    """
    Write input files for calculations on a list of atoms objects
    Args:
        cells (Atoms): List of atoms objects.
        calculator (calculator): Calculator to run calculation.
        workdir (str/Path): working directory
    """
    workdirs = [Path(workdir) / f'{ii:05d}' for ii, _ in enumerate(cells)]

    for cell, wdir in zip(cells, workdirs):
        with cwd(wdir, mkdir=True):
            cell.set_calculator(calculator)
            try:
                calculator.write_input(cell)
            except AttributeError:
                print("Calculator has no input just attaching to cell")
    return cells, workdirs


# def compute_forces(cells, calculator, workdir, trajectory=None):
#     """
#     Compute forces in given list of atoms objects
#     Args:
#         cells: list of atoms objects
#         calculator: ase.calculator for calculating the forces
#         workdir: the working directory to compute forces in
#
#     Returns:
#         list of forces for each atoms object in cells
#     """
#
#     if trajectory:
#         force_sets = return_trajectory(trajectory)
#
#     force_sets = []
#     for ii, cell in enumerate(cells):
#         folder_with_disp = Path(workdir) / f'disp-{ii:03d}'
#         folder_with_disp.mkdir(parents=True, exist_ok=True)
#         cell.write(folder_with_disp / 'geometry.in')
#         cell = calculate(cell, calculator, folder_with_disp)
#         force = cell.get_forces()
#         force_sets.append(force)
#     return force_sets
#
# def compute_forces_socketio(cells, calculator, port, workdir, traj_file,
#                             log_file):
#     """
#     Compute forces in given list of atoms objects utilizing the SocketIOCalculator
#     Args:
#         cells: list of atoms objects
#         calculator: ase.calculator for calculating the forces, needs to support socketio
#         port: the port to be used for socket communication
#         workdir: the workingdirectory for storing input and output files
#         traj_file: file to store ase.Trajectory() object in. Basically a list of atom objects
#         log_file: file to pipe status messages to
#
#     Returns:
#         list of forces for each atoms object in cells
#
#     """
#     force_sets = []
#     workdir.mkdir(exist_ok=True)
#     with ExitStack() as stack, cwd(workdir):
#         calc = stack.enter_context(
#             SocketIOCalculator(calculator, log=log_file.open('w'), port=port))
#         traj = stack.enter_context(Trajectory(str(traj_file), mode='a'))
#         atoms = cells[0].copy()
#         atoms.calc = calc
#         for _, cell in enumerate(cells):
#             atoms.positions = cell.positions
#             atoms.calc.calculate(atoms, system_changes=['positions'])
#             force_sets.append(atoms.get_forces())
#             traj.write(atoms)
#     return force_sets
