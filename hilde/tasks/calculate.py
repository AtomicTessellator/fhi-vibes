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


def return_trajectory(cells, calculator, trajectory, force=False):
    """Return either a trajectory of computed atoms objects or an empty one

    Returns:
        (traj, is_calculated): the trajectory and wether it was already
        calculated"""

    # Case 1: Trajectory does not yet exist
    if not Path(trajectory).exists():
        Path(trajectory).parent.mkdir(parents=True)
        traj = Trajectory(str(trajectory), mode='a')
        return traj, False

    # Case 2: Overwrite trajectory
    if force:
        Path(trajectory).rename(str(trajectory) + '.bak')
        traj = Trajectory(str(trajectory), mode='a')
        return traj, False

    # Case 3: Read trajectory and try to return if calculations seem plausible
    try:
        traj = [pAtoms(ase_atoms=atoms) for atoms in
                Trajectory(str(trajectory), mode='r')]

        atoms_coincide = all(cell == atoms for cell, atoms in zip(cells, traj))
        calcs_coincide = all((p1 == p2 for p1, p2 in
                              zip_longest(calculator.parameters,
                                          atoms.calc.parameters))
                             for atoms in traj)
        if atoms_coincide and calcs_coincide:
            return traj, True
        # else
        message = (f'atoms objects and calculators in trajectory ' +
                   f'{trajectory} do not coincide. Compute explicitly.')
        warn(message=message)

    except Exception as inst:
        print(inst)
        exit(f'probably the Trajectory in {trajectory} is empty')

    Path(trajectory).rename(str(trajectory) + '.bak')
    traj = Trajectory(str(trajectory), mode='a')
    return traj, False


def calculate_multiple(cells,
                       calculator,
                       workdir,
                       trajectory=None,
                       workdirs=None,
                       force=False):
    """Calculate several atoms object sharing the same calculator.

    Args:
        cells (Atoms): List of atoms objects.
        calculator (calculator): Calculator to run calculation.
        workdir (str/Path): working directory
        workdirs: list of str/paths: specify the directory to calculate each cell in
        trajectory (Trajectory): store (and read) information from trajectory
        read (bool): Read from trajectory if True.

    Returns:
        list: Atoms objects with calculation performed.

    """

    if trajectory:
        traj_file = Path(workdir) / trajectory
        traj, is_calculated = return_trajectory(cells, calculator,
                                                traj_file, force)
        if is_calculated and not force:
            return traj
    if workdirs is None:
        workdirs = [Path(workdir) / f'{ii:05d}' for ii, _ in enumerate(cells)]
    else:
        assert len(cells) == len(workdirs)

    cells_calculated = []
    for cell, wdir in zip(cells, workdirs):
        cell = calculate(cell, calculator, wdir)
        cells_calculated.append(cell)
        if trajectory:
            traj.write(cell)

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
            except:
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
