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

def cells_and_workdirs(cells, base_dir):
    """ generate tuples of atoms object and workingdirectory path """
    for ii, cell in enumerate(cells):
        workdir = Path(base_dir) / f'{ii:05d}'
        yield cell, workdir


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
    warn('dont compare stupid things -> replace me with hashes')
    try:
        calcs_coincide = all(p1 == p2 for p1, p2 in
                             zip_longest(calc1.parameters,
                                         calc2.parameters))
    except AttributeError:
        calcs_coincide = False
    return calcs_coincide


def return_from_trajectory(cells, calculator, trajectory,
                           ordered=True):
    """ return pre computed atoms objects from trajectory """
    cells_computed = []
    if ordered:
        for cell in cells:
            for atoms in trajectory:
                if ((atoms == cell) and calculators_coincide(
                    calculator, atoms.calc)):
                    cells_computed.append(atoms)
            else:
                cells_computed.append(cell)
    else:
        raise Exception('non ordered read from trajectory not yet implemented.')
    return cells_computed


def calculate_multiple(cells, calculator, workdir,
                       trajectory=None,
                       trajectory_to=None):
    """Calculate several atoms object sharing the same calculator.

    Args:
        cells (Atoms): List of atoms objects.
        calculator (calculator): Calculator to run calculation.
        workdir (str/Path): working directory
        trajectory (Trajectory): trajectory for caching
        trajectory_to (Trajectory): store information to trajectory
        read (bool): Read from trajectory if True.

    Returns:
        list: Atoms objects with calculation performed.

    """

    # If trajectory_from is given, try to read pre-computed atoms from it
    cells_pre_calculated = cells
    if trajectory is not None and trajectory_to is None:
        trajectory_to = trajectory

    if trajectory is not None:
        traj_file = os.path.join(workdir, trajectory)
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

    cells_calculated = []
    for cell, wdir in cells_and_workdirs(cells_pre_calculated, workdir):
        if cell is None:
            cells_calculated.append(cell)
            continue
        if cell.calc is None:
            cell = calculate(cell, calculator, wdir)
        cells_calculated.append(cell)

        if trajectory_to is not None:
            with Trajectory(traj_file, mode='a') as traj:
                traj.write(cell)

    return cells_calculated


def calculate_multiple_socketio(cells, calculator, workdir, port,
                                trajectory='socketio.traj',
                                trajectory_to=None,
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

    # If trajectory_from is given, try to read pre-computed atoms from it
    cells_pre_calculated = cells
    if trajectory is not None and trajectory_to is None:
        trajectory_to = trajectory

    if trajectory is not None:
        traj_file = os.path.join(workdir, trajectory)
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


    cells_calculated = []
    workdir.mkdir(exist_ok=True)
    with ExitStack() as stack, cwd(workdir, debug=True):
        calc = stack.enter_context(
            SocketIOCalculator(calculator, log=log_file.open('w'), port=port))
        atoms = cells[0].copy()
        atoms.calc = calc
        for _, cell in enumerate(cells_pre_calculated):

            if cell is None:
                cells_calculated.append(cell)
                continue

            if cell.calc is None:
                atoms.positions = cell.positions
                atoms.calc.calculate(atoms, system_changes=['positions'])
                computed_cell = atoms.copy()
                # restore the 'original' calculator
                computed_cell.calc = calculator
                computed_cell.calc.results = atoms.calc.results
            else:
                computed_cell = cell

            cells_calculated.append(computed_cell)

            if trajectory_to is not None:
                with Trajectory(traj_file, mode='a') as traj:
                    traj.write(computed_cell)

    return cells_calculated


def setup_multiple(cells, calculator, workdir):
    """
    Write input files for calculations on a list of atoms objects
    Args:
        cells (Atoms): List of atoms objects.
        calculator (calculator): Calculator to run calculation.
        workdir (str/Path): working directory
    """
    workdirs = []
    for cell, wdir in cells_and_workdirs(cells, workdir):
        workdirs.append(wdir)
        with cwd(wdir, mkdir=True):
            cell.set_calculator(calculator)
            try:
                calculator.write_input(cell)
            except AttributeError:
                print("Calculator has no input just attaching to cell")
    return cells, workdirs
