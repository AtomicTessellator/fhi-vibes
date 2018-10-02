from hilde.parsers import read_aims_output
from hilde.helpers.paths import cwd
from ase.calculators.socketio import SocketIOCalculator
from ase.io import Trajectory
from contextlib import ExitStack

def compute_forces(cells, calculator, workdir):
    """
    Compute forces in given list of atoms objects
    Args:
        cells: list of atoms objects
        calculator: ase.calculator for calculating the forces
        workdir: the working directory to compute forces in

    Returns:
        list of forces for each atoms object in cells
    """
    force_sets = []
    workdir.mkdir(exist_ok=True)
    for ii, cell in enumerate(cells):
        folder_with_disp = workdir / f'disp-{ii:03d}'
        folder_with_disp.mkdir(parents=True, exist_ok=True)
        try:
            force = read_aims_output(folder_with_disp / 'aims.out')[0].get_forces()
        except (FileNotFoundError, IndexError):
            cell.calc = calculator
            with cwd(folder_with_disp):
                try:
                    cell.calc.calculate(cell)
                except Error as inst:
                    print(inst)
            force = read_aims_output(folder_with_disp / 'aims.out')[0].get_forces()
        force_sets.append(force)
    return force_sets

def compute_forces_socketio(cells, calculator, port, workdir, traj_file, log_file):
    """
    Compute forces in given list of atoms objects utilizing the SocketIOCalculator
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
    force_sets = []
    tmp_dir.mkdir(exist_ok=True)
    with ExitStack() as stack, cwd(workdir):
        calc = stack.enter_context(SocketIOCalculator(calculator, log=log_file.open('w'), port=port))
        traj = stack.enter_context(Trajectory(str(traj_file), mode='a'))
        atoms = cells[0].copy()
        atoms.calc = calc
        for ii, cell in enumerate(cells):
            try:
                atoms.positions = cell.positions
                atoms.calc.calculate(atoms, system_changes=['positions'])
                force_sets.append(atoms.get_forces())
                traj.write(atoms)
            except Error as inst:
                print(inst)
    return force_sets
