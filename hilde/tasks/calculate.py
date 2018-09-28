from hilde.parsers import read_aims_output
from hilde.helpers.paths import cwd

def compute_forces(scs, calculator, workdir):
    force_sets = []
    workdir.mkdir(exist_ok=True)
    for ii, cell in enumerate(scs):
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


