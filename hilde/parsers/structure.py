from hilde.structure import Cell
from ase.io import read as ase_read

# Parse geometry.in file
def read_aims(fname, symprec=symprec, sorted = False):
    from .structure import Cell
    latvecs = []
    positions = []
    scaled_positions = []
    symbols = []
    pbc = False
    constraints_pos = []
    constrains_lv   = []

    with open(fname,'r') as f:
        for line in f:
            if line.strip().startswith('#'):
                pass
            if line.strip().startswith('lattice_vector'):
                latvecs.append([float(el) for el in line.strip().split()[1:4]])
                pbc = True
            if line.strip().startswith('atom '):
                positions.append([float(el) for el in line.strip().split()[1:4]])
                symbols.append(line.strip().split()[4])
            if line.strip().startswith('atom_frac'):
                scaled_positions.append([float(el) for el in line.strip().split()[1:4]])
                symbols.append(line.strip().split()[4])

    kwargs = {
        'symbols': symbols,
        'cell': latvecs,
        'pbc': pbc
    }

    if positions:
        kwargs['positions'] = positions
    elif scaled_positions:
        kwargs['scaled_positions'] = scaled_positions
    else:
        exit(f'** Please specify atomic positions in {fname}.')

    #
    # Create cell object from this
    cell = Cell(symprec=symprec, **kwargs)

    if sorted:
        cell.sort_positions()

    return cell

def read_aims_output(fname):
    return ase_read(fname, ':', 'aims-output')
