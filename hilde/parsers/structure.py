from hilde.structure import pAtoms
from hilde.konstanten.symmetry import symprec
from ase.io import read as ase_read

# Parse geometry.in file
def read_structure(fname, symprec=symprec, format='aims'):
    return pAtoms(ase_read(fname, 0, format), symprec=symprec)


def read_aims(fname, symprec=symprec):
    print('** Please use hilde.parsers.read_structure instead of .read_aims')
    return read_structure(fname, symprec, format)

# def read_aims(fname, symprec=symprec, sorted = False):
#     from .structure import pAtoms
#     latvecs = []
#     positions = []
#     scaled_positions = []
#     symbols = []
#     pbc = False
#     constraints_pos = []
#     constrains_lv   = []
#
#     with open(fname,'r') as f:
#         for line in f:
#             if line.strip().startswith('#'):
#                 pass
#             if line.strip().startswith('lattice_vector'):
#                 latvecs.append([float(el) for el in line.strip().split()[1:4]])
#                 pbc = True
#             if line.strip().startswith('atom '):
#                 positions.append([float(el) for el in line.strip().split()[1:4]])
#                 symbols.append(line.strip().split()[4])
#             if line.strip().startswith('atom_frac'):
#                 scaled_positions.append([float(el) for el in line.strip().split()[1:4]])
#                 symbols.append(line.strip().split()[4])
#
#     kwargs = {
#         'symbols': symbols,
#         'cell': latvecs,
#         'pbc': pbc
#     }
#
#     if positions:
#         kwargs['positions'] = positions
#     elif scaled_positions:
#         kwargs['scaled_positions'] = scaled_positions
#     else:
#         exit(f'** Please specify atomic positions in {fname}.')
#
#     #
#     # Create cell object from this
#    cell = pAtoms(symprec=symprec, **kwargs)
#
#    if sorted:
#        cell.sort_positions()
#
#    return cell

def read_output(fname, format='aims-output'):
    """ Right now this is just wrapper for ase.io.read(file, ':', 'aims-output')"""
    return ase_read(fname, ':', format)


def read_aims_output(fname):
    """ Right now this is just wrapper for ase.io.read(file, ':', 'aims-output')"""
    print('** Please use hilde.parsers.read_output instead of read_aims_output')
    return ase_read(fname, ':', 'aims-output')


def read_lammps_output(fname):
    """ Right now this is just wrapper for ase.io.read(file, ':', 'aims-output')"""
    print('** Please use hilde.parsers.read_output instead of ' +
          'read_lammps_output')
    return ase_read(fname, ':', 'lammps')
