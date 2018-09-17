from playground.konstanten.io import n_geom_digits
from playground.konstanten.symmetry import symprec
from scipy.linalg import norm
import datetime

def get_aims_string(cell, decorated=True, scaled=True, velocities=False):
    """ print the string that is geometry.in """
    if decorated:
        string  = '#=====================================================\n'
        string += f'# libflo:  geometry.in \n'
        #string += '#   Material: {:s}\n'.format(cell.get_chemical_formula())
        string += f"#   Date:    {datetime.datetime.now().isoformat(' ', timespec='seconds')}\n"
        string += '#=====================================================\n'
        string += f'#   Material:          {cell.get_chemical_formula()}\n'
        string += f'#   No. atoms:         {cell.n_atoms}\n'
        if cell.spacegroup is not None:
            string += '#   Spacegroup:        {:d}\n'.format(cell.spacegroup.number)
            string += '#   Wyckoff positions: {:s}\n'.format(' '.join([ss for ss in cell.spacegroup.wyckoffs]))
        if any(cell.pbc):
            string += '#   Unit cell volume:  {:f} AA^3\n'.format(cell.get_volume())
        if hasattr(cell, 'tags'):
            for ii, tag in enumerate(cell.tags):
                string += f'#   Tag {ii+1:2d}:            {tag}\n'
        # Supercell
        if hasattr(cell, 'smatrix'):
            string += f'#   Supercell matrix:  {cell.smatrix.flatten()}\n'

        # string += '\n'
    else:
        string  = ''
    # Write lattice
    # if decorated and cell.pbc:
    #     string += '# Lattice:\n'
    #
    # Order lattice by lengths (doesn't do anything right now)
    if any(cell.pbc):
        lengths = [norm(lv) for lv in cell.get_cell()]
    else:
        lengths = [norm(pos) for pos in cell.get_positions()]
    lv_args = range(len(lengths)) #np.argsort(lengths)

    # for latvec, constraint in zip(latvecs[lv_args], cell.constraints_lv[lv_args]):
    if any(cell.pbc):
        latvecs = cell.get_cell()
        for latvec, constraint in zip(latvecs, cell.constraints_lv):

            if decorated:
                string += f'  lattice_vector '
            else:
                string += f'lattice_vector '
            #
            string += f'{latvec[0]: .{n_geom_digits}e} '
            string += f'{latvec[1]: .{n_geom_digits}e} '
            string += f'{latvec[2]: .{n_geom_digits}e}\n'
            if constraint:
                string += 'constrain_relaxation .true.\n'

    # if not pbc: direct positions!
    else:
        scaled = False
    #
    # Write (preferably) scaled positions
    symbols = cell.get_chemical_symbols()
    if scaled:
        # if decorated:
        #     string += '\n# Scaled positions:\n'
        #
        positions = cell.get_scaled_positions()
        atompos   = 'atom_frac'
    else:
        # if decorated:
        #     string += '\n# Cartesian positions:\n'
        #
        positions = cell.get_positions()
        atompos   = 'atom'
    #
    if velocities:
        vels = cell.get_velocities()
    #
    for ii, (pos, sym) in enumerate(zip(positions, symbols)):
        if decorated:
            string += f'  {atompos:9s}  '
        else:
            string += f'{atompos:s}  '
        #
        string += f'{pos[lv_args[0]]: .{n_geom_digits}e} '
        string += f'{pos[lv_args[1]]: .{n_geom_digits}e} '
        string += f'{pos[lv_args[2]]: .{n_geom_digits}e}  '
        string += f'{sym:3s}\n'
        if velocities:
            vel = vels[ii]
            if decorated:
                string += '    velocity '
            else:
                string += 'velocity '
            string += f'{vel[lv_args[0]]: .{n_geom_digits}e} '
            string += f'{vel[lv_args[1]]: .{n_geom_digits}e} '
            string += f'{vel[lv_args[2]]: .{n_geom_digits}e}\n'
    #
    return string

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

def inform(cell, dft=False, fname=None, verbosity=0):
    # Structure info:
    print(f'\nGeometry info for:')
    print(f'  input geometry:    {cell.sysname}')
    if fname: print(f'  from:              {fname}')
    print(f'  Symmetry prec.:    {cell.symprec}')
    print(f'  Number of atoms:   {cell.n_atoms}')
    print(f"  Species:           {', '.join(cell.get_unique_symbols()[0])}")
    print(f'  Periodicity:       {cell.pbc}')
    print(f'')

    if cell.spacegroup:
        sds = cell.spacegroup
        print(f'  Spacegroup:          {sds.international} ({sds.number})')
        print(f'  Wyckoff positions:   {sds.wyckoffs}')
        print(f'  Equivalent atoms:    {sds.equivalent_atoms}')
        print(f'  Standard lattice:  ')
        for vec in sds.spglib_std_lattice:
            print(f'    {vec}')

    # lengths and angles
    la = cell.get_cell_lengths_and_angles()
    print('\nCell lengths and angles [\u212B, °]:')
    print('  a, b, c: {}'.format(' '.join([f'{l:15.8f}' for l in la[:3]])))
    print('  \u03B1, \u03B2, \u03B3: {}'.format(' '.join([f'{l:15.8f}' for l in la[3:]])))
    print(f'  Volume:  {cell.get_volume():15.8f} \u212B**3')

    for (ii, tag) in enumerate(tag for tag in cell.tags if tag):
        print(f'  Tag {ii+1:2}:                {tag}')

    if hasattr(cell, 'dft_calculation') and dft:
        print(f'  DFT calc performed:    yes')
        if verbosity == 1:
            dft_dict = cell.dft_calculation.__dict__
            print(f'Full Summary of DFT Results:')
            for key in dft_dict.keys():
                print(f'  {key:25s}: {dft_dict[key]}')
        #
        print(f'Compact summary of DFT Results:')
        print(f'  Forces computed:       {len(cell.get_forces()) > 0}')
        print(f'  Stress computed:       {cell.get_stress() is not None}')
        if cell.get_stress() is not None:
            print(f'  Pressure:              {cell.get_pressure():12.4e} eV')
        print(f'  Total energy:          {cell.get_total_energy():12.4e} eV')
        if cell.get_chemical_potential():
            print(f'  Fermi level:           {cell.get_chemical_potential():12.4e} eV')
        else:
            print(f'  Fermi level:           None')
        print(f'  HOMO-LUMO gap:         {cell.get_homo_lumo_gap():12.4e} eV')