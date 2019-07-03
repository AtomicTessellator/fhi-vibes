"""read, format, and write structures and inform about them"""

import datetime
import numpy as np
from scipy.linalg import norm
from hilde.konstanten import v_unit
from hilde.konstanten.io import n_geom_digits
from hilde.konstanten.symmetry import symprec
from hilde.helpers.numerics import clean_matrix
from hilde.helpers import talk
from hilde.helpers.brillouinzone import get_special_points
from hilde.structure.misc import get_sysname
from hilde.spglib.wrapper import get_symmetry_dataset


def get_aims_string(cell, decorated=True, scaled=None, velocities=False, wrap=True):
    """print the string that is geometry.in

    Parameters
    ----------
    cell: ase.atoms.Atoms
        The cell to convert to geometry.in
    decorated: bool
        If True add header to the geoemtry.in string
    scaled: bool
        If True use scaled positions
    velocities: bool
        If True include velocities
    wrap: bool
        If True wrap the scaled positions

    Returns
    -------
    string: str
        string representation of geometry.in
    """
    if scaled is None:
        if "supercell" in cell.tags:
            scaled = False
        else:
            scaled = True

    if decorated:
        sds = get_symmetry_dataset(cell, symprec=symprec)
        string = "#=====================================================\n"
        string += f"# libflo:  geometry.in \n"
        # string += '#   Material: {:s}\n'.format(cell.get_chemical_formula())
        string += f"#   Date:    {datetime.datetime.now().isoformat(' ', timespec='seconds')}\n"
        string += "#=====================================================\n"
        string += f"#   Material:          {cell.get_chemical_formula()}\n"
        string += f"#   No. atoms:         {cell.n_atoms}\n"
        string += f"#   Spacegroup:        {sds.number:d}\n"
        # string += (f'#   Wyckoff positions: ' +
        #             ', '.join(f'{c}*{w}' for (w, c) in sds.wyckoffs_unique) +
        #             '\n')
        if any(cell.pbc):
            string += f"#   Unit cell volume:  {cell.get_volume():f} AA^3\n"
        if hasattr(cell, "tags"):
            for ii, tag in enumerate(cell.tags):
                string += f"#   Tag {ii+1:2d}:            {tag}\n"
        # Supercell
        if hasattr(cell, "smatrix"):
            string += f"#   Supercell matrix:  {cell.smatrix.flatten()}\n"

        # string += '\n'
    else:
        string = ""
    # Write lattice
    # if decorated and cell.pbc:
    #     string += '# Lattice:\n'
    #
    # Order lattice by lengths (doesn't do anything right now)
    if any(cell.pbc):
        lengths = [norm(lv) for lv in cell.get_cell()]
    else:
        lengths = [norm(pos) for pos in cell.get_positions()]
    lv_args = range(len(lengths))  # np.argsort(lengths)

    # for latvec, constraint in zip(latvecs[lv_args], cell.constraints_lv[lv_args]):
    if any(cell.pbc):
        latvecs = clean_matrix(cell.get_cell())
        for latvec, constraint in zip(latvecs, cell.constraints_lv):

            if decorated:
                string += f"  lattice_vector "
            else:
                string += f"lattice_vector "
            #
            string += f"{latvec[0]: .{n_geom_digits}e} "
            string += f"{latvec[1]: .{n_geom_digits}e} "
            string += f"{latvec[2]: .{n_geom_digits}e}\n"
            if constraint:
                string += "constrain_relaxation .true.\n"

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
        positions = clean_matrix(cell.get_scaled_positions(wrap=wrap))
        atompos = "atom_frac"
    else:
        # if decorated:
        #     string += '\n# Cartesian positions:\n'
        #
        positions = clean_matrix(cell.get_positions())
        atompos = "atom"
    #
    if velocities:
        vels = cell.get_velocities()
    #
    for ii, (pos, sym) in enumerate(zip(positions, symbols)):
        if decorated:
            string += f"  {atompos:9s}  "
        else:
            string += f"{atompos:s}  "
        #
        string += f"{pos[lv_args[0]]: .{n_geom_digits}e} "
        string += f"{pos[lv_args[1]]: .{n_geom_digits}e} "
        string += f"{pos[lv_args[2]]: .{n_geom_digits}e}  "
        string += f"{sym:3s}\n"
        if velocities:
            vel = vels[ii]
            if decorated:
                string += "    velocity "
            else:
                string += "velocity "
            string += f"{vel[lv_args[0]]: .{n_geom_digits}e} "
            string += f"{vel[lv_args[1]]: .{n_geom_digits}e} "
            string += f"{vel[lv_args[2]]: .{n_geom_digits}e}\n"
    #
    return string


def inform(cell, fname=None, verbosity=1, symprec=symprec):
    """geometry information

    Parameters
    ----------
    cell: ase.atoms.Atoms
        The cell to convert to geometry.in
    fname: str
        Path to the geometry.in file
    verbosity: int
        How much information to print to the screen
    symprec: float
        Tolerance for determining the symmetry and space group of a material
    """
    unique_symbols, multiplicity = np.unique(cell.symbols, return_counts=True)
    # Structure info:
    talk(f"\nGeometry info for:")
    talk(f"  input geometry:    {get_sysname(cell)}")
    if fname:
        talk(f"  from:              {fname}")
    talk(f"  Symmetry prec.:    {symprec}")
    talk(f"  Number of atoms:   {len(cell)}")

    msg = ", ".join([f"{s} ({m})" for (s, m) in zip(unique_symbols, multiplicity)])
    talk(f"  Species:           {msg}")
    talk(f"  Periodicity:       {cell.pbc}")
    if verbosity > 0 and any(cell.pbc):
        talk(f"  Lattice:  ")
        for vec in cell.cell:
            talk(f"    {vec}")

    print("")

    if symprec is not None:
        sds = get_symmetry_dataset(cell, symprec=symprec)

        talk(f"  Spacegroup:          {sds.international} ({sds.number})")
        if sds.number > 1:
            msg = "  Wyckoff positions:   "
            talk(msg + ", ".join(f"{c}*{w}" for (w, c) in sds.wyckoffs_unique))
            msg = "  Equivalent atoms:    "
            talk(msg + ", ".join(f"{c}*{a}" for (a, c) in sds.equivalent_atoms_unique))

        if verbosity > 1:
            talk(f"  Standard lattice:  ")
            for vec in sds.std_lattice:
                talk(f"    {vec}")

        if verbosity > 1:
            talk(f"  Special k points:")
            for key, val in get_special_points(cell).items():
                talk(f"    {key}: {val}")

    # Info
    for (key, val) in cell.info.items():
        talk(f"  {key:10s}: {val}")

    # lengths and angles
    if verbosity > 0:
        la = cell.get_cell_lengths_and_angles()
        talk("\nCell lengths and angles [\u212B, °]:")
        talk("  a, b, c: {}".format(" ".join([f"{l:11.4f}" for l in la[:3]])))
        angles = "  \u03B1, \u03B2, \u03B3: "
        values = "{}".format(" ".join([f"{l:11.4f}" for l in la[3:]]))
        talk(angles + values)
        talk(f"  Volume:  {cell.get_volume():11.4f} \u212B**3")

        if cell.get_velocities() is not None:
            v = cell.get_momenta().sum(axis=0) / v_unit / cell.get_masses().sum()
            talk(f"\n Net velocity: {v} \u212B/ps")
