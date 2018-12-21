import datetime
import numpy as np
from scipy.linalg import norm
from hilde.konstanten.io import n_geom_digits
from hilde.konstanten.symmetry import symprec
from hilde.helpers.numerics import clean_matrix
from hilde.structure.misc import get_sysname
from hilde.spglib.wrapper import get_symmetry_dataset


def get_aims_string(cell, decorated=True, scaled=None, velocities=False, wrap=True):
    """ print the string that is geometry.in """
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


def inform(cell, dft=False, fname=None, verbosity=0, symprec=symprec):
    # Structure info:
    print(f"\nGeometry info for:")
    print(f"  input geometry:    {get_sysname(cell)}")
    if fname:
        print(f"  from:              {fname}")
    print(f"  Symmetry prec.:    {symprec}")
    print(f"  Number of atoms:   {len(cell)}")
    print(f"  Species:           {', '.join(np.unique(cell.symbols))}")
    print(f"  Periodicity:       {cell.pbc}")
    if any(cell.pbc):
        print(f"  Lattice:  ")
        for vec in cell.cell:
            print(f"    {vec}")

    print(f"")

    if symprec is not None:
        sds = get_symmetry_dataset(cell, symprec=symprec)

        print(f"  Spacegroup:          {sds.international} ({sds.number})")
        print(
            f"  Wyckoff positions:   "
            + ", ".join(f"{c}*{w}" for (w, c) in sds.wyckoffs_unique)
        )
        print(
            f"  Equivalent atoms:    "
            + ", ".join(f"{c}*{a}" for (a, c) in sds.equivalent_atoms_unique)
        )
        print(f"  Standard lattice:  ")
        for vec in sds.std_lattice:
            print(f"    {vec}")

    # Info
    for ii, (key, val) in enumerate(cell.info.items()):
        print(f"  {key:10s}: {val}")

    # lengths and angles
    la = cell.get_cell_lengths_and_angles()
    print("\nCell lengths and angles [\u212B, °]:")
    print("  a, b, c: {}".format(" ".join([f"{l:11.4f}" for l in la[:3]])))
    print(
        "  \u03B1, \u03B2, \u03B3: {}".format(" ".join([f"{l:11.4f}" for l in la[3:]]))
    )
    print(f"  Volume:  {cell.get_volume():11.4f} \u212B**3")

    if hasattr(cell, "dft_calculation") and dft:
        print(f"  DFT calc performed:    yes")
        if verbosity == 1:
            dft_dict = cell.dft_calculation.__dict__
            print(f"Full Summary of DFT Results:")
            for key in dft_dict.keys():
                print(f"  {key:25s}: {dft_dict[key]}")
        #
        print(f"Compact summary of DFT Results:")
        print(f"  Forces computed:       {len(cell.get_forces()) > 0}")
        print(f"  Stress computed:       {cell.get_stress() is not None}")
        if cell.get_stress() is not None:
            print(f"  Pressure:              {cell.get_pressure():12.4e} eV")
        print(f"  Total energy:          {cell.get_total_energy():12.4e} eV")
        if cell.get_chemical_potential():
            print(f"  Fermi level:           {cell.get_chemical_potential():12.4e} eV")
        else:
            print(f"  Fermi level:           None")
        print(f"  HOMO-LUMO gap:         {cell.get_homo_lumo_gap():12.4e} eV")
