""" hilde quality of life """
import numpy as np

from ase.io import read
from phonopy.file_IO import parse_FORCE_CONSTANTS
from phonopy.structure.atoms import PhonopyAtoms

from hilde.helpers import Timer, progressbar
from hilde.structure.convert import to_Atoms
from hilde.helpers.fileformats import last_from_yaml
from hilde.helpers.converters import input2dict
from hilde.phonopy._defaults import displacement_id_str
from hilde.spglib.wrapper import get_symmetry_dataset
from hilde.structure.convert import to_Atoms


def last_calculation_id(trajectory):
    """ return the id of the last computed supercell """
    disp_id = -1

    try:
        dct = last_from_yaml(trajectory)
        disp_id = dct["info"][displacement_id_str]
    except (FileNotFoundError, KeyError):
        pass

    return disp_id


def to_phonopy_atoms(atoms):
    """ convert ase.Atoms to PhonopyAtoms """
    phonopy_atoms = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        masses=atoms.get_masses(),
        positions=atoms.get_positions(wrap=True),
    )
    return phonopy_atoms


def enumerate_displacements(cells, info_str=displacement_id_str):
    """ Assign a displacemt id to every atoms obect in cells.

    Parameters:
        cells (list): atoms objects created by, e.g., phonopy
        info_str (str): how to name the child

    Returns:
        list: cells with id attached to atoms.info (inplace)

    """
    for nn, scell in enumerate(cells):
        if scell is None:
            continue
        scell.info[info_str] = nn


def get_supercells_with_displacements(phonon):
    """ Create a phonopy object and supercells etc. """

    supercell = to_Atoms(
        phonon.get_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": phonon.get_supercell_matrix().T.flatten().tolist(),
        },
    )

    scells = phonon.get_supercells_with_displacements()

    supercells_with_disps = [to_Atoms(cell) for cell in scells]

    enumerate_displacements(supercells_with_disps)

    return phonon, supercell, supercells_with_disps


def metadata2dict(phonon, calculator):
    """ convert metadata information to plain dict """

    atoms = to_Atoms(phonon.get_primitive())

    prim_data = input2dict(atoms)

    phonon_dict = {
        "version": phonon.get_version(),
        "primitive": prim_data["atoms"],
        "supercell_matrix": phonon.get_supercell_matrix().T.astype(int).tolist(),
        "symprec": float(phonon.get_symmetry().get_symmetry_tolerance()),
        "displacement_dataset": phonon.get_displacement_dataset(),
    }

    try:
        displacements = phonon.get_displacements()
        phonon_dict.update({"displacements": displacements})
    except AttributeError:
        pass

    supercell = to_Atoms(phonon.get_supercell())
    supercell_data = input2dict(supercell, calculator)

    return {str(phonon.__class__.__name__): phonon_dict, **supercell_data}


def get_force_constants_from_trajectory(
    trajectory, supercell=None, reduce_fc=False, two_dim=False
):
    """
    Remaps the phonopy force constants into an fc matrix for a new structure
    Parameters:
        trajectory (hilde.Trajectory): phonopy trajectory
        supercell (ase.Atoms): Atoms Object to map force constants onto
        reduce_fc (bool): return in [N_prim, N_sc, 3, 3]  shape
        two_dim (bool): return in [3*N_sc, 3*N_sc] shape
    Returns:
        (np.ndarray): new force constant matrix
    """
    from hilde.phonopy.postprocess import postprocess

    if reduce_fc and two_dim:
        raise IOError("Only one of reduce_fc and two_dim can be True")

    phonon = postprocess(trajectory)

    if supercell is None:
        supercell = to_Atoms(phonon.get_supercell())
        if reduce_fc:
            return phonon.get_force_constants()

    return remap_force_constants(
        phonon.get_force_constants(),
        to_Atoms(phonon.get_primitive()),
        to_Atoms(phonon.get_supercell()),
        supercell,
        reduce_fc,
        two_dim,
    )


def remap_force_constants(
    force_constants,
    primitive,
    supercell,
    new_supercell=None,
    reduce_fc=False,
    two_dim=False,
    tol=1e-5,
    eps=1e-13,
):
    """remap force constants [N_prim, N_sc, 3, 3] to [N_sc, N_sc, 3, 3]

    Parameters:
        force_constants (np.ndarray): force constants in [N_prim, N_sc, 3, 3] shape
        primitive (ase.Atoms): primitive cell for reference
        supercell (ase.Atoms): supercell for reference
        new_supercell (ase.Atoms, optional): supercell to map to (default)
        reduce_fc (bool): return in [N_prim, N_sc, 3, 3]  shape
        two_dim (bool): return in [3*N_sc, 3*N_sc] shape
        tol (float): tolerance to discern pairs
        eps (float): finite zero

    Returns:
        newforce_constants (np.ndarray)

    """
    from hilde.spglib.wrapper import get_symmetry_dataset  # , standardize_cell

    timer = Timer("remap force constants")

    if new_supercell is None:
        new_supercell = supercell.copy()

    primitive.wrap(eps=tol)
    supercell.wrap(eps=tol)

    # prim = standardize_cell(supercell, to_primitve=True, no_idealize=True)
    # prim.set_scaled_positions(prim.get_scaled_positions(wrap=True))
    # primitive.set_scaled_positions(primitive.get_scaled_positions(wrap=True))
    # supercell.set_scaled_positions(supercell.get_scaled_positions(wrap=True))

    # pos_diff = np.abs(primitive.get_scaled_positions() - prim.get_scaled_positions())
    # pos_diff = pos_diff.flatten()
    # pos_diff -= np.floor(pos_diff + eps)
    # fail_cell = np.max(np.abs(primitive.cell - prim.cell).flatten()) > 1000 * eps
    # fail_pos = np.max(pos_diff) > 1000 * eps
    # if fail_cell or fail_pos:
    #     msg = "primitive cell of the supercell and given primitive cell NOT equal"
    #     raise IOError(msg)

    n_sc = len(supercell)
    n_sc_new = len(new_supercell)

    sds = get_symmetry_dataset(new_supercell)
    map2prim = sds.mapping_to_primitive
    uc_index = np.unique(map2prim)

    sc_r = np.zeros((force_constants.shape[0], force_constants.shape[1], 3))
    for aa, a1 in enumerate(primitive):
        diff = supercell.positions - a1.position
        p2s = np.where(np.sum(np.abs(diff), axis=1) < tol)[0][0]
        sc_r[aa] = supercell.get_distances(p2s, range(n_sc), mic=True, vector=True)

    ref_struct_pos = new_supercell.get_scaled_positions(wrap=True)

    fc_out = np.zeros((n_sc_new, n_sc_new, 3, 3))
    for a1, r0 in progressbar(enumerate(new_supercell.positions)):
        uc_index = map2prim[a1]
        for sc_a2, sc_r2 in enumerate(sc_r[uc_index]):
            r_pair = r0 + sc_r2
            sc_temp = new_supercell.get_cell(complete=True)
            r_pair = np.linalg.solve(sc_temp.T, r_pair.T).T % 1.0
            for a2 in range(n_sc_new):
                r_diff = np.abs(r_pair - ref_struct_pos[a2])
                # Integer value is the equivalent of 0.0
                r_diff -= np.floor(r_diff + eps)
                if np.sum(np.abs(r_diff)) < tol:
                    fc_out[a1, a2, :, :] += force_constants[uc_index, sc_a2, :, :]

    timer()

    if two_dim:
        return fc_out.swapaxes(1, 2).reshape(2 * (3 * fc_out.shape[1],))

    if reduce_fc:
        return reduce_force_constants(fc_out, map2prim)

    return fc_out


def reduce_force_constants(fc_full, map2prim):
    """reduce force constants from [N_sc, N_sc, 3, 3] to [N_prim, N_sc, 3, 3]"""
    uc_index = np.unique(map2prim)
    fc_out = np.zeros((len(uc_index), fc_full.shape[1], 3, 3))
    for ii, _ in enumerate(uc_index):
        fc_out[ii, :, :, :] = fc_full[uc_index, :, :, :]

    return fc_out


def parse_phonopy_force_constants(
    uc_filename="geometry.primitive",
    sc_filename="geometry.supercell",
    fc_filename="FORCE_CONSTANTS",
    two_dim=True,
    eps=1e-13,
    tol=1e-5,
    format="aims",
):
    """parse phonopy FORCE_CONSTANTS file and return as 2D array

    Parameters:
        uc_filename (str/Path): primitive unit cell
        sc_filename (str/Path): supercell
        fc_filename (str/Path): phonopy forceconstant file to parse
        two_dim (bool): return in [3*N_sc, 3*N_sc] shape
        eps (float): finite zero
        tol (float): tolerance to discern pairs

    Returns:
            force_constant (np.ndarray(dtype=float)):
                Force constants in (3*N_sc, 3*N_sc) shape

    """

    if "poscar" in uc_filename.lower():
        format = "vasp"

    uc = read(uc_filename, format=format)
    sc = read(sc_filename, format=format)

    fc = parse_FORCE_CONSTANTS(fc_filename)

    fc = remap_force_constants(fc, uc, sc, two_dim=two_dim, tol=tol, eps=eps)

    return fc
