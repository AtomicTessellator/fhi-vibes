""" vibes quality of life """
from pathlib import Path

from ase import Atoms
from ase.io import read
from phonopy.file_IO import parse_FORCE_CONSTANTS, read_force_constants_hdf5
from phonopy.structure.atoms import PhonopyAtoms

from vibes.helpers.converters import input2dict
from vibes.helpers.fileformats import last_from_yaml
from vibes.helpers.force_constants import remap_force_constants
from vibes.phonopy._defaults import displacement_id_str
from vibes.structure.convert import to_Atoms


_prefix = "phonopy.utils"


def last_calculation_id(trajectory_file):
    """Return the id of the last computed supercell

    Parameters
    ----------
    trajectory_file: str or Path
        Path to the trajectory file with calculated supercells

    Returns
    -------
    disp_id: int
        The id of the last computed supercell
    """
    disp_id = -1

    try:
        dct = last_from_yaml(trajectory_file, allow_empty=True)
        if dct is not None:
            disp_id = dct["info"][displacement_id_str]
    except (FileNotFoundError):
        pass

    return disp_id


def to_phonopy_atoms(atoms):
    """Convert ase.atoms.Atoms to PhonopyAtoms

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Atoms to convert

    Returns
    -------
    phonopy_atoms: PhonopyAtoms
        The PhonopyAtoms for the same structure as atoms
    """
    phonopy_atoms = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        masses=atoms.get_masses(),
        positions=atoms.get_positions(wrap=True),
    )
    return phonopy_atoms


def enumerate_displacements(cells, info_str=displacement_id_str):
    """Assign a displacemt id to every atoms obect in cells.

    Parameters
    ----------
    cells: list
        Atoms objects created by, e.g., phonopy
    info_str: str
        How to name the cell

    Returns
    -------
    list
        cells with id attached to atoms.info (inplace)
    """
    for nn, scell in enumerate(cells):
        if scell is None:
            continue
        scell.info[info_str] = nn


def get_supercells_with_displacements(phonon):
    """ Create a phonopy object and supercells etc.

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The phonopy object with displacement_dataset

    Returns
    -------
    phonon: phonopy.Phonopy
        The phonopy object with displacement_dataset, and displaced supercells
    supercell: ase.atoms.Atoms
        The undisplaced supercell
    supercells_with_disps: list of ase.atoms.Atoms
        All of the supercells with displacements
    """

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
    """Convert metadata information to plain dict

    Parameters
    ----------
    phonon: phonopy.Phonopy
        Phonopy Object to get metadata for
    calculator: ase.calculators.calulator.Calculator
        The calculator for the force calculation

    Returns
    -------
    dict
        Metadata as a plain dict with the following items

        atoms: dict
            Dictionary representation of the supercell
        calculator: dict
            Dictionary representation of calculator
        Phonopy/Phono3py: dict
            Phonopy/Phono3py metadata with the following items

            version: str
                Version string for the object
            primitive: dict
                dictionary representation of the primitive cell
            supercell_matrix: list
                supercell matrix as a list
            symprec: float
                Tolerance for determining the symmetry/space group of the primitive cell
            displacement_dataset: dict
                The displacement dataset for the phonon calculation
    """

    atoms = to_Atoms(phonon.get_unitcell())

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
    trajectory_file, supercell=None, reduce_fc=False, two_dim=False
):
    """Remaps the phonopy force constants into an fc matrix for a new structure

    Parameters
    ----------
    trajectory_file: vibes.Trajectory
        phonopy trajectory
    supercell: ase.atoms.Atoms
        Atoms Object to map force constants onto
    reduce_fc: bool
        return in [N_prim, N_sc, 3, 3]  shape
    two_dim: bool
        return in [3*N_sc, 3*N_sc] shape

    Returns
    -------
    np.ndarray
         new force constant matrix

    Raises
    ------
    IOError
        If both reduce_fc and two_dim are True
    """
    from vibes.phonopy.postprocess import postprocess

    if reduce_fc and two_dim:
        raise IOError("Only one of reduce_fc and two_dim can be True")

    phonon = postprocess(trajectory_file)

    if supercell is None:
        supercell = to_Atoms(phonon.get_supercell())
        if reduce_fc:
            return phonon.get_force_constants()

    return remap_force_constants(
        phonon.get_force_constants(),
        to_Atoms(phonon.get_unitcell()),
        to_Atoms(phonon.get_supercell()),
        supercell,
        reduce_fc,
        two_dim,
    )


def parse_phonopy_force_constants(
    fc_file="FORCE_CONSTANTS",
    primitive="geometry.primitive",
    supercell="geometry.supercell",
    fortran=True,
    symmetrize=True,
    two_dim=True,
    eps=1e-5,
    tol=1e-4,
    format="aims",
):
    """parse phonopy FORCE_CONSTANTS file and return as 2D array

    Args:
        fc_file (Pathlike): phonopy forceconstant file to parse
        primitive (Atoms or Pathlike): either unitcell as Atoms or where to find it
        supercell (Atoms or Pathlike): either supercell as Atoms or where to find it
        symmetrize (bool): make force constants symmetric
        two_dim (bool): return in [3*N_sc, 3*N_sc] shape
        eps: finite zero
        tol: tolerance to discern pairs and wrap
        format: File format for the input geometries

    Returns:
        Force constants in (3*N_sc, 3*N_sc) shape
    """
    if "hdf5" in str(fc_file):
        fc = read_force_constants_hdf5(fc_file)
    else:
        fc = parse_FORCE_CONSTANTS(fc_file)

    if fc.shape[0] == fc.shape[1]:
        if two_dim:
            fc = fc.swapaxes(1, 2).reshape(3 * fc.shape[1], 3 * fc.shape[1])

        return fc

    if two_dim:
        if isinstance(primitive, Atoms):
            uc = primitive
        elif Path(primitive).exists():
            uc = read(primitive, format=format)
        else:
            raise RuntimeError("primitive cell missing")

        if isinstance(supercell, Atoms):
            sc = supercell
        elif Path(supercell).exists():
            sc = read(supercell, format=format)
        else:
            raise RuntimeError("supercell missing")

        fc = remap_force_constants(
            fc,
            uc,
            sc,
            two_dim=two_dim,
            tol=tol,
            eps=eps,
            fortran=fortran,
            symmetrize=symmetrize,
        )

    return fc
