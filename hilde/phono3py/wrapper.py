"""
A leightweight wrapper for Phono3py
"""

from collections import namedtuple
import numpy as np
from phono3py.phonon3 import Phono3py
from hilde import konstanten as const
from hilde.helpers.supercell import make_cubic_supercell
from hilde.phonopy import enumerate_displacements, metadata2dict
from hilde.structure.convert import to_Atoms, to_phonopy_atoms
from hilde.helpers.attribute_dict import AttributeDict as adict

defaults = adict(
    {
        "displacement": 0.03,
        "cutoff_pair_distance": 100.0,
        "symprec": 1e-5,
        "q_mesh": [11, 11, 11],
        "log_level": 2,
    }
)


def prepare_phono3py(
    atoms,
    fc2_supercell_matrix=None,
    fc3_supercell_matrix=None,
    q_mesh=defaults.q_mesh,
    fc2=None,
    fc3=None,
    disp=defaults.displacement,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    symmetrize_fc3q=False,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
):
    """
    Prepare a Phono3py object
    Args:
        atoms (ASE Atoms Object): Structure to calculate phonons
        fc2_supercell_matrix (np.ndarray): second order supercell matrix
        fc3_supercell_matrix (np.ndarray): third order supercell matrix
        q_mesh (np.ndarray): q-mesh size
        fc2 (np.ndarray): second order force constants
        fc3 (np.ndarray): third order force constants
        disp (float): finite displacement value
        cutoff_pair_distance (float): Sets all interatomic interactions to zero for atoms farther apart than this number
        symmetrize_fc3 (bool): if True symmeterize third order force constants
        symprec (float): Tolerance factor detecting the space group
        log_level (int): log level for the phonopy object
    Returns: The Phono3py object
    """
    if fc3_supercell_matrix is None:
        raise ValueError("Please define an fc3_supercell_matrix")
    ph_atoms = to_phonopy_atoms(atoms, wrap=True)

    phonon3 = Phono3py(
        ph_atoms,
        supercell_matrix=np.transpose(fc3_supercell_matrix),
        phonon_supercell_matrix=np.transpose(fc2_supercell_matrix),
        mesh=q_mesh,
        symprec=symprec,
        is_symmetry=True,
        symmetrize_fc3q=symmetrize_fc3q,
        frequency_factor_to_THz=const.omega_to_THz,
        log_level=log_level,
    )

    phonon3.generate_displacements(
        distance=disp, cutoff_pair_distance=cutoff_pair_distance
    )

    if fc2 is not None:
        phonon3.set_fc2(fc2)
    if fc3 is not None:
        phonon3.set_fc3(fc3)

    return phonon3

def preprocess_fireworks(
    atoms,
    calc,
    kpt_density=None,
    supercell_matrix=None,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    up_kpoint_from_pc=False,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
):
    """ wrap the Phono3py preprocess for workflow """
    if supercell_matrix is None:
        raise InputError("The supercell_matrix must be defined")

    # Phonopy preprocess
    phonon3, _, supercell, _, scs = phono3py_preprocess(
        atoms=atoms,
        fc3_supercell_matrix=supercell_matrix,
        n_atoms_in_sc_3=n_atoms_in_sc,
        disp=displacement,
        cutoff_pair_distance=cutoff_pair_distance,
        symprec=symprec,
    )
    # make sure forces are computed (aims only)
    if calc.name == "aims":
        calc.parameters["compute_forces"] = True

    if kpt_density is not None:
        update_k_grid(supercell, calc, kpt_density)

    metadata = metadata2dict(atoms, calc, phonon3)

    to_run_scs = []
    for sc in scs:
        if sc:
            sc.calc = calc
            to_run_scs.append(sc)
    return calc, supercell, to_run_scs, phonon3, metadata

def preprocess(
    atoms,
    fc2_supercell_matrix=None,
    n_atoms_in_sc_2=None,
    fc3_supercell_matrix=None,
    n_atoms_in_sc_3=None,
    q_mesh=defaults.q_mesh,
    disp=defaults.displacement,
    cutoff_pair_distance=defaults.cutoff_pair_distance,
    symprec=defaults.symprec,
    log_level=defaults.log_level,
):
    """
    Set up a Phono3py object and generate all the necessary supercells
    Args:
        atoms (ASE Atoms Object): Structure to calculate phonons
        fc2_supercell_matrix (np.ndarray): second order supercell matrix
        fc3_supercell_matrix (np.ndarray): third order supercell matrix
        q_mesh (np.ndarray): q-mesh size
        disp (float): finite displacement value
        cutoff_pair_distance (float): Sets all interatomic interactions to zero for atoms farther apart than this number
        symprec (float): Tolerance factor detecting the space group
        log_level (int): log level for the phonopy object
    Returns:
        The Phono3py object, second order supercell, third order supercell, list of second order supercells with displacements, list of third order supercells with displacements
    """
    if fc3_supercell_matrix is None and n_atoms_in_sc_3:
        _, fc3_supercell_matrix = make_cubic_supercell(atoms, n_atoms_in_sc_3)
    elif fc3_supercell_matrix is None:
        raise InputError("Either fc3_supercell_matrix or natoms_in_sc must be specified")
    if fc2_supercell_matrix is None and n_atoms_in_sc_2:
        _, fc2_supercell_matrix = make_cubic_supercell(atoms, n_atoms_in_sc_2)
    elif fc2_supercell_matrix is None:
        fc2_supercell_matrix = fc3_supercell_matrix
    phonon3 = prepare_phono3py(
        atoms,
        fc2_supercell_matrix=np.array(fc2_supercell_matrix).reshape(3, 3),
        fc3_supercell_matrix=np.array(fc3_supercell_matrix).reshape(3, 3),
        q_mesh=q_mesh,
        disp=disp,
        cutoff_pair_distance=cutoff_pair_distance,
        symprec=symprec,
        log_level=log_level,
    )

    # phonpoy supercells
    fc2_supercell = to_Atoms(
        phonon3.get_phonon_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": np.asarray(fc2_supercell_matrix).flatten().tolist(),
        },
    )

    fc3_supercell = to_Atoms(
        phonon3.get_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": np.asarray(fc3_supercell_matrix).flatten().tolist(),
        },
    )

    scells = phonon3.get_phonon_supercells_with_displacements()
    fc2_supercells_with_disps = [to_Atoms(cell) for cell in scells]

    scells = phonon3.get_supercells_with_displacements()
    fc3_supercells_with_disps = [to_Atoms(cell) for cell in scells]

    enumerate_displacements([*fc2_supercells_with_disps, *fc3_supercells_with_disps])
    cols = (
        "phonon3",
        "fc2_supercell",
        "fc3_supercell",
        "fc2_supercells_with_displacements",
        "fc3_supercells_with_displacements",
    )
    pp = namedtuple("phono3py_preprocess", " ".join(cols))

    return pp(
        phonon3,
        fc2_supercell,
        fc3_supercell,
        fc2_supercells_with_disps,
        fc3_supercells_with_disps,
    )


def get_forces(supercells_computed):
    """ Return force_sets taking care of supercells that were not computed
    because of cutoff. """

    zero_force = np.zeros([len(supercells_computed[0]), 3])
    force_sets = []
    for scell in supercells_computed:
        if scell is None:
            force_sets.append(zero_force)
        else:
            force_sets.append(scell.get_forces())

    if len(force_sets) != len(supercells_computed):
        print(
            "len(force_sets), len(supercells):",
            len(force_sets),
            len(supercells_computed),
        )
        raise RuntimeError("Number of computed supercells incorrect.")

    return force_sets
