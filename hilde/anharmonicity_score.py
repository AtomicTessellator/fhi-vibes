''' Functions to calculate the coefficient of determination for a set of thermally displaced atoms'''
import numpy as np
from hilde.trajectory import reader
from hilde.harmonic_analysis.displacements import get_dR
from hilde.helpers.converters import dict2atoms
from hilde.phonopy.postprocess import postprocess as postprocess_ph
from hilde.structure.convert import to_phonopy_atoms
from hilde.spglib.wrapper import get_symmetry_dataset

def get_harmonic_forces(sc, ref_structure, force_constants):
    return -force_constants @ get_dR(sc, ref_structure).flatten()

def get_force_sets_from_trajectory(
    trajectory="statistical_sampling/trajectory.son",
    force_constants=None,
    ref_structure=None,
    temperature=None
):
    # Get the displaced atoms
    scs_displaced, meta = reader(trajectory, True)

    if ref_structure is None:
        ref_structure = dict2atoms(meta["atoms"])

    if force_constants is None:
        force_constants = np.array(meta["force_constants"])

    r2 = dict()
    forces_dft = dict()
    forces_harmonic = dict()
    for sc in scs_displaced:
        if temperature and sc.info["temperature"] != temperature:
            continue
        if sc.info["temperature"] not in forces_dft:
            forces_dft[sc.info['temperature']] = list()
            forces_harmonic[sc.info['temperature']] = list()
        forces_dft[sc.info['temperature']] += list(sc.get_forces().flatten())
        forces_harmonic[sc.info['temperature']] += list(get_harmonic_forces(sc, ref_structure, force_constants))
    if temperature is None:
        return forces_dft, forces_harmonic
    return forces_dft[temperature], forces_harmonic[temperature]

def reshape_forces(forces_dft, forces_harmonic):
    forces_dft = forces_dft.copy().reshape(-1,3)
    forces_harmonic = forces_harmonic.copy().reshape(-1,3)
    return forces_dft, forces_harmonic

def R2(f_data, f_model):
    '''Calculates the coefficient of determination between data points f_data and model f_model'''
    f_data_mean = np.mean(f_data, axis=0)
    Sres = (f_data - f_model) @ (f_data - f_model)
    Stot = (f_data - f_data_mean) @ (f_data - f_data_mean)
    return 1 - Sres / Stot

def get_r2(forces_dft, forces_harmonic):
    '''Calculates the coefficient of determination with forces_dft as f_data and harmonic forces as f_model'''
    return R2(forces_dft, forces_harmonic)

def per_atom_r2(
    ref_structure,
    forces_dft,
    forces_harmonic,
    by_equivilent_atom=True,
    by_atom_species=False,
):
    forces_dft, forces_harmonic = reshape_forces(forces_dft, forces_harmonic)
    sds = get_symmetry_dataset(ref_structure)

    forces_by_atom = []

    if by_equivilent_atom:
        compare_list = sds.equivalent_atoms
        unique_atoms, multiplicity  = np.unique(equiv_atoms, return_counts=1)
        for at in unique_atoms:
            forces_by_atom.append([[], []])
    elif by_atom_species:
        compare_list = ref_structure.symbols.copy()
        unique_atoms, multiplicity  = np.unique(ref_structure.symbols, return_counts=1)
        for at in unique_atoms:
            forces_by_atom.append([[], []])
    else:
        compare_list = range(len(ref_structure))
        unique_atoms = compare_list.copy()
        for at in unique_atoms:
            forces_by_atom.append([[], []])

    for ii, atom in enumerate(compare_list):
        atom_ind = np.where(unique_atoms==atom)[0][0]
        forces_by_atom[atom_ind][0] += forces_dft[ii]
        forces_by_atom[atom_ind][1] += forces_harmonic[ii]

    r2_by_atoms = []

    for forces in forces_by_type:
        r2_by_atoms.append(R2(np.array(forces[0]).flatten(), np.array(forces[1]).flatten()))

    return unique_atoms, r2_by_atoms

def get_deviation(forces_dft, forces_harmonic, ref_structure):
    forces_dft, forces_harmonic = reshape_forces(forces_dft, forces_harmonic)

    vol = ref_structure.get_volume()

    diff = np.abs(forces_dft - forces_harmonic)
    dev_vol_nromalized = np.linalg.norm(diff/vol, axis=1)

    for ii in range(3):
        diff[:,ii] /= np.linalg.norm(forces_harmonic, axis=1)

    dev_force_normalized = diff.copy()

    return dev_vol_nromalized, dev_force_normalized

def get_rmse(forces_dft, forces_harmonic):
    rmse = np.sqrt(np.sum((forces_dft - forces_harmonic)**2.0)/forces_harmonic.shape[0])
    return rmse

def get_normalized_rmse(forces_dft, forces_harmonic, ref_structure):
    vol = ref_structure.get_volume()
    rmse_vol_normalized = get_rmse(forces_dft, forces_harmonic) / vol

    forces_dft, forces_harmonic = reshape_forces(forces_dft, forces_harmonic)
    force_harmonic_mag = np.linalg.norm(forces_harmonic, axis=1)
    forces_dft      /= force_harmonic_mag
    forces_harmonic /= force_harmonic_mag
    rmse_force_normalized = get_rmse(forces_dft.flatten(), forces_harmonic.flatten())

    return rmse_vol_normalized, rmse_force_normalized

def get_anharmonicity_score(
    trajectory="statistical_sampling/trajectory.son",
    force_constants=None,
    ref_structure=None,
    r2_file="r2.dat",
    **kwargs
):
    '''
    Calculates the R^2 anharmonicity parameter for a given structure
    Args:
        trajectory (str): File name for the trajectory of statistically sampled supercell displacements
        phonon_file (str): File name for the phonopy trajectory used to generate the force constants used to generate the structures in trajectory
        r2_file(str): File name to output the temperature and r^2 value
    Returns (tuple(float, float)): Temperature, Coefficient of Determination between actual forces and harmonic ones
    '''
    r2 = []
    forces_dft, forces_harmonic = get_force_sets_from_trajectory(
        trajectory,
        force_constants,
        ref_structure,
    )
    for temp in sorted(forces_dft.keys()):
        r2.append([temp, get_r2(np.array(forces_dft[temp]), np.array(forces_harmonic[temp]))])
    np.savetxt(r2_file, r2, header="Temperature, R^2")
    return r2
