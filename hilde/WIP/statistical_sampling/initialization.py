''' Initialize thermally displaced cells for statistical sampling '''
import ase
from ase import units as u
from ase.md.velocitydistribution import PhononHarmonics

import numpy as np

from phonopy.file_IO import parse_FORCE_CONSTANTS

from hilde import konstanten as const
from hilde.helpers.converters import input2dict
from hilde.helpers.supercell import make_supercell, find_cubic_cell
from hilde.phonopy.postprocess import postprocess as postprocess_ph
from hilde.phonopy.utils import get_force_constants_from_trajectory
from hilde.phonopy.wrapper import preprocess as preprocess_ph
from hilde.structure.convert import to_Atoms

def prepare_phonon_harmonic_sampling(
    atoms,
    force_constants,
    temperature,
    n_samples=1,
    deterministic=True,
    rng=np.random,
    **kwargs
):
    '''
    Generates a list of displaced supercells based on a thermal excitation of phonons
    Parameters:
        atoms (ase.atoms.Atoms): Non-displaced supercell
        force_constants(np.ndarray(shape=(3*len(atoms), 3*len(atoms)))): Force constant matrix for atoms
        temperature (float): Temperature to populate the phonons modes at
        n_samples(int): number of samples to generate
        deterministic(bool): If True then displace atoms with +/- the amplitude according to PRB 94, 075125
    Returns (list of ase.atoms.Atoms): The thermally displaced supercells
    '''
    thermally_disp_cells = []
    for ii in range(n_samples):
        td_cell = atoms.copy()
        PhononHarmonics(
            td_cell,
            force_constants,
            quantum=False,
            temp=temperature * u.kB,
            rng=rng,
            plus_minus=deterministic,
            failfast=True,
        )
        thermally_disp_cells.append(td_cell)
    return thermally_disp_cells

def preprocess(
    atoms,
    phonon_file,
    temperatures=None,
    debye_temp_fact=None,
    n_samples=1,
    deterministic=True,
    rng_seed=None,
    **kwargs,
):
    '''
    Sets ups the statistical sampling
    Parameters:
        atoms (ase.atoms.Atoms): structure to perform statistical sampling on
        phonon_file (str): String to the phonopy trajectory
        temperatures (list (floats)): List of temperatures to excite the phonons to
        debye_temp_fact (list(floats)): List of factors to multiply the debye temperature by to populate temperatures
        n_samples (int): number of samples to calculate
        deterministic(bool): If True then displace atoms with +/- the amplitude according to PRB 94, 075125
        rng_seed (int): Seed for the random number generator used by PhononHarmonics
    Additional arguments in kwargs:
        supercell_matrix (np.ndarray(int)): Supercell matrix used to create the supercell from atoms
    Returns: (list(dicts)): A list of thermally displaced supercells to calculate the forces on
    '''
    if temperatures is None:
        temperatures = list()

    # Set up supercell and Force Constants
    if "supercell_matrix" in kwargs:
        sc = make_supercell(atoms, np.array(kwargs["supercell_matrix"]).reshape(3,3))
    else:
        sc = None
    if phonon_file:
        phonon = postprocess_ph(phonon_file, write_files=False, calculate_full_force_constants=False)
        if sc is None:
            sc = to_Atoms(phonon.get_supercell())
        force_constants = get_force_constants_from_trajectory(phonon_file, sc, two_dim=True)
    else:
        phonon = None

    if sc is None:
        sc = atoms.copy()
    assert(force_constants.shape[0] == force_constants.shape[1] == len(sc)*3 )

    # If using Debye temperature calculate it
    if debye_temp_fact is not None:
        if phonon is None:
            raise IOError("Debye Temperature must be calculated with phonopy, please add phonon_file")
        phonon.set_mesh([51, 51, 51])
        phonon.set_total_DOS(freq_pitch=0.01)
        phonon.set_Debye_frequency()
        debye_temp = phonon.get_Debye_frequency() * const.THzToEv / const.kB
        temperatures += [tt * debye_temp for tt in debye_temp_fact]
    elif temperatures is None:
        raise IOError("temperatures must be given to do harmonic analysis")

    if "deterministic" in kwargs:
        deterministic = kwargs["deterministic"]

    # Generate metadata
    metadata = {
        "ase_vesrsion": ase.__version__,
        "n_samples_per_temperature": n_samples,
        "temperatures": temperatures,
        "deterministic": deterministic,
        "force_constants": force_constants,
        "supercell": input2dict(sc),
        "primitive": input2dict(atoms),
        **input2dict(sc),
    }
    if not deterministic:
        if rng_seed is None:
            rng_seed = np.random.randint(2**32-1)
        elif not isinstance(rng_seed, int):
            rng_seed = int(rng_seed)
        metadata["rng_seed"] = rng_seed
    rng = np.random.RandomState(rng_seed)

    # Set up thermally displaced supercells
    td_cells = list()
    for temp in temperatures:
        sc.info["temperature"] = temp
        td_cells += prepare_phonon_harmonic_sampling(
            sc,
            force_constants,
            temp,
            n_samples,
            deterministic,
            rng,
        )
    return td_cells, metadata
