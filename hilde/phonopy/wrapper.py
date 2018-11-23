"""
A leightweight wrapper for Phonopy()
"""

from collections import namedtuple
import json
from pathlib import Path
import numpy as np
from phonopy import Phonopy
from hilde import konstanten as const
from hilde.helpers import brillouinzone as bz
from hilde.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy, to_dict
from hilde.phonopy import enumerate_displacements, displacement_id_str
from hilde.structure.convert import to_Atoms, to_phonopy_atoms
from hilde.helpers.maths import get_3x3_matrix
from hilde.spglib.wrapper import map_unique_to_atoms

default_disp = 0.01
default_q_mesh = [35, 35, 35]


def prepare_phonopy(
    atoms, supercell_matrix, fc2=None, disp=default_disp, symprec=1e-5, trigonal=False
):
    """Create a Phonpopy object

    Args:
        atoms (atoms): reference structure (unit cell)
        supercell_matrix (list): supercell matrix
        fc2 (ndarray): second order force constants (optional)
        disp (float): finite displacemt in angstrom (default .01)
        symprec (float): symmetry precision
        trigonal (bool): wether the system is trigonal

    Returns:
        Phonopy(): a phonopy object

    """

    ph_atoms = to_phonopy_atoms(atoms, wrap=True)

    supercell_matrix = get_3x3_matrix(supercell_matrix)

    phonon = Phonopy(
        ph_atoms,
        supercell_matrix=supercell_matrix,
        symprec=symprec,
        is_symmetry=True,
        factor=const.eV_to_THz,
    )

    phonon.generate_displacements(
        distance=disp, is_plusminus="auto", is_diagonal=True, is_trigonal=trigonal
    )

    if fc2 is not None:
        phonon.set_force_constants(fc2)

    return phonon


def preprocess(
    atoms, supercell_matrix, disp=default_disp, symprec=1e-5, trigonal=False, **kwargs
):
    """
    Creates a phonopy object from given input
    Args:
        atoms: atoms object that represents the (primitive) unit cell
        supercell_matrix: supercell matrix
        disp: displacement for the finite displacemt

    Returns:
        namedtuple with the phonon object, the supercell
        and the supercells_with_displacements as ase.atoms
    """

    phonon = prepare_phonopy(
        atoms, supercell_matrix, disp=disp, symprec=symprec, trigonal=trigonal
    )

    supercell = to_Atoms(
        phonon.get_supercell(),
        info={
            "supercell": True,
            "supercell_matrix": np.asarray(supercell_matrix).flatten().tolist(),
        },
    )

    scells = phonon.get_supercells_with_displacements()

    supercells_with_disps = [to_Atoms(cell) for cell in scells]

    enumerate_displacements(supercells_with_disps)

    pp = namedtuple(
        "phonopy_preprocess", "phonon supercell supercells_with_displacements"
    )

    return pp(phonon, supercell, supercells_with_disps)


def get_force_constants(phonon, force_sets=None):
    """
    fkdev: is this necessary?
    Take a Phonopy object and produce force constants from the given forces
    """
    n_atoms = phonon.get_supercell().get_number_of_atoms()

    phonon.produce_force_constants(force_sets)

    force_constants = phonon.get_force_constants()

    if force_constants is not None:
        # convert forces from (N, N, 3, 3) to (3*N, 3*N)
        force_constants = (
            phonon.get_force_constants().swapaxes(1, 2).reshape(2 * (3 * n_atoms,))
        )
        return force_constants
    # else
    raise ValueError("Force constants not yet created, specify force_sets.")


def _postprocess_init(phonon, force_sets=None):
    """
    Make sure that force_constants are present before the actual postprocess is
    performed.
    Args:
        phonon: pre-processed phonon object
        force_constants: computed force_constants (optional)
        force_sets: computed forces (optional)

    Returns:

    """
    if phonon.get_force_constants() is None:
        if force_sets is not None:
            phonon.produce_force_constants(force_sets)
        else:
            exit("** Cannot run postprocess, force_sets have not been provided.")


def get_dos(
    phonon,
    q_mesh=default_q_mesh,
    freq_min=0,
    freq_max="auto",
    freq_pitch=0.1,
    tetrahedron_method=True,
    write=False,
    filename="total_dos.dat",
    force_sets=None,
):
    """ Compute the DOS (and save to file) """

    _postprocess_init(phonon, force_sets)

    phonon.set_mesh(q_mesh)

    if freq_max == "auto":
        freq_max = phonon.get_mesh()[2].max() * 1.05

    phonon.set_total_DOS(
        freq_min=freq_min,
        freq_max=freq_max,
        freq_pitch=freq_pitch,
        tetrahedron_method=tetrahedron_method,
    )

    if write:
        phonon.write_total_DOS()
        Path("total_dos.dat").rename(filename)

    return phonon.get_total_DOS()


def get_bandstructure(phonon, paths=None, force_sets=None):
    """
    Compute bandstructure for given path

    Returns:
        (qpoints, distances, frequencies, eigenvectors, labels)

    """

    _postprocess_init(phonon, force_sets)

    bands, labels = bz.get_bands_and_labels(phonon.primitive, paths)

    phonon.set_band_structure(bands)

    return (*phonon.get_band_structure(), labels)


def plot_bandstructure(phonon, file="bandstructure.pdf", paths=None, force_sets=None):
    """ Plot bandstructure for given path and save to file """

    *_, labels = get_bandstructure(phonon, paths, force_sets)

    plt = phonon.plot_band_structure(labels=labels)

    plt.savefig(file)


def plot_bandstructure_and_dos(
    phonon, q_mesh=default_q_mesh, partial=False, file="bands_and_dos.pdf"
):
    """ Plot bandstructure and PDOS """

    *_, labels = get_bandstructure(phonon)

    if partial:
        phonon.set_mesh(q_mesh, is_eigenvectors=True, is_mesh_symmetry=False)
        phonon.set_partial_DOS(tetrahedron_method=True)
        pdos_indices = map_unique_to_atoms(phonon.get_primitive())
    else:
        phonon.set_mesh(q_mesh)
        phonon.set_total_DOS(tetrahedron_method=True)
        pdos_indices = None

    plt = phonon.plot_band_structure_and_dos(pdos_indices=pdos_indices, labels=labels)
    plt.savefig(file)


def summarize_bandstructure(phonon, fp_file=None):
    from hilde.konstanten.einheiten import THz_to_cm

    get_bandstructure(phonon)

    qpts = np.array(phonon.band_structure.qpoints).reshape(-1, 3)

    freq = np.array(phonon.band_structure.frequencies).reshape(qpts.shape[0], -1)

    gamma_freq = freq[np.where((qpts == np.zeros(3)).all(-1))[0][0]]
    max_freq = np.max(freq.flatten())

    if fp_file:
        print(f"Saving the fingerprint to {fp_file}")
        fp = get_phonon_bs_fingerprint_phononpy(phonon, binning=False)
        fp_dict = to_dict(fp)
        for key,val in fp_dict.items():
            fp_dict[key] = val.tolist()
        with open(fp_file, 'w') as outfile:
            json.dump(fp_dict, outfile, indent=4)

    mf = max_freq
    mf_cm = mf * THz_to_cm
    print(f"The maximum frequency is: {mf:.3f} THz ({mf_cm:.3f} cm^-1)")
    print(f"The frequencies at the gamma point are:")
    print(f"              THz |        cm^-1")
    for ii, freq in enumerate(gamma_freq):
        print(f"{ii+1:3d}: {freq:-12.5f} | {freq*THz_to_cm:-12.5f}")
    return gamma_freq, max_freq


if __name__ == "__main__":
    main()
