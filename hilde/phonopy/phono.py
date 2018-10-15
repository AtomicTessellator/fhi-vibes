"""
A leightweight wrapper for Phonopy()
"""

from collections import namedtuple
from pathlib import Path
from phonopy import Phonopy
from hilde import konstanten as const
from hilde.structure import pAtoms
from hilde.helpers import brillouinzone as bz


def preprocess(atoms, supercell_matrix, disp=0.01, symprec=1e-5, trigonal=False):
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

    ph_atoms = atoms.to_phonopy_atoms()

    phonon = Phonopy(ph_atoms,
                     supercell_matrix=supercell_matrix,
                     symprec=symprec,
                     is_symmetry=True,
                     factor=const.eV_to_THz)

    phonon.generate_displacements(distance=disp,
                                  is_plusminus='auto',
                                  is_diagonal=True,
                                  is_trigonal=trigonal)

    supercell = pAtoms(phonopy_atoms=phonon.get_supercell(),
                       tags=['supercell',
                             ('smatrix', list(supercell_matrix.T.flatten()))])

    supercells_with_disps = [pAtoms(phonopy_atoms=disp, symprec=None,
                                    tags=['supercell',
                                          ('smatrix',
                                           list(supercell_matrix.T.flatten()))])
                             for disp in phonon.get_supercells_with_displacements()]

    pp = namedtuple('phonopy_preprocess', 'phonon supercell supercells_with_displacements')

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
        force_constants = phonon.get_force_constants().swapaxes(1, 2).reshape(2*(3*n_atoms, ))
        return force_constants
    # else
    print('**Force constants not yet created, please specify force_sets.')


def _postprocess_init(phonon,
                     force_sets=None):
    """
    Make sure that force_constants are present before the actual postprocess is performed.
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
            exit('** Cannot run postprocess, force_sets have not been provided.')


def get_dos(phonon,
            q_mesh=[10, 10, 10],
            freq_min=0,
            freq_max='auto',
            freq_pitch=.1,
            tetrahedron_method=True,
            write=False,
            filename='total_dos.dat',
            force_sets=None):
    """
    Compute the DOS (and save to file)
    Args:
        phonon: Phonopy object
        q_mesh: q mesh to evaluate D(q)
        freq_min: freq. to start with [THz]
        freq_max: freq. to stop with [THz]
        freq_pitch: freq. step [THz]
        tetrahedron_method: use the tetrahedron method
        write: should a file be written?
        filename: file that total DOS is written to
        force_sets: give force_sets if force_constants are missing

    Returns:
        tuple: (frequencies, DOS)

    """

    _postprocess_init(phonon, force_sets)

    phonon.set_mesh(q_mesh)

    if freq_max == 'auto':
        freq_max = phonon.get_mesh()[2].max() * 1.05

    phonon.set_total_DOS(freq_min=freq_min,
                         freq_max=freq_max,
                         freq_pitch=freq_pitch,
                         tetrahedron_method=tetrahedron_method)

    if write:
        phonon.write_total_DOS()
        Path('total_dos.dat').rename(filename)

    return phonon.get_total_DOS()


def get_bandstructure(phonon,
                      paths=None,
                      force_sets=None):
    """
    Compute bandstructure for given path

    Args:
        phonon: Phonopy object
        path: path in the Brillouin zone
        force_sets: (optional)

    Returns:
        (qpoints, distances, frequencies, eigenvectors)

    """

    _postprocess_init(phonon, force_sets)


    bands, labels = bz.get_bands_and_labels(phonon.primitive,
                                            paths)

    phonon.set_band_structure(bands)

    return (*phonon.get_band_structure(), labels)
