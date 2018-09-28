"""
A leightweight wrapper for Phonopy()
"""

import os
import numpy as np
from hilde import konstanten as const
from hilde.structure.convert import ASE_to_phonopy_atoms, phonopy_to_ASE_atoms
from phonopy import Phonopy
from collections import namedtuple
from pathlib import Path

def preprocess(atoms, supercell_matrix, disp=0.01, symprec=1e-5, trigonal=False):
    """
    Creates a phonopy object from given input
    Args:
        atoms: atoms object that represents the (primitive) unit cell
        smatrix: supercell matrix
        disp: displacement for the finite displacemt

    Returns:
        namedtuple with the phonon object, the supercell
        and the supercells_with_displacements as ase.atoms
    """

    ph_atoms = ASE_to_phonopy_atoms(atoms)

    phonon = Phonopy(ph_atoms,
                     supercell_matrix = supercell_matrix,
                     symprec          = symprec,
                     is_symmetry      = True,
                     factor           = const.eV_to_THz)

    phonon.generate_displacements(distance     = disp,
                                  is_plusminus ='auto',
                                  is_diagonal  = True,
                                  is_trigonal  = trigonal)

    supercell = phonopy_to_ASE_atoms(phonon.get_supercell())
    supercells_with_disps = [phonopy_to_ASE_atoms(disp)
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

    fc = phonon.get_force_constants()

    if fc is not None:
        # convert forces from (N, N, 3, 3) to (3*N, 3*N)
        force_constants = phonon.get_force_constants().swapaxes(1, 2).reshape(2*(3*n_atoms, ))
        return force_constants
    else:
        print('**Force constants not yet created, please specify force_sets.')
        return None

def postprocess_init(phonon,
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
            freq_max=25,
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

    postprocess_init(phonon, force_sets)

    phonon.set_mesh(q_mesh)
    phonon.set_total_DOS(freq_min=freq_min,
                         freq_max=freq_max,
                         freq_pitch=freq_pitch,
                         tetrahedron_method=tetrahedron_method)

    if write:
        phonon.write_total_DOS()
        Path('total_dos.dat').rename(filename)

    return phonon.get_total_DOS()

def get_bandstructure(phonon,
                      path,
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

    postprocess_init(phonon, force_sets)

    phonon.set_band_structure(path)

    return phonon.get_band_structure()
