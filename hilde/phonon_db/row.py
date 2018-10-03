from random import randint
import numpy as np
from ase.db.row import AtomsRow, atoms2dict
from ase.atoms import Atoms
from phonopy import Phonopy
from hilde.materials_fp.MaterialsFingerprints import get_phonon_bs_fingerprint_phononpy, get_phonon_dos_fingerprint_phononpy
from hilde.structure.structure import pAtoms
from hilde.structure.convert import to_phonopy_atoms
def phonopy_to_ASE_atoms(structure):
    ASE_atoms= Atoms(
        symbols   = structure.get_chemical_symbols(),
        cell      = structure.get_cell(),
        masses    = structure.get_masses(),
        positions = structure.get_positions(),
        pbc       = True
    )
    return pAtoms(ase_atoms=ASE_atoms)

# def ASE_to_phonopy_atoms(structure):
#     phonopy_atoms= PhonopyAtoms(
#         symbols   = structure.get_chemical_symbols(),
#         cell      = structure.get_cell(),
#         masses    = structure.get_masses(),
#         positions = structure.get_positions(wrap=True))
#     return phonopy_atoms

def phonon2dict(phonon):
    dct = atoms2dict(phonopy_to_ASE_atoms(phonon.get_primitive()))
    if phonon.get_supercell_matrix() is not None:
        dct['supercell_matrix'] = list(phonon.get_supercell_matrix().flatten())
    if phonon.get_force_constants() is not None:
        dct['force_constants'] = phonon.get_force_constants()
    if phonon._total_dos is not None:
        dct['phonon_dos_fp'] = get_phonon_dos_fingerprint_phononpy(phonon)
    if phonon.band_structure is not None:
        dct['qpoints'] = {}
        for ii in range( len(phonon.band_structure.qpoints) ):
            if list(phonon.band_structure.qpoints[ii][ 0]) not in dct['qpoints'].values():
                dct['qpoints'][phonon.band_structure.distances[ii][ 0]] = list(phonon.band_structure.qpoints[ii][ 0])
            if list(phonon.band_structure.qpoints[ii][-1]) not in dct['qpoints'].values():
                dct['qpoints'][phonon.band_structure.distances[ii][-1]] = list(phonon.band_structure.qpoints[ii][-1])
        dct['phonon_bs_fp'] = get_phonon_bs_fingerprint_phononpy(phonon, dct['qpoints'])
    return dct

class PhononRow(AtomsRow):
    def __init__(self, dct):
        if isinstance(dct, dict):
            dct = dct.copy()
        else:
            dct = phonon2dict(dct)
        assert 'numbers' in dct
        assert 'cell' in dct
        self._constraints = dct.pop('constraints', [])
        self._constrained_forces = None
        self._data = dct.pop('data', {})
        kvp = dct.pop('key_value_pairs', {}) # If the dictionary has additional keys that are not default add them here
        self._keys = list(kvp.keys())
        self.__dict__.update(kvp)
        self.__dict__.update(dct)

    def to_phonon(self, attach_calculator=False, add_additional_information=False):
        phonon = Phonopy(to_phonopy_atoms(pAtoms(ase_atoms=self.toatoms())),
                supercell_matrix = np.array(self.supercell_matrix).reshape(3,3),
                symprec          = 1e-5,
                is_symmetry      = True,
                factor           = 15.633302,
                log_level        = 0)
        if("force_constants" in self):
            if(len(self.force_constants.shape) < 4):
                self.force_constants.reshape( (int(np.sqrt(len(self.force_constants))/3), int(np.sqrt(len(self.force_constants))/3),3,3) )
            phonon.set_force_constants(self.force_constants)
        return phonon
