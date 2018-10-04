from random import randint
import numpy as np
from ase.db.row import AtomsRow, atoms2dict
from ase.atoms import Atoms
from phonopy import Phonopy
from hilde.materials_fp.MaterialsFingerprints import get_phonon_bs_fingerprint_phononpy, get_phonon_dos_fingerprint_phononpy
from hilde.structure.structure import pAtoms
from hilde.structure.convert import to_phonopy_atoms

def phonon2dict(phonon):
    dct = atoms2dict(pAtoms(phonopy_atoms=phonon.get_primitive()))
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
    if phonon.thermal_properties is not None:
        dct['thermal_prop_ZPE'] = phonon.thermal_properties.zero_point_energy
        dct['thermal_prop_high_T_S'] = phonon.thermal_properties.high_T_entropy
        dct['thermal_prop_T' ], dct['thermal_prop_A' ], dct['thermal_prop_S' ], dct['thermal_prop_Cv'] = phonon.get_thermal_properties()
    return dct

class PhononRow(AtomsRow):
    def __init__(self, dct):
        if isinstance(dct, dict):
            dct = dct.copy()
            if type(dct['supercell_matrix']) is not list and type(dct['supercell_matrix']) is not str:
                dct['supercell_matrix'] = list(dct['supercell_matrix'].flatten())
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
        print(type(self.supercell_matrix))
        phonon = Phonopy(to_phonopy_atoms(pAtoms(ase_atoms=self.toatoms())),
                supercell_matrix = np.array(self.supercell_matrix).reshape(3,3),
                symprec          = 1e-5,
                is_symmetry      = True,
                factor           = 15.633302,
                log_level        = 0)
        if "force_constants" in self:
            if(len(self.force_constants.shape) < 4):
                self.force_constants.reshape( (int(np.sqrt(len(self.force_constants))/3), int(np.sqrt(len(self.force_constants))/3),3,3) )
            phonon.set_force_constants(self.force_constants)
        if "thermal_prop_T" in self:
            phonon.set_thermal_properties(temperatures=self.thermal_prop_T)
        return phonon

    def thermal_heat_capacity_v(self, T):
        return self.thermal_prop_Cv[np.where( self.thermal_prop_T == T)[0] ][0]

    def thermal_entropy(self, T):
        return self.thermal_prop_S[np.where( self.thermal_prop_T == T)[0] ][0]

    def thermal_free_energy(self, T):
        return self.thermal_prop_A[np.where( self.thermal_prop_T == T)[0] ][0]

