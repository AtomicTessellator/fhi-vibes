"""
A file defining a row object for the phonon database
"""
import numpy as np
from ase.db.row import AtomsRow, atoms2dict
from phonopy import Phonopy

from hilde.helpers.utility_functions import reshape_fc
from hilde.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy
from hilde.materials_fp.material_fingerprint import get_phonon_dos_fingerprint_phononpy
from hilde.materials_fp.material_fingerprint import to_dict
from hilde.structure.structure import pAtoms
from hilde.structure.convert import to_phonopy_atoms

def phonon2dict(phonon, to_mongo=False):
    '''
    Converts a phonopy object to a dictionary
    Args:
        phonon: the phonopy object to be converted
        to_mongo: True if it is being sent to a mongo database
    Returns:
        dct: the dictionary representation of phonon
    '''
    dct = atoms2dict(pAtoms(phonopy_atoms=phonon.get_primitive()))
    if phonon.get_supercell_matrix() is not None:
        dct['supercell_matrix'] = list(phonon.get_supercell_matrix().flatten())
        try:
            dct['natoms_in_sc'] = len(phonon.get_supercell().symbols)
        except:
            dct['natoms_in_sc'] = len(dct['numbers'])
    if phonon.get_force_constants() is not None:
        dct['force_constants'] = phonon.get_force_constants()
    if phonon.mesh is not None:
        dct['qmesh'] = phonon.mesh.mesh_numbers
    if phonon._total_dos is not None:
        dct['phonon_dos_fp'] = to_dict(get_phonon_dos_fingerprint_phononpy(phonon))
    if phonon.band_structure is not None:
        dct['qpoints'] = {}
        for ii, q_pt in enumerate(phonon.band_structure.qpoints):
            if list(q_pt[0]) not in dct['qpoints'].values():
                dct['qpoints'][phonon.band_structure.distances[ii][0]] = list(q_pt[0])
            if list(q_pt[-1]) not in dct['qpoints'].values():
                dct['qpoints'][phonon.band_structure.distances[ii][-1]] = list(q_pt[-1])
        dct['phonon_bs_fp'] = to_dict(get_phonon_bs_fingerprint_phononpy(phonon, dct['qpoints']), to_mongo)
        if to_mongo:
            q_pt = {}
            for pt in dct['qpoints']:
                q_pt[str(pt)] = dct['qpoints'][pt]
            dct['qpoints'] = q_pt
    if phonon.thermal_properties is not None:
        dct['tp_ZPE'] = phonon.thermal_properties.zero_point_energy
        dct['tp_high_T_S'] = phonon.thermal_properties.high_T_entropy
        dct['tp_T'], dct['tp_A'], dct['tp_S'], dct['tp_Cv'] = phonon.get_thermal_properties()
    return dct

class PhononRow(AtomsRow):
    '''
    Class that is largely based off of the ASE AtomsRow object but expanded for phonopy
    '''
    def __init__(self, dct):
        '''
        Constructor for the PhononRow.
        Args:
            dct: a phonopy object or a dict
                representation of the phonopy object to be added to the database
        '''
        if isinstance(dct, dict):
            dct = dct.copy()
            if 'supercell_matrix' in dct and not isinstance(dct['supercell_matrix'], list) and not isinstance(dct['supercell_matrix'], str):
                dct['supercell_matrix'] = list(dct['supercell_matrix'].flatten())
            elif 'supercell_matrix' in dct and isinstance(dct['supercell_matrix'], list) and isinstance(dct['supercell_matrix'][0], str):
                for i, sc_el in enumerate(dct['supercell_matrix']):
                    dct['supercell_matrix'][i] = int(sc_el)
        else:
            dct = phonon2dict(dct)
        assert 'numbers' in dct
        assert 'cell' in dct
        self._constraints = dct.pop('constraints', [])
        self._constrained_forces = None
        self._data = dct.pop('data', {})
        # If the dictionary has additional keys that are not default add them here
        kvp = dct.pop('key_value_pairs', {})
        self._keys = list(kvp.keys())
        self.__dict__.update(kvp)
        self.__dict__.update(dct)

    def to_phonon(self):
        '''
        Converts the row back into a phonopy object
        Returns:
            phonon: The phonopy object the PhononRow represents
        '''
        phonon = Phonopy(to_phonopy_atoms(pAtoms(ase_atoms=self.toatoms())),
                         supercell_matrix=np.array(self.supercell_matrix).reshape(3, 3),
                         symprec=1e-5,
                         is_symmetry=True,
                         factor=15.633302,
                         log_level=0
                        )
        if "force_constants" in self:
            self.force_constants = np.asarray(self.force_constants)
            if len(self.force_constants.shape) < 4:
                self.force_constants = reshape_fc(self.force_constants)
            phonon.set_force_constants(self.force_constants)
        if "qmesh" in self and self.qmesh is not None:
            phonon.set_mesh(self.qmesh)
        if "tp_T" in self and self.tp_T is not None:
            phonon.set_thermal_properties(temperatures=self.tp_T)
        return phonon

    def thermal_heat_capacity_v(self, T):
        '''
        Gets the Cv of the material at a given temperature
        Args:
            T: float
                The temperature
        Returns:
            Cv : float
                the heat_capacity_v at temperature T
        '''
        return self.tp_Cv[np.where(self.tp_T == T)[0]][0]

    def thermal_entropy(self, T):
        '''
        Gets the entropy of the material at a given temperature
        Args:
            T: float
                The temperature
        Returns:
            S: float
                The entropy at temperature T
        '''
        return self.tp_S[np.where(self.tp_T == T)[0]][0]

    def thermal_free_energy(self, T):
        '''
        Gets the Hemholtz free energy of the material at a given temperature
        Args:
            T: float
                The temperature
        Returns:
            A: float
                The Hemholtz free energy at temperature T
        '''
        return self.tp_A[np.where(self.tp_T == T)[0]][0]
