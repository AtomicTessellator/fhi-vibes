"""
A file defining a row object for the phonon database
"""
import numpy as np
from ase.db.row import AtomsRow, atoms2dict
from ase.io.jsonio import decode
from phonopy import Phonopy
from phono3py.phonon3 import Phono3py

from hilde.helpers.utility_functions import reshape_fc_2, reshape_fc_3
from hilde import konstanten as const
from hilde.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy
from hilde.materials_fp.material_fingerprint import get_phonon_dos_fingerprint_phononpy
from hilde.materials_fp.material_fingerprint import to_dict
from hilde.structure.structure import pAtoms
from hilde.structure.convert import to_phonopy_atoms


def phonon_to_dict(phonon, to_mongo=False):
    """
    Converts a phonopy object to a dictionary
    Args:
        phonon: the phonopy object to be converted
        to_mongo: True if it is being sent to a mongo database
    Returns:
        dct: the dictionary representation of phonon
    """
    dct = atoms2dict(pAtoms(phonopy_atoms=phonon.get_primitive()))
    dct["symprec"] = phonon._symprec
    if phonon.get_supercell_matrix() is not None:
        dct["sc_matrix_2"] = list(np.array(phonon.get_supercell_matrix()).flatten())
        try:
            dct["natoms_in_sc_2"] = len(phonon.get_supercell().symbols)
        except:
            dct["natoms_in_sc_2"] = len(dct["numbers"])
    if phonon.get_force_constants() is not None:
        dct["fc_2"] = phonon.get_force_constants()
    if phonon.mesh is not None:
        dct["qmesh"] = phonon.mesh.mesh_numbers
    if phonon._total_dos is not None:
        dct["phonon_dos_fp"] = to_dict(get_phonon_dos_fingerprint_phononpy(phonon))

    if phonon.band_structure is not None:
        dct["qpoints"] = {}
        for ii, q_pt in enumerate(phonon.band_structure.qpoints):
            if list(q_pt[0]) not in dct["qpoints"].values():
                dct["qpoints"][phonon.band_structure.distances[ii][0]] = list(q_pt[0])
            if list(q_pt[-1]) not in dct["qpoints"].values():
                dct["qpoints"][phonon.band_structure.distances[ii][-1]] = list(q_pt[-1])
        dct["phonon_bs_fp"] = to_dict(
            get_phonon_bs_fingerprint_phononpy(phonon, dct["qpoints"]), to_mongo
        )
        if to_mongo:
            q_pt = {}
            for pt in dct["qpoints"]:
                q_pt[str(pt)] = dct["qpoints"][pt]
            dct["qpoints"] = q_pt
    if phonon.thermal_properties is not None:
        dct["tp_ZPE"] = phonon.thermal_properties.zero_point_energy
        dct["tp_high_T_S"] = phonon.thermal_properties.high_T_entropy
        dct["tp_T"], dct["tp_A"], dct["tp_S"], dct["tp_Cv"] = phonon.get_thermal_properties()
    return dct


def phonon3_to_dict(phonon3, phonon=None, to_mongo=False):
    """
    Converts a phonopy object to a dictionary
    Args:
        phonon: the phonopy object to be converted
        to_mongo: True if it is being sent to a mongo database
    Returns:
        dct: the dictionary representation of phonon
    """
    dct = atoms2dict(pAtoms(phonopy_atoms=phonon3.get_primitive()))
    dct["symprec"] = phonon3._symprec
    if phonon is None:
        phonon = Phonopy(
            phonon3.get_unitcell(),
            supercell_matrix=phonon3._phonon_supercell_matrix,
            symprec=phonon3._symprec,
            is_symmetry=phonon3._is_symmetry,
            factor=phonon3._frequency_factor_to_THz,
            log_level=phonon3._log_level,
        )
        if phonon3.get_fc2():
            phonon.set_force_constants(phonon3.get_fc2())
    if phonon3.get_supercell_matrix() is not None:
        dct["sc_matrix_3"] = list(np.array(phonon3.get_supercell_matrix()).flatten())
        try:
            dct["natoms_in_sc_3"] = len(phonon3.get_supercell().symbols)
        except:
            dct["natoms_in_sc_3"] = len(dct["numbers"])
    if phonon3.get_fc3() is not None:
        dct["fc_3"] = phonon3.get_fc3()

    if phonon3.get_thermal_conductivity():
        dct["tp_T"] = phonon3.get_thermal_conductivity().get_temperatures()
        dct["tp_kappa"] = phonon3.get_thermal_conductivity().get_kappa()
        dct["qmesh"] = phonon3.get_thermal_conductivity().get_mesh_numbers()
        phonon.set_mesh(dct["qmesh"])
        phonon.set_thermal_properties(temperatures=dct["tp_T"])

    for key, val in phonon_to_dict(phonon).items():
        dct[key] = val
    return dct


class PhononRow(AtomsRow):
    """
    Class that is largely based off of the ASE AtomsRow object but expanded for phonopy
    """

    def __init__(self, dct=None, phonon3=None, phonon=None):
        """
        Constructor for the PhononRow.
        Args:
            dct: a phonopy object or a dict
                representation of the phonopy object to be added to the database
        """
        if dct:
            dct = dct.copy()
            if "sc_matrix_2" in dct and isinstance(dct["sc_matrix_2"], np.ndarray):
                dct["sc_matrix_2"] = list(dct["sc_matrix_2"].flatten())
            elif (
                "sc_matrix_2" in dct
                and isinstance(dct["sc_matrix_2"], list)
                and isinstance(dct["sc_matrix_2"][0], str)
            ):
                for i, sc_el in enumerate(dct["sc_matrix_2"]):
                    dct["sc_matrix_2"][i] = int(sc_el)
            if "sc_matrix_3" in dct and isinstance(dct["sc_matrix_3"], np.ndarray):
                dct["sc_matrix_3"] = list(dct["sc_matrix_3"].flatten())
            elif (
                "sc_matrix_3" in dct
                and isinstance(dct["sc_matrix_3"], list)
                and isinstance(dct["sc_matrix_3"][0], str)
            ):
                for i, sc_el in enumerate(dct["sc_matrix_3"]):
                    dct["sc_matrix_3"][i] = int(sc_el)
        elif phonon3:
            dct = phonon3_to_dict(phonon3, phonon)
        elif phonon:
            dct = phonon_to_dict(phonon)
        else:
            raise AttributeError("dct, phonon3 or phonon must be defined")
        assert "numbers" in dct
        assert "cell" in dct
        self._constraints = dct.pop("constraints", [])
        self._constrained_forces = None
        self._data = dct.pop("data", {})
        # If the dictionary has additional keys that are not default add them here
        kvp = dct.pop("key_value_pairs", {})
        self._keys = list(kvp.keys())
        self.__dict__.update(kvp)
        self.__dict__.update(dct)

    def to_phonon(self):
        """
        Converts the row back into a phonopy object
        Returns:
            phonon: The phonopy object the PhononRow represents
        """
        phonon = Phonopy(
            to_phonopy_atoms(pAtoms(ase_atoms=self.toatoms())),
            supercell_matrix=np.array(self.sc_matrix_2).reshape(3, 3),
            symprec=self.symprec,
            is_symmetry=True,
            factor=15.633_302,
            log_level=0,
        )
        if "fc_2" in self:
            self.fc_2 = np.asarray(self.fc_2)
            if len(self.fc_2.shape) < 4:
                self.fc_2 = reshape_fc_2(self.fc_2)
            phonon.set_force_constants(self.fc_2)
        if "qmesh" in self and self.qmesh is not None:
            phonon.set_mesh(self.qmesh)
            if "tp_T" in self and self.tp_T is not None:
                phonon.set_thermal_properties(temperatures=self.tp_T)
        return phonon

    def to_phonon3(self, mesh=None):
        """
        Converts the row back into a phono3py object
        Returns:
            phonon3: The phono3py object the PhononRow represents
        """
        if mesh is None and "qmesh" in self and self.qmesh is not None:
            phonon3 = Phono3py(
                to_phonopy_atoms(pAtoms(ase_atoms=self.toatoms())),
                supercell_matrix=np.array(self.sc_matrix_3).reshape(3, 3),
                phonon_supercell_matrix=np.array(self.sc_matrix_2).reshape(3, 3),
                symprec=self.symprec,
                is_symmetry=True,
                frequency_factor_to_THz=const.eV_to_THz,
                log_level=0,
                mesh=self.qmesh,
            )
        elif mesh is None:
            phonon3 = Phono3py(
                to_phonopy_atoms(pAtoms(ase_atoms=self.toatoms())),
                supercell_matrix=np.array(self.sc_matrix_3).reshape(3, 3),
                phonon_supercell_matrix=np.array(self.sc_matrix_2).reshape(3, 3),
                symprec=self.symprec,
                is_symmetry=True,
                frequency_factor_to_THz=const.eV_to_THz,
                log_level=0,
            )
        else:
            phonon3 = Phono3py(
                to_phonopy_atoms(pAtoms(ase_atoms=self.toatoms())),
                supercell_matrix=np.array(self.sc_matrix_3).reshape(3, 3),
                phonon_supercell_matrix=np.array(self.sc_matrix_2).reshape(3, 3),
                symprec=self.symprec,
                is_symmetry=True,
                frequency_factor_to_THz=const.eV_to_THz,
                log_level=0,
                mesh=mesh,
            )
        if "fc_2" in self:
            self.fc_2 = np.asarray(self.fc_2)
            if len(self.fc_2.shape) < 4:
                self.fc_2 = reshape_fc_2(self.fc_2)
            phonon3.set_fc2(self.fc_2)
        if "fc_3" in self:
            self.fc_3 = np.asarray(self.fc_3)
            if len(self.fc_3.shape) < 6:
                self.fc_3 = reshape_fc_3(self.fc_3)
            phonon3.set_fc3(self.fc_3)
        return phonon3

    def thermal_heat_capacity_v(self, T):
        """
        Gets the Cv of the material at a given temperature
        Args:
            T: float
                The temperature
        Returns:
            Cv : float
                the heat_capacity_v at temperature T
        """
        return self.tp_Cv[np.where(self.tp_T == T)[0]][0]

    def thermal_entropy(self, T):
        """
        Gets the entropy of the material at a given temperature
        Args:
            T: float
                The temperature
        Returns:
            S: float
                The entropy at temperature T
        """
        return self.tp_S[np.where(self.tp_T == T)[0]][0]

    def thermal_free_energy(self, T):
        """
        Gets the Hemholtz free energy of the material at a given temperature
        Args:
            T: float
                The temperature
        Returns:
            A: float
                The Hemholtz free energy at temperature T
        """
        return self.tp_A[np.where(self.tp_T == T)[0]][0]

    def thermal_conductivity(self, T):
        """
        Gets the thermal conductivity of the material at a given temperature
        Args:
            T: float
                The temperature
        Returns:
            A: float
                The Hemholtz free energy at temperature T
        """
        return self.tp_kappa[np.where(self.tp_T == T)[0]][0]
