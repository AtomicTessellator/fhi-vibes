"""
A file defining a row object for the phonon database
"""
import copy

import numpy as np

from ase.db.row import AtomsRow

from vibes.spglib.wrapper import get_spacegroup

from vibes import konstanten as const
from vibes.ase.db.dict_converters import atoms2dict

# from vibes.phonon_db.ase_converters import dict2atoms
from vibes.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy
from vibes.materials_fp.material_fingerprint import get_phonon_dos_fingerprint_phononpy
from vibes.materials_fp.material_fingerprint import to_dict
from vibes.structure.convert import to_Atoms, to_phonopy_atoms


def standardize_sc_matrix(sc_matrix):
    """convert sc_matrix into a list of ints"""
    return list(np.array(sc_matrix, dtype=int).T.flatten())


def phonon_to_dict(phonon, to_mongo=False, add_fc=False):
    """Converts a phonopy object to a dictionary

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The Phonopy object to be converted
    to_mongo: bool
        If True then it is being sent to a mongo database

    Returns
    -------
    dct: dict
        the dictionary representation of phonon
    """
    dct = atoms2dict(to_Atoms(phonon.get_primitive()))
    dct["key_value_pairs"] = {"symprec": phonon._symprec}
    if phonon.get_supercell_matrix() is not None:
        dct["sc_matrix_2"] = standardize_sc_matrix(phonon.get_supercell_matrix())
        dct["natoms_in_sc_2"] = len(dct["numbers"]) * int(
            round(np.linalg.det(phonon.get_supercell_matrix()))
        )
    if add_fc:
        dct["_fc_2"] = np.array(phonon.get_force_constants())
    else:
        displacement_dataset = copy.deepcopy(phonon._displacement_dataset)
        dct["force_2"] = []
        for disp1 in displacement_dataset["first_atoms"]:
            if "forces" in disp1:
                dct["force_2"].append(disp1.pop("forces"))

        dct["force_2"] = np.array(dct["force_2"])

    dct["displacement_dataset_2"] = displacement_dataset

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
                q_pt[str(pt).replace(".", "_")] = dct["qpoints"][pt]
            dct["qpoints"] = q_pt
    if phonon.thermal_properties is not None:
        dct["tp_ZPE"] = phonon.thermal_properties.zero_point_energy
        dct["tp_high_T_S"] = phonon.thermal_properties.high_T_entropy
        dct["tp_T"], dct["tp_A"], dct["tp_S"], dct[
            "tp_Cv"
        ] = phonon.get_thermal_properties()

    return dct


def phonon3_to_dict(phonon3, store_second_order=False, to_mongo=False):
    """Converts a phonopy object to a dictionary

    Parameters
    ----------
    phonon3: phono3py.phonon3.Phono3py
        The Phono3py object to be converted
    store_second_order: bool
        If True store the second order properties of the phonopy object
    to_mongo: bool
        If True then it is being sent to a mongo database

    Returns
    -------
    dct: dict
        the dictionary representation of phonon3
    """
    dct = atoms2dict(to_Atoms(phonon3.get_primitive()))
    dct["symprec"] = phonon3._symprec

    if phonon3.get_supercell_matrix() is not None:
        dct["sc_matrix_3"] = standardize_sc_matrix(phonon3.get_supercell_matrix())
        dct["natoms_in_sc_3"] = len(dct["numbers"]) * int(
            round(np.linalg.det(phonon3.get_supercell_matrix()))
        )

    if store_second_order and phonon3._phonon_supercell_matrix is not None:
        dct["sc_matrix_2"] = standardize_sc_matrix(phonon3._phonon_supercell_matrix)
        dct["natoms_in_sc_2"] = len(dct["numbers"]) * int(
            round(np.linalg.det(phonon3._phonon_supercell_matrix))
        )

    dct["force_3"] = []
    get_forces = True
    displacement_dataset = copy.deepcopy(phonon3._displacement_dataset)
    for disp1 in displacement_dataset["first_atoms"]:
        if "forces" in disp1:
            if disp1["forces"].shape[0] == dct["natoms_in_sc_3"]:
                dct["force_3"].append(disp1.pop("forces"))
            else:
                get_forces = False
                break

    if get_forces:
        for ii, disp1 in enumerate(displacement_dataset["first_atoms"]):
            for disp2 in disp1["second_atoms"]:
                if "delta_forces" in disp2:
                    dct["force_3"].append(
                        disp2.pop("delta_forces") + dct["force_3"][ii]
                    )
        dct["force_3"] = np.array(dct["force_3"])
    else:
        print("Warning not storing forces because of an inconsistent number of atoms")
        dct["force_3"] = np.ndarray(0)

    dct["force_3"] = np.array(dct["force_3"])
    dct["displacement_dataset_3"] = displacement_dataset

    if phonon3.get_thermal_conductivity():
        dct["tp_T"] = phonon3.get_thermal_conductivity().get_temperatures()
        dct["tp_kappa"] = phonon3.get_thermal_conductivity().get_kappa()[0]
        dct["qmesh"] = phonon3.get_thermal_conductivity().get_mesh_numbers()

    return dct


class PhononRow(AtomsRow):
    """ASE AtomsRow object but expanded for phonopy"""

    def __init__(self, dct, store_second_order=False):
        """Constructor for the PhononRow.

        Parameters
        ----------
        dct: dict
            A dictionary representation of the PhononRow
        Raises
        ------
        AttributeError
            If dct, phonon3, and phonon are all None
        AssertionError
            If dct does not have numbers OR
            If dct does not have cell
        """
        from phonopy import Phonopy
        from ase.atoms import Atoms

        if isinstance(dct, dict):
            if "sc_matrix_2" in dct:
                dct["sc_matrix_2"] = standardize_sc_matrix(dct["sc_matrix_2"])
            if dct.get("sc_matrix_3", None):
                dct["sc_matrix_3"] = standardize_sc_matrix(dct["sc_matrix_3"])
            if "force_2" in dct:
                dct["key_value_pairs"]["has_fc2"] = True
            if "force_3" in dct:
                dct["key_value_pairs"]["has_fc3"] = True
        elif isinstance(dct, Phonopy):
            dct = phonon_to_dict(dct)
            dct["key_value_pairs"]["has_fc2"] = True
        elif isinstance(dct, Atoms):
            dct = atoms2dict(dct)
            dct["key_value_pairs"] = {}
        else:
            try:
                from phono3py.phonon3 import Phono3py

                if isinstance(dct, Phono3py):
                    dct = phonon3_to_dict(dct)
                    dct["key_value_pairs"]["has_fc3"] = True
                    if store_second_order:
                        dct["key_value_pairs"]["has_fc2"] = True
            except ModuleNotFoundError:
                pass

        atoms = AtomsRow(dct).toatoms()
        sg = get_spacegroup(atoms)
        if sg is None:
            dct["key_value_pairs"]["space_group"] = -1
        else:
            dct["key_value_pairs"]["space_group"] = int(sg.split("(")[-1].split(")")[0])

        super(PhononRow, self).__init__(dct)

        self.clean_displacement_dataset()

    def clean_displacement_dataset(self):
        """Cleans the displacement dataset

        Storing in the database can convert some numbers to strings, fix this problem
        """
        if "displacement_dataset_2" in self and self.displacement_dataset_2 is not None:
            for disp1 in self.displacement_dataset_2["first_atoms"]:
                disp1["number"] = int(disp1["number"])
        if "displacement_dataset_3" in self and self.displacement_dataset_3 is not None:
            for disp1 in self.displacement_dataset_3["first_atoms"]:
                disp1["number"] = int(disp1["number"])
                for disp2 in disp1["second_atoms"]:
                    disp2["number"] = int(disp2["number"])

    @property
    def fc_2(self):
        """The second order force constants"""
        if "_fc_2" in self:
            return self._fc_2
        else:
            self.__dict__["_fc_2"] = self.get_fc_2()
            return self._fc_2

    def get_fc_2(self):
        """Calculate the second order force constants from data stored in the row

        Returns
        -------
        np.ndarray
            Second order force constants
        """
        from phonopy import Phonopy

        phonon = Phonopy(
            to_phonopy_atoms(self.toatoms()),
            supercell_matrix=np.array(self.sc_matrix_2).reshape(3, 3),
            symprec=self.get("symprec", 1e-5),
            is_symmetry=True,
            factor=const.omega_to_THz,
            log_level=0,
        )
        phonon.set_displacement_dataset(self.displacement_dataset_2)
        if "force_2" in self and len(self.force_2) > 0:
            phonon.produce_force_constants(
                self.force_2, calculate_full_force_constants=False
            )
        return phonon.get_force_constants()

    @property
    def fc_3(self):
        """The second order force constants"""
        if "_fc_3" in self:
            return self._fc_3
        else:
            self.__dict__["_fc_3"] = self.get_fc_3()
            return self._fc_3

    def get_fc_3(self):
        """Calculate the third order force constants from data stored in the row

        Returns
        -------
        np.ndarray
            Third order force constants
        """
        try:
            from phono3py.phonon3 import Phono3py
        except ModuleNotFoundError:
            return None

        phonon3 = Phono3py(
            to_phonopy_atoms(self.toatoms()),
            supercell_matrix=np.array(self.sc_matrix_3).reshape(3, 3),
            phonon_supercell_matrix=np.array(self.sc_matrix_2).reshape(3, 3),
            symprec=self.symprec,
            is_symmetry=True,
            frequency_factor_to_THz=const.omega_to_THz,
            log_level=0,
            mesh=self.qmesh,
        )
        phonon3._phonon_displacement_dataset = self.displacement_dataset_2.copy()
        phonon3.set_displacement_dataset(self.displacement_dataset_3)

        if "force_3" in self and len(self.force_3) > 0:
            phonon3.produce_fc3(self.force_3)
        return phonon3.get_fc3()

    @property
    def spacegroup(self):
        if self.key_value_pairs["space_group"] > 0:
            return self.key_value_pairs["space_group"]
        else:
            return None

    def to_phonon(self):
        """Converts the row back into a phonopy object

        Returns
        -------
        phonon: phonopy.Phonopy
            The phonopy object the PhononRow represents
        """
        from phonopy import Phonopy

        phonon = Phonopy(
            to_phonopy_atoms(self.toatoms()),
            supercell_matrix=np.array(self.sc_matrix_2).reshape(3, 3).transpose(),
            symprec=self.symprec,
            is_symmetry=True,
            factor=const.omega_to_THz,
        )
        phonon.set_displacement_dataset(self.displacement_dataset_2)
        phonon.set_force_constants(self.fc_2)

        if "qmesh" in self and self.qmesh is not None:
            phonon.set_mesh(self.qmesh)
            if "tp_T" in self and self.tp_T is not None:
                phonon.set_thermal_properties(temperatures=self.tp_T)
        return phonon

    def to_phonon3(self, mesh=None):
        """Converts the row back into a phono3py object

        Returns
        -------
        phonon3: Phonoepy Object
            The phono3py object the PhononRow represents
        """
        from phono3py.phonon3 import Phono3py

        phonon3 = Phono3py(
            to_phonopy_atoms(self.toatoms()),
            supercell_matrix=np.array(self.sc_matrix_3).reshape(3, 3).transpose(),
            phonon_supercell_matrix=np.array(self.sc_matrix_2)
            .reshape(3, 3)
            .transpose(),
            symprec=self.symprec,
            is_symmetry=True,
            frequency_factor_to_THz=const.omega_to_THz,
            log_level=0,
            mesh=self.qmesh if "qmesh" in self else None,
        )
        phonon3._phonon_displacement_dataset = self.displacement_dataset_2.copy()
        phonon3.set_displacement_dataset(self.displacement_dataset_3)

        if "force_2" in self and len(self.force_2) > 0:
            phonon3.produce_fc2(
                self.force_2, displacement_dataset=self.displacement_dataset_2
            )
            self.__dict__["_fc_2"] = phonon3.get_fc2()

        if "force_3" in self and len(self.force_3) > 0:
            phonon3.produce_fc3(
                self.force_3, displacement_dataset=self.displacement_dataset_3
            )
            self.__dict__["_fc_3"] = phonon3.get_fc3()

        if mesh is None and "qmesh" in self and self.qmesh is not None:
            phonon3._mesh = np.array(self.qmesh, dtype="intc")
        elif mesh is not None:
            phonon3._mesh = np.array(mesh, dtype="intc")

        if "tp_T" in self:
            phonon3.run_thermal_conductivity(temperatures=self.tp_T, write_kappa=True)
        return phonon3

    def thermal_heat_capacity_v(self, T):
        """Gets the Cv of the material at a given temperature

        Parameters
        ----------
        T: float
            The temperature

        Returns
        -------
        Cv : float
            the heat_capacity_v at temperature T
        """
        return self.tp_Cv[np.where(self.tp_T == T)[0]][0]

    def thermal_entropy(self, T):
        """Gets the entropy of the material at a given temperature

        Parameters
        ----------
        T: float
            The temperature

        Returns
        -------
        S: float
            The entropy at temperature T
        """
        return self.tp_S[np.where(self.tp_T == T)[0]][0]

    def thermal_free_energy(self, T):
        """Gets the Hemholtz free energy of the material at a given temperature

        Parameters
        ----------
        T: float
            The temperature

        Returns
        -------
        A: float
            The Hemholtz free energy at temperature T
        """
        return self.tp_A[np.where(self.tp_T == T)[0]][0]

    def thermal_conductivity(self, T):
        """Gets the thermal conductivity of the material at a given temperature

        Parameters
        ----------
        T: float
            The temperature

        Returns
        -------
        A: float
            The Hemholtz free energy at temperature T
        """
        return self.tp_kappa[np.where(self.tp_T == T)[0]][0]
