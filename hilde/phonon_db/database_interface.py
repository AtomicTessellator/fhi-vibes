from ase.atoms import Atoms
from ase.symbols import symbols2numbers

from phonopy import Phonopy
from phono3py.phonon3 import Phono3py

import numpy as np

from hilde import konstanten as const
from hilde.helpers.converters import atoms2dict, dict2atoms, dict2results
from hilde.helpers.hash import hash_atoms, hash_atoms_and_calc, hash_traj
from hilde.helpers.warnings import warn
from hilde.phonon_db.phonon_db import connect
from hilde.phonon_db.row import phonon_to_dict
from hilde.phonon_db.row import phonon3_to_dict
from hilde.phono3py.wrapper import produce_fc3
from hilde.structure.convert import to_Atoms_db, to_Atoms, to_phonopy_atoms

from hilde.trajectory.reader import reader as traj_reader

results_keys = [
    "force_2",
    "force_3",
    "displacement_dataset_2",
    "displacement_dataset_3",
    "qmesh",
    "phonon_dos_fp",
    "qpoints",
    "phonon_bs_fp",
    "tp_ZPE",
    "tp_high_T_S",
    "tp_T",
    "tp_A",
    "tp_S",
    "tp_Cv",
    "tp_kappa",
    "fmax",
    "smax",
]

def obj2dict(obj):
    if isinstance(obj, str):
        calculated_atoms, metadata = traj_reader(obj, True)
        if "Phonopy" in metadata:
            ph =  phonon = Phonopy(
                to_phonopy_atoms(dict2results(metadata["Phonopy"]["primitive"])),
                supercell_matrix=np.array(metadata["Phonopy"]["supercell_matrix"]).reshape(3,3),
                is_symmetry=True,
                symprec=metadata["Phonopy"]["symprec"],
                factor=const.omega_to_THz,
            )
            ph.set_displacement_dataset(metadata["Phonopy"]["displacement_dataset"])
            ph.set_forces([at.get_forces() for at in calculated_atoms])
            return phonon_to_dict(ph)
        elif "Phono3py" in metadata:
            ph3 = Phono3py(
                to_phonopy_atoms(dict2results(metadata["Phono3py"]["primitive"])),
                supercell_matrix=np.array(metadata["Phono3py"]["supercell_matrix"]).reshape(3,3),
                is_symmetry=True,
                frequency_factor_to_THz=const.omega_to_THz,
            )
            ph3 = produce_fc3(
                ph3,
                metadata,
                np.array([at.get_forces() for at in calculated_atoms]),
            )
            return phonon3_to_dict(ph3)
        else:
            raise TypeError("Trajectory file can't be added to database")
    elif isinstance(obj, dict):
        return obj.copy()
    elif isinstance(obj, Atoms):
        return atoms2dict(obj)
    elif isinstance(obj, Phonopy):
        return phonon_to_dict(obj)
    else:
        try:
            from phono3py.phonon3 import Phono3py
            if isinstance(obj, Phono3py):
                return phonon3_to_dict(obj)
        except:
            raise InputError("obj has to be a dict, Phonopy, Phono3py, or Atoms object")

def from_database(db_path, identifier=None, obj=None, selection=[], get_id=False, get_atoms=False, get_phonon=False, get_phonon3=False, **kwargs):
    db = connect(db_path)
    if identifier:
        selection = [(identifier[0], "=", identifier[1])]
    else:
        for key, val in kwargs.items():
            selection.append((key, "=", val))

        dct = obj2dict(obj) if obj is not None else {}

        for key, val in dct.items():
            if key not in results_keys:
                selection.append(key, "=", val)

    row = db.get(selection)

    to_ret = []
    if get_id:
        to_ret.append(row.id)
    if get_atoms:
        to_ret.append(row.toatoms())
    if get_phonon:
        to_ret.append(row.to_phonon())
    if get_phonon3:
        to_ret.append(row.to_phonon3())

    if len(to_ret) > 1:
        to_ret = tuple(to_ret)
    else:
        to_ret = to_ret[0]
    return to_ret

def to_database(db_path, obj, calc=None, key_val_pairs=None):
    selection_no_sc = []
    selection = []
    dct = obj2dict(obj)
    hashes = {}
    if not key_val_pairs:
        key_val_pairs = {}
    if isinstance(obj, str):
        calculated_atoms, metadata = traj_reader(obj, True)
        if "Phonopy" in metadata:
            at_dict = metadata["Phonopy"]['primitive'].copy()
            for key, val in metadata['calculator'].items():
                at_dict[key] = val
            at_dict["numbers"] = symbols2numbers(at_dict['symbols'])
            at_dict['pbc'] = [True, True, True]
            at_dict['calculator'] = at_dict['calculator'].lower()
            atoms = dict2atoms(at_dict)
            typ = "ph"
            selection.append(("sc_matrix_2", "=", dct["sc_matrix_2"]))
            key_val_pairs["displacement"] = metadata["Phonopy"]['displacement_dataset']['first_atoms'][0]['displacement'][0]
        elif "Phono3py" in metadata:
            at_dict = metadata["Phono3py"]['primitive'].copy()
            for key, val in metadata['calculator'].items():
                at_dict[key] = val
            at_dict["numbers"] = symbols2numbers(at_dict['symbols'])
            at_dict['pbc'] = [True, True, True]
            at_dict['calculator'] = at_dict['calculator'].lower()
            atoms = dict2atoms(at_dict)
            typ = "ph3"
            selection.append(("sc_matrix_3", "=", dct["sc_matrix_3"]))
            key_val_pairs["pair_distance_cutoff"] = metadata['displacement_dataset']['cutoff_distance']
            key_val_pairs["displacement"] = metadata["Phono3py"]['displacement_dataset']['first_atoms'][0]['displacement'][0]
        else:
            raise TypeError("Trajectory file can't be added to database")

        hashes['traj_hash'], hashes[typ + '_hash'] = hash_traj(calculated_atoms, metadata, True)
        hashes['cell_hash'] = hash_atoms_and_calc(atoms)[0]

        selection_no_sc.append(("cell_hash", '=', hashes["cell_hash"]))

        np.round(atoms.cell, 15)
        np.round(atoms.positions, 15)
        if "masses" in atoms.__dict__:
            del(atoms.__dict__["masses"])
    elif isinstance(obj, dict):
        atoms = dict2atoms(dct)
    elif isinstance(obj, Atoms):
        atoms = obj.copy()
    elif isinstance(obj, Phonopy):
        atoms = to_Atoms_db(obj.get_primitive())

        hashes["cell_hash"] = hash_atoms_and_calc(to_Atoms(obj.get_primitive()))[0]

        selection_no_sc.append(("cell_hash", '=', hashes["cell_hash"]))
        selection.append(("sc_matrix_2", "=", dct["sc_matrix_2"]))

        key_val_pairs["displacement"] = obj._displacement_dataset['first_atoms'][0]['displacement'][0]
    else:
        try:
            from phono3py.phonon3 import Phono3py
            if isinstance(obj, Phono3py):
                atoms = to_Atoms_db(obj.get_primitive())

                hashes["cell_hash"] = hash_atoms_and_calc(to_Atoms(obj.get_primitive()))[0]

                selection_no_sc.append(("cell_hash", '=', hashes["cell_hash"]))
                selection.append(("sc_matrix_3", "=", dct["sc_matrix_3"]))

                key_val_pairs["displacement"] = obj._displacement_dataset['first_atoms'][0]['displacement'][0]
                key_val_pairs["pair_distance_cutoff"] = obj._displacement_dataset['cutoff_distance']
        except:
            pass
    atoms.set_calculator(calc)
    if "calc_hash" in key_val_pairs or "atoms_hash" in key_val_pairs:
        print(key_val_pairs)
        warn("Replacing the atoms and calc hashes")
    hashes['atoms_hash'], hashes['calc_hash'] = hash_atoms_and_calc(atoms)

    for key, val in atoms2dict(atoms).items():
        dct[key] = val

    selection_no_sc.append(("atoms_hash", "=", hashes["atoms_hash"]))
    selection_no_sc.append(("calc_hash", "=", hashes["calc_hash"]))

    if "symprec" in key_val_pairs.keys():
        selection_no_sc.append(("symprec", "=", key_val_pairs["symprec"]))

    for sel in selection_no_sc:
        selection.append(sel)

    key_val_pairs.update(hashes)
    db = connect(db_path)
    try:
        try:
            rows = list(db.select(selection=selection))
            rows_no_sc = list(db.select(selection=selection_no_sc))
            for row in rows_no_sc:
                row_sc_mat_2 = False
                if "sc_matrix_2" not in row.__dict__ or row.sc_matrix_2 is None:
                    row_sc_mat_2 = True

                row_sc_mat_3 = False
                if "sc_matrix_3" not in row.__dict__ or row.sc_matrix_3 is None:
                    row_sc_mat_3 = True

                if (row_sc_mat_2 and "sc_matrix_2" in dct) or (row_sc_mat_3 and "sc_matrix_3" in dct):
                    rows.append(row)
            if not rows:
                raise KeyError
            for row in rows:
                db.update(
                    row.id,
                    dct=dct,
                    has_fc2=("forces_2" in dct or "forces_2" in row),
                    has_fc3=("forces_3" in dct or "forces_3" in row),
                    store_second_order=True,
                    **key_val_pairs,
                )
        except KeyError:
            db.write(
                dct,
                has_fc2=("forces_2" in dct),
                has_fc3=("forces_3" in dct),
                **key_val_pairs,
            )

    except ValueError:
        print(f"Error in adding the object to the database {db_path}")
    return hashes

def update_phonon_db(db_path, atoms, phonon, calc_type="calc", symprec=1e-5, **kwargs):
    """
    Adds a phonon dictionary to a database defined by db_path
    Args:
        db_path: str
            path to the database
        atoms: ASE Atoms, dict
            Initial atoms object for the phonopy calculations
        phonon: dict, Phonopy Object, Phono3Py Object, or ASE Atoms object
            A dictionary representation of the phonopy object to be added to the database
        calc_type: str
            keyword describing the calculation
        symprec: float
            symprec value used for making the phonopy object
        kwargs: dict
            keys to add to the database
    """
    print(f"Adding phonon calculations to the database {db_path}")
    if isinstance(atoms, dict):
        atoms = dict2atoms(atoms)
    atoms_hash, calc_hash = hash_atoms_and_calc(atoms)
    if isinstance(phonon, Atoms):
        phonon = atoms2dict(phonon)
    elif isinstance(phonon, Phonopy):
        phonon = phonon_to_dict(phonon)
    else:
        try:
            from phono3py.phonon3 import Phono3py
            if isinstance(phonon, Phono3py):
                phonon = phonon3_to_dict(phonon)
        except:
            pass
    try:
        db = connect(db_path)
        selection = [
            ("symprec", "=", symprec),
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("calc_type", "=", calc_type),
        ]
        selection_no_sc_mat = selection.copy()
        if (kwargs is not None) and ("original_atoms_hash" in kwargs):
            selection.append(("original_atoms_hash", "=", kwargs["original_atoms_hash"]))
        if (kwargs is not None) and ("sc_matrix_2" in phonon):
            selection.append(("sc_matrix_2", "=", phonon["sc_matrix_2"]))
            del(kwargs["sc_matrix_2"])
        if (kwargs is not None) and ("sc_matrix_3" in phonon):
            selection.append(("sc_matrix_3", "=", phonon["sc_matrix_3"]))
            del(kwargs["sc_matrix_3"])
        try:
            rows = list(db.select(selection=selection))
            rows_no_sc = list(db.select(selection=selection_no_sc_mat))
            for row in rows_no_sc:
                if (row.sc_matrix_2 is None and "sc_matrix_2" in phonon) != (row.sc_matrix_3 is None and "sc_matrix_3" in phonon):
                    rows.append(row)
            if not rows:
                raise KeyError
            for row in rows:
                db.update(
                    row.id,
                    dct=phonon,
                    has_fc2=("forces_2" in phonon or "forces_2" in row),
                    has_fc3=("forces_3" in phonon or "forces_3" in row),
                    calc_type=calc_type,
                    **kwargs,
                )
        except KeyError:
            db.write(
                phonon,
                symprec=symprec,
                atoms_hash=atoms_hash,
                calc_hash=calc_hash,
                has_fc2=("forces_2" in phonon),
                has_fc3=("forces_3" in phonon),
                calc_type=calc_type,
                **kwargs,
            )
    except ValueError:
        print(f"Error in adding the object to the database {db_path}")
