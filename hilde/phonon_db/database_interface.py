from ase.atoms import Atoms
from phonopy import Phonopy

from hilde.helpers.converters import atoms2dict, dict2atoms
from hilde.helpers.hash import hash_atoms, hash_atoms_and_calc
from hilde.helpers.warnings import warn
from hilde.phonon_db.phonon_db import connect
from hilde.phonon_db.row import phonon_to_dict
from hilde.phonon_db.row import phonon3_to_dict
from hilde.structure.convert import to_Atoms_db as to_Atoms

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
    if isinstance(obj, dict):
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

def from_database(db_path, obj=None, selection=[], get_id=False, get_atoms=False, get_phonon=False, get_phonon3=False, **kwargs):
    db = connect(db_path)

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

def to_database(db_path, obj, calc, key_val_pairs={}):
    selection_no_sc = []
    selection = []
    dct = obj2dict(obj)

    if isinstance(obj, dict):
        atoms = dict2atoms(dct)
    elif isinstance(obj, Atoms):
        atoms = obj.copy()
    elif isinstance(obj, Phonopy):
        atoms = to_Atoms(obj.get_primitive())
        selection.append(("sc_matrix_2", "=", dct["sc_matrix_2"]))
        key_val_pairs["displacement"] = obj._displacement_dataset['first_atoms'][0]['displacement'][0]
    else:
        try:
            from phono3py.phonon3 import Phono3py
            if isinstance(obj, Phono3py):
                atoms = to_Atoms(obj.get_primitive())
                selection.append(("sc_matrix_3", "=", dct["sc_matrix_3"]))
                key_val_pairs["displacement"] = obj._displacement_dataset['first_atoms'][0]['displacement'][0]
                key_val_pairs["pair_distance_cutoff"] = obj._displacement_dataset['cutoff_distance']
        except:
            pass
    atoms.set_calculator(calc)
    atoms_hash, calc_hash = hash_atoms_and_calc(atoms)

    if "calc_hash" in key_val_pairs or "atoms_hash" in key_val_pairs:
        warn("Replacing the atoms and calc hashes")
    key_val_pairs['atoms_hash'] = atoms_hash
    key_val_pairs['calc_hash'] = calc_hash
    for key, val in atoms2dict(atoms).items():
        dct[key] = val
    if "symprec" in key_val_pairs.keys():
        selection_no_sc.append(("symprec", "=", key_val_pairs["symprec"]))

    for sel in selection_no_sc:
        selection.append(sel)

    db = connect(db_path)
    try:
        try:
            rows = list(db.select(selection=selection))
            rows_no_sc = list(db.select(selection=selection_no_sc))
            for row in rows_no_sc:
                if (row.sc_matrix_2 is None and "sc_matrix_2" in dct) or (row.sc_matrix_3 is None and "sc_matrix_3" in dct):
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
