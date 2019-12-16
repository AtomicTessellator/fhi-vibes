"""Interface Functions for the database"""

from ase.atoms import Atoms
from phonopy import Phonopy

from vibes.helpers.hash import hash_atoms_and_calc, hash_dict, hash_traj
from vibes.helpers.warnings import warn
from vibes.helpers.converters import dict2atoms as traj_dict2atoms
from vibes.ase.db.dict_converters import dict2atoms, atoms2dict
from vibes.phonon_db.phonon_db import connect
from vibes.phonon_db.row import phonon_to_dict
from vibes.phonopy.postprocess import postprocess as ph_postprocess
from vibes.structure.convert import to_Atoms_db, to_Atoms

from vibes.trajectory import reader

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


def traj_to_database(db_path, traj, ret_all_hashes=False):
    """Processes a trajectory file and adds it to the database

    Parameters
    ----------
    db_path: str
        path to the database
    traj: str
        The trajectory's path
    ret_all_hashes: bool
        If True return all hashes

    Returns
    -------
    str
        The hash of the dictionary representation of phonon
    dict
        The dictionary of all the hashes

    Raises
    ------
    IOError
        If trajectory does not have Phonopy or Phono3py metadata
    """
    calc_atoms, metadata = reader(traj, True, verbose=False)

    if "Phonopy" in metadata:
        phonon = ph_postprocess(traj, pickle_file=None)
        at_dict = metadata["Phonopy"]["primitive"].copy()
    # elif "Phono3py" in metadata:
    #     phonon = ph3_postprocess(traj, pickle_file=None)
    #     from vibes.phono3py.postprocess import postprocess as ph3_postprocess

    #     at_dict = metadata["Phono3py"]["primitive"].copy()
    else:
        raise IOError("This trajectory can't be processed")

    calc = traj_dict2atoms(at_dict, metadata["calculator"], False).get_calculator()

    return to_database(
        db_path,
        phonon,
        calc=calc,
        ret_all_hashes=ret_all_hashes,
        traj_hash=hash_traj(calc_atoms, metadata, True),
    )


def to_database(
    db_path, phonon, calc=None, key_val_pairs=None, ret_all_hashes=False, traj_hash=None
):
    """Adds a Phonopy, Phono3py or ASE Atoms object to the database

    Parameters
    ----------
    db_path: str
        Path to the database
    phonon: ase.atoms.Atoms, phonopy.Phonopy, or phono3py.phonon3.Phono3py
        Object to be added to the database
    calc: ase.calculators.calulator.Calculator
        Calculator parameters to add to the Database
    key_val_pairs: dict
        Additional key_val_pairs to add to the database
    ret_all_hashes: bool
        if True return a dict of all hashes

    Returns
    -------
    str
        The hash of the dictionary representation of phonon
    hashes: dict
        The dictionary of all the hashes

    Raises
    ------
    IOError
        If Phono3py is not installed OR
        If phonon is not a dict, Atoms, Phonopy, or Phono3py object
    """
    selection_no_sc = []
    selection = []
    dct = obj2dict(phonon)
    hashes = {}

    if not key_val_pairs:
        key_val_pairs = {}
    if traj_hash:
        hashes["traj_hash"] = traj_hash[0]
        hashes["meta_hash"] = traj_hash[1]
    if isinstance(phonon, dict):
        atoms = dict2atoms(dct)
    elif isinstance(phonon, Atoms):
        atoms = phonon.copy()
    elif isinstance(phonon, Phonopy):
        atoms = to_Atoms_db(phonon.get_primitive())

        hashes["cell_hash"] = hash_atoms_and_calc(to_Atoms(phonon.get_primitive()))[0]

        selection_no_sc.append(("cell_hash", "=", hashes["cell_hash"]))
        selection.append(("sc_matrix_2", "=", dct["sc_matrix_2"]))

        key_val_pairs["displacement"] = phonon._displacement_dataset["first_atoms"][0][
            "displacement"
        ][0]
    else:
        # try:
        #     from phono3py.phonon3 import Phono3py
        #     from vibes.phonon_db.row import phonon3_to_dict

        #     if isinstance(phonon, Phono3py):
        #         atoms = to_Atoms_db(phonon.get_primitive())

        #         hashes["cell_hash"] = hash_atoms_and_calc(
        #             to_Atoms(phonon.get_primitive())
        #         )[0]

        #         selection_no_sc.append(("cell_hash", "=", hashes["cell_hash"]))
        #         selection.append(("sc_matrix_3", "=", dct["sc_matrix_3"]))

        #         key_val_pairs["displacement"] = phonon._displacement_dataset[
        #             "first_atoms"
        #         ][0]["displacement"][0]
        #         key_val_pairs["pair_distance_cutoff"] = phonon._displacement_dataset[
        #             "cutoff_distance"
        #         ]
        # except:
        raise ValueError("Object is not defined")

    atoms.set_calculator(calc)

    if "calc_hash" in key_val_pairs or "atoms_hash" in key_val_pairs:
        warn("Replacing the atoms and calc hashes")

    hashes["atoms_hash"], hashes["calc_hash"] = hash_atoms_and_calc(atoms)

    for key, val in atoms2dict(atoms).items():
        dct[key] = val
    selection_no_sc.append(("atoms_hash", "=", hashes["atoms_hash"]))
    selection_no_sc.append(("calc_hash", "=", hashes["calc_hash"]))

    if "symprec" in key_val_pairs.keys():
        selection_no_sc.append(("symprec", "=", key_val_pairs["symprec"]))

    for sel in selection_no_sc:
        selection.append(sel)

    if "displacement" in key_val_pairs and "pair_distance_cutoff" in key_val_pairs:
        hashes["ph3_hash"] = hash_dict(dct, ignore_keys=["unique_id"])
    elif "displacement" in key_val_pairs:
        hashes["ph_hash"] = hash_dict(dct, ignore_keys=["unique_id"])
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

                if (row_sc_mat_2 and "sc_matrix_2" in dct) or (
                    row_sc_mat_3 and "sc_matrix_3" in dct
                ):
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
    if ret_all_hashes:
        return hashes
    return hash_dict(dct)


def obj2dict(obj):
    """Converts a Phonopy, Phono3py, or ase.atoms.Atoms to a dict

    Parameters
    ----------
    obj: phonopy.Phonopy, phono3py.phonon3.Phono3py, or ase.atoms.Atoms
        object to be converted to a dict

    Returns
    -------
    dict
        The dictionary representation of the obj

    Raises
    ------
    IOError
        If obj is not a dict, phonopy.Phonopy, ase.atoms.Atoms, or phono3py.phonon3.Phono3py object OR
        If obj is not a dict, phonopy.Phonopy, or ase.atoms.Atoms and Phono3py is not installed
    """
    if isinstance(obj, dict):
        return obj.copy()
    elif isinstance(obj, Atoms):
        return atoms2dict(obj)
    elif isinstance(obj, Phonopy):
        return phonon_to_dict(obj)
    else:
        # try:
        #     from phono3py.phonon3 import Phono3py

        #     if isinstance(obj, Phono3py):
        #         return phonon3_to_dict(obj)
        # except ImportError:
        #     raise IOError(
        #         "Phono3py is not installed, obj has to be a dict, Phonopy, or Atoms object"
        #     )
        raise ValueError("obj has to be a dict, Phonopy, or Atoms object")


def from_database(
    db_path,
    selection=None,
    get_id=False,
    get_atoms=False,
    get_phonon=False,
    get_phonon3=False,
    **kwargs,
):
    """Pulls an object from the database

    Parameters
    ----------
    db_path: str
        Path to the database
    selection: list
        A list of tuples used as selection parameters
    get_id: bool
        If True return row id
    get_atoms:(ool
        If True return the ASE Atoms object
    get_phonon: bool
        If True return the Phonopy Object
    get_phonon3: bool
        If True return the Phono3py Object
    kwargs: dict
        Optional kwargs to add to selection, if "phonopy_hash" or "phono3py_hash" in kwargs, then get_phonon or get_phonon3 is automatically set to True. with the rest set to False

    Returns
    -------
    to_ret: tuple
        All the desired outputs
    """
    db = connect(db_path)
    if selection is None:
        selection = []
    for key, val in kwargs.items():
        selection.append((key, "=", val))
    if "ph_hash" in kwargs:
        get_phonon = True
        get_atoms = False
        get_phonon3 = False
    if "ph3_hash" in kwargs:
        get_phonon3 = True
        get_atoms = False
        get_phonon = False

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
