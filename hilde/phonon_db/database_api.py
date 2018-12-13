from ase.atoms import Atoms
from phonopy import Phonopy

from hilde.helpers.converters import atoms2dict, dict2atoms
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonon_db.phonon_db import connect
from hilde.phonon_db.row import phonon_to_dict
from hilde.phonon_db.row import phonon3_to_dict

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
