''' Define FireTasks for electronic structure calculations '''
import numpy as np
from fireworks import FWAction, Firework, PyTask
from ase.db.row import atoms2dict, AtomsRow

from hilde.helpers.hash import hash_atoms
from hilde.phonon_db.phonon_db import connect
from hilde.phonon_db.row import phonon2dict
from hilde.phonopy import phono as ph
from hilde.structure.structure import pAtoms
from hilde.tasks.calculate import setup_multiple
from hilde.tasks import fireworks_tasks as fwt

def initialize_phonopy(atoms_ideal, smatrix, workdir):
    '''
    A wrapper function to initialize all phonopy calculations and add the new
    FireWorks as detours in the workflow
    Args:
        atoms_ideal: dict generated from atoms2dict
            The dictionary representation of the atoms or pAtoms object of the ndisplaced atoms
        smatrix: list
            The supercell matrix of the calculation
        workdir: str
            The base work directory for the calculation
    Returns: FWAction
        A FWAction that will add the single point force calculations to the workflow
        as detours (Adds child FireWorks to the one calling this function and transfers
        its current children to the new FireWorks)
    '''
    atoms = pAtoms(AtomsRow(dict(atoms_ideal)).toatoms(attach_calculator=True))
    smatrix = np.array(smatrix).reshape(3, 3)
    _, _, supercells_with_disps = ph.preprocess(atoms, smatrix.T)
    supercells_with_disps, workdirs = setup_multiple(supercells_with_disps,
                                                     atoms.get_calculator(),
                                                     workdir)
    atom_dicts = [atoms2dict(cell) for cell in supercells_with_disps]
    return FWAction(update_spec={"atom_dicts": atom_dicts, "workdirs": workdirs})


def analyze_phonopy(atoms_ideal, smatrix, db_path, calc_atoms):
    '''
    A wrapper function to calculate 2nd order force constants with phonopy where
    the displacement cells were calculated in its parent FireWork and adds the
    results to a phonon_db
    Args:
        atoms_ideal: dict generated from atoms2dict
            The dictionary representation of the atoms or pAtoms object of the undisplayed atoms
        smatrix: list
            The supercell matrix of the calculation
        db_path: str
            Path to the phonon_db where these calculations would be included
        calc_atoms: str
            The key in the spec where the displaced atoms with calculated forces are stored
    Returns: FWAction
        A FWAction that will store the calculated phonopy objects in the MongoDB for FireWorks
    '''
    atoms = pAtoms(AtomsRow(dict(atoms_ideal)).toatoms(attach_calculator=True))
    disp_cells = [pAtoms(AtomsRow(dict(ca)).toatoms(attach_calculator=True)) for ca in calc_atoms]
    smatrix = np.array(smatrix).reshape(3, 3)
    phonon, _, _ = ph.preprocess(atoms, smatrix.T)
    disp_cells = sorted(disp_cells, key=lambda x: x.calc_id)
    phonon.set_forces([cell.get_forces() for cell in disp_cells])
    atoms_hash, calc_hash = hash_atoms(atoms)
    phonon.produce_force_constants()
    try:
        db = connect(db_path)
        try:
            rows = list(db.select(selection=[("supercell_matrix", "=", smatrix),
                                             ("atoms_hash", "=", atoms_hash),
                                             ("calc_hash", "=", calc_hash)
                                            ]))
            if not rows:
                raise KeyError
            for row in rows:
                db.update(row.id, phonon=phonon, has_fc=not phonon.get_force_constants() is None)
        except KeyError:
            db.write(phonon,
                     atoms_hash=atoms_hash,
                     calc_hash=calc_hash,
                     has_fc=(phonon.get_force_constants() is not None))
    except ValueError:
        print(f"Fireworker could not access the database {db_path}")
    return FWAction(stored_data={'phonopy_calc': phonon2dict(phonon)})
