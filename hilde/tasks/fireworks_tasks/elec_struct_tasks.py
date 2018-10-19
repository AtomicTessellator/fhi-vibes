''' Define FireTasks for electronic structure calculations '''
from fireworks import FWAction
from ase.db.row import atoms2dict, AtomsRow
from hilde.tasks.calculate import calculate

def calculate_sp(atoms_dict, workdir, out_spec):
    '''
    A wrapper function for calculate to work within FireWorks
    Args:
        atoms_dict: dict generated from atoms2dict
            The dictionary representation of the atoms or pAtoms object
        workdir: str
            The base work directory for the calculation
        out_spec: str
            A key in the analysis' FireWork's spec to add the calculated
            results
    Returns: FWAction
        A FWAction that will modify the spec of the analysis FireWork to
        include the calculated results (pushes it to the end of a list)
    '''
    atoms = AtomsRow(atoms_dict).toatoms(attach_calculator=True)
    temp_atoms = calculate(atoms, atoms.get_calculator(), workdir)
    return FWAction(mod_spec=[{'_push': {out_spec: atoms2dict(temp_atoms)}}])
