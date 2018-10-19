''' Define FireTasks for electronic structure calculations '''
from fireworks import FWAction, PyTask, Firework
from hilde.structure.structure import patoms2dict, dict2patoms
from hilde.tasks.calculate import calculate
from hilde.tasks import fireworks_tasks as fwt

def calculate_sp(atoms_dict, workdir, out_spec):
    '''
    A wrapper function for calculate to work within FireWorks
    Args:
        atoms_dict: dict generated from patoms2dict
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
    atoms = dict2patoms(atoms_dict)
    temp_atoms = calculate(atoms, atoms.get_calculator(), workdir)
    return FWAction(mod_spec=[{'_push': {out_spec: patoms2dict(temp_atoms)}}])

def calculate_multiple_sp(atom_dicts, workdirs):
    '''
    A wrapper function that generate FireWorks for a set of atoms and associated work
    directories
    Args:
        atom_dicts: list of dicts
            A list of dictionary representing atoms objects for the calculation
        workdirs: list of str
            A list of the paths to perform the calculations
        Returns: FWAction
            A FWAction that will add the single point force calculations to the workflow
            as detours (Adds child FireWorks to the one calling this function and transfers
            its current children to the new FireWorks)
    '''
    atom_list = [dict2patoms(ad) for ad in atom_dicts]
    firework_detours = []
    for i, cell in enumerate(atom_list):
        cell.set_calc_id(i)
        task = PyTask({"func": fwt.calculate_py_task,
                       "args": [patoms2dict(cell), workdirs[i], "calc_atoms"]})
        firework_detours.append(Firework(task))
    return FWAction(detours=firework_detours)
