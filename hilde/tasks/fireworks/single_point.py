""" Define FireTasks for electronic structure calculations """
from fireworks import FWAction, PyTask, Firework
from hilde.helpers.converters import atoms2dict, dict2atoms
from hilde.tasks.calculate import calculate as calc_hilde
from hilde.tasks import fireworks as fw

module_name = __name__


def calculate_none(out_spec):
    return FWAction(mod_spec=[{"_push": {out_spec: None}}])


def calculate(workdir, out_spec, atoms_dict, calc=None):
    """
    A wrapper function for calculate to work within FireWorks
    Args:
        atoms_dict: dict generated from atoms2dict
            The dictionary representation of the atoms or Atoms object
        workdir: str
            The base work directory for the calculation
        out_spec: str
            A key in the analysis' FireWork's spec to add the calculated
            results
        calc: dict
            Dictionary of the ase calculator
    Returns: FWAction
        A FWAction that will modify the spec of the analysis FireWork to
        include the calculated results (pushes it to the end of a list)
    """
    print(f"Beginning single point calc in {workdir}")
    if calc:
        for key, val in calc.items():
            atoms_dict[key] = val
    atoms = dict2atoms(atoms_dict)
    temp_atoms = calc_hilde(atoms, atoms.get_calculator(), workdir)
    return FWAction(mod_spec=[{"_push": {out_spec: atoms2dict(temp_atoms)}}])


def calculate_multiple(workdirs, atom_dicts, calculator=None, out_spec="calc_atoms", spec_qad=None):
    """
    A wrapper function that generate FireWorks for a set of atoms and associated work
    directories
    Args:
        atom_dicts: list of dicts
            A list of dictionary representing atoms objects for the calculation
        workdirs: list of str
            A list of the paths to perform the calculations
        calc_mods: dict
            A dictionary describing all modifications needed for the calculator
        spec_qad: dict
            Updated spec for job submission
        Returns: FWAction
            A FWAction that will add the single point force calculations to the workflow
            as detours (Adds child FireWorks to the one calling this function and transfers
            its current children to the new FireWorks)
    """
    print(f"Setting up calculations in {workdirs}, and appending their Fireworks to the workflow")
    firework_detours = []
    if spec_qad is None:
        spec_qad = {}
    for i, cell in enumerate(atom_dicts):
        if cell is None:
            task = PyTask({"func": fw.calculate_none.name, "args": [out_spec]})
        else:
            task = PyTask(
                {"func": fw.calculate.name, "args": [workdirs[i], out_spec, cell, calculator]}
            )
        firework_detours.append(Firework(task, name=f"calc_{i}", spec=spec_qad))
    return FWAction(detours=firework_detours)


calculate.name = f"{module_name}.{calculate.__name__}"
calculate_none.name = f"{module_name}.{calculate_none.__name__}"
calculate_multiple.name = f"{module_name}.{calculate_multiple.__name__}"
