"""FWAction generators for relaxations"""

from fireworks import FWAction

from vibes.fireworks.workflows.firework_generator import generate_firework, time2str
from vibes.helpers.converters import dict2atoms
from vibes.helpers.k_grid import k2d

from vibes.fireworks.tasks.postprocess.relax import check_aims


def check_aims_complete(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """A function that checks if a relaxation is converged (if outputs is True)

    either stores the relaxed structure in the MongoDB or
    appends another Firework as its child to restart the relaxation

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The original atoms at the start of this job
    calc: ase.calculators.calulator.Calculator
        The original calculator
    outputs: ase.atoms.Atoms
        The geometry of the final relaxation step
    func: str
        Path to function that performs the MD like operation
    func_fw_out: str
        Path to this function
    func_kwargs: dict
        keyword arguments for func
    func_fw_kwargs: dict
        Keyword arguments for fw_out function
    fw_settings: dict
        FireWorks specific settings

    Returns
    -------
    FWAction
        The correct action (restart or updated spec) if convergence is reached

    Raises
    ------
    RuntimeError
        If the FHI-Aims calculation fails
    """
    completed, calc_number, new_atoms_dict, walltime = check_aims(
        atoms, calc, outputs, **func_kwargs
    )

    update_spec = {}
    if completed:
        update_spec = {}
        if "out_spec_atoms" in fw_settings:
            update_spec[fw_settings["out_spec_atoms"]] = new_atoms_dict
        if "out_spec_calc" in fw_settings:
            update_spec[fw_settings["out_spec_calc"]] = calc

        return FWAction(update_spec=update_spec)

    update_spec = {}
    if "in_spec_atoms" in fw_settings:
        update_spec[fw_settings["in_spec_atoms"]] = new_atoms_dict
    if "in_spec_calc" in fw_settings:
        update_spec[fw_settings["in_spec_calc"]] = calc
    update_spec["kgrid"] = k2d(
        dict2atoms(new_atoms_dict), calc["calculator_parameters"]["k_grid"]
    )
    calc.parameters["walltime"] = walltime
    if fw_settings and "spec" in fw_settings and "_queueadapter" in fw_settings["spec"]:
        fw_settings["spec"]["_queueadapter"]["walltime"] = time2str(walltime)
    func_kwargs["walltime"] = walltime
    del calc["results"]
    fw_settings["fw_name"] = fw_settings["fw_base_name"] + str(calc_number)
    fw_settings["spec"].update(update_spec)
    fw_settings["from_db"] = False

    if "to_launchpad" in fw_settings and fw_settings["to_launchpad"]:
        fw_settings["to_launchpad"] = False

    fw = generate_firework(
        func=func,
        func_fw_out=func_fw_out,
        func_kwargs=func_kwargs,
        atoms=new_atoms_dict,
        calc=calc,
        func_fw_out_kwargs=func_fw_kwargs,
        fw_settings=fw_settings,
    )
    return FWAction(detours=[fw], update_spec=update_spec)
