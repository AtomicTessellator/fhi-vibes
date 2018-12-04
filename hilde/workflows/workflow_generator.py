from hilde.helpers.hash import hash_atoms
from hilde.helpers.k_grid import update_k_grid
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.templates.aims import setup_aims


def get_step_fw(step_settings, atoms=None):
    if "geometry" in step_settings:
        atoms, calc = setup_aims(settings=step_settings)
    else:
        calc = setup_aims(settings=step_settings)

    # print(calc.parameters, step_settings)
    update_k_grid(atoms, calc, step_settings.control_kpt.density)
    atoms.set_calculator(calc)
    atoms_hash, calc_hash = hash_atoms(atoms)
    fw_list = []
    if (
        "from_db" not in step_settings.fw_settings
        or not step_settings.fw_settings.from_db
    ):
        at = atoms
        cl = calc
    else:
        at = step_settings.fw_settings.in_spec_atoms
        cl = step_settings.fw_settings.in_spec_calc

    if "fw_spec" in step_settings:
        step_settings.fw_settings["spec"] = dict(step_settings.fw__spec)
    else:
        step_settings.fw_settings["spec"] = {}
    if "fw_spec_qadapter" in step_settings:
        step_settings.fw_settings["spec"]["_queueadapter"] = dict(
            step_settings.fw_spec_qadapter
        )

    if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
        db_kwargs = {
            "db_path": step_settings.db_storage.db_path,
            "calc_type": f"relaxation_{step_settings.basisset.type}",
            "symprec": 1e-5,
            "original_atom_hash": atoms_hash,
        }
    else:
        db_kwargs = {}
    if "fw_name" in step_settings:
        step_settings["fw_name"] += (
            atoms.symbols.get_chemical_formula() + "_" + hash_atoms(atoms)[0]
        )
    if "relaxation" in step_settings:
        fw_list.append(
            generate_firework(
                "hilde.relaxation.bfgs.relax",
                "hilde.tasks.fireworks.fw_action_outs.check_relaxation_complete",
                step_settings.relaxation,
                at,
                cl,
                func_fw_out_kwargs=db_kwargs,
                atoms_calc_from_spec=step_settings.fw_settings.from_db,
                fw_settings=step_settings.fw_settings,
                update_calc_settings=step_settings.control,
            )
        )
    elif "aims_relaxation" in step_settings:
        db_kwargs["relax_step"] = 0
        fw_list.append(
            generate_firework(
                "hilde.tasks.calculate.calculate",
                "hilde.tasks.fireworks.fw_action_outs.check_aims_relaxation_complete",
                step_settings.aims_relaxation,
                at,
                cl,
                func_fw_out_kwargs=db_kwargs,
                atoms_calc_from_spec=step_settings.fw_settings.from_db,
                fw_settings=step_settings.fw_settings,
                update_calc_settings=step_settings.control,
            )
        )
    elif "kgrid_opt" in step_settings:
        fw_list.append(
            generate_firework(
                "hilde.auto_tune_parameters.k_grid.converge_kgrid.converge_kgrid",
                "hilde.tasks.fireworks.fw_action_outs.check_kgrid_opt_completion",
                step_settings.kgrid_opt,
                at,
                cl,
                fw_settings=step_settings.fw_settings,
                update_calc_settings=None,
            )
        )
    elif "phonopy" in step_settings:
        if step_settings.fw_settings["serial"]:
            if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
                step_settings.phonopy["db_path"] = step_settings.db_storage.db_path
                step_settings.phonopy["original_atom_hash"] = atoms_hash
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.workflow.run",
                    "hilde.tasks.fireworks.fw_action_outs.serial_phonopy_continue",
                    dict(step_settings.phonopy),
                    at,
                    cl,
                    atoms_calc_from_spec=step_settings.fw_settings.from_db,
                    fw_settings=step_settings.fw_settings,
                    update_calc_settings=step_settings.control,
                )
            )
        else:
            kwargs_init = {"supercell_matrix": step_settings.phonopy.supercell_matrix}
            kwargs_init_fw_out = {"workdir": step_settings.phonopy.workdir}
            if "displacement" in step_settings.phonopy:
                kwargs_init["displacement"] = step_settings.phonopy.displacement
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.workflow.preprocess_fireworks",
                    "hilde.tasks.fireworks.fw_action_outs.add_phonon_force_calcs",
                    kwargs_init,
                    at,
                    cl,
                    func_fw_out_kwargs=kwargs_init_fw_out,
                    atoms_calc_from_spec=step_settings.fw_settings.from_db,
                    fw_settings=step_settings.fw_settings,
                    update_calc_settings=step_settings.control,
                )
            )
            kwargs = {"fireworks": True}
            if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
                kwargs["db_path"] = step_settings.db_storage.db_path
                kwargs["original_atom_hash"] = atoms_hash

            anal_keys = [
                "trajectory",
                "analysis_workdir",
                "force_constant_file",
                "displacement",
                "workdir",
            ]
            for key in anal_keys:
                if key in step_settings.phonopy:
                    kwargs[key] = step_settings.phonopy[key]
            if "analysis_workdir" in kwargs:
                kwargs["workdir"] = kwargs["analysis_workdir"]
                del (kwargs["analysis_workdir"])
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.postprocess.postprocess",
                    "hilde.tasks.fireworks.fw_action_outs.fireworks_no_mods_gen_function",
                    args=[],
                    inputs=["phonon", step_settings.fw_settings["mod_spec_add"]],
                    func_kwargs=kwargs,
                    fw_settings=step_settings.fw_settings,
                )
            )
    elif "phono3py" in step_settings:
        if step_settings.fw_settings["serial"]:
            if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
                step_settings.phono3py["db_path"] = step_settings.db_storage.db_path
                step_settings.phono3py["original_atom_hash"] = atoms_hash
            fw_list.append(
                generate_firework(
                    "hilde.phono3py.workflow.run",
                    "hilde.tasks.fireworks.fw_action_outs.serial_phonopy_continue",
                    dict(step_settings.phono3py),
                    at,
                    cl,
                    atoms_calc_from_spec=step_settings.fw_settings.from_db,
                    fw_settings=step_settings.fw_settings,
                    update_calc_settings=step_settings.control,
                )
            )
        else:
            kwargs_init = {
                "supercell_matrix": step_settings.phono3py.supercell_matrix,
            }
            kwargs_init_fw_out = {"workdir": step_settings.phono3py.workdir}
            if "displacement" in step_settings.phono3py:
                kwargs_init["displacement"] = step_settings.phono3py.displacement
            fw_list.append(
                generate_firework(
                    "hilde.phono3py.workflow.preprocess_fireworks",
                    "hilde.tasks.fireworks.fw_action_outs.add_phonon_force_calcs",
                    kwargs_init,
                    at,
                    cl,
                    func_fw_out_kwargs=kwargs_init_fw_out,
                    atoms_calc_from_spec=step_settings.fw_settings.from_db,
                    fw_settings=step_settings.fw_settings,
                    update_calc_settings=step_settings.control,
                )
            )
            kwargs = {"fireworks": True}
            if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
                kwargs["db_path"] = step_settings.db_storage.db_path
                kwargs["original_atom_hash"] = atoms_hash

            anal_keys = [
                "trajectory",
                "analysis_workdir",
                "force_constant_file",
                "displacement",
            ]
            for key in anal_keys:
                if key in step_settings.phono3py:
                    kwargs[key] = step_settings.phono3py[key]
            if "analysis_workdir" in kwargs:
                kwargs["workdir"] = kwargs["analysis_workdir"]
                del (kwargs["analysis_workdir"])
            fw_list.append(
                generate_firework(
                    "hilde.phono3py.postprocess.postprocess",
                    "hilde.tasks.fireworks.fw_action_outs.fireworks_no_mods_gen_function",
                    args=[],
                    inputs=["phonon", step_settings.fw_settings["mod_spec_add"]],
                    func_kwargs=kwargs,
                    fw_settings=step_settings.fw_settings,
                )
            )
    else:
        raise ValueError("Type not defiend")

    return fw_list
