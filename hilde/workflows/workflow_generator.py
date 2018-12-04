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

    if "from_db" not in fw_settings or not fw_settings.from_db:
        at = atoms
        cl = calc
    else:
        at = fw_settings.in_spec_atoms
        cl = fw_settings.in_spec_calc

    fw_settings = step_settings.fw_settings
    if "fw_spec" in step_settings:
        fw_settings["spec"] = dict(step_settings.fw__spec)
    else:
        fw_settings["spec"] = {}
    if "fw_spec_qadapter" in step_settings:
        fw_settings["spec"]["_queueadapter"] = dict(
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
    task_spec_list = []
    if "relaxation" in step_settings:
        task_spec_list.append(
            TaskSpec(
                "hilde.relaxation.bfgs.relax",
                "hilde.tasks.fireworks.fw_action_outs.check_relaxation_complete",
                True,
                step_settings.relaxation,
                func_fw_out_kwargs=db_kwargs,
            )
        )
    elif "aims_relaxation" in step_settings:
        db_kwargs["relax_step"] = 0
        task_spec_list.append(
            TaskSpec(
                "hilde.tasks.calculate.calculate",
                "hilde.tasks.fireworks.fw_action_outs.check_aims_relaxation_complete",
                True,
                step_settings.aims_relaxation,
                func_fw_out_kwargs=db_kwargs,
            )
        )
    elif "kgrid_opt" in step_settings:
        task_spec_list.append(
            TaskSpec(
                "hilde.auto_tune_parameters.k_grid.converge_kgrid.converge_kgrid",
                "hilde.tasks.fireworks.fw_action_outs.check_kgrid_opt_completion",
                True,
                step_settings.kgrid_opt,
            )
        )
    elif "phonopy" in step_settings or "phono3py" in step_settings:
        if fw_settings["serial"]:
            if "phonopy" in step_settings and "phono3py" in step_settings:
                raise ArgumentError("For a serial calculation phonopy and phono3py must be separate steps")
            if "phonopy" in step_settings:
                func = "hilde.phonopy.workflow.run"
                func_kwargs = step_settings.phonopy
            else:
                func = "hilde.phono3py.workflow.run"
                func_kwargs = step_settings.phono3py
            if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
                step_settings.phonopy["db_path"] = step_settings.db_storage.db_path
                step_settings.phonopy["original_atom_hash"] = atoms_hash
            task_spec_list.append(
                TaskSpec(
                    func,
                    "hilde.tasks.fireworks.fw_action_outs.check_kgrid_opt_completion",
                    True,
                    func_kwargs,
                )
            )
        else:
            if "phonopy" in step_settings:
                kwargs_init = {"supercell_matrix": step_settings.phonopy.supercell_matrix}
                kwargs_init_fw_out = {"workdir": step_settings.phonopy.workdir}
                if "displacement" in step_settings.phonopy:
                    kwargs_init["displacement"] = step_settings.phonopy.displacement
                task_spec_list.append(
                    TaskSpec(
                        "hilde.phonopy.workflow.preprocess_fireworks",
                        "hilde.tasks.fireworks.fw_action_outs.add_phonon_force_calcs",
                        kwargs_init,
                        func_fw_out_kwargs=kwargs_init_fw_out,
                    )
                )
            if "phono3py" in step_settings:
                kwargs_init = {"supercell_matrix": step_settings.phono3py.supercell_matrix}
                kwargs_init_fw_out = {"workdir": step_settings.phono3py.workdir}
                if "displacement" in step_settings.phono3py:
                    kwargs_init["displacement"] = step_settings.phono3py.displacement
                task_spec_list.append(
                    TaskSpec(
                        "hilde.phono3py.workflow.preprocess_fireworks",
                        "hilde.tasks.fireworks.fw_action_outs.add_phonon_force_calcs",
                        kwargs_init,
                        func_fw_out_kwargs=kwargs_init_fw_out,
                    )
                )
    else:
        raise ValueError("Type not defiend")

    fw_list = [
        generate_firework(
            task_spec_list,
            at,
            cl,
            atoms_calc_from_spec=step_settings fw_settings.from_db,
            fw_settings=fw_settings
            update_calc_settings=step_settings.control,
        )
    ]
    task_spec_list = []
    if ("phonopy" in step_settings or "phono3py" in step_settings) and not fw_settings["serial"]:
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
        if "phonopy" in step_settings:
            func_analysis = "hilde.phonopy.postprocess.postprocess"
            for key in anal_keys:
                if key in step_settings.phonopy:
                    kwargs[key] = step_settings.phonopy[key]
            if "analysis_workdir" in kwargs:
                kwargs["workdir"] = kwargs["analysis_workdir"]
                del (kwargs["analysis_workdir"])
            task_spec_list.append(
                TaskSpec(
                    "hilde.phono3py.postprocess.postprocess",
                    "hilde.tasks.fireworks.fw_action_outs.fireworks_no_mods_gen_function",
                    False,
                    inputs=["phonon", fw_settings["mod_spec_add"]],
                    func_kwargs=kwargs,
                )
            )

        if "phono3py" in step_settings:
            func_analysis = "hilde.phono3py.postprocess.postprocess"
            for key in anal_keys:
                if key in step_settings.phono3py:
                    kwargs[key] = step_settings.phono3py[key]
            if "analysis_workdir" in kwargs:
                kwargs["workdir"] = kwargs["analysis_workdir"]
                del (kwargs["analysis_workdir"])
            task_spec_list.append(
                TaskSpec(
                    "hilde.phono3py.postprocess.postprocess",
                    "hilde.tasks.fireworks.fw_action_outs.fireworks_no_mods_gen_function",
                    False,
                    inputs=["phonon", fw_settings["mod_spec_add"]],
                    func_kwargs=kwargs,
                )
            )
        fw_list.append(
            generate_firework(
                task_spec_list,
                at,
                cl,
                fw_settings=fw_settings
            )
        )
    return fw_list


