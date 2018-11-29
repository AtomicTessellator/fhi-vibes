from hilde.helpers.hash import hash_atoms
from hilde.helpers.k_grid import update_k_grid
from hilde.phonopy.workflow import phonopy
from hilde.relaxation.bfgs import relax
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
    if "from_db" not in step_settings.fw_settings or not step_settings.fw_settings.from_db:
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
        step_settings.fw_settings["spec"]["_queueadapter"] = dict(step_settings.fw_spec_qadapter)

    if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
        db_kwargs = {
            "db_path": step_settings.db_storage.db_path,
            "calc_type": f"relaxation_{step_settings.basisset.type}",
            "symprec": 1e-5,
            "original_atom_hash": atoms_hash,
        }
    else:
        db_kwargs = {}

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
                step_settings.phonopy[
                    "db_path"
                ] = step_settings.db_storage.db_path
                step_settings.phonopy["original_atom_hash"] = atoms_hash
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.workflow.phonopy",
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
            kwargs_init = {
                "supercell_matrix": step_settings.phonopy.supercell_matrix
            }
            kwargs_init_fw_out = {"workdir": step_settings.phonopy.workdir}
            if "displacement" in step_settings.phonopy:
                kwargs_init["displacement"] = step_settings.phonopy.displacement
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.workflow.initialize_phonopy_attach_calc",
                    "hilde.tasks.fireworks.fw_action_outs.add_phonopy_force_calcs",
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
    else:
        raise ValueError("Type not defiend")

    return fw_list


def run_step(config_file, hilde_defaults_config_file, atoms):
    settings = Settings(hilde_defaults_config_file)
    step_settings = Settings(config_file)
    # reset default calculator values
    for key, val in step_settings.control.items():
        if val is None and key in settings.control:
            del (settings.control[key])
        elif val is not None:
            settings.control[key] = val
    if "basisset" in step_settings:
        step_settings.control["basisset_type"] = step_settings.basisset.type
        settings.basisset.type = step_settings.basisset.type
    # Construct calculator
    calc = setup_aims(settings=settings)
    update_k_grid(atoms, calc, step_settings.control_kpt.density)
    atoms.set_calculator(calc)
    atoms_hash, calc_hash = hash_atoms(atoms)
    if step_settings.calculation_step.type == "relaxation":
        relax(
            atoms,
            calc,
            socketio_port=settings.socketio.port,
            **step_settings.function_kwargs,
        )
    elif step_settings.calculation_step.type == "phonons":
        if "db_path" in step_settings.db_storage:
            step_settings.function_kwargs["db_path"] = step_settings.db_storage.db_path
            step_settings.function_kwargs["original_atom_hash"] = atoms_hash
        phonopy(
            atoms,
            calc,
            socketio_port=settings.socketio.port,
            **step_settings.function_kwargs,
        )
    else:
        raise ValueError("Type not defiend")
