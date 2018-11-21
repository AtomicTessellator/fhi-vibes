from hilde.helpers.k_grid import update_k_grid
from hilde.helpers.hash import hash_atoms
from hilde.phonopy.workflow import phonopy
from hilde.relaxation.bfgs import relax
from hilde.settings import Settings
from hilde.tasks.fireworks.general_py_task import generate_firework
from hilde.templates.aims import setup_aims

def get_step_fw(config_file, hilde_defaults_config_file, atoms):
    settings = Settings(hilde_defaults_config_file)
    step_settings = Settings(config_file)
    if "basisset" in step_settings:
        step_settings.control["basisset_type"] = step_settings.basisset.type
    calc = setup_aims(settings=settings)
    update_k_grid(atoms, calc, step_settings.control_kpt.density)
    atoms.set_calculator(calc)
    atoms_hash, calc_hash = hash_atoms(atoms)
    fw_list = []
    if not step_settings.fw.from_db:
        at = atoms
        cl = calc
    else:
        at = step_settings.fw.in_spec_atoms
        cl = step_settings.fw.in_spec_calc
    if step_settings.calculation_step.type == "relaxation":
        if "db_path" in step_settings.db_storage:
            db_kwargs = {
                "db_path": step_settings.db_storage.db_path,
                "calc_type": f"relaxation_{settings.basisset.type}",
                "symprec": 1e-5,
                "original_atom_hash": atoms_hash,
            }
        else:
            db_kwargs = {}
        fw_list.append(
            generate_firework(
                "hilde.relaxation.bfgs.relax",
                "hilde.tasks.fireworks.fw_action_outs.cont_md_out_fw_action",
                step_settings.function_kwargs,
                at,
                cl,
                func_fw_out_kwargs=db_kwargs,
                atoms_calc_from_spec= step_settings.fw.from_db,
                fw_settings=step_settings.fw,
                update_calc_settings=step_settings.control,
            )
        )
    elif step_settings.calculation_step.type == "phonons":
        if step_settings.fw['serial']:
            if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
                step_settings.function_kwargs["db_path"] = step_settings.db_storage.db_path
                step_settings.function_kwargs["original_atom_hash"] = atoms_hash
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.workflow.phonopy",
                    "hilde.tasks.fireworks.fw_action_outs.return_null_atoms",
                    dict(step_settings.function_kwargs),
                    at,
                    cl,
                    atoms_calc_from_spec= step_settings.fw.from_db,
                    fw_settings=step_settings.fw,
                    update_calc_settings=step_settings.control,
                )
            )
        else:
            kwargs_init = {"supercell_matrix": step_settings.function_kwargs.supercell_matrix}
            kwargs_init_fw_out = {"workdir": step_settings.function_kwargs.workdir}
            if "displacement" in step_settings.function_kwargs:
                kwargs_init["displacement"] = step_settings.function_kwargs.displacement
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.workflow.initialize_phonopy_attach_calc",
                    "hilde.tasks.fireworks.fw_action_outs.fw_out_initialize_phonopy",
                    kwargs_init,
                    at,
                    cl,
                    func_fw_out_kwargs=kwargs_init_fw_out,
                    atoms_calc_from_spec= step_settings.fw.from_db,
                    fw_settings=step_settings.fw,
                    update_calc_settings=step_settings.control,
                )
            )
            kwargs = {"fireworks": True}
            if "db_storage" in step_settings and "db_path" in step_settings.db_storage:
                kwargs["db_path"] = step_settings.db_storage.db_path
                kwargs["original_atom_hash"] = atoms_hash

            anal_keys = ["trajectory", "workdir", "force_constant_file", "displacement"]
            for key in anal_keys:
                if key in step_settings.function_kwargs:
                    kwargs[key] = step_settings.function_kwargs[key]
            fw_list.append(
                generate_firework(
                    "hilde.phonopy.postprocess.postprocess",
                    "hilde.tasks.fireworks.fw_action_outs.return_null_general",
                    args=[],
                    inputs=["phonon", step_settings.fw["mod_spec_add"]],
                    func_kwargs=kwargs,
                    fw_settings=step_settings.fw,
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
            del(settings.control[key])
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
            **step_settings.function_kwargs
        )
    elif step_settings.calculation_step.type == "phonons":
        if "db_path" in step_settings.db_storage:
            step_settings.function_kwargs["db_path"] = step_settings.db_storage.db_path
            step_settings.function_kwargs["original_atom_hash"] = atoms_hash
        phonopy(atoms, calc, socketio_port=settings.socketio.port, **step_settings.function_kwargs)
    else:
        raise ValueError("Type not defiend")
