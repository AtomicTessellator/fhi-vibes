


def serial_phonopy_continue(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that checks if a set of force calculations are completed, and if not adds another Firework to the Workflow to continue the calculation
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    update_spec = {}
    at = atoms
    cl = calc
    fw_settings['from_db'] = False
    if "kpoint_density_spec" in fw_settings:
        del(fw_settings["kpoint_density_spec"])

    if outputs:
        if converged:
            return FWAction(update_spec=update_spec)
        else:
            new_func_kwargs = func_kwargs.copy()
            if "supercell_matrix" in new_func_kwargs:
                sc_matrix = np.array(new_func_kwargs.pop("supercell_matrix")).reshape(3,3)
                new_func_kwargs["n_atoms_in_sc"] = int(np.linalg.det(sc_matrix)) + 50
            elif "n_atoms_in_sc" in new_func_kwargs:
                new_func_kwargs["n_atoms_in_sc"] += 50
            else:
                raise AttributeError("either supercell_matrix or n_atoms_in_sc must be defined")
            new_func_kwargs["workdir"] += "/n_atoms_in_sc_" + str(new_func_kwargs["n_atoms_in_sc"]) + "/"
            fw = generate_firework(
                func=func,
                func_fw_out=func_fw_out,
                func_kwargs=new_func_kwargs,
                atoms=at,
                calc=cl,
                func_fw_out_kwargs=func_fw_kwargs,
                fw_settings=fw_settings,
            )
            return FWAction(detours[fw], update_spec)
    fw = generate_firework(
        func=func,
        func_fw_out=func_fw_out,
        func_kwargs=func_kwargs,
        atoms=at,
        calc=cl,
        func_fw_out_kwargs=func_fw_kwargs,
        fw_settings=fw_settings,
    )
    return FWAction(detours=[fw],update_spec=update_spec)

def post_bootsrap(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    detours = []
    update_spec = {}
    fw_settings = fw_settings.copy()

    if "phononpy" in outputs:
        ph_fw_set = fw_settings.copy()
        ph_outputs = outputs["phononpy"]
        ph_settings = func_kwargs["phonon_settings"].copy()
        update_spec["metadata_ph"] = ph_outputs["metadata"]
        ph_fw_set["metadata_spec"] = "metadata_ph"
        calc_dict = calc2dict(ph_outputs["calculator"])
        if ph_settings["serial"]:
            update_spec["ph_calculated_atoms"] = [atoms2dict(at) for at in ph_outputs["atoms_to_calculate"]]
            update_spec["ph_calculator"] = calc_dict
            ph_fw_set["calc_atoms_spec"] = "ph_calculated_atoms"
            ph_fw_set["calc_spec"] = "ph3_calculator"

            calc_kwargs = {}
            calc_keys = ["trajectory", "workdir", "backup_folder"]
            for key in calc_keys:
                if key in ph_settings:
                    calc_kwargs[key] = ph_settings[key]
            fw = generate_firework(
                func="hilde.tasks.calculate.calculate_socket",
                func_fw_out="hilde.tasks.fireworks.fw_action_outs.serial_calc",
                func_kwargs=calc_kwargs,
                atoms_calc_from_spec=False,
                inputs=["ph_calculated_atoms", "ph_calculator", "metadata_ph"],
                fw_settings=ph_fw_set,
            )
            detour.append(fw)
        else:
            fw_settings["mod_spec_add"] = "ph_forces"
            detours = add_to_detours(detours, ph_settings, atoms, ph_outputs["atoms_to_calculate"], calc_dict, fw_settings, "ph")

    if "phonon3py" in outputs:
        ph3_fw_set = fw_settings.copy()
        ph3_outputs = outputs["phonon3py"]
        ph3_settings = func_kwargs["phonon3_settings"].copy()
        update_spec["metadata_ph3"] = ph3_outputs["metadata"]
        ph3_fw_set["metadata_spec"] = "metadata_ph3"
        calc_dict = calc2dict(ph3_outputs["calculator"])
        if ph_settings["serial"]:
            update_spec["ph3_calculated_atoms"] = [atoms2dict(at) for at in ph3_outputs["atoms_to_calculate"]]
            update_spec["ph3_calculator"] = calc_dict
            ph3_fw_set["calc_atoms_spec"] = "ph_calculated_atoms"
            ph3_fw_set["calc_spec"] = "ph3_calculator"

            calc_kwargs = {}
            calc_keys = ["trajectory", "workdir", "backup_folder"]
            for key in calc_keys:
                if key in ph3_settings:
                    calc_kwargs[key] = ph3_settings[key]
            fw = generate_firework(
                func="hilde.tasks.calculate.calculate_socket",
                func_fw_out="hilde.tasks.fireworks.fw_action_outs.serial_calc",
                func_kwargs=calc_kwargs,
                atoms_calc_from_spec=False,
                inputs=["ph3_calculated_atoms", "ph3_calculator", "metadata_ph3"],
                fw_settings=ph3_fw_set,
            )
            detour.append(fw)
        else:
            fw_settings["mod_spec_add"] = "ph3_forces"
            detours = add_to_detours(detours, ph3_settings, atoms, ph3_outputs["atoms_to_calculate"], calc_dict, fw_settings, "ph3")

    return FWAction(update_spec=update_spec, detours=detours)

def add_phonon_force_calcs(
    atoms, calc, outputs, func, func_fw_out, func_kwargs, func_fw_kwargs, fw_settings
):
    """
    A function that initializes and adds all phonopy force calculations to the FireWorks Workflow,
    and adds the Phonopy Object to the spec
    Args:
        atoms: ASE Atoms object
            The original atoms at the start of this job
        calc: ASE Calculator object
            The original calculator
        outputs:
            The outputs from the function (assumes to be a single bool output)
        func: str
            Path to function that performs the MD like operation
        func_fw_out: str
            Path to this function
        func_kwargs: dict
            keyword arguments for func
        fw_settings: dict
            FireWorks specific settings
    """
    detours = []
    update_spec = {}
    fw_settings = fw_settings.copy()
    if "spec" in fw_settings and "kpoint_density_spec" in fw_settings:
        fw_settings["spec"][fw_settings["kpoint_density_spec"]] = k2d(dict2atoms(atoms), calc["calculator_parameters"]["k_grid"])
    elif "kpoint_density_spec" in fw_settings:
        fw_settings["spec"] = {fw_settings["kpoint_density_spec"]: k2d(dict2atoms(atoms), calc["calculator_parameters"]["k_grid"])}
    if outputs[0]:
        update_spec["metadata_ph"] = outputs[0][4]
        calc_dict = calc2dict(outputs[0][0])
        fw_settings["mod_spec_add"] = "ph_forces"
        detours = add_to_detours(detours, func_fw_kwargs["phonopy_settings"], atoms, outputs[0][2], calc_dict, fw_settings, "ph")

    if outputs[1]:
        update_spec["metadata_ph3"] = outputs[1][4]
        calc_dict = calc2dict(outputs[1][0])
        fw_settings["mod_spec_add"] = "ph3_forces"
        detours = add_to_detours(detours, func_fw_kwargs["phono3py_settings"], atoms, outputs[1][2], calc_dict, fw_settings, "ph3")

    return FWAction(update_spec=update_spec, detours=detours)

def add_to_detours(detours, func_fw_kwargs, atoms, atoms_list, calc_dict, fw_settings, prefix):
    for i, sc in enumerate(atoms_list):
        if not sc:
            continue
        fw_settings=fw_settings.copy()
        sc.info["displacement_id"] = i
        sc_dict = atoms2dict(sc)
        for key, val in calc_dict.items():
            sc_dict[key] = val
        calc_kwargs = {"workdir": func_fw_kwargs["workdir"] + f"/{i:05d}"}
        fw_settings["fw_name"] = prefix + f"forces_{Symbols(atoms['numbers']).get_chemical_formula()}_{i}"
        detours.append(
            generate_firework(
                func="hilde.tasks.calculate.calculate",
                func_fw_out="hilde.tasks.fireworks.fw_action_outs.mod_spec_add",
                func_kwargs=calc_kwargs,
                atoms=sc_dict,
                calc=calc_dict,
                atoms_calc_from_spec=False,
                fw_settings=fw_settings,
            )
        )
    return detours
