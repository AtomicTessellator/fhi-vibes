from fireworks import FWAction, PyTask, Firework
from hilde.helpers.k_grid import update_k_grid_calc_dict
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonon_db.phonon_db import connect
from hilde.parsers.structure import read_structure
from hilde.helpers.converters import atoms2dict, dict2atoms, calc2dict
from hilde.tasks import fireworks as fw

module_name = __name__


def mod_calc(param_key, calc_spec, calc, new_val, atoms=None, spec_key=None):
    """
    Function to modify a calculator within the MongoDB
    Args:
        param_key (str): key in the calculator dictionary to change
        calc (dict): a dict representing an ASE Calculator
        new_val: the new value calc[param_key] should be updated to
        atoms (dict): A dict representing an ASE Atoms object
        spec_key (str): The key in the MongoDB to update the new_val (used to pass the param down the Workflow)
    """
    if param_key is "command":
        calc[param_key] = new_val
    elif param_key == "basisset_type":
        sd = calc["calculator_parameters"]["species_dir"].split("/")
        sd[-1] = val
        calc["calculator_parameters"]["species_dir"] = "/".join(sd)
    elif param_key == "k_grid_density":
        recipcell = np.linalg.pinv(atoms["cell"]).transpose()
        update_k_grid_calc_dict(calc, recipcell, atoms["pbc"], new_val)
    else:
        calc["calculator_parameters"][param_key] = new_val
    up_spec = {calc_spec: calc}
    if spec_key:
        up_spec[spec_key] = new_val
    return FWAction(update_spec=up_spec)

def update_calc(calc_dict, key, val):
    if key == "command":
        calc_dict[key] = val
    elif key == "basisset_type":
        sd = calc_dict["calculator_parameters"]["species_dir"].split("/")
        sd[-1] = val
        calc_dict["calculator_parameters"]["species_dir"] = "/".join(sd)
    else:
        if val is None and key in calc_dict["calculator_parameters"]:
            del (calc_dict["calculator_parameters"][key])
        elif val is not None:
            calc_dict["calculator_parameters"][key] = val
    return calc_dict

def update_calc_in_db(calc_spec, update_calc_params, calc):
    """
    Updates a calculator in the MongoDB with a new set of parameters
    calc_spec (str): spec to store the new calculator
    update_calc_params (dict): A dictionary describing the new parameters to update the calc with
    calc (dict): A dict representing an ASE Calculator
    """
    del_key_list = []
    safe_key_list = ["k_grid", "species_dir", "command"]
    for key in calc["calculator_parameters"].keys():
        if key not in update_calc_params and key not in safe_key_list:
            del_key_list.append(key)
    for key in del_key_list:
        del (calc["calculator_parameters"][key])
    for key, val in update_calc_params.items():
        calc = update_calc(calc, key, val)
    return FWAction(update_spec={calc_spec: calc})


def add_result_to_spec(result_key, spec_key, atoms_calc):
    return FWAction(update_spec={spec_key: atoms_calc["results"][result_key]})


def transfer_spec(key, val):
    return FWAction(update_spec={key: val})


def check_convergence(
    prop_spec,
    atoms_spec,
    final_atoms_spec,
    loss_function,
    mutate_function,
    workdir,
    atoms_dicts,
    criteria=0.01,
):
    print(f"Checking convergence criteria, threshold value is {criteria}")
    loss_toks = loss_function.rsplit(".", 1)
    mut_toks = mutate_function.rsplit(".", 1)
    if len(loss_toks) == 2:
        modname, funcname = loss_toks
        mod = __import__(modname, globals(), locals(), [str(funcname)], 0)
        loss_func = getattr(mod, funcname)
    else:
        # Handle built in functions.
        loss_func = getattr(builtins, loss_toks[0])

    if len(mut_toks) == 2:
        modname, funcname = mut_toks
        mod = __import__(modname, globals(), locals(), [str(funcname)], 0)
        mut_func = getattr(mod, funcname)
    else:
        # Handle built in functions.
        mut_func = getattr(builtins, mut_toks[0])

    if loss_func(atoms_dicts[-1], atoms_dicts[-2], criteria=criteria):
        print("Convergence criteria met")
        _, cur_val = mut_func(atoms_dicts[-2])
        return FWAction(
            update_spec={prop_spec: cur_val, final_atoms_spec: atoms_dicts[-1]}
        )
    print("Convergence criteria not met. Mutating system and running another iteration")

    new_atoms, _ = mut_func(atoms_dicts[-1])
    new_workdir = (
        "/".join(workdir.split("/")[:-1]) + f'/{int(workdir.split("/")[-1])+1:05d}'
    )
    check_conv_args = [
        prop_spec,
        atoms_spec,
        final_atoms_spec,
        loss_function,
        mutate_function,
        new_workdir,
    ]
    task_list = [
        PyTask(
            {"func": fw.calculate.name, "args": [new_workdir, atoms_spec, new_atoms]}
        )
    ]
    task_list.append(
        PyTask(
            {
                "func": fw.check_convergence.name,
                "args": check_conv_args,
                "inputs": [atoms_spec],
                "kwargs": {"criteria": criteria},
            }
        )
    )
    new_firework = Firework(
        task_list, name="converging", spec={atoms_spec: atoms_dicts}
    )
    return FWAction(detours=[new_firework])


def get_relaxed_structure(new_struct_fname, out_atoms_spec, cur_atoms):
    try:
        new_atoms = read_structure(new_struct_fname)
        new_atoms.sym_block = cur_atoms["sym_block"]
    except:
        print("WARNING: new structure not found, using current atoms")
        new_atoms = dict2atoms(cur_atoms)
    return FWAction(update_spec={out_atoms_spec: atoms2dict(new_atoms)})


def add_phonon_to_db(
    db_path, atoms_ideal, phonon_dict, calc_type="calc", symprec=1e-5, **kwargs
):
    """
    Adds a phonon dictionary to a database defined by db_path
    Args:
        phonon: dict
            A dictionary representation of the phonopy object to be added to the database
        atoms_ideal: dict generated from atoms2dict
            The dictionary representation of the atoms or Atoms object of the undisplayed atoms
        db_path: str
            String to the database path
    """
    print(f"Adding phonon calculations to the database {db_path}")
    atoms = dict2atoms(atoms_ideal)
    atoms_hash, calc_hash = hash_atoms_and_calc(atoms)
    try:
        db = connect(db_path)
        selection = [
            ("symprec", "=", symprec),
            ("atoms_hash", "=", atoms_hash),
            ("calc_hash", "=", calc_hash),
            ("calc_type", "=", calc_type),
        ]
        if (kwargs is not None) and ("original_atoms_hash" in kwargs):
            selection.append(
                ("original_atoms_hash", "=", kwargs["original_atoms_hash"])
            )
        if (kwargs is not None) and ("supercell_matrix" in phonon_dict):
            selection.append(("supercell_matrix", "=", phonon_dict["supercell_matrix"]))
        try:
            rows = list(db.select(selection=selection))
            if not rows:
                raise KeyError
            for row in rows:
                db.update(
                    row.id,
                    dct=phonon_dict,
                    has_fc=("fc_2" in phonon_dict),
                    calc_type=calc_type,
                    **kwargs,
                )
        except KeyError:
            db.write(
                phonon_dict,
                symprec=symprec,
                atoms_hash=atoms_hash,
                calc_hash=calc_hash,
                has_fc2=("fc_2" in phonon_dict),
                has_fc3=("fc_3" in phonon_dict),
                calc_type=calc_type,
                **kwargs,
            )
    except ValueError:
        print(f"Fireworker could not access the database {db_path}")
    return FWAction(update_spec={"phonon_dict": {}})


mod_calc.name = f"{module_name}.{mod_calc.__name__}"
update_calc_in_db.name = f"{module_name}.{update_calc_in_db.__name__}"
add_result_to_spec.name = f"{module_name}.{add_result_to_spec.__name__}"
transfer_spec.name = f"{module_name}.{transfer_spec.__name__}"
check_convergence.name = f"{module_name}.{check_convergence.__name__}"
get_relaxed_structure.name = f"{module_name}.{get_relaxed_structure.__name__}"
add_phonon_to_db.name = f"{module_name}.{add_phonon_to_db.__name__}"
