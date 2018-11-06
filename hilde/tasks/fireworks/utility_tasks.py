from fireworks import FWAction, PyTask, Firework
from hilde.helpers.hash import hash_atoms
from hilde.phonon_db.phonon_db import connect
from hilde.parsers.structure import read_structure
from hilde.structure.structure import patoms2dict, dict2patoms
from hilde.tasks import fireworks as fw
module_name = __name__

def mod_calc(param_key, calc, new_val, spec_key=None):
    if param_key is "command":
        calc[param_key] = new_val
    else:
        calc["calculator_parameters"][param_key] = new_val
    up_spec = {"calculator" : calc}
    if spec_key:
        up_spec[spec_key] = new_val
    return FWAction(update_spec=up_spec)

def add_result_to_spec(result_key, spec_key, atoms_calc):
    return FWAction(update_spec={spec_key:atoms_calc["results"][result_key]})

def transfer_spec(key, val):
    return FWAction(update_spec={key: val})

def check_convergence(prop_spec,
                      atoms_spec,
                      final_atoms_spec,
                      loss_function,
                      mutate_function,
                      workdir,
                      atoms_dicts,
                      criteria=0.01):
    loss_toks = loss_function.rsplit('.', 1)
    mut_toks = mutate_function.rsplit('.', 1)
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
        _, cur_val = mut_func(atoms_dicts[-2])
        return FWAction(update_spec={prop_spec: cur_val, final_atoms_spec: atoms_dicts[-1]})

    new_atoms, _ = mut_func(atoms_dicts[-1])
    new_workdir = "/".join(workdir.split("/")[:-1]) + f'/{int(workdir.split("/")[-1])+1:05d}'
    check_conv_args = [prop_spec,
                       atoms_spec,
                       final_atoms_spec,
                       loss_function,
                       mutate_function,
                       new_workdir]
    task_list = [PyTask({"func": fw.calculate.name,
                         "args": [new_workdir, atoms_spec, new_atoms]})]
    task_list.append(PyTask({"func": fw.check_convergence.name,
                             "args": check_conv_args,
                             "inputs": [atoms_spec],
                             "kwargs": {"criteria": criteria}}))
    new_firework = Firework(task_list, name="converging", spec={atoms_spec: atoms_dicts})
    return FWAction(detours=[new_firework])

def get_relaxed_structure(new_struct_fname, out_atoms_spec, cur_atoms):
    try:
        new_atoms = read_structure(new_struct_fname)
        new_atoms.sym_block = cur_atoms['sym_block']
    except:
        print("WARNING: new structure not found, using current atoms")
        new_atoms = dict2patoms(cur_atoms)
    return FWAction(update_spec={out_atoms_spec: patoms2dict(new_atoms)})

def add_phonon_to_db(db_path, atoms_ideal, phonon_dict,symprec=1e-5, **kwargs):
    """
    Adds a phonon dictionary to a database defined by db_path
    Args:
        phonon: dict
            A dictionary representation of the phonopy object to be added to the database
        atoms_ideal: dict generated from patoms2dict
            The dictionary representation of the atoms or pAtoms object of the undisplayed atoms
        db_path: str
            String to the database path
    """
    atoms = dict2patoms(atoms_ideal)
    atoms_hash, calc_hash = hash_atoms(atoms)
    try:
        db = connect(db_path)
        selection = [("supercell_matrix", "=", phonon_dict["supercell_matrix"]),
                     ("symprec", "=", symprec),
                     ("atoms_hash", "=", atoms_hash),
                     ("calc_hash", "=", calc_hash)]
        if **kwargs is not None and "original_atoms_hash" in **kwargs:
            selection.append(original_atoms_hash)
        try:
            rows = list(db.select(selection=selection))
            if not rows:
                raise KeyError
            for row in rows:
                db.update(row.id, phonon=phonon_dict, has_fc=("force_constants" in phonon_dict), **kwargs)
        except KeyError:
            db.write(phonon_dict,
                     symprec=symprec,
                     atoms_hash=atoms_hash,
                     calc_hash=calc_hash,
                     has_fc=("force_constants" in phonon_dict),
                     **kwargs)
    except ValueError:
        print(f"Fireworker could not access the database {db_path}")
    return FWAction(stored_data={'phonopy_calc': phonon_dict})

mod_calc.name = f'{module_name}.{mod_calc.__name__}'
add_result_to_spec.name = f'{module_name}.{add_result_to_spec.__name__}'
transfer_spec.name = f'{module_name}.{transfer_spec.__name__}'
check_convergence.name = f'{module_name}.{check_convergence.__name__}'
get_relaxed_structure.name = f'{module_name}.{get_relaxed_structure.__name__}'
add_phonon_to_db.name = f'{module_name}.{add_phonon_to_db.__name__}'
