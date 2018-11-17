from ase.atoms import Atoms

from fireworks import Firework, PyTask

from hilde.helpers.converters import atoms2dict, dict2atoms, calc2dict
from hilde.tasks import fireworks as fw

module_name = __name__

def func_to_fire_works(func, func_fw_out, func_kwargs, atoms, calc, fw_settings, update_calc_settings=None):
    '''
    '''
    if not isinstance(atoms, dict):
        atoms_dict = atoms2dict(atoms)
    else:
        atoms_dict = atoms
    if not isinstance(calc, dict):
        calc_dict = calc2dict(calc)
    else:
        calc_dict = calc
    if update_calc_settings:
        for key, val in update_calc_settings.items():
            if val is None and key in calc_dict:
                    del(calc_dict[key])
            elif key in calc_dict:
                calc_dict[key] = val

    if "fw_name" not in fw_settings:
        fw_name = None
        fw_settings["fw_base_name"] = ""
    elif "fw_base_name" not in fw_settings:
        fw_settings["fw_base_name"] = fw_settings["fw_name"]

    ft_args = [
        func,
        func_fw_out,
        dict(func_kwargs),
        atoms_dict,
        calc_dict,
        dict(fw_settings),
    ]
    task = PyTask(
        {
            "func": fw.general_ase_calc_fxn_as_pytask.name,
            "args": ft_args,
        }
    )
    return Firework(task, name=fw_settings["fw_name"])

def get_func(func_path):
    toks = func_path.rsplit('.', 1)
    if len(toks) == 2:
        modname, funcname = toks
        mod = __import__(modname, globals(), locals(), [str(funcname)], 0)
        return getattr(mod, funcname)
    # Handle built in functions.
    return getattr(builtins, toks[0])

def general_ase_calc_fxn_as_pytask(func_path, func_fw_out_path, func_kwargs,  atoms_dict, calc, fw_settings=None):
    if fw_settings is None:
        fw_settings = {}
    func = get_func(func_path)
    func_fw_out = get_func(func_fw_out_path)

    for key, val in calc.items():
        atoms_dict[key] = val
    atoms = dict2atoms(atoms_dict)

    outputs = func(atoms, atoms.calc, **func_kwargs)

    return func_fw_out(atoms_dict, calc, outputs, func_path, func_fw_out_path, func_kwargs, fw_settings)

general_ase_calc_fxn_as_pytask.name = f"{module_name}.{general_ase_calc_fxn_as_pytask.__name__}"
