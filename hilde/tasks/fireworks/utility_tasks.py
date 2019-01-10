from fireworks import FWAction, PyTask, Firework

import numpy as np

from hilde.helpers.k_grid import update_k_grid_calc_dict
from hilde.helpers.hash import hash_atoms_and_calc
from hilde.phonon_db.phonon_db import connect
# from hilde.parsers.structure import read_structure
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
    ''' update the calculator dictionary '''
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
    del_key_list = [
        "relax_geometry",
        "relax_unit_cell",
        "use_sym",
    ]
    for key in del_key_list:
        if key in calc["calculator_parameters"]:
            del (calc["calculator_parameters"][key])
    for key, val in update_calc_params.items():
        calc = update_calc(calc, key, val)
    return FWAction(update_spec={calc_spec: calc})


mod_calc.name = f"{module_name}.{mod_calc.__name__}"
update_calc_in_db.name = f"{module_name}.{update_calc_in_db.__name__}"
