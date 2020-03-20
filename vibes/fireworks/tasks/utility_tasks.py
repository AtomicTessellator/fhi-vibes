"""Utility functions used in for HiLDe"""
import numpy as np
from fireworks import FWAction

from vibes.helpers.k_grid import update_k_grid_calc_dict


def mod_calc(param_key, calc_spec, calculator, val, atoms=None, spec_key=None):
    """Function to modify a calculator within the MongoDB

    Parameters
    ----------
    param_key: str
        key in the calculator dictionary to change
    calc_spec: str
        key for the calculator spec
    calculator: dict
        a dict representing an ASE Calculator
    val: the
        w value calculator[param_key] should be updated to
    atoms: dict
        A dict representing an ASE Atoms object
    spec_key: str
        The key in the MongoDB to update the val

    Returns
    -------
    FWAction
        An FWAction that modifies the calculator inside the spec
    """
    if param_key == "command":
        calculator[param_key] = val
    elif param_key == "basisset_type":
        sd = calculator["calculator_parameters"]["species_dir"].split("/")
        sd[-1] = val
        calculator["calculator_parameters"]["species_dir"] = "/".join(sd)
    elif param_key == "k_grid_density":
        recipcell = np.linalg.pinv(atoms["cell"]).transpose()
        update_k_grid_calc_dict(calculator, recipcell, val)
    else:
        calculator["calculator_parameters"][param_key] = val
    up_spec = {calc_spec: calculator}
    if spec_key:
        up_spec[spec_key] = val
    return FWAction(update_spec=up_spec)


def update_calc(calculator_dict, key, val):
    """Update the calculator dictionary

    Parameters
    ----------
    calculator_dict: dict
        The dictionary representation of the ASE Calculator
    key: str
        The key string of the parameter to be changed
    val: The
        dated value associated with the key string

    Returns
    -------
    dict
        The updated clac_dict
    """
    if key == "command":
        calculator_dict[key] = val
    elif key == "basisset_type":
        sd = calculator_dict["calculator_parameters"]["species_dir"].split("/")
        sd[-1] = val
        calculator_dict["calculator_parameters"]["species_dir"] = "/".join(sd)
    elif key == "use_pimd_wrapper" and isinstance(val, int):
        calculator_dict["calculator_parameters"][key] = ("localhost", val)
    else:
        if val is None and key in calculator_dict["calculator_parameters"]:
            del calculator_dict["calculator_parameters"][key]
        elif val is not None:
            calculator_dict["calculator_parameters"][key] = val
    return calculator_dict


def update_calc_in_db(calc_spec, update_calc_params, calculator):
    """Updates a calculator in the MongoDB with a new set of parameters

    Parameters
    ----------
    calc_spec: str
        spec to store the new calculator
    update_calc_params: dict
        A dictionary describing the new parameters to update the calculator with
    calculator: dict
        A dict representing an ASE Calculator

    Returns
    -------
    FWAction
        An FWAction that updates the calculator in the spec
    """
    del_key_list = ["relax_geometry", "relax_unit_cell"]
    for key in del_key_list:
        if key in calculator["calculator_parameters"]:
            del calculator["calculator_parameters"][key]
    for key, val in update_calc_params.items():
        calculator = update_calc(calculator, key, val)
    return FWAction(update_spec={calc_spec: calculator})
