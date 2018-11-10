""" tools for conerting atoms objects to json representations """


import json
from pathlib import Path
from ase.db.row import atoms2dict
from ase.io.jsonio import MyEncoder
from ase.calculators.calculator import all_properties


def get_json(obj):
    "Return json representation of obj"
    return json.dumps(obj, cls=MyEncoder, sort_keys=True)


def atoms2json(
    atoms, ignore_results=False, ignore_keys=["unique_id"], ignore_calc_params=[]
):
    """ return json representation of atoms and calculator objects.
        possibility to remove certain keys from the atoms dictionary, e.g. for hashing
    """

    # dictionary contains all the information in atoms object
    atomsdict = atoms2dict(atoms)

    # remove unwanted keys from atomsdict
    for name in ignore_keys:
        atomsdict.pop(name)

    calcdict = {}

    # move physical properties from atomsdict to calcdict if they are wanted
    for key in all_properties:
        if key in atomsdict:
            value = atomsdict.pop(key)
            if not ignore_results:
                calcdict[key] = value

    # clean calculator entries
    if "calculator_parameters" in atomsdict:
        calculator_params = atomsdict["calculator_parameters"]
        for name in [key for key in calculator_params  if key in ignore_calc_params]:
            calculator_params.pop(name)

        if "species_dir" in calculator_params:
            calculator_params["species_dir"] = Path(
                calculator_params["species_dir"]
            ).parts[-1]

    for name in ["calculator", "calculator_parameters"]:
        if name in atomsdict:
            calcdict[name] = atomsdict.pop(name)

    return get_json(atomsdict), get_json(calcdict)
