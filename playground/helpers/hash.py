from hashlib import sha1 as hash_sha
from ase.db.row import atoms2dict
from ase.io.jsonio import encode, MyEncoder
from ase.calculators.calculator import all_properties
from pathlib import Path
from configparser import ConfigParser
import json

def hashfunc(string, default=''):
    """ Wrap the sha hash function and check for empty objects """
    if (string == ''
        or string == '[]'
        or string == '{}'
        or string == 'None'):
        string = ''
    return hash_sha(string.encode('utf8'))

def atoms2json(atoms, ignore_atoms=[], ignore_calc=[]):
    atomsdict = atoms2dict(atoms)
    for name in ['unique_id'] + ignore_atoms:
        atomsdict.pop(name)
    calcdict = {}

    # clean calculator entries
    if 'calculator_parameters' in atomsdict:
        cp = atomsdict['calculator_parameters']
        for name in [k for k in cp if k in ignore_calc]:
            cp.pop(name)
        if 'species_dir' in cp:
            cp['species_dir'] = Path(cp['species_dir']).parts[-1]

    for name in all_properties:
        atomsdict.pop(name, None)
    for name in ['calculator', 'calculator_parameters']:
        if name in atomsdict:
            calcdict[name] = atomsdict.pop(name)
    getjson = lambda obj: json.dumps(obj, cls=MyEncoder, sort_keys=True)
    return getjson(atomsdict), getjson(calcdict)


def hash_atoms(atoms, ignore_atoms=[], ignore_calc=[], ignore_file=None):
    """ Hash atoms and calculator object, with possible ignores"""
    if ignore_file is not None:
        fil = Path(ignore_file)
        if fil.exists():
            cg = ConfigParser()
            cg.read(fil)
            ignore_atoms = [k for k in cg['atoms'] if not cg.getboolean('atoms', k)]
            ignore_calc = [k for k in cg['calculator'] if not cg.getboolean('calculator', k)]

    atomsjson, calcjson = atoms2json(atoms, ignore_atoms, ignore_calc)
    atomshash = hashfunc(atomsjson).hexdigest()
    calchash = hashfunc(calcjson).hexdigest()
    return atomshash, calchash
