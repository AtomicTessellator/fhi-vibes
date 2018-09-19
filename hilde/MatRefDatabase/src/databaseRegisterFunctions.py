# Import json for json string based encodings
import json
# Import modules from ASE to input/output atom objects
from ase import Atoms
from ase import Atom
from ase.constraints import dict2constraint
from ase.calculators.calculator import get_calculator, all_properties
from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import chemical_symbols, atomic_masses
from ase.io.jsonio import decode
from ase.utils import formula_metal, basestring

# Helpers for registering atoms/Fingerprints
from helpers.hash import atoms2json
from MaterialsFingerprints import MaterialsFingerprint
from ase.io.jsonio import decode

import numpy as np

atomSymbolList= [ 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
def convert_fingerprint(s):
    sSplit = s.split(b";")
    fingerprint = {}
    nbins = int(sSplit[0])
    isB = bool(sSplit[1])
    isE = bool(sSplit[2])
    ii = 3;
    while ii < len(sSplit):
        pt = sSplit[ii].decode("ascii")
        ii += 1
        fingerprint[pt] = np.ndarray(shape=(nbins,2), dtype=float)
        for jj in range(nbins):
            fingerprint[pt][jj,:] = list(map(float,sSplit[ii:ii+2]))
            ii += 2
    return MaterialsFingerprint(isE, isB, nbins=nbins, fp=fingerprint)

def adapt_fingerprint(fp):
    frmt = ("%i;%r;%r" % (fp.nbins, fp.isB, fp.isElec))
    for pt in fp.fingerprint:
        frmt += (";%s" % (pt))
        for ff in fp.fingerprint[pt]:
            frmt += ";%f;%f" % (ff[0], ff[1])
    return frmt.encode()

def adapt_atom(at):
    atomStr, calcStr = atoms2json(at)
    if(calcStr != '{}'):
        atomStr = atomStr[:-1] + ", " + calcStr[1:]
    return atomStr.encode()

def convert_atom(s):
    atomsDict = decode(s.decode())
    return dict2atoms(atomsDict)

def dict2atoms(dct):
    atomList = []
    if(len(dct["numbers"]) != len(dct['positions']) ):
        raise ValueError("Length of the numbers and positions for the atom dictionary must be the same")
    for aa in range(len(dct['numbers'])):
        atomList.append( Atom( symbol=atomSymbolList[dct['numbers'][aa]-1], position=dct['positions'][aa]  ) )
    atoms = Atoms(symbols=atomList)
    # atoms.numbers   = np.array(dct['numbers'])
    # atoms.positions = dct['positions']
    if 'pbc' in dct:
        atoms.pbc = dct['pbc']
        atoms.cell = dct['cell']
    if 'initial_magmoms' in dct:
        atoms.set_initial_magnetic_moments(dct['initial_magmoms'])
    if 'initial_charges' in dct:
        atoms.set_initial_charges(dct['initial_charges'])
    if 'masses' in dct:
        atoms.set_masses(dct['masses'])
    if 'tags' in dct:
        atoms.set_tags(dct['tags'])
    if 'momenta' in dct:
        atoms.set_momenta(dct['momenta'])
    if 'constraints' in dct:
        atoms.constraints = [dict2constraint(dt) for dt in dct["constraints"]]
    if 'calculator' in dct:
        atoms.set_calculator( get_calculator(dct['calculator'])() )
        # atoms.calc = get_calculator(dct['calculator'])
        atoms.calc.parameters = {}
        for dt in dct['calculator_parameters']:
            atoms.calc.parameters[dt] = dct["calculator_parameters"][dt]
        if 'energy' in dct:
            atoms.calc.results['energy'] = dct['energy']
        if 'forces' in dct:
            atoms.calc.results['forces'] = dct['forces']
        if 'stress' in dct:
            atoms.calc.results['stress'] = dct['stress']
        if 'dipole' in dct:
            atoms.calc.results['dipole'] = dct['dipole']
        if 'charges' in dct:
            atoms.calc.results['charges'] = dct['charges']
        if 'magmom' in dct:
            atoms.calc.results['magmom'] = dct['magmom']
        if 'magmoms' in dct:
            atoms.calc.results['magmoms'] = dct['magmoms']
        if 'free_energy' in dct:
            atoms.calc.results['free_energy'] = dct['free_energy']
    return atoms