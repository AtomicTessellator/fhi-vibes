from ase import Atoms
from phonopy.structure.atoms import PhonopyAtoms

def ASE_to_phonopy_atoms(structure):
    phonopy_atoms= PhonopyAtoms(
        symbols   = structure.get_chemical_symbols(),
        cell      = structure.get_cell(),
        masses    = structure.get_masses(),
        positions = structure.get_positions(wrap=True))
    return phonopy_atoms

def phonopy_to_ASE_atoms(structure):
    ASE_atoms= Atoms(
        symbols   = structure.get_chemical_symbols(),
        cell      = structure.get_cell(),
        masses    = structure.get_masses(),
        positions = structure.get_positions(),
        pbc       = True
    )
    return ASE_atoms
