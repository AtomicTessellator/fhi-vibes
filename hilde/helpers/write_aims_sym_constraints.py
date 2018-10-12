'''
A set of functionalities to write a symmetry blocks for aims calculations
'''
import os
import re
from pathlib import Path
from phonopy import Phonopy
import numpy as np

from hilde.structure.structure import pAtoms
from hilde.parsers.structure import read_aims

def write_sym_default(zero_pos, **kwargs):
    '''
    The default function for writing symmetry parameters
    Args:
        out_geo: file object
            The file that is being written to
        zero_pos: list or np.ndarray of floats (size = 3)
            The position of the atom if all symmetry parameters are 0
    Raises: TypeError:
        This function is a default that will be used for an atoms input struct,
        but it should never actually be called
    '''
    raise TypeError("Default write function is being used.")

def write_sym_xyz(zero_pos, op=["", "", ""], **kwargs):
    '''
    Writes a symmetry line in the out_geo file for an atom with symmetry defined in
    Cartesian coordinates
    Args:
        out_geo: file object
            The file that is being written to
        zero_pos: list or np.ndarray of floats (size = 3)
            The position of the atom if all symmetry parameters are 0
        op: list or np.ndarray of str (size=3)
            The list of symmetry operators for the (x, y, z) coordinates,
            an operator is {+,-,*}symmetry_parameter
    '''
    return(f"\nsymmetry_frac {zero_pos[0]:.5f}" + (f"{op[0]}" if op[0] != '' else "") +
                         f", {zero_pos[1]:.5f}" + (f"{op[1]}" if op[1] != '' else "") +
                         f", {zero_pos[2]:.5f}" + (f"{op[2]}" if op[2] != '' else ""))

def write_sym_ijk(zero_pos,
                  op=["", "", ""],
                  coord_i=0,
                  use_xyz=[True, True, True],
                  strict=True, **kwargs):
    '''
    Writes a symmetry line in the out_geo file for an atom with symmetry defined from an arbitrary
    axial vector i
    Args:
        out_geo: file object
            The file that is being written to
        zero_pos: list or np.ndarray of floats (size = 3)
            The position of the atom if all symmetry parameters are 0
        op: list or np.ndarray of str (size=3)
            The list of symmetry operators for the (i, j, k) coordinates, an operator is
            {+,-,*}symmetry_parameter
        coord_i: int
            The axis that is defined as i must be between 0 and 2
        use_xyz: list or np.ndarray of bools (size=3)
            A list of bools in the (x, y, z) axis where if false then symmetry operators do
            not act on that coordinate
        strict: bool
            If True: if coord_i is not in range(0,3) then raise error, else: coord_i = coord_i % 3
    '''
    if(coord_i < 0 or coord_i > 2):
        if not strict:
            print("WARNING: No coord_i value was not between 0 and 2, \
                setting coord_i to coord_i % 3")
            coord_i %= 3
        else:
            raise ValueError("i coordinate is outside the acceptable range of 0 to 2")
    else:
        return f"\nsymmetry_frac " + \
               f"{zero_pos[0]:.5f}" if(op[(0-coord_i)%3] == "" or not use_xyz[0]) else \
               op[(0-coord_i)%3] + ", " + \
               f"{zero_pos[1]:.5f}" if(op[(1-coord_i)%3] == "" or not use_xyz[1]) else \
               op[(1-coord_i)%3] + ", " + \
               f"{zero_pos[2]:.5f}" if(op[(2-coord_i)%3] == "" or not use_xyz[2]) else \
               op[(2-coord_i)%3]

def write_sym_fixed(zero_pos, **kwargs):
    '''
    Writes a symmetry line in the out_geo file for an atom with where the fractional coordinates
    are fixed at zero_pos
    Args:
        out_geo: file object
            The file that is being written to
        zero_pos: list or np.ndarray of floats (size = 3)
            The position of the atom if all symmetry parameters are 0
    Raises: TypeError:
        This function is a default that will be used for an atoms input struct, but it should never
        actually be called
    '''
    return f"\nsymmetry_frac {zero_pos[0]:.5f}, {zero_pos[1]:.5f}, {zero_pos[2]:.5f}"

def distort_atom_default(pos, **kwargs):
    '''
    The default function for writing an atomic position that copies the position from an input file
    Args:
        out_geo: file object
            The file that is being written to
        atom: str
            The symbol of the atom
        pos: list or np.ndarray of floats (size = 3)
            The position of the atom
    '''
    return pos

def distort_atom_xyz(pos, distortion=[0.0, 0.0, 0.0], **kwargs):
    '''
    The default function for writing an atomic position that copies the position from an
    input file
    Args:
        out_geo: file object
            The file that is being written to
        atom: str
            The symbol of the atom
        pos: list or np.ndarray of floats (size = 3)
            The position of the atom
        distortion: list or np.ndarray of floats (size = 3)
            A list describing how to distort the atoms in the (x, y, z) directions
    '''
    return np.array([pos[ii] + distortion[ii] for ii in range(3)])

def distort_atom_ijk(pos, distortion=[0.0, 0.0, 0.0], coord_i=0, **kwargs):
    '''
    The default function for writing an atomic position that copies the position from an
    input file
    Args:
        out_geo: file object
            The file that is being written to
        atom: str
            The symbol of the atom
        pos: list or np.ndarray of floats (size = 3)
            The position of the atom
        distortion: list or np.ndarray of floats (size = 3)
            A list describing how to distort the atoms in the (i, j, k) directions
        coord_i: int
            The axis that is defined as i must be between 0 and 2
    '''
    return np.array([pos[ii]+distortion[(ii-coord_i)%3] for ii in range(3)])

def get_coord_i_default(coords, **kwargs):
    '''
    The function for finding the i coordinate for a vector, that only is used in the
    Cartesian coordinate system (default behavior)
    Args:
        coords: list or np.ndarray of floats (size = 3)
            The vector that you want to find the i coordinate of
    Returns: int
        -1 since for xyz objects the coord_i value should not be used
    '''
    return 0

def get_coord_i_from_value(coords, val=0.0, error_tolerance=1e-10, **kwargs):
    '''
    The function for finding the i coordinate for a vector, that looks to find a
    certain value in the vector
    Args:
        coords: list or np.ndarray of floats (size = 3)
            The vector that you want to find the i coordinate of
        val: float
            The value in the vector that sets the i coordinate
        error_tolerance: float
            The allowed error in finding the correct value (used in cases of
            distorted structures)
    Returns: int
        The coordinate whose value matches val within the error tolerance
    '''
    # Check where the atomic position matches the input variable within some buffered region
    check_coords = np.where(np.abs(coords - np.floor(coords) - val) < error_tolerance)[0]
    assert len(check_coords) == 1
    return check_coords[0]

def get_type_default(coords, distortion):
    '''
    Gets the atom type for a given atom species
    Args:
        coords: list or np.ndarray of floats (size=3)
            The atomic coordinates
        distortion: list or np.ndarray of floats (size=3)
            The distortion of the current atom type
    Returns:
        type_num: int
            The type of the atom (default is always 0)
    '''
    return 0

class AtomInputs:
    '''
    Storage class for atom type parameters
    name: str
        {2 digit atomic number}_{atomic symbol}
    zero_position: list or np.ndarray of str (size 3)
        The position of the atom in all symmetry parameters are 0
        if @xx, @yy, @zz use atomic cooridnates
    sym_ops: list or np.ndarray of strs (size=3)
        A list of base symmetry operators for the atom type
    type_num: int
        An optional type parameter if there are multiple types of the same atom within
        a cell
    num_in_cell: int
        Total number of atoms with a given name and type num in the cell
    distort_atom: function
        A function that can be used to write a line in the atomic position portion of
        the geometry.in file
    write_sym: function
        A function that can be used to write a line in the geometry.in.sym file
        (only symmetry parameters)
    get_coord_i: function
        A function that can get the i coordinate of an atomic position vector
    get_type: function
        A function that can get the type number of the atom
    distortion_atomlist or np.ndarray of floats (size=3)
        A vector describing how to distort the zero_positions for each atom type when
        writing the atomic positions
    distortion_sym: list or np.ndarray of floats (size=3)
        A vector describing how to distort the zero_positions for each atom type when
        writing symmetry parameters
    use_xyz: list or np.ndarray of bools(size=3)
        A list of bools in the (x, y, z) axis where if false then symmetry operators do
        not act on that coordinate (used for i,j,k representationss)
    '''
    __slots__ = ["name",
        "zero_position",
        "sym_params",
        "type_num",
        "num_in_cell",
        "distort_atom",
        "write_sym",
        "get_coord_i",
        "get_type",
        "distortion_atom",
        "distortion_sym",
        "use_xyz",
        "coord_i_val",
        "coord_i_et"]
    def __init__(self,
                 name,
                 zero_position=["@xx", "@xx", "@xx"],
                 sym_params=["", "", ""],
                 type_num=0,
                 num_in_cell=0,
                 distort_atom=distort_atom_default,
                 write_sym=write_sym_default,
                 get_coord_i=get_coord_i_default,
                 get_type=get_type_default,
                 distortion_atom=[0.0, 0.0, 0.0],
                 distortion_sym=[0.0, 0.0, 0.0],
                 use_xyz=[True, True, True],
                 coord_i_val=0.0,
                 coord_i_et=1e-10):
        self.name = name
        self.zero_position = np.array([str(ii) for ii in zero_position])
        self.sym_params = np.array(sym_params)
        self.type_num = type_num
        self.num_in_cell = num_in_cell
        self.distort_atom = distort_atom
        self.write_sym = write_sym
        self.get_coord_i = get_coord_i
        self.get_type = get_type
        self.distortion_atom = np.array(distortion_atom)
        self.distortion_sym = np.array(distortion_sym)
        self.use_xyz = use_xyz
        self.coord_i_val = coord_i_val
        self.coord_i_et = coord_i_et

def write_sym_constraints_geo(in_file,
                              out_file,
                              out_file_sym,
                              atom_list,
                              symmetry_lv,
                              lat_param_list,
                              atom_param_list=[],
                              smatrix=[1, 0, 0, 0, 1, 0, 0, 0, 1],
                              **kwargs):
    '''
    Take in a base geometry.in file and generate a distorted geometry.in and geometry.in.sym
    file describing the symmetry parameters, and appends geometry.in.sym to geometry.in
    in_file: str
        File name for the input geometry.in file (should be generated by ASE)
    out_file: str
        File name for the distorted geometry.in file
    out_file_sym: str
        File name for the geometry.in.sym file describing the symmetry operations for
        the system
    atom_list: list of AtomInputs
        A list describing all atom types in the system
    symmetry_lv: np.ndarray of str (size=3,3)
        An array storing the lattice symmetry operators for the cell
    lat_param_list: list of str
        A list storing all lattice symmetry parameters in lat_param_list
    use_neg: bool
        If true allow negative atoms
    '''
    # Convert geometry file into a standardized version
    unitcell = read_aims(in_file)
    unitcell_scaled_pos = unitcell.get_scaled_positions()
    phonon = Phonopy(unitcell.to_phonopy_atoms(wrap=False),
                     supercell_matrix = smatrix,
                     symprec          = 1e-5,
                     is_symmetry      = True,
                     factor           = 15.633302,
                     log_level        = 0)
    supercell = pAtoms(phonopy_atoms=phonon.get_supercell())
    out_dir = Path(out_file).parent
    os.makedirs(out_dir, exist_ok=True)
    scaled_positions = supercell.get_scaled_positions()
    u2u_map = phonon.supercell.get_unitcell_to_unitcell_map()
    s2u_map = phonon.supercell.get_supercell_to_unitcell_map()
    sc_inds2uc_incd = [u2u_map[s2u_map[aa]] for aa in range(supercell.n_atoms)]
    coord_prims = unitcell_scaled_pos[sc_inds2uc_incd]
    for aa in range(len(scaled_positions)):
        for atom in atom_list:
            if atom.name[3:] == supercell.symbols[aa] and atom.get_type(coord_prims[aa], atom.distortion_sym) == atom.type_num:
                atom.num_in_cell += 1
                break
    # Generate the atom_param_list by parsing each atom types operator list
    for atom in atom_list:
        param_list = []
        for pp in atom.sym_params:
            param_list += re.split(r"\+|-|\*", pp)
        param_list = list(filter(''.__ne__, param_list))
        for param in param_list:
            if bool(re.search("@nn", param)):
                for ii in range(atom.num_in_cell):
                    if param.replace("@nn", str(ii)) not in atom_param_list:
                        atom_param_list.append(param.replace("@nn", str(ii)))
            elif param != "" and param not in atom_param_list:
                atom_param_list.append(param)
    # Create the geometry.in.sym file
    supercell.symmetry_block.append("\n# The Symmetry Paramters for the cell")
    supercell.symmetry_block.append("\n# format: symmetry_n_params [n n_lv n_fracpos]")
    supercell.symmetry_block.append(f"\nsymmetry_n_params {len(atom_param_list) + len(lat_param_list)} {len(lat_param_list)} {len(atom_param_list)}")
    supercell.symmetry_block.append("\nsymmetry_params")
    for param in lat_param_list+atom_param_list:
        supercell.symmetry_block.append(" " + param)
    # Write Lattice Parameters
    supercell.symmetry_block.append(f"\nsymmetry_lv {symmetry_lv[0][0]}, {symmetry_lv[0][1]}, {symmetry_lv[0][2]}")
    supercell.symmetry_block.append(f"\nsymmetry_lv {symmetry_lv[1][0]}, {symmetry_lv[1][1]}, {symmetry_lv[1][2]}")
    supercell.symmetry_block.append(f"\nsymmetry_lv {symmetry_lv[2][0]}, {symmetry_lv[2][1]}, {symmetry_lv[2][2]}")
    # Write atomic positions using atom specific functions (defined in atom_list)
    for aa in range(len(scaled_positions)):
        diff = supercell.positions[aa] - unitcell.positions[sc_inds2uc_incd[aa]]
        for atom in atom_list:
            if (supercell.symbols[aa] == atom.name[3:]) and (atom.get_type(coord_prims[aa], atom.distortion_sym) == atom.type_num):
                atom.num_in_cell -= 1
                # zero_pos = np.zeros(3)
                op = [atom.sym_params[ii].replace("@nn", str(atom.num_in_cell) ) for ii in range(3)]
                coord_i = atom.get_coord_i(coord_prims[aa],
                                           val=atom.coord_i_val,
                                           error_tolerance=atom.coord_i_et)
                zero_pos = np.array([float(atom.zero_position[(ii-coord_i)%3]) if atom.zero_position[(ii-coord_i)%3] != '@xx' else scaled_positions[aa][ii] for ii in range(3)])
                if np.any(atom.zero_position != "@xx"):
                    inds = (np.where(atom.zero_position != "@xx")[0]+coord_i)%3
                    zero_pos[inds] = np.linalg.solve( supercell.get_cell().T, (np.dot(zero_pos, unitcell.get_cell()) + diff).T).T[inds]
                    zero_pos[inds] -= np.floor(zero_pos[inds])
                supercell.symmetry_block.append(atom.write_sym(zero_pos,
                                                               op=op,
                                                               coord_i =coord_i,
                                                               use_xyz=atom.use_xyz,
                                                               strict=True,
                                                               **kwargs))
                supercell.positions[aa] = np.dot(atom.distort_atom(scaled_positions[aa],
                                                                   distortion=atom.distortion_atom,
                                                                   coord_i=coord_i),
                                                 supercell.get_cell())
    supercell.write(filename=out_file)
    with open(out_file_sym, 'w') as out_sym:
        for line in supercell.symmetry_block:
            out_sym.write(line)
