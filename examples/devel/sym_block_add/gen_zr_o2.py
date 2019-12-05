'''
Example of how to use the symmetry block generating function in helpers
'''
from argparse import ArgumentParser as argpars
from pathlib import Path
import numpy as np

from vibes.helpers.write_aims_sym_constraints import AtomInputs, distort_atom_default, write_sym_constraints_geo, write_sym_xyz, get_coord_i_default

parser = argpars(description='What to do')
parser.add_argument('--in_file', type=str, default="", help='input geometry file')
parser.add_argument('--out_dir', type=str, default="", help='base output directory')
parser.add_argument('--smatrix',
                    type=int, nargs='+',
                    default=[1, 0, 0, 0, 1, 0, 0, 0, 1],
                    help='supercell lattice vectors')
args = parser.parse_args()

def get_zr_o2_type(coords, sym_param):
    '''
    Determines the type number for ZrO2
    Parameters:
        coords: np.ndarray(dtype=float, shape=3)
            Coordinates of the atom
        sym_param: np.ndarray(dtype=float, shape=3)
            current values of the symmetry operators in the (x,y,z) direction
    '''
    if np.max(np.abs(coords - sym_param)) < 1e-4:
        return 0
    if np.max(np.abs(np.array([coords[0]+sym_param[0]-1.0, coords[1]+sym_param[1]-1.0, coords[2]-sym_param[2]-0.5]))) < 1e-4:
        return 1
    if np.max(np.abs(np.array([coords[0]-sym_param[0]+0.5, coords[1]+sym_param[1]-1.0, coords[2]-sym_param[2]]))) < 1e-4:
        return 2
    if np.max(np.abs(np.array([coords[0]+sym_param[0]-1.5, coords[1]-sym_param[1], coords[2]-sym_param[2]-0.5]))) < 1e-4:
        return 3
    return -1

# Type A
# Atom lists
atom_list = []
atom_list.append(AtomInputs("40_Zr", [0.0, 0.0, 0.0], ["+xZr@nn", "+yZr@nn", "+zZr@nn"], 0, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.530, 0.267, 0.356]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("40_Zr", [0.0, 0.0, 0.5], ["-xZr@nn", "-yZr@nn", "+zZr@nn"], 1, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.530, 0.267, 0.356]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("40_Zr", [0.5, 0.0, 0.0], ["+xZr@nn", "-yZr@nn", "+zZr@nn"], 2, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.530, 0.267, 0.356]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("40_Zr", [0.5, 0.0, 0.5], ["-xZr@nn", "+yZr@nn", "+zZr@nn"], 3, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.530, 0.267, 0.356]), [True, True, True], 0.0, 1e-10))

atom_list.append(AtomInputs("08_O", [0.0, 0.0, 0.0], ["+xOI@nn", "+yOI@nn", "+zOI@nn"], 0, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.771, 0.537, 0.106]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("08_O", [0.0, 0.0, 0.5], ["-xOI@nn", "-yOI@nn", "+zOI@nn"], 1, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.771, 0.537, 0.106]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("08_O", [0.5, 0.0, 0.0], ["+xOI@nn", "-yOI@nn", "+zOI@nn"], 2, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.771, 0.537, 0.106]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("08_O", [0.5, 0.0, 0.5], ["-xOI@nn", "+yOI@nn", "+zOI@nn"], 3, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.771, 0.537, 0.106]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("08_O", [0.0, 0.0, 0.0], ["+xOII@nn", "+yOII@nn", "+zOII@nn"], 0, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.639, 0.068, 0.000]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("08_O", [0.0, 0.0, 0.5], ["-xOII@nn", "-yOII@nn", "+zOII@nn"], 1, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.639, 0.068, 0.000]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("08_O", [0.5, 0.0, 0.0], ["+xOII@nn", "-yOII@nn", "+zOII@nn"], 2, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.639, 0.068, 0.000]), [True, True, True], 0.0, 1e-10))
atom_list.append(AtomInputs("08_O", [0.5, 0.0, 0.5], ["-xOII@nn", "+yOII@nn", "+zOII@nn"], 3, 0, distort_atom_default, write_sym_xyz, get_coord_i_default, get_zr_o2_type, np.array([0.0, 0.0, 0.0]), np.array([0.639, 0.068, 0.000]), [True, True, True], 0.0, 1e-10))

symmetry_lv = np.array(["a", "0.0", "0.0", "0.0", "b", "0.0", "0.0", "0.0", "c"]).reshape(3, 3)
lat_param_list = ["a", "b", "c"]
direc = args.out_dir + f"/sc_{args.smatrix[0]}_{args.smatrix[1]}_{args.smatrix[2]}_{args.smatrix[3]}_{args.smatrix[4]}_{args.smatrix[5]}_{args.smatrix[6]}_{args.smatrix[7]}_{args.smatrix[8]}"
args.smatrix = np.array(args.smatrix).reshape(3, 3)
write_sym_constraints_geo(args.in_file,
                          Path(direc+"/geometry.in"),
                          Path(direc+"/geometry.in.sym"),
                          atom_list,
                          symmetry_lv,
                          lat_param_list,
                          args.smatrix)
