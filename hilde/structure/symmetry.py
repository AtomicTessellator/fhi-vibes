from subprocess import Popen, PIPE, STDOUT
import json
from sys import exit
import numpy as np
from numpy import sin, cos, sqrt
import scipy.linalg as la
import spglib as spg
from hilde.helpers.cell import cell_to_cellpar, reciprocal_lattice
from hilde.helpers.maths import clean_matrix
from hilde.konstanten.numerics import loose_tol, wrap_tol, eps
from hilde.konstanten.symmetry import symprec
from .io import inform
from .misc import generate_lattice


class SymmetryOperation:
    """ Class to represent a symmetry operation"""
    def __init__(self, rotation, translation):
        self.rotation = rotation
        self.translation = translation

class Spacegroup:
    def __init__(self, atoms, symprec=None, mode=0, devel=False):
        """Should be similar to the symmetry dataset from spglib

        Keys (P = property):
            number (P)
            hall_number (P)
            international (P)
            hall (P)
            choice (P)
            transformation_matrix (P)
            origin_shift (P)
            rotations
            translations
            wyckoffs (P)
            equivalent_atoms (P)
            mapping_to_primitive (P)
            std_lattice
            std_types
            std_positions
            std_mapping_to_primitive
            pointgroup

        So: wrap the dictionaries obtained from spglib or aflow
        """
        # Initialize datasets
        self.aflow_sgdata = None
        self.aflow_edata = None
        self.aflow_dataset = None
        self._spglib_dataset = None

        # Attributes: moved to properties
        self.aflow_std_lattice = None
        self.aflow_std_positions = None
        self.aflow_spglib_matchlist = None
        self._site_symmetries = None
        # fmap representations of symmetry operations
        self.fmaps = None
        # set up
        self.atoms = atoms
        self.symprec = symprec
        self.setup(mode=mode, devel=devel)

    def setup(self, mode=0, devel=False):
        """
        :param mode: 0 for spglib dataset
                     1 for aflow + spglib
                     2 for full aflow analysis + spglib + sanity check
        """

        # Always perform spglib (if not done yet)
        if self._spglib_dataset is None:
            self._spglib_dataset = spg.get_symmetry_dataset(
                self.atoms.to_spglib_cell(), symprec=self.symprec)

            self.lat = self.atoms.get_cell()
            self.ilat = la.inv(self.lat)

        if mode == 0:
            # set mode
            self.mode = 0
        #
        elif mode == 1:
            self.mode = 1
            if self.aflow_sgdata is None:
                self._set_aflow_sgdata()
            if self.aflow_dataset is None:
                self._set_aflow_dataset()

            # Compare spglib and aflow
            aflow_spacegroup = int(self.aflow_sgdata['space_group_number'])
            if aflow_spacegroup != self.number:
                print(f'*Spacegroup number from spglib and aflow do not coincide!')
                print(f'  spglib: {self.number}')
                print(f'  aflow:  {aflow_spacegroup}')
                exit('Break.')
        #
        elif mode > 1 :
            if mode > 2:
                print(f'*Currently implemented symmetry modes: 0(spg), 1(aflow+spg),'
                      f'2(full aflow + spg.')
                print(f'*  fall back from {mode} to 1')

            self.mode = 2

            if self.aflow_sgdata is None:
                self._set_aflow_sgdata()
            if self.aflow_dataset is None:
                self._set_aflow_dataset()
            if self.aflow_edata is None:
                self._set_aflow_edata()
            #
            # Compare spglib and aflow
            aflow_spacegroup = int(self.aflow_edata['space_group_number'])
            tolerance = float(self.aflow_edata['wyccar']['title'].split('sym_eps:')[1].split()[0])
            if aflow_spacegroup != self.number:
                print(f'*Spacegroup number from spglib and aflow do not coincide!')
                print(f'  spglib: {self.number}')
                print(f'  aflow:  {aflow_spacegroup}')
                print(f'Aflow with self-consistent tolerance has been used.')
                print(f'  Tolerance:   {tolerance}')
                exit('Break.')
        #
        # if aflow_dataset has been computed, compare the symmetry operations:
        if self.aflow_dataset:
            fgroup = self.aflow_dataset['fgroup']
            aflow_frac_rotations    = np.array([np.array(el['Uf'])   for el in fgroup])
            aflow_frac_translations = np.array([np.array(el['ftau']) for el in fgroup])
            aflow_cart_rotations    = np.array([np.array(el['Uc'])   for el in fgroup])
            aflow_cart_translations = np.array([np.array(el['ctau']) for el in fgroup])

            #
            # Compare aflow and spglib fractional transformations:
            matchlist = []
            for ii in range(len(self.frac_rotations)):
                spglib_frot = self.frac_rotations[ii]
                # remove full lattice vectors
                spglib_ftau = (self.frac_translations[ii] + wrap_tol) % 1 - wrap_tol
                for jj in range(len(aflow_frac_rotations)):
                    aflow_frot = aflow_frac_rotations[jj]
                    aflow_ftau = (aflow_frac_translations[jj] + wrap_tol) % 1 - wrap_tol

                    # fkdev: must it be loose_tol? 1e-6 issue with aflowsym
                    if (la.norm(spglib_frot - aflow_frot) < loose_tol):
                        if (la.norm(spglib_ftau - aflow_ftau) < loose_tol):

                            # if found, append to the matchlist and leave the loop
                            matchlist.append(jj)
                            break
                        elif devel:
                            print(f'Translation does NOT match for pair ({ii}, {jj})')
                            print('spglib')
                            print(spglib_frot, spglib_ftau)
                            print('aflow')
                            print(aflow_frot, aflow_ftau)
                            print(self.frac_translations[ii])
                            print(aflow_frac_translations[jj] % 1)
                            print(f'norms: {la.norm(spglib_frot - aflow_frot)}  ' +
                                  f' {la.norm(spglib_ftau - aflow_ftau)}')
            #
            # Have all symmetry ops been found?
            if self.n_symops != len(matchlist):
                print(f'*Number of symops:                              {self.n_symops}')
                print(f'*Symops found by comparing spglib and aflowsym: {len(matchlist)}')
                exit('*Error in Spacegroup.setup() ')

    # Properties
    @property
    def number(self):
        return self._spglib_dataset['number']

    @property
    def hall_number(self):
        return self._spglib_dataset['hall_number']

    @property
    def international(self):
        return self._spglib_dataset['international']

    @property
    def hall(self):
        return self._spglib_dataset['hall']

    @property
    def choice(self):
        return self._spglib_dataset['choice']

    @property
    def transformation_matrix(self):
        return self._spglib_dataset['transformation_matrix']

    @property
    def origin_shift(self):
        return self._spglib_dataset['origin_shift']

    @property
    def frac_rotations(self):
        return self._spglib_dataset['rotations']

    @property
    def rotations(self):
        return [clean_matrix(self.lat.T @ frot @ self.ilat.T) for frot in self.frac_rotations]

    @property
    def frac_translations(self):
        """ return cleaned and wrapped fractional translations """
        frac_translations = self._spglib_dataset['translations']
        frac_translations += wrap_tol
        frac_translations = (frac_translations % 1 % 1 - wrap_tol)
        return clean_matrix(frac_translations)

    @property
    def translations(self):
        return clean_matrix([self.lat.T @ ft for ft in self.frac_translations])

    @property
    def wyckoffs(self):
        return self._spglib_dataset['wyckoffs']

    @property
    def wyckoffs_unique(self):
        uwcks, count = np.unique(self.wyckoffs, return_counts=True)
        return [(w, c) for (w, c) in zip(uwcks, count)]

    @property
    def equivalent_atoms(self):
        return self._spglib_dataset['equivalent_atoms']

    @property
    def equivalent_atoms_unique(self):
        ats, count = np.unique(self.equivalent_atoms, return_counts=True)
        return zip(ats, count)

    @property
    def mapping_to_primitive(self):
        return self._spglib_dataset['mapping_to_primitive']

    @property
    def spglib_std_lattice(self):
        return self._spglib_dataset['std_lattice']

    @property
    def spglib_std_positions(self):
        return self._spglib_dataset['std_positions']

    @property
    def std_mapping_to_primitive(self):
        return self._spglib_dataset['std_mapping_to_primitive']

    @property
    def pointgroup(self):
        return self._spglib_dataset['pointgroup']

    @property
    def n_symops(self):
        return len(self.frac_rotations)

    def inform(self):
        inform(self.atoms, self)

    # Symmetry elements from spglib
    def get_site_symmetries(self):
        if self._site_symmetries is not None:
            if len(self.site_symmetries) == self.atoms.get_n_atoms():
                return self._site_symmetries
        #
        # if not yet calculated, do it
        fpos = self.atoms.get_scaled_positions()
        frot = self.frac_rotations
        ftrl = self.frac_translations
        lattice = self.atoms.get_cell()
        all_site_symmetries = []

        for iat, pos in enumerate(fpos):
            site_symmetries = []
            for ii, (fr, ft) in enumerate(zip(frot, ftrl)):
                rot_pos = fr @ pos + ft
                diff = pos - rot_pos
                diff -= np.rint(diff)
                diff = np.dot(diff, lattice)
                if la.norm(diff) < self.symprec:
                    site_symmetries.append(ii)
            all_site_symmetries.append(site_symmetries)

        self._site_symmetries = np.array(all_site_symmetries, dtype=int)
        return self._site_symmetries

    @property
    def symbol(self, format='Herman_Mauguin'):
        return self.get_symbol(format='Herman_Mauguin')

    def get_symbol(self, format='Herman_Mauguin'):
        self.setup(0)
        if format.lower() == 'herman_mauguin':
            return self.international
        elif format.lower() == 'hall':
            self.setup(1)
            return self.aflow_sgdata['space_group_Hall']
        elif format.lower() == 'schoenflies':
            self.setup(1)
            return self.aflow_sgdata['space_group_Schoenflies']
        else:
            print(f'*Available spacegroup symbols:')
            print(f'*  Herman_Maugin, Hall, Schoenflies')
            print(f'*  return Herman_Mauguin')
            return self.aflow_sgdata['space_group_Hermann_Mauguin']
    #
    # Wrap the aflow binary calls (maybe move to dedicated module at some point)
    def _set_aflow_sgdata(self):
        struct = self.atoms.get_string(decorated=False).encode()
        asym = Popen(['aflow', '--sgdata', '--screen_only', '--print=json'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        jsonin = asym.communicate(input=struct)[0]
        result = json.loads(jsonin)
        self.aflow_sgdata = result

    def _set_aflow_edata(self):
        struct = self.atoms.get_string(decorated=False).encode()
        asym = Popen(['aflow', '--edata', '--print=json'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        jsonin = asym.communicate(input=struct)[0]
        result = json.loads(jsonin)
        self.aflow_edata = result

    def _set_aflow_dataset(self):
        struct = self.atoms.get_string(decorated=False).encode()
        asym = Popen(['aflow', '--fullsym', '--screen_only', '--print=json'], shell=False, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        jsonin = asym.communicate(input=struct)[0]
        result = json.loads(jsonin)

        #> fkdev
        # print( [el['ftau'] for el in result['fgroup']])
        # with open('jsonin.json', 'bw') as f:
        #     f.write(jsonin)
        # exit()
        self.aflow_dataset = result

    def get_std_cell(self, typ='prim'):
        from hilde.structure import pAtoms
        self.setup(mode=2)
        #
        if 'prim' in typ.lower():
            struct = self.aflow_edata['standard_primitive_structure']
        elif 'conv' in typ.lower():
            struct = self.aflow_edata['standard_conventional_structure']
        else:
            exit(f'Cell type {typ} unknown (must be prim or conv)')
        #
        scale = struct['scale']
        atoms_list = struct['atoms']
        symbols = [atom['name'] for atom in atoms_list]

        # Lattice and positions
        lattice = scale * np.array(struct['lattice'])
        if struct['coordinates_type'] == 'direct':
            scaled_positions = [atom['position'] for atom in atoms_list]
            positions = None
        else:
            positions = [atom['position'] for atom in atoms_list]
            scaled_positions = None
        #
        multiplicity = struct['number_each_type']

        return lattice, positions, scaled_positions, symbols
    #
    def get_primitive_cell(self):
        self.setup(mode=2)
        return self.get_std_cell(typ='prim')
    #
    def get_conventional_cell(self):
        self.setup(mode=2)
        return self.get_std_cell(typ='conv')
    #
    def get_conventional_smatrix(self):
        self.setup(mode=2)
        """ C = S.P
        <=> S = C.P^-1 """
        plat = self.aflow_edata['standard_primitive_structure']['lattice']
        clat = self.aflow_edata['standard_conventional_structure']['lattice']

        smatrix = np.array(np.rint(clat @ la.inv(plat)), dtype=int)
        return smatrix

    def get_crystal_system(self):
        """
        Get the crystal system for the structure, e.g., (triclinic,
        orthorhombic, cubic, etc.).

        Returns:
            (str): Cell system for structure.
        """
        n = self.number

        f = lambda i, j: i <= n <= j
        cs = {"triclinic": (1, 2), "monoclinic": (3, 15),
              "orthorhombic": (16, 74), "tetragonal": (75, 142),
              "trigonal": (143, 167), "hexagonal": (168, 194),
              "cubic": (195, 230)}

        crystal_sytem = None

        for k, v in cs.items():
            if f(*v):
                crystal_sytem = k
                break
        return crystal_sytem

    def get_lattice_type(self):
        """
        Get the lattice for the structure, e.g., (triclinic,
        orthorhombic, cubic, etc.).This is the same than the
        crystal system with the exception of the hexagonal/rhombohedral
        lattice

        Returns:
            (str): Lattice type for structure.
        """
        n = self.number

        lattice_type = self.get_crystal_system()
        if n in [146, 148, 155, 160, 161, 166, 167]:
            return "rhombohedral"
        elif lattice_type == "trigonal":
            return "hexagonal"
        else:
            return lattice_type

    def is_cubic(self):
        """ is the structure cubic? """
        if self.number in range(195, 231):
            return True
        else:
            return False

    def get_kpath(self):
       lattice_type = self.get_lattice_type()
       spg_symbol = self.international
       # conv = self.get_conventional_standardized()
       # prim = self.get_primitive_standardized()

       if lattice_type == "cubic":
           if "P" in spg_symbol:
               name = "CUB"
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'X': np.array([0.0, 0.5, 0.0]),
                          'R': np.array([0.5, 0.5, 0.5]),
                          'M': np.array([0.5, 0.5, 0.0])}
               path = [["\\Gamma", "X", "M", "\\Gamma", "R", "X"], ["M", "R"]]
           elif "F" in spg_symbol:
               name = "FCC"
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'K': np.array([3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]),
                          'L': np.array([0.5, 0.5, 0.5]),
                          'U': np.array([5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]),
                          'W': np.array([0.5, 1.0 / 4.0, 3.0 / 4.0]),
                          'X': np.array([0.5, 0.0, 0.5])}
               path = [["\\Gamma", "X", "W", "K",
                        "\\Gamma", "L", "U", "W", "L", "K"], ["U", "X"]]
           elif "I" in spg_symbol:
               name = "BCC"
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'H': np.array([0.5, -0.5, 0.5]),
                          'P': np.array([0.25, 0.25, 0.25]),
                          'N': np.array([0.0, 0.0, 0.5])}
               path = [["\\Gamma", "H", "N", "\\Gamma", "P", "H"], ["P", "N"]]
           else:
               warn("Unexpected value for spg_symbol: %s" % spg_symbol)
       # end cubic k path
       elif lattice_type == "tetragonal":
           if "P" in spg_symbol:
               name = "TET"
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'A': np.array([0.5, 0.5, 0.5]),
                          'M': np.array([0.5, 0.5, 0.0]),
                          'R': np.array([0.0, 0.5, 0.5]),
                          'X': np.array([0.0, 0.5, 0.0]),
                          'Z': np.array([0.0, 0.0, 0.5])}
               path = [["\\Gamma", "X", "M", "\\Gamma", "Z", "R", "A", "Z"], ["X", "R"],
                       ["M", "A"]]
           elif "I" in spg_symbol:
               conv_landang = cell_to_cellpar(self.get_conventional_cell())
               a = conv_landang[0]
               c = conv_landang[2]
               if c < a:
                   name = "BCT1"
                   eta = (1 + c ** 2 / a ** 2) / 4.0
                   kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                              'M': np.array([-0.5, 0.5, 0.5]),
                              'N': np.array([0.0, 0.5, 0.0]),
                              'P': np.array([0.25, 0.25, 0.25]),
                              'X': np.array([0.0, 0.0, 0.5]),
                              'Z': np.array([eta, eta, -eta]),
                              'Z_1': np.array([-eta, 1 - eta, eta])}
                   path = [["\\Gamma", "X", "M", "\\Gamma", "Z", "P", "N", "Z_1", "M"],
                           ["X", "P"]]
               else:
                   name = "BCT2"
                   eta = (1 + a ** 2 / c ** 2) / 4.0
                   zeta = a ** 2 / (2 * c ** 2)
                   kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                              'N': np.array([0.0, 0.5, 0.0]),
                              'P': np.array([0.25, 0.25, 0.25]),
                              '\\Sigma': np.array([-eta, eta, eta]),
                              '\\Sigma_1': np.array([eta, 1 - eta, -eta]),
                              'X': np.array([0.0, 0.0, 0.5]),
                              'Y': np.array([-zeta, zeta, 0.5]),
                              'Y_1': np.array([0.5, 0.5, -zeta]),
                              'Z': np.array([0.5, 0.5, -0.5])}
                   path = [["\\Gamma", "X", "Y", "\\Sigma", "\\Gamma", "Z",
                            "\\Sigma_1", "N", "P", "Y_1", "Z"], ["X", "P"]]
           else:
               warn("Unexpected value for spg_symbol: %s" % spg_symbol)
       # end tetragonal k path
       elif lattice_type == "orthorhombic":
           conv_landang = cell_to_cellpar(self.get_conventional_cell())
           a = conv_landang[0]
           b = conv_landang[1]
           c = conv_landang[2]

           if "P" in spg_symbol:
               name = "ORC"
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'R': np.array([0.5, 0.5, 0.5]),
                          'S': np.array([0.5, 0.5, 0.0]),
                          'T': np.array([0.0, 0.5, 0.5]),
                          'U': np.array([0.5, 0.0, 0.5]),
                          'X': np.array([0.5, 0.0, 0.0]),
                          'Y': np.array([0.0, 0.5, 0.0]),
                          'Z': np.array([0.0, 0.0, 0.5])}
               path = [["\\Gamma", "X", "S", "Y", "\\Gamma",
                        "Z", "U", "R", "T", "Z"], ["Y", "T"], ["U", "X"], ["S", "R"]]
           elif "F" in spg_symbol:
               if 1 / a ** 2 > 1 / b ** 2 + 1 / c ** 2:
                   name = "ORCF1"
                   zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
                   eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4

                   kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                              'A': np.array([0.5, 0.5 + zeta, zeta]),
                              'A_1': np.array([0.5, 0.5 - zeta, 1 - zeta]),
                              'L': np.array([0.5, 0.5, 0.5]),
                              'T': np.array([1, 0.5, 0.5]),
                              'X': np.array([0.0, eta, eta]),
                              'X_1': np.array([1, 1 - eta, 1 - eta]),
                              'Y': np.array([0.5, 0.0, 0.5]),
                              'Z': np.array([0.5, 0.5, 0.0])}
                   path = [["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
                           ["T", "X_1"], ["X", "A", "Z"], ["L", "\\Gamma"]]
               elif 1 / a ** 2 < 1 / b ** 2 + 1 / c ** 2:
                   name = "ORCF2"
                   phi = (1 + c ** 2 / b ** 2 - c ** 2 / a ** 2) / 4
                   eta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
                   delta = (1 + b ** 2 / a ** 2 - b ** 2 / c ** 2) / 4
                   kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                              'C': np.array([0.5, 0.5 - eta, 1 - eta]),
                              'C_1': np.array([0.5, 0.5 + eta, eta]),
                              'D': np.array([0.5 - delta, 0.5, 1 - delta]),
                              'D_1': np.array([0.5 + delta, 0.5, delta]),
                              'L': np.array([0.5, 0.5, 0.5]),
                              'H': np.array([1 - phi, 0.5 - phi, 0.5]),
                              'H_1': np.array([phi, 0.5 + phi, 0.5]),
                              'X': np.array([0.0, 0.5, 0.5]),
                              'Y': np.array([0.5, 0.0, 0.5]),
                              'Z': np.array([0.5, 0.5, 0.0])}
                   path = [["\\Gamma", "Y", "C", "D", "X", "\\Gamma",
                            "Z", "D_1", "H", "C"], ["C_1", "Z"], ["X", "H_1"], ["H", "Y"],
                           ["L", "\\Gamma"]]
               else:
                   name = "ORCF3"
                   zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
                   eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4
                   kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                              'A': np.array([0.5, 0.5 + zeta, zeta]),
                              'A_1': np.array([0.5, 0.5 - zeta, 1 - zeta]),
                              'L': np.array([0.5, 0.5, 0.5]),
                              'T': np.array([1, 0.5, 0.5]),
                              'X': np.array([0.0, eta, eta]),
                              'X_1': np.array([1, 1 - eta, 1 - eta]),
                              'Y': np.array([0.5, 0.0, 0.5]),
                              'Z': np.array([0.5, 0.5, 0.0])}
                   path = [["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
                           ["X", "A", "Z"], ["L", "\\Gamma"]]
           elif "I" in spg_symbol:
               name = "ORCI"
               zeta = (1 + a ** 2 / c ** 2) / 4
               eta = (1 + b ** 2 / c ** 2) / 4
               delta = (b ** 2 - a ** 2) / (4 * c ** 2)
               mu = (a ** 2 + b ** 2) / (4 * c ** 2)
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'L': np.array([-mu, mu, 0.5 - delta]),
                          'L_1': np.array([mu, -mu, 0.5 + delta]),
                          'L_2': np.array([0.5 - delta, 0.5 + delta, -mu]),
                          'R': np.array([0.0, 0.5, 0.0]),
                          'S': np.array([0.5, 0.0, 0.0]),
                          'T': np.array([0.0, 0.0, 0.5]),
                          'W': np.array([0.25, 0.25, 0.25]),
                          'X': np.array([-zeta, zeta, zeta]),
                          'X_1': np.array([zeta, 1 - zeta, -zeta]),
                          'Y': np.array([eta, -eta, eta]),
                          'Y_1': np.array([1 - eta, eta, -eta]),
                          'Z': np.array([0.5, 0.5, -0.5])}
               path = [["\\Gamma", "X", "L", "T", "W", "R", "X_1", "Z",
                        "\\Gamma", "Y", "S", "W"], ["L_1", "Y"], ["Y_1", "Z"]]
           elif "C" in spg_symbol or "A" in spg_symbol:
               name = "ORCC"
               zeta = (1 + a ** 2 / b ** 2) / 4
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'A': np.array([zeta, zeta, 0.5]),
                          'A_1': np.array([-zeta, 1 - zeta, 0.5]),
                          'R': np.array([0.0, 0.5, 0.5]),
                          'S': np.array([0.0, 0.5, 0.0]),
                          'T': np.array([-0.5, 0.5, 0.5]),
                          'X': np.array([zeta, zeta, 0.0]),
                          'X_1': np.array([-zeta, 1 - zeta, 0.0]),
                          'Y': np.array([-0.5, 0.5, 0]),
                          'Z': np.array([0.0, 0.0, 0.5])}
               path = [["\\Gamma", "X", "S", "R", "A", "Z",
                        "\\Gamma", "Y", "X_1", "A_1", "T", "Y"], ["Z", "T"]]
           else:
               warn("Unexpected value for spg_symbol: %s" % spg_symbol)
       # end orthorhombic k path
       elif lattice_type == "hexagonal":
           name = "HEX"
           kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                      'A': np.array([0.0, 0.0, 0.5]),
                      'H': np.array([1.0 / 3.0, 1.0 / 3.0, 0.5]),
                      'K': np.array([1.0 / 3.0, 1.0 / 3.0, 0.0]),
                      'L': np.array([0.5, 0.0, 0.5]),
                      'M': np.array([0.5, 0.0, 0.0])}
           path = [["\\Gamma", "M", "K", "\\Gamma", "A", "L", "H", "A"], ["L", "M"],
                   ["K", "H"]]

       elif lattice_type == "rhombohedral":
           prim_landang = cell_to_cellpar(self.get_primitive_cell())
           alpha = prim_landang[3] * pi / 180
           if alpha < 90:
               name = "RHL1"
               eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
               nu = 3.0 / 4.0 - eta / 2.0
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'B': np.array([eta, 0.5, 1.0 - eta]),
                          'B_1': np.array([1.0 / 2.0, 1.0 - eta, eta - 1.0]),
                          'F': np.array([0.5, 0.5, 0.0]),
                          'L': np.array([0.5, 0.0, 0.0]),
                          'L_1': np.array([0.0, 0.0, -0.5]),
                          'P': np.array([eta, nu, nu]),
                          'P_1': np.array([1.0 - nu, 1.0 - nu, 1.0 - eta]),
                          'P_2': np.array([nu, nu, eta - 1.0]),
                          'Q': np.array([1.0 - nu, nu, 0.0]),
                          'X': np.array([nu, 0.0, -nu]),
                          'Z': np.array([0.5, 0.5, 0.5])}
               path = [["\\Gamma", "L", "B_1"], ["B", "Z", "\\Gamma", "X"],
                       ["Q", "F", "P_1", "Z"], ["L", "P"]]
           else:
               name = "RHL2"
               eta = 1 / (2 * tan(alpha / 2.0) ** 2)
               nu = 3.0 / 4.0 - eta / 2.0
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'F': np.array([0.5, -0.5, 0.0]),
                          'L': np.array([0.5, 0.0, 0.0]),
                          'P': np.array([1 - nu, -nu, 1 - nu]),
                          'P_1': np.array([nu, nu - 1.0, nu - 1.0]),
                          'Q': np.array([eta, eta, eta]),
                          'Q_1': np.array([1.0 - eta, -eta, -eta]),
                          'Z': np.array([0.5, -0.5, 0.5])}
               path = [["\\Gamma", "P", "Z", "Q", "\\Gamma",
                        "F", "P_1", "Q_1", "L", "Z"]]
       # end rhombohedral k path
       elif lattice_type == "monoclinic":
           conv_landang = cell_to_cellpar(self.get_conventional_cell())
           a, b, c = conv_landang[0:3]
           alpha = conv_landang[3] * pi / 180
           #beta = self._conv.lattice.lengths_and_angles[1][1]

           if "P" in spg_symbol:
               name = "MCL"
               eta = (1 - b * cos(alpha) / c) / (2 * sin(alpha) ** 2)
               nu = 0.5 - eta * c * cos(alpha) / b
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'A': np.array([0.5, 0.5, 0.0]),
                          'C': np.array([0.0, 0.5, 0.5]),
                          'D': np.array([0.5, 0.0, 0.5]),
                          'D_1': np.array([0.5, 0.5, -0.5]),
                          'E': np.array([0.5, 0.5, 0.5]),
                          'H': np.array([0.0, eta, 1.0 - nu]),
                          'H_1': np.array([0.0, 1.0 - eta, nu]),
                          'H_2': np.array([0.0, eta, -nu]),
                          'M': np.array([0.5, eta, 1.0 - nu]),
                          'M_1': np.array([0.5, 1 - eta, nu]),
                          'M_2': np.array([0.5, 1 - eta, nu]),
                          'X': np.array([0.0, 0.5, 0.0]),
                          'Y': np.array([0.0, 0.0, 0.5]),
                          'Y_1': np.array([0.0, 0.0, -0.5]),
                          'Z': np.array([0.5, 0.0, 0.0])}
               path = [["\\Gamma", "Y", "H", "C", "E", "M_1", "A", "X", "H_1"],
                       ["M", "D", "Z"], ["Y", "D"]]
           elif "C" in spg_symbol:
               rec_prim_cell = reciprocal_lattice(self.get_primitive_cell())
               rec_prim_landang = cell_to_cellpar(rec_prim_cell)
               kgamma = rec_prim_landang[5]
               if kgamma > 90:
                   name = "MCLC1"
                   zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
                   eta = 0.5 + 2 * zeta * c * cos(alpha) / b
                   psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
                   phi = psi + (0.75 - psi) * b * cos(alpha) / c
                   kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                              'N': np.array([0.5, 0.0, 0.0]),
                              'N_1': np.array([0.0, -0.5, 0.0]),
                              'F': np.array([1 - zeta, 1 - zeta, 1 - eta]),
                              'F_1': np.array([zeta, zeta, eta]),
                              'F_2': np.array([-zeta, -zeta, 1 - eta]),
                              #'F_3': np.array([1 - zeta, -zeta, 1 - eta]),
                              'I': np.array([phi, 1 - phi, 0.5]),
                              'I_1': np.array([1 - phi, phi - 1, 0.5]),
                              'L': np.array([0.5, 0.5, 0.5]),
                              'M': np.array([0.5, 0.0, 0.5]),
                              'X': np.array([1 - psi, psi - 1, 0.0]),
                              'X_1': np.array([psi, 1 - psi, 0.0]),
                              'X_2': np.array([psi - 1, -psi, 0.0]),
                              'Y': np.array([0.5, 0.5, 0.0]),
                              'Y_1': np.array([-0.5, -0.5, 0.0]),
                              'Z': np.array([0.0, 0.0, 0.5])}
                   path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "F_1"],
                           ["Y", "X_1"], ["X", "\\Gamma", "N"], ["M", "\\Gamma"]]
               if kgamma == 90:
                   name = "MCLC2"
                   zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
                   eta = 0.5 + 2 * zeta * c * cos(alpha) / b
                   psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
                   phi = psi + (0.75 - psi) * b * cos(alpha) / c
                   kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                              'N': np.array([0.5, 0.0, 0.0]),
                              'N_1': np.array([0.0, -0.5, 0.0]),
                              'F': np.array([1 - zeta, 1 - zeta, 1 - eta]),
                              'F_1': np.array([zeta, zeta, eta]),
                              'F_2': np.array([-zeta, -zeta, 1 - eta]),
                              'F_3': np.array([1 - zeta, -zeta, 1 - eta]),
                              'I': np.array([phi, 1 - phi, 0.5]),
                              'I_1': np.array([1 - phi, phi - 1, 0.5]),
                              'L': np.array([0.5, 0.5, 0.5]),
                              'M': np.array([0.5, 0.0, 0.5]),
                              'X': np.array([1 - psi, psi - 1, 0.0]),
                              'X_1': np.array([psi, 1 - psi, 0.0]),
                              'X_2': np.array([psi - 1, -psi, 0.0]),
                              'Y': np.array([0.5, 0.5, 0.0]),
                              'Y_1': np.array([-0.5, -0.5, 0.0]),
                              'Z': np.array([0.0, 0.0, 0.5])}
                   path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "F_1"],
                           ["N", "\\Gamma", "M"]]
               if kgamma < 90:
                   if b * cos(alpha * pi / 180) / c\
                           + b ** 2 * sin(alpha) ** 2 / a ** 2 < 1:
                       name = "MCLC3"
                       mu = (1 + b ** 2 / a ** 2) / 4.0
                       delta = b * c * cos(alpha) / (2 * a ** 2)
                       zeta = mu - 0.25 + (1 - b * cos(alpha) / c)\
                           / (4 * sin(alpha) ** 2)
                       eta = 0.5 + 2 * zexamples/fireworks/fireworks_tutorial/FireWorks_Tutorials.ipynbeta * c * cos(alpha) / b
                       phi = 1 + zeta - 2 * mu
                       psi = eta - 2 * delta
                       kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                                  'F': np.array([1 - phi, 1 - phi, 1 - psi]),
                                  'F_1': np.array([phi, phi - 1, psi]),
                                  'F_2': np.array([1 - phi, -phi, 1 - psi]),
                                  'H': np.array([zeta, zeta, eta]),
                                  'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                                  'H_2': np.array([-zeta, -zeta, 1 - eta]),
                                  'I': np.array([0.5, -0.5, 0.5]),
                                  'M': np.array([0.5, 0.0, 0.5]),
                                  'N': np.array([0.5, 0.0, 0.0]),
                                  'N_1': np.array([0.0, -0.5, 0.0]),
                                  'X': np.array([0.5, -0.5, 0.0]),
                                  'Y': np.array([mu, mu, delta]),
                                  'Y_1': np.array([1 - mu, -mu, -delta]),
                                  'Y_2': np.array([-mu, -mu, -delta]),
                                  'Y_3': np.array([mu, mu - 1, delta]),
                                  'Z': np.array([0.0, 0.0, 0.5])}
                       path = [["\\Gamma", "Y", "F", "H", "Z", "I", "F_1"],
                               ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
                   if b * cos(alpha * pi / 180) / c \
                           + b ** 2 * sin(alpha) ** 2 / a ** 2 == 1:
                       name = "MCLC4"
                       mu = (1 + b ** 2 / a ** 2) / 4.0
                       delta = b * c * cos(alpha) / (2 * a ** 2)
                       zeta = mu - 0.25 + (1 - b * cos(alpha) / c)\
                           / (4 * sin(alpha) ** 2)
                       eta = 0.5 + 2 * zeta * c * cos(alpha) / b
                       phi = 1 + zeta - 2 * mu
                       psi = eta - 2 * delta
                       kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                                  'F': np.array([1 - phi, 1 - phi, 1 - psi]),
                                  'F_1': np.array([phi, phi - 1, psi]),
                                  'F_2': np.array([1 - phi, -phi, 1 - psi]),
                                  'H': np.array([zeta, zeta, eta]),
                                  'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                                  'H_2': np.array([-zeta, -zeta, 1 - eta]),
                                  'I': np.array([0.5, -0.5, 0.5]),
                                  'M': np.array([0.5, 0.0, 0.5]),
                                  'N': np.array([0.5, 0.0, 0.0]),
                                  'N_1': np.array([0.0, -0.5, 0.0]),
                                  'X': np.array([0.5, -0.5, 0.0]),
                                  'Y': np.array([mu, mu, delta]),
                                  'Y_1': np.array([1 - mu, -mu, -delta]),
                                  'Y_2': np.array([-mu, -mu, -delta]),
                                  'Y_3': np.array([mu, mu - 1, delta]),
                                  'Z': np.array([0.0, 0.0, 0.5])}
                       path = [["\\Gamma", "Y", "F", "H", "Z", "I"],
                               ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
                   if b * cos(alpha * pi / 180) / c \
                           + b ** 2 * sin(alpha) ** 2 / a ** 2 > 1:
                       name = "MCLC5"
                       zeta = (b ** 2 / a ** 2 + (1 - b * cos(alpha) / c)
                               / sin(alpha) ** 2) / 4
                       eta = 0.5 + 2 * zeta * c * cos(alpha) / b
                       mu = eta / 2 + b ** 2 / (4 * a ** 2) \
                           - b * c * cos(alpha) / (2 * a ** 2)
                       nu = 2 * mu - zeta
                       rho = 1 - zeta * a ** 2 / b ** 2
                       omega = (4 * nu - 1 - b ** 2 * sin(alpha) ** 2 / a ** 2)\
                           * c / (2 * b * cos(alpha))
                       delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25
                       kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                                  'F': np.array([nu, nu, omega]),
                                  'F_1': np.array([1 - nu, 1 - nu, 1 - omega]),
                                  'F_2': np.array([nu, nu - 1, omega]),
                                  'H': np.array([zeta, zeta, eta]),
                                  'H_1': np.array([1 - zeta, -zeta, 1 - eta]),
                                  'H_2': np.array([-zeta, -zeta, 1 - eta]),
                                  'I': np.array([rho, 1 - rho, 0.5]),
                                  'I_1': np.array([1 - rho, rho - 1, 0.5]),
                                  'L': np.array([0.5, 0.5, 0.5]),
                                  'M': np.array([0.5, 0.0, 0.5]),
                                  'N': np.array([0.5, 0.0, 0.0]),
                                  'N_1': np.array([0.0, -0.5, 0.0]),
                                  'X': np.array([0.5, -0.5, 0.0]),
                                  'Y': np.array([mu, mu, delta]),
                                  'Y_1': np.array([1 - mu, -mu, -delta]),
                                  'Y_2': np.array([-mu, -mu, -delta]),
                                  'Y_3': np.array([mu, mu - 1, delta]),
                                  'Z': np.array([0.0, 0.0, 0.5])}
                       path = [["\\Gamma", "Y", "F", "L", "I"], ["I_1", "Z", "H", "F_1"],
                               ["H_1", "Y_1", "X", "\\Gamma", "N"], ["M", "\\Gamma"]]
           else:
               warn("Unexpected value for spg_symbol: %s" % spg_symbol)
       # end monoclinic k path
       elif lattice_type == "triclinic":
           rec_prim_cell = reciprocal_lattice(self.get_primitive_cell())
           rec_prim_landang = cell_to_cellpar(rec_prim_cell)
           kalpha = rec_prim_landang[3]
           kbeta = rec_prim_landang[4]
           kgamma = rec_prim_landang[5]

           if kalpha > 90 and kbeta > 90 and kgamma >= 90:
               name = "TRI1a"
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'L': np.array([0.5, 0.5, 0.0]),
                          'M': np.array([0.0, 0.5, 0.5]),
                          'N': np.array([0.5, 0.0, 0.5]),
                          'R': np.array([0.5, 0.5, 0.5]),
                          'X': np.array([0.5, 0.0, 0.0]),
                          'Y': np.array([0.0, 0.5, 0.0]),
                          'Z': np.array([0.0, 0.0, 0.5])}
               path = [["X", "\\Gamma", "Y"], ["L", "\\Gamma", "Z"],
                       ["N", "\\Gamma", "M"], ["R", "\\Gamma"]]
           if kalpha < 90 and kbeta < 90 and kgamma <= 90:
               name = "TRI1b"
               kpoints = {'\\Gamma': np.array([0.0, 0.0, 0.0]),
                          'L': np.array([0.5, -0.5, 0.0]),
                          'M': np.array([0.0, 0.0, 0.5]),
                          'N': np.array([-0.5, -0.5, 0.5]),
                          'R': np.array([0.0, -0.5, 0.5]),
                          'X': np.array([0.0, -0.5, 0.0]),
                          'Y': np.array([0.5, 0.0, 0.0]),
                          'Z': np.array([-0.5, 0.0, 0.5])}
               path = [["X", "\\Gamma", "Y"], ["L", "\\Gamma", "Z"],
                       ["N", "\\Gamma", "M"], ["R", "\\Gamma"]]
       # end triclinic k path
       else:
           warn("Unknown lattice type %s" % lattice_type)

       return {'name': name, 'kpoints': kpoints, 'path': path}
    # end get_kpath

    def refine(self, primitive=True):
        from hilde.structure import pAtoms
        lattice, scaled_positions, numbers = spg.standardize_cell(
            self.atoms, to_primitive=primitive, no_idealize=0,
            symprec=self.symprec)
        refined_cell = pAtoms(cell=lattice, scaled_positions=scaled_positions,
                              numbers=numbers, pbc=True, symprec=self.symprec)
        #
        refined_cell.wrap()
        return refined_cell

    def get_conventional_standardized(self, international_monoclinic=False):
        """ [pymatgen, adapted]
        Gives a structure with a conventional cell according to certain
        standards. The standards are defined in Setyawan, W., & Curtarolo,
        S. (2010). High-throughput electronic band structure calculations:
        Challenges and tools. Computational Materials Science,
        49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
        They basically enforce as much as possible
        norm(a1)<norm(a2)<norm(a3)

        Returns:
            The structure in a conventional standardized cell
        """
        from hilde.structure import pAtoms
        symprec = self.symprec
        struct = self.refine(primitive=False)
        latt_type = self.get_lattice_type()
        latt = struct.cell
        landang = struct.get_cell_lengths_and_angles()
        lengths = landang[0:3]
        sorted_lengths = sorted(lengths)
        sorted_dic = sorted([{'vec': latt[i],
                              'length': lengths[i],
                              'orig_index': i} for i in [0, 1, 2]],
                            key=lambda k: k['length'])

        if latt_type in ("orthorhombic", "cubic"):
            # you want to keep the c axis where it is
            # to keep the C- settings
            transf = np.zeros(shape=(3, 3))
            if self.symbol.startswith("C"):
                transf[2] = [0, 0, 1]
                a, b = sorted(lengths[:2])
                sorted_dic = sorted([{'vec': latt[i],
                                      'length': lengths[i],
                                      'orig_index': i} for i in [0, 1]],
                                    key=lambda k: k['length'])
                for i in range(2):
                    transf[i][sorted_dic[i]['orig_index']] = 1
                c = lengths[2]
            else:
                for i in range(len(sorted_dic)):
                    transf[i][sorted_dic[i]['orig_index']] = 1
                a, b, c = sorted_lengths
            latt = generate_lattice(a, b, c, lattice_type='orthorhombic')

        elif latt_type == "tetragonal":
            # find the "a" vectors
            # it is basically the vector repeated two times
            transf = np.zeros(shape=(3, 3))
            a, b, c = sorted_lengths
            for d in range(len(sorted_dic)):
                transf[d][sorted_dic[d]['orig_index']] = 1

            if abs(b - c) < loose_tol:
                a, c = c, a
                transf = np.dot([[0, 0, 1], [0, 1, 0], [1, 0, 0]], transf)
            latt = generate_lattice(a, c=c, lattice_type='tetragonal')

        elif latt_type in ("hexagonal", "rhombohedral"):
            # for the conventional cell representation,
            # we allways show the rhombohedral lattices as hexagonal

            # check first if we have the refined structure shows a rhombohedral
            # cell
            # if so, make a supercell
            a, b, c = lengths
            if np.all(np.abs([a - b, c - b, a - c]) < 0.001):
                superstruct = struct.make_supercell([[1, -1, 0], [0, 1, -1], [1, 1, 1]])
                a, b, c = sorted(superstruct.get_cell_lengths_and_angles()[0:3])

            if abs(b - c) < 0.001:
                a, c = c, a
            new_matrix = [[a / 2, -a * sqrt(3) / 2, 0],
                          [a / 2, a * sqrt(3) / 2, 0],
                          [0, 0, c]]
            latt = new_matrix
            transf = np.eye(3, 3)

        elif latt_type == 'monoclinic':
            # You want to keep the c axis where it is to keep the C- settings
            # a vector is crystal axis if alpha < 90 (Curtarolo definition)
            if self.symbol.startswith('C'):
                transf = np.zeros(shape=(3, 3))
                transf[2] = [0, 0, 1]
                sorted_dic = sorted([{'vec': latt[i],
                                      'length': lengths[i],
                                      'orig_index': i} for i in [0, 1]],
                                    key=lambda k: k['length'])
                a = sorted_dic[0]['length']
                b = sorted_dic[1]['length']
                c = lengths[2]
                new_matrix = None
                for t in itertools.permutations(list(range(2)), 2):
                    m = latt
                    landang = cell_to_cellpar([m[t[0]], m[t[1]], m[2]])
                    if landang[3] > 90:
                        # if the angle is > 90 we invert a and b to get
                        # an angle < 90
                        landang = cell_to_cellpar([-m[t[0]], -m[t[1]], m[2]])
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][2] = 1
                        a, b, c = landang[0:3]
                        alpha = pi * landang[3] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]
                        continue

                    elif landang[3] < 90:
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][2] = 1
                        a, b, c = landang[0:3]
                        alpha = pi * landang[3] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]

                if new_matrix is None:
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [[a, 0, 0],
                                  [0, b, 0],
                                  [0, 0, c]]
                    transf = np.zeros(shape=(3, 3))
                    for c in range(len(sorted_dic)):
                        transf[c][sorted_dic[c]['orig_index']] = 1
            #if not C-setting
            else:
                # try all permutations of the axis
                # keep the ones with the non-90 angle=alpha
                # and b<c
                new_matrix = None
                for t in itertools.permutations(list(range(3)), 3):
                    m = latt
                    landang = cell_to_cellpar([m[t[0]], m[t[1]], m[t[2]]])
                    if landang[3] > 90 and landang[1] < landang[2]:
                        landang = cell_to_cellpar([-m[t[0]], -m[t[1]], m[t[2]]])
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][t[2]] = 1
                        a, b, c = landang[0:3]
                        alpha = pi * landang[3] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]
                        continue
                    elif landang[3] < 90 and landang[1] < landang[2]:
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][t[2]] = 1
                        a, b, c = landang[0:3]
                        alpha = pi * landang[3] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]
                if new_matrix is None:
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [[sorted_lengths[0], 0, 0],
                                  [0, sorted_lengths[1], 0],
                                  [0, 0, sorted_lengths[2]]]
                    transf = np.zeros(shape=(3, 3))
                    for c in range(len(sorted_dic)):
                        transf[c][sorted_dic[c]['orig_index']] = 1

            if international_monoclinic:
                # The above code makes alpha the non-right angle.
                # The following will convert to proper international convention
                # that beta is the non-right angle.
                op = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
                transf = np.dot(op, transf)
                new_matrix = np.dot(op, new_matrix)
                beta = cell_to_cellpar(new_matrix)[4]
                if beta < 90:
                    op = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
                    transf = np.dot(op, transf)
                    new_matrix = np.dot(op, new_matrix)

            latt = new_matrix
        # end monoclinic

        elif latt_type == "triclinic":
            #we use a LLL Minkowski-like reduction for the triclinic cells
            reduced = struct.get_reduced_structure("LLL")
            cart_coord = struct.get_positions()
            struct = pAtoms(cell=reduced, positions=cart_coord,
                            numbers=struct.numbers, pbc=True,
                            symprec=self.symprec)

            a, b, c = landang[0:3]
            alpha, beta, gamma = [pi * i / 180
                                  for i in landang[3:6]]
            new_matrix = None
            test_matrix = [[a, 0, 0],
                          [b * cos(gamma), b * sin(gamma), 0.0],
                          [c * cos(beta),
                           c * (cos(alpha) - cos(beta) * cos(gamma)) /
                           sin(gamma),
                           c * sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                         - cos(beta) ** 2
                                         + 2 * cos(alpha) * cos(beta)
                                         * cos(gamma)) / sin(gamma)]]

            def is_all_acute_or_obtuse(m):
                recp_cell = reciprocal_lattice(m)
                recp_angles = cell_to_cellpar(recp_cell)[3:6]
                return np.all(recp_angles <= 90) or np.all(recp_angles > 90)

            if is_all_acute_or_obtuse(test_matrix):
                transf = np.eye(3)
                new_matrix = test_matrix

            test_matrix = [[-a, 0, 0],
                           [b * cos(gamma), b * sin(gamma), 0.0],
                           [-c * cos(beta),
                            -c * (cos(alpha) - cos(beta) * cos(gamma)) /
                            sin(gamma),
                            -c * sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                           - cos(beta) ** 2
                                           + 2 * cos(alpha) * cos(beta)
                                           * cos(gamma)) / sin(gamma)]]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0],
                          [0, 1, 0],
                          [0, 0, -1]]
                new_matrix = test_matrix

            test_matrix = [[-a, 0, 0],
                           [-b * cos(gamma), -b * sin(gamma), 0.0],
                           [c * cos(beta),
                            c * (cos(alpha) - cos(beta) * cos(gamma)) /
                            sin(gamma),
                            c * sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                          - cos(beta) ** 2
                                          + 2 * cos(alpha) * cos(beta)
                                          * cos(gamma)) / sin(gamma)]]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0],
                          [0, -1, 0],
                          [0, 0, 1]]
                new_matrix = test_matrix

            test_matrix = [[a, 0, 0],
                           [-b * cos(gamma), -b * sin(gamma), 0.0],
                           [-c * cos(beta),
                            -c * (cos(alpha) - cos(beta) * cos(gamma)) /
                            sin(gamma),
                            -c * sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                           - cos(beta) ** 2
                                           + 2 * cos(alpha) * cos(beta)
                                           * cos(gamma)) / sin(gamma)]]
            if is_all_acute_or_obtuse(test_matrix):
                transf = [[1, 0, 0],
                          [0, -1, 0],
                          [0, 0, -1]]
                new_matrix = test_matrix

            latt = new_matrix
        # end triclinic

        new_coords = np.dot(transf,
                            np.transpose(struct.get_scaled_positions())).T
        new_atoms = pAtoms(cell=latt,scaled_positions=new_coords,
                           numbers=struct.numbers, pbc=True,
                           symprec=self.symprec)
        return new_atoms
    # end get_conventional_standardized

    def get_primitive_standardized(self, international_monoclinic=False):
        """ [pymatgen, adapted]
        Gives a structure with a primitive cell according to certain standards
        the standards are defined in Setyawan, W., & Curtarolo, S. (2010).
        High-throughput electronic band structure calculations:
        Challenges and tools. Computational Materials Science,
        49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010

        Returns:
            The structure in a primitive standardized cell
        """
        from hilde.structure import pAtoms
        conv = self.get_conventional_standardized(
                international_monoclinic=international_monoclinic)
        lattice = self.get_lattice_type()

        if "P" in self.symbol or lattice == "hexagonal":
            return conv

        if lattice == "rhombohedral":
            # check if the conventional representation is hexagonal or
            # rhombohedral
            landang = conv.get_cell_lengths_and_angles()
            lengths = landang[0:3]
            angles = landang[3:6]
            if abs(lengths[0]-lengths[2]) < 0.0001:
                transf = np.eye
            else:
                transf = np.array([[-1, 1, 1], [2, 1, 1], [-1, -2, 1]],
                                  dtype=np.float) / 3
        elif "I" in self.symbol:
            transf = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]],
                              dtype=np.float) / 2
        elif "F" in self.symbol:
            transf = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]],
                              dtype=np.float) / 2
        elif "C" in self.symbol or "A" in self.symbol:
            if self.get_crystal_system() == "monoclinic":
                transf = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]],
                                  dtype=np.float) / 2
            else:
                transf = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]],
                                  dtype=np.float) / 2
        else:
            transf = np.eye(3)

        # Define two helper functions
        # fractional distance between two coordinates
        pbc_diff = lambda x,y: np.subtract(x,y) - np.round(np.subtract(x,y))
        # return true if arguments are periodic images of each other (just for scaled pos.)
        is_periodic_image = lambda x,y: np.allclose(pbc_diff(x,y), [0,0,0], atol=1e-8)

        new_scaled_positions = []
        new_numbers = []
        latt = np.dot(transf, conv.cell)
        for pos, num in zip(conv.get_positions(), conv.numbers):
            new_sc_pos = np.dot(np.linalg.inv(np.transpose(latt)),pos)
            new_sc_pos = np.mod(new_sc_pos,1)   # to unit cell
            if not any( [is_periodic_image(new_sc_pos,n_p) for n_p in new_scaled_positions] ):
                new_scaled_positions.append(new_sc_pos)
                new_numbers.append(num)

        new_atoms = pAtoms(cell=latt, scaled_positions=new_scaled_positions,
                           numbers=new_numbers, pbc=True,
                           symprec=self.symprec)

        if lattice == "rhombohedral":
            landang = new_atoms.get_cell_lengths_and_angles()
            a = landang[0]
            alpha = pi * landang[3] / 180
            latt = [[a * cos(alpha / 2), -a * sin(alpha / 2), 0],
                    [a * cos(alpha / 2), a * sin(alpha / 2), 0],
                    [a * cos(alpha) / cos(alpha / 2), 0,
                     a * sqrt(1 - (cos(alpha) ** 2 / (cos(alpha / 2) ** 2)))]]
            new_scaled_positions = []
            new_numbers = []
            for sc_pos, num in zip(new_atoms.get_scaled_positions(), new_atoms.numbers):
                pos = np.dot(np.transpose(new_atoms.cell),sc_pos)
                new_sc_pos = np.dot(np.linalg.inv(np.transpose(latt)),pos)
                new_sc_pos = np.mod(new_sc_pos,1)   # to unit cell
                if not any( [is_periodic_image(new_sc_pos,n_p) for n_p in new_scaled_positions] ):
                    new_scaled_positions.append(new_sc_pos)
                    new_numbers.append(num)
            new_atoms = pAtoms(cell=latt,scaled_positions=new_scaled_positions,
                               numbers=new_numbers,
                               pbc=True, symprec=self.symprec)
        return new_atoms
    # end get_primitive_standardized

def get_spacegroup(atoms, symprec=symprec, mode=0):
    return Spacegroup(atoms, symprec=symprec, mode=mode)
