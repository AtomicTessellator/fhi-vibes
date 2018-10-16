import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path

from ase.atoms import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import get_cellinfo, special_paths, bandpath
from ase.io import write, read

from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive.fitting import Optimizer
from hiphive.structure_generation import generate_mc_rattled_structures
from hiphive.fitting import CrossValidationEstimator as CVE

from hilde.helpers.supercell import find_cubic_cell, make_supercell
from hilde.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy, get_phonon_dos_fingerprint_phononpy, scalar_product
from hilde.parsers import read_structure
from hilde.phonopy import phono as ph
from hilde.structure.structure import pAtoms
from hilde.templates.lammps import setup_lammps_si
from hilde.tasks.calculate import compute_forces

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from random import randint

def get_smatrix(atoms, n_target=64):
    """ Return the supercell matrix for atoms with target size """
    target_size = n_target / len(atoms)
    return find_cubic_cell(cell=atoms.cell, target_size=target_size)

def setup_workdir(atoms, smatrix):
    """ Set up a working directory """
    vol = atoms.get_volume()
    workdir = Path('./sc_{}{}{}_{}{}{}_{}{}{}_{:.3f}_lammps'.format(
        *smatrix.flatten(), vol)).absolute()
    workdir.mkdir(parents=True, exist_ok=True)
    return workdir

def get_bs(phonon, filename):
    unitcell = pAtoms(phonopy_atoms=phonon.unitcell)
    cell_info = get_cellinfo(unitcell.cell)
    qpaths = special_paths[cell_info.lattice].split(",")
    bands = []
    for path in qpaths:
        for ii, _ in enumerate(path[:-1]):
            bands.append(bandpath([cell_info.special_points[path[ii]], cell_info.special_points[path[ii+1]]], unitcell.cell)[0])
    phonon.set_band_structure(bands)
    plt = phonon.plot_band_structure(labels=qpaths)
    plt.ylabel('Frequency [THz]')
    plt.savefig(filename)
    return get_phonon_bs_fingerprint_phononpy(phonon, cell_info.special_points, False)

def get_thermal(phonon):
    phonon.set_mesh(mesh, is_eigenvectors=True, is_mesh_symmetry=False)
    phonon.set_thermal_displacements(temperatures=temperatures)
    _, msds = phonon.get_thermal_displacements()
    msds = np.sum(msds, axis=1)  # sum up the MSD over x,y,z
    phonon.set_thermal_properties(temperatures=temperatures)
    _, free_ener, _,_ = phonon.get_thermal_properties()
    for temperature, msd, A in zip(temperatures, msds, free_ener):
        print('T = {:4d} K    MSD = {:.5f} A**2, A = {:.5f} KJ/mol'.format(temperature, msd, A))

# parameters
structures_fname = 'structures/rattled_structures.extxyz'
number_of_structures = 5
rattle_std = 0.03
minimum_distance=2.4
# setup
atoms_ideal = read_structure('si.in')

smatrix = get_smatrix(atoms_ideal, 64)
ph_calc = ph.preprocess(atoms_ideal, smatrix.T) + (setup_workdir(atoms_ideal, smatrix=smatrix),)
calc = setup_lammps_si(ph_calc[3])
try:
    fcp = ForceConstantPotential.read("si.fcp")
except:
    hiphive_supercell = Atoms(ph_calc[1])
    structures = generate_mc_rattled_structures(hiphive_supercell, number_of_structures, rattle_std, minimum_distance, seed=randint(0,2**32-1))
    for structure in structures:
        structure.set_calculator(calc)
        forces = structure.get_forces()
        displacements = structure.positions - hiphive_supercell.get_positions()
        structure.new_array('displacements', displacements)
        structure.new_array('forces', forces)
        structure.positions = hiphive_supercell.get_positions()
        structure.calc = None

    # set up cluster space
    cutoffs = [5.0, 5.0, 5.0, 5.0]
    cs = ClusterSpace(structures[0], cutoffs)

    # ... and structure container
    sc = StructureContainer(cs)
    for structure in structures:
        sc.add_structure(structure)

    # train model
    # opt = Optimizer(sc.get_fit_data())
    opt = Optimizer(sc.get_fit_data(),
                    fit_method='ardr',
                    seed=randint(0,2**32-1)
                   )
    opt.train()

    # construct force constant potential
    fcp = ForceConstantPotential(cs, opt.parameters)
    fcp.write('si.fcp')
    # print(fcp)

# parameters
mesh = [32, 32, 32]  # q-point mesh for MSD calculation
temperatures = [300, 600, 900, 1200]  # temperatures for evaluating MSD

ph, sc, scs, workdir = ph_calc
fcs = fcp.get_force_constants(sc)
print("phonopy")
calc = setup_lammps_si(workdir)
ph.produce_force_constants(compute_forces(cells=scs, calculator=calc, workdir=workdir))
fc_phonopy = ph.get_force_constants()
get_thermal(ph)
bs_phonopy = get_bs(ph, "bs_phonopy.pdf")
print("hiphive")
ph.set_force_constants(fcs.get_fc_array(order=2))
get_thermal(ph)
bs_hiphive = get_bs(ph, "bs_hiphive.pdf")
print(f"The frequency error between the phonopy and hiphive band structure at the Gamma point is: {np.linalg.norm(bs_hiphive[0][0]-bs_phonopy[0][0])/np.linalg.norm(bs_phonopy[0][0])}")
print(f"The FC2 error between the phonopy and hiphive force constant matrix  is: {np.linalg.norm(fc_phonopy-ph.get_force_constants())/np.linalg.norm(fc_phonopy)}")

