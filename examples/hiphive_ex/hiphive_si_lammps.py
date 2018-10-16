""" Example of hiphive calculation for silicon in comparison to phonopy """

from pathlib import Path
from random import randint
import numpy as np

from hilde.helpers.supercell import make_cubic_supercell
from hilde.materials_fp import get_phonon_bs_fingerprint_phononpy
from hilde.parsers import read_structure
from hilde.phonopy import phono as ph
from hilde.templates.lammps import setup_lammps_si
from hilde.tasks.calculate import calculate_multiple
from hilde.helpers.brillouinzone import get_bands_and_labels
from hilde.helpers.geometry import inscribed_sphere_in_box

from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive.fitting import Optimizer
from hiphive.structure_generation import generate_mc_rattled_structures


def setup_workdir(atoms, smatrix):
    """ Set up a working directory """
    vol = atoms.get_volume()
    workdir = Path('./sc_{}{}{}_{}{}{}_{}{}{}_{:.3f}_lammps'.format(
        *smatrix.flatten(), vol)).absolute()
    workdir.mkdir(parents=True, exist_ok=True)
    return workdir

def get_bs(phonon, filename):
    bands, labels = get_bands_and_labels(phonon.primitive)
    phonon.set_band_structure(bands)
    plt = phonon.plot_band_structure(labels=labels)
    plt.ylim([0, 18])
    plt.ylabel('Frequency [THz]')
    plt.savefig(filename)
    return get_phonon_bs_fingerprint_phononpy(phonon, binning=False)

def get_thermal(phonon):
    phonon.set_mesh(mesh, is_eigenvectors=True, is_mesh_symmetry=False)
    phonon.set_thermal_displacements(temperatures=temperatures)
    _, msds = phonon.get_thermal_displacements()
    msds = np.sum(msds, axis=1)  # sum up the MSD over x,y,z
    phonon.set_thermal_properties(temperatures=temperatures)
    _, free_ener, _, _ = phonon.get_thermal_properties()
    for temperature, msd, A in zip(temperatures, msds, free_ener):
        print('T = {:4d} K    MSD = {:.5f} A**2, A = {:.5f} KJ/mol'.format(temperature, msd, A))


def get_structures_with_displacements(supercell, structures):
    """ Compute displacements wrt. supercell and set them as new array"""
    s_with_disps = []
    for structure in structures:
        struct = supercell.copy()
        forces = structure.get_forces()
        displacements = structure.positions - supercell.positions
        struct.new_array('displacements', displacements)
        struct.new_array('forces', forces)
        s_with_disps.append(struct)
    return s_with_disps


def get_fcp(supercell, calculator,
            cutoffs,
            number_of_structures,
            rattle_std,
            minimum_distance,
            fit_method='ardr',
            seed=randint(0, 2**32 - 1),
            filename="si.fcp",
            force=False):
    """ Get a ForceConstantPotential for the material """
    try:
        if force:
            raise
        fcp = ForceConstantPotential.read(filename)
    except:
        structures = generate_mc_rattled_structures(supercell,
                                                    number_of_structures,
                                                    rattle_std,
                                                    minimum_distance,
                                                    seed=seed)
        # set up cluster space
        cs = ClusterSpace(supercell, cutoffs)

        # compute forces
        structures = calculate_multiple(structures,
                                        calculator,
                                        workdir='hiphive',
                                        trajectory='hiphive.traj')
        # generate displacements
        structures = get_structures_with_displacements(supercell,
                                                       structures)

        # ... and structure container
        sc = StructureContainer(cs)
        for structure in structures:
            sc.add_structure(structure)

        # train model
        # opt = Optimizer(sc.get_fit_data())
        opt = Optimizer(sc.get_fit_data(),
                        fit_method=fit_method,
                        seed=seed
                       )
        opt.train()

        # construct force constant potential
        fcp = ForceConstantPotential(cs, opt.parameters)
        fcp.write('si.fcp')

    return fcp


# setup
atoms_ideal = read_structure('si.in')
_, smatrix = make_cubic_supercell(atoms_ideal, 64)

# explicit is better than implicit
phonon, supercell, supercells_with_disps = ph.preprocess(atoms_ideal, smatrix.T)
workdir = setup_workdir(atoms_ideal, smatrix=smatrix)

calc = setup_lammps_si(workdir)

# force constant potential
cutoff_max = inscribed_sphere_in_box(supercell.cell) - .01
fcp_params = {
    'number_of_structures': 5,
    'rattle_std': 0.01,
    'minimum_distance': 2.4,
    'cutoffs': [cutoff_max, 3.0],
    'force': False
}

fcp = get_fcp(supercell, calc, **fcp_params)
fcs = fcp.get_force_constants(supercell)

# parameters
mesh = [32, 32, 32]  # q-point mesh for MSD calculation
temperatures = [300, 600, 900, 1200]  # temperatures for evaluating MSD

scs_calculated = calculate_multiple(supercells_with_disps,
                                    calc,
                                    workdir=workdir,
                                    trajectory='lammps.traj')
phonon.produce_force_constants([sc.get_forces() for sc in scs_calculated])
fc_phonopy = phonon.get_force_constants()

print("phonopy")
get_thermal(phonon)
bs_phonopy = get_bs(phonon, "bs_phonopy.pdf")

print("hiphive")
phonon.set_force_constants(fcs.get_fc_array(order=2))
get_thermal(phonon)
bs_hiphive = get_bs(phonon, "bs_hiphive.pdf")

diff = lambda x, y, z: np.linalg.norm(x - y) / np.linalg.norm(z)
freq_diff = diff(bs_hiphive[0][0], bs_phonopy[0][0], bs_phonopy[0][0])
fc2_diff = diff(fc_phonopy, phonon.get_force_constants(), fc_phonopy)

print("The frequency error between the phonopy and hiphive band structure at " +
      f"the Gamma point is: {freq_diff:.2e}")
print("The FC2 error between the phonopy and hiphive force constant matrix is:"
      + f"{fc2_diff:.2e}")
