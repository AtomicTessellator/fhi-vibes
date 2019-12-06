""" Example of hiphive calculation for silicon in comparison to phonopy """

from pathlib import Path
import numpy as np

from ase.build import bulk
from ase.calculators.emt import EMT

from vibes.structure import pAtoms
from vibes.helpers.supercell import make_cubic_supercell
from vibes.phonopy import wrapper as ph
from vibes.tasks.calculate import calculate_multiple
from vibes.helpers.geometry import inscribed_sphere_in_box

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
            seed=42):
    """ Get a ForceConstantPotential for the material """
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
    fcp.write('al.fcp')

    return fcp


# setup
atoms_ideal = pAtoms(ase_atoms=bulk('Al'))
_, smatrix = make_cubic_supercell(atoms_ideal, 28)

# explicit is better than implicit
phonon, supercell, supercells_with_disps = ph.preprocess(atoms_ideal, smatrix.T)
workdir = setup_workdir(atoms_ideal, smatrix=smatrix)

calc = EMT()

# force constant potential
cutoff_max = inscribed_sphere_in_box(supercell.cell) - .01
fcp_params = {
    'number_of_structures': 5,
    'rattle_std': 0.01,
    'minimum_distance': 2.4,
    'cutoffs': [cutoff_max]
    }

fcp = get_fcp(supercell, calc, **fcp_params)
fcs = fcp.get_force_constants(supercell)

# parameters
mesh = 3 * [5]  # q-point mesh for MSD calculation
temperatures = [300, 900]  # temperatures for evaluating MSD

scs_calculated = calculate_multiple(supercells_with_disps,
                                    calc,
                                    workdir=workdir,
                                    trajectory='lammps.traj')

phonon.produce_force_constants([sc.get_forces() for sc in scs_calculated])
fc_phonopy = phonon.get_force_constants()


phonon = ph.prepare_phonopy(atoms_ideal, smatrix.T, fc2=fc_phonopy)
phonon.set_mesh(mesh)
_, _, frequencies, _ = phonon.get_mesh()
freq_max_phonopy = frequencies.max()

phonon = ph.prepare_phonopy(atoms_ideal, smatrix.T,
                            fc2=fcs.get_fc_array(order=2))
phonon.set_mesh(mesh)
_, _, frequencies, _ = phonon.get_mesh()
freq_max_hiphive = frequencies.max()

diff = lambda x, y, z: np.linalg.norm(x - y) / np.linalg.norm(z)
freq_diff = diff(freq_max_phonopy, freq_max_hiphive, 1)
fc2_diff = diff(fc_phonopy, phonon.get_force_constants(), fc_phonopy)

print("The largest frequency difference between phonopy and hiphive " +
      f"on a {mesh} grid is: {freq_diff:.2e} THz")
print("The FC2 error between the phonopy and hiphive force constant matrix is:"
      + f"{fc2_diff:.2e}")

assert abs(freq_diff) < .1
