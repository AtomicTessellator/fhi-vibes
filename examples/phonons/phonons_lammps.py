""" Example on how to run a phonopy calculation with lammps as calculator """

import os
from time import time
from pathlib import Path
import matplotlib.pyplot as plt

from hilde.parsers import read_aims
from hilde.helpers.supercell import find_cubic_cell, make_supercell
from hilde.phonopy import phono as ph
from hilde.tasks.calculate import compute_forces
from hilde.templates.lammps import setup_lammps_si

def get_smatrix(atoms, n_target=64):
    """ Return the supercell matrix for atoms with target size """
    target_size = n_target / len(atoms)
    return find_cubic_cell(cell=atoms.cell, target_size=target_size)

def setup_workdir(atoms, smatrix):
    """ Set up a working directory """
    vol = atoms.get_volume()
    workdir = Path('./{}_{}{}{}_{}{}{}_{}{}{}_{:.3f}_lammps'.format(
        atoms.sysname, *smatrix.flatten(), vol)).absolute()
    workdir.mkdir(exist_ok=True)
    return workdir


def plot_dos_and_bandstructure(phonon, workdir):
    """ Plot DOS and bandstructure in the working directory """
    dos = ph.get_dos(phonon,
                     q_mesh=[25, 25, 25],
                     freq_max='auto')
    plt.plot(dos[0], dos[1])
    plt.savefig(workdir / f'dos.pdf')

    *_, labels = ph.get_bandstructure(phonon)

    phonon.plot_band_structure(labels=labels)
    plt.savefig(workdir / 'bandstructure.pdf')


def main():
    """ Main function to run the example calculation """
    atoms = read_aims('si.in')

    smatrix = get_smatrix(atoms)
    print(make_supercell(atoms, smatrix).cell)

    workdir = setup_workdir(atoms, smatrix=smatrix)

    lammps = setup_lammps_si(workdir)

    # set the phonon object
    phonon, _, scs = ph.preprocess(atoms, smatrix.T)
    print(f'{len(scs)} supercells created.')

    tmp_dir = workdir

    stime = time()
    force_sets = compute_forces(cells=scs,
                                calculator=lammps,
                                workdir=tmp_dir)
    timing = time() - stime
    print(f'.. done in {timing:.2f}s')

    phonon.produce_force_constants(force_sets)

    plot_dos_and_bandstructure(phonon, workdir)

main()
