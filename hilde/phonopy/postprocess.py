""" Provide a full highlevel phonopy workflow """
from pathlib import Path

import numpy as np

from hilde.helpers.converters import dict2results
from hilde.helpers import Timer
from hilde.phonopy.wrapper import prepare_phonopy, get_force_constants
from hilde.trajectory import reader
from hilde.helpers.pickle import psave
from hilde.io import write


def postprocess(
    trajectory="phonopy/trajectory.yaml",
    pickle_file="phonon.pick",
    write_files=True,
    born_charges_file=None,
    silent=False,
    **kwargs,
):
    """ Phonopy postprocess """

    timer = Timer()
    trajectory = Path(trajectory)
    if not silent:
        print("Start phonopy postprocess:")

    calculated_atoms, metadata = reader(trajectory, True)
    for disp in metadata["Phonopy"]["displacement_dataset"]["first_atoms"]:
        disp["number"] = int(disp["number"])
    primitive = dict2results(metadata["Phonopy"]["primitive"])
    supercell = dict2results(metadata["atoms"])
    supercell_matrix = metadata["Phonopy"]["supercell_matrix"]
    supercell.info = {"supercell_matrix": str(supercell_matrix)}
    symprec = metadata["Phonopy"]["symprec"]

    phonon = prepare_phonopy(primitive, supercell_matrix, symprec=symprec)
    phonon._displacement_dataset = metadata["Phonopy"]["displacement_dataset"].copy()

    force_sets = [atoms.get_forces() for atoms in calculated_atoms]

    phonon.produce_force_constants(force_sets)

    # born charges?
    if born_charges_file:
        from phonopy.file_IO import get_born_parameters

        prim = phonon.get_primitive()
        psym = phonon.get_primitive_symmetry()
        if not silent:
            print(f".. read born effective charges from {born_charges_file}")
        nac_params = get_born_parameters(open(born_charges_file), prim, psym)
        phonon.set_nac_params(nac_params)

    # save pickled phonopy object
    if pickle_file and write_files:
        fname = trajectory.parent / pickle_file
        psave(phonon, fname)
        if not silent:
            print(f".. Pickled phonopy object written to {fname}")

    if write_files:
        # Save the supercell
        fname = "geometry.in.supercell"
        write(supercell, fname)
        if not silent:
            print(f".. Supercell written to {fname}")

        force_constants = get_force_constants(phonon)
        fname = "force_constants.dat"
        np.savetxt(fname, force_constants)
        if not silent:
            print(f".. Force constants saved to {fname}.")
    if not silent:
        timer("done")

    return phonon


def extract_results(
    phonon,
    plot_bandstructure=True,
    plot_dos=False,
    plot_pdos=False,
    tdep=False,
    tdep_reduce_fc=True,
    force_constants_file="FORCE_CONSTANTS",
):
    """ Extract results from phonopy object and present them.
        With `tdep=True`, the necessary input files for TDEP's
          `convert_phonopy_to_forceconstant`
        are written. """
    from hilde.phonopy.wrapper import plot_bandstructure, plot_bandstructure_and_dos
    from hilde.structure.convert import to_Atoms_db
    from phonopy.file_IO import write_FORCE_CONSTANTS

    plot_bandstructure(phonon, file="bandstructure.pdf")
    if plot_dos:
        plot_bandstructure_and_dos(phonon, file="bands_and_dos.pdf")
    if plot_pdos:
        plot_bandstructure_and_dos(phonon, partial=True, file="bands_and_pdos.pdf")

    primitive = to_Atoms_db(phonon.get_primitive())
    supercell = to_Atoms_db(phonon.get_supercell())

    if tdep:
        write_settings = {"format": "vasp", "direct": True, "vasp5": True}
        fnames = {"primitive": "infile.ucposcar", "supercell": "infile.ssposcar"}

        # reproduce reduces force constants
        if tdep_reduce_fc:
            phonon.produce_force_constants(calculate_full_force_constants=False)

        write_FORCE_CONSTANTS(
            phonon.get_force_constants(),
            filename=force_constants_file,
            p2s_map=phonon.get_primitive().get_primitive_to_supercell_map(),
        )

        print(f"Reduced force constants saved to {force_constants_file}.")

    else:
        write_settings = {"format": "aims", "scaled": True}
        fnames = {
            "primitive": "geometry.in.primitive",
            "supercell": "geometry.in.supercell",
        }

    fname = fnames["primitive"]
    primitive.write(fname, **write_settings)
    print(f"Primitive cell written to {fname}")

    fname = fnames["supercell"]
    supercell.write(fname, **write_settings)
    print(f"Supercell cell written to {fname}")

    # save as force_constants.dat
    if tdep_reduce_fc:
        phonon.produce_force_constants()
    n_atoms = phonon.get_supercell().get_number_of_atoms()

    force_constants = (
        phonon.get_force_constants().swapaxes(1, 2).reshape(2 * (3 * n_atoms,))
    )

    fname = "force_constants.dat"
    np.savetxt(fname, force_constants)
    print(f"Full force constants as numpy matrix written to {fname}.")
