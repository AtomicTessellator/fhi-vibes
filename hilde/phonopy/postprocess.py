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
    **kwargs,
):
    """ Phonopy postprocess """

    timer = Timer()
    print("Start phonopy postprocess:")
    trajectory = Path(trajectory)

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
    if pickle_file and write_files:
        fname = trajectory.parent / pickle_file
        psave(phonon, fname)
        print(f".. Pickled phonopy object written to {fname}")

    if write_files:
        # Save the supercell
        fname = "geometry.in.supercell"
        write(supercell, fname)
        print(f".. Supercell written to {fname}")

        force_constants = get_force_constants(phonon)
        fname = "force_constants.dat"
        np.savetxt(fname, force_constants)
        print(f".. Force constants saved to {fname}.")

    timer("done")

    return phonon


def extract_results(
    phonon,
    plot_bandstructure=True,
    plot_dos=False,
    plot_pdos=False,
    tdep=False,
    force_constants_file="FORCE_CONSTANTS",
):
    """ Extract results from phonopy object and present them.
        With `tdep=True`, the necessary input files for TDEP's
          `convert_phonopy_to_forceconstant`
        are written. """
    from hilde.phonopy.wrapper import plot_bandstructure, plot_bandstructure_and_dos
    from hilde.structure.convert import to_Atoms
    from phonopy.file_IO import write_FORCE_CONSTANTS

    plot_bandstructure(phonon, file="bandstructure.pdf")
    if plot_dos:
        plot_bandstructure_and_dos(phonon, file="bands_and_dos.pdf")
    if plot_pdos:
        plot_bandstructure_and_dos(phonon, partial=True, file="bands_and_pdos.pdf")

    # print structures:
    primitive = to_Atoms(phonon.get_primitive())
    supercell = to_Atoms(phonon.get_supercell())

    if tdep:
        write_settings = {"format": "vasp", "direct": True, "vasp5": True}
        fnames = {"primitive": "infile.ucposcar", "supercell": "infile.ssposcar"}
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

    # reproduce reduces force constants
    phonon.produce_force_constants(calculate_full_force_constants=False)

    write_FORCE_CONSTANTS(
        phonon.get_force_constants(),
        filename=args.force_constants_file,
        p2s_map=phonon.get_primitive().get_primitive_to_supercell_map(),
    )

    print(f"Force constants saved to {args.force_constants_file}.")
