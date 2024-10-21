"""Provide a full highlevel phono3py workflow"""

from pathlib import Path

import numpy as np

from vibes.filenames import filenames
from vibes.helpers import Timer, warn
from vibes.helpers.converters import dict2atoms
from vibes.helpers.paths import cwd
from vibes.io import write
from vibes.phono3py.wrapper import prepare_phono3py
from vibes.phonopy import displacement_id_str
from vibes.structure.convert import to_Atoms
from vibes.trajectory import reader


def postprocess(
    trajectory=filenames.trajectory,
    workdir=".",
    verbose=True,
    **kwargs,
):
    """
    Phonopy postprocess

    Args:
    ----
        trajectory: The trajectory file to process
        workdir: The working directory where trajectory is stored
        verbose: be verbose

    Returns:
    -------
        phono3py.Phono3py: The Phono3py object with the force constants calculated

    """
    timer = Timer("Start phono3py postprocess:")

    trajectory = Path(workdir) / trajectory

    calculated_atoms, metadata = reader(trajectory, get_metadata=True)

    atoms_ids = [atoms.info[displacement_id_str] for atoms in calculated_atoms]
    for nn in range(atoms_ids[-1]):
        if nn not in atoms_ids:
            calculated_atoms.insert(nn, None)
        elif calculated_atoms[nn].info[displacement_id_str] != nn:
            warn(f"Displacement ids are not in order. Inspect {trajectory}!", level=2)

    excluded_displacements = []
    for disp in metadata["Phono3py"]["displacement_dataset"]["first_atoms"]:
        disp["number"] = int(disp["number"])
        for d in disp["second_atoms"]:
            d["number"] = int(d["number"])
            if not d["included"]:
                excluded_displacements.append(int(d["id"]) - 1)

    primitive = dict2atoms(metadata["Phono3py"]["primitive"])
    supercell = dict2atoms(metadata["atoms"])
    supercell_matrix = metadata["Phono3py"]["supercell_matrix"]
    supercell.info = {"supercell_matrix": str(supercell_matrix)}
    symprec = metadata["Phono3py"]["symprec"]

    phonon3 = prepare_phono3py(primitive, supercell_matrix, symprec=symprec)
    phonon3.dataset = metadata["Phono3py"]["displacement_dataset"].copy()

    scs = phonon3.supercells_with_displacements

    n_sc = len(scs)
    n_calc = len(calculated_atoms)
    for disp_id in range(n_calc, n_sc):
        if disp_id not in excluded_displacements:
            msg = f"No. of supercells {n_sc} != no. of calculated atoms: {n_calc}"
            raise RuntimeError(msg)
        n_calc += 1
        calculated_atoms.append(None)

    force_sets = []
    for pa, a in zip(scs, calculated_atoms):
        if pa is None:
            force_sets.append(np.zeros([len(supercell), 3]))
            continue
        force_sets.append(a.get_forces())

    phonon3.forces = np.array(force_sets)

    phonon3.produce_fc2()
    phonon3.produce_fc3()

    if verbose:
        timer("done")

    return phonon3


def extract_results(phonon3, output_dir="output"):
    from phono3py import file_IO as io

    from .wrapper import phono3py_save

    primitive = phonon3.unitcell
    supercell = phonon3.supercell
    p2s_map = phonon3.primitive.p2s_map

    dds = phonon3.dataset

    fc2 = phonon3.fc2
    fc3 = phonon3.fc3

    with cwd(output_dir, mkdir=True):
        p = to_Atoms(primitive)
        write(p, filenames.primitive)
        s = to_Atoms(supercell)
        write(s, filenames.supercell)

        io.write_disp_fc3_yaml(dds, supercell)

        io.write_fc2_to_hdf5(fc2, p2s_map=p2s_map)
        io.write_fc3_to_hdf5(fc3, p2s_map=p2s_map)

        # save yaml
        phono3py_save(phonon3)
