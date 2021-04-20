import json
from pathlib import Path

import numpy as np
import xarray
from ase import Atoms
from hiphive import (
    ClusterSpace,
    ForceConstantPotential,
    StructureContainer,
    enforce_rotational_sum_rules,
)
from hiphive.fitting import Optimizer
from phonopy import Phonopy
from phonopy.file_IO import write_FORCE_CONSTANTS
from vibes.helpers.geometry import inscribed_sphere_in_box
from vibes.phonopy.wrapper import plot_bandstructure_and_dos


def get_fit_structures(reference_atoms, dataset, stride=100):
    """get the fit structures from the  dataset using samples separated by STRIDE"""

    fit_structures = []
    atoms_ideal = reference_atoms.copy()

    forces = dataset.forces.data[::stride]
    displacements = dataset.displacements.data[::stride]

    enlarge = False
    if len(atoms_ideal) == 8 * forces.shape[1]:
        enlarge = True

    for f, d in zip(forces, displacements):
        structure = atoms_ideal.copy()

        if enlarge:
            f = np.tile(f, (2 ** 3, 1))
            d = np.tile(d, (2 ** 3, 1))

        structure.arrays["forces"] = f
        structure.arrays["displacements"] = d

        structure.positions = atoms_ideal.positions
        fit_structures.append(structure)

    return fit_structures


def main(
    file: Path,
    sum_rules: bool = False,
    outfile: Path = "FORCE_CONSTANTS_hiphive",
    plot: bool = False,
    stride: int = 100,
    verbose: bool = False,
):
    """fit hiphive force_constants

    Args:
        file: trajectory dataset file
        sum_rules: enfore rotational sum rules
        outfile: where the fitted FC are stored
        stride: use every STRIDE steps for fit
        verbose: be verbose

    """

    ds = xarray.open_dataset(file)

    primitive = Atoms(**json.loads(ds.attrs["atoms_primitive"]))
    supercell = Atoms(**json.loads(ds.attrs["atoms_supercell"]))
    supercell.set_tags(np.arange(len(supercell)))

    supercell_matrix = (
        (supercell.cell @ primitive.cell.reciprocal().T).round().astype(int)
    )

    # build phonopy object and sanity check
    phonon = Phonopy(primitive, supercell_matrix=supercell_matrix.T)
    n0, n1 = len(supercell), len(phonon.supercell)
    c0, c1 = supercell.get_cell(), phonon.supercell.get_cell()
    assert n0 == n1, (n0, n1)
    assert np.allclose(c0, c1), (c0, c1)

    supercell_fit = supercell

    # if enlarge_cutoff:
    #     supercell_fit = supercell.repeat(2)
    # else:
    #     supercell_fit = supercell

    cutoffs = [0.99 * inscribed_sphere_in_box(supercell_fit.cell)]  # max. cutoff

    cs = ClusterSpace(primitive, cutoffs=cutoffs)

    if verbose:
        print(cs)

    fit_structures = get_fit_structures(supercell_fit, ds, stride=stride)

    sc = StructureContainer(cs)

    sc.delete_all_structures()
    for ii, s in enumerate(fit_structures):
        print(f"\nAdd structure {ii:4d} / {len(fit_structures)}")
        sc.add_structure(s)

    if verbose:
        print(sc)

    opt = Optimizer(sc.get_fit_data())
    opt.train()

    print(opt)

    parameters = opt.parameters
    if sum_rules:
        print(".. enforce rotational sum rules")
        sum_rules = ["Huang", "Born-Huang"]
        parameters = enforce_rotational_sum_rules(cs, parameters, sum_rules)
        outfile += "_sum_rules"

    fcp = ForceConstantPotential(cs, parameters)

    fcs = fcp.get_force_constants(supercell_fit)
    fc = fcs.get_fc_array(2)

    # # reduce from enlarged supercell -> this doesn't really work
    # new_fc = np.zeros((len(supercell), len(supercell), 3, 3))
    # for i, j in itertools.product(range(len(supercell_fit)), repeat=2):
    #     ti = supercell_fit[i].tag
    #     tj = supercell_fit[j].tag
    #     new_fc[ti, tj] += fc[i, j]

    # if enlarge_cutoff:
    #     new_fc /= 2 ** 3

    # reduced fc
    print(f".. write fc to {outfile}")
    p2s_map = phonon.primitive.get_primitive_to_supercell_map()
    write_FORCE_CONSTANTS(fc[p2s_map, :, :, :], filename=outfile)

    if plot:
        phonon.set_force_constants(fc)
        file = "bandstructure_" + outfile + ".pdf"
        plot_bandstructure_and_dos(phonon, q_mesh=3 * (31,), file=file)


if __name__ == "__main__":
    import typer

    typer.run(main)
