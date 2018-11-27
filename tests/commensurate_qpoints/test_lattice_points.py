from pathlib import Path
from hilde.io import read

from hilde.helpers.supercell.lattice_points import (
    get_lattice_points,
    get_commensurate_q_points,
)

materials = ['gan_simple'] #["gan", "gao", "si"]

for material in materials:

    print(f"test {material}")

    primitive = read(Path(material) / "geometry.in")
    supercell = read(Path(material) / "geometry.in.supercell")

    lattice_points = get_lattice_points(primitive, supercell, verbose=True)
    inv_lattice_points = get_commensurate_q_points(primitive, supercell, verbose=True)

    for (ilp, L) in ((ilp, L) for ilp in inv_lattice_points for L in supercell.cell):
        print(abs((ilp @ L + .01) % 1 - .01) < 1e-5)

