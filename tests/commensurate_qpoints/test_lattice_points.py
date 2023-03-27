from pathlib import Path

from vibes.helpers.lattice_points import get_commensurate_q_points, get_lattice_points
from vibes.io import read


tolerance = 1e-5
parent = Path(__file__).parent
materials = ["si", "gan", "gao"]


def check_q_points(q_points, superlattice, tol=tolerance):
    return all(
        (abs((ilp @ L + tol) % 1 - tol) < tol for ilp in q_points for L in superlattice)
    )


for material in materials:

    print(f"Test {material}")

    primitive = read(parent / material / "geometry.in")
    supercell = read(parent / material / "geometry.in.supercell")

    print("\nReal space lattice points")
    lattice_points, _ = get_lattice_points(primitive.cell, supercell.cell,)

    print("\nMomentum space lattice points")
    inv_lattice_points = get_commensurate_q_points(primitive.cell, supercell.cell,)

    check = check_q_points(inv_lattice_points, supercell.cell)
    print(f"\n-> All q_points orthogonal to all lattice vectors: {check}")

    assert check, f"*** failed for {material}"

    for (ilp, L) in ((ilp, L) for ilp in inv_lattice_points for L in supercell.cell):
        if not abs((ilp @ L + 0.001) % 1 - 0.001) < 1e-9:
            print(f"ERROR: q: {ilp} L: {L} q.L: {ilp @ L}")

    print("\n")
