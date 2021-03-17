from pathlib import Path

from ase.io import read
from phonopy.file_IO import write_FORCE_CONSTANTS

from vibes.hiphive import enforce_rotational_sum_rules
from vibes.io import parse_force_constants


def main(
    fc_file: Path = "FORCE_CONSTANTS",
    primitive_file: Path = "geometry.in.primitive",
    supercell_file: Path = "geometry.in.supercell",
    format: str = "aims",
):
    """enforce sum rules on force constants

    Args:
        fc_file: trajectory dataset file


    """
    primitive = read(primitive_file, format=format)
    supercell = read(supercell_file, format=format)

    fc = parse_force_constants(fc_file)

    kw = {"primitive": primitive, "supercell": supercell}
    new_fc = enforce_rotational_sum_rules(fc, **kw)

    outfile = fc_file + "_sum_rules"
    print(f".. write to {outfile}")
    write_FORCE_CONSTANTS(new_fc, filename=outfile)


if __name__ == "__main__":
    import typer

    typer.run(main)
