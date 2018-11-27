"""
Script to initialize positions and velocities with force constants.
Similar to canonical_sampling from TDEP.
"""
from argparse import ArgumentParser as argpars
import numpy as np
from ase.io import read
from hilde.structure.io import inform
from hilde.konstanten.einheiten import eV_to_THz
from hilde.molecular_dynamics.velocitydistribution import PhononHarmonics


def get_frequencies(atoms, force_constants):
    masses = atoms.get_masses()
    # Build dynamical matrix
    rminv = (masses ** -0.5).repeat(3)
    dynamical_matrix = force_constants * rminv[:, None] * rminv[None, :]

    # Solve eigenvalue problem to compute phonon spectrum and eigenvectors
    w2_s, _ = np.linalg.eigh(dynamical_matrix)

    return w2_s * (eV_to_THz) ** 2


def main():
    parser = argpars(description="Read geometry create supercell")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument("temperature", type=int)
    parser.add_argument("-fc", "--force_constants", default="force_constants.dat")
    parser.add_argument("-n", type=int, default=1, help="no. of samples")
    parser.add_argument("--quantum", action="store_true")
    parser.add_argument("--ignore_negative", action="store_false")
    parser.add_argument("--format", default="aims")
    args = parser.parse_args()

    atoms = read(args.geom, format=args.format)
    inform(atoms)

    force_constants = np.loadtxt(args.force_constants)

    # Check dyn. matrix
    w2 = get_frequencies(atoms, force_constants)
    print("The first 10 frequencies:")
    for ii, freq in enumerate(w2[:10]):
        print(f" {ii + 1:4d}: {np.sign(freq) * np.sqrt(abs(freq))}")

    print("Highest 6 frequencies")
    for ii, freq in enumerate(w2[-6:]):
        print(f" {len(w2) - ii:4d}: {np.sign(freq) * np.sqrt(abs(freq))}")

    for ii in range(args.n):
        sample = atoms.copy()

        PhononHarmonics(
            sample,
            force_constants,
            args.temperature,
            quantum=args.quantum,
            failfast=args.ignore_negative,
        )

        sample.write(
            f"{args.geom}.{ii:03d}",
            velocities=True,
            info_str=["created from force constants", f"T = {args.temperature} K"],
            format=args.format,
        )

        print(f"Temperature in sample: {sample.get_temperature():.3f}K")


if __name__ == "__main__":
    main()