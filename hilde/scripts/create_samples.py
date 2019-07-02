"""
Script to initialize positions and velocities with force constants.
Similar to canonical_sampling from TDEP.
"""
from argparse import ArgumentParser as argpars
import numpy as np

from ase.io import read
from ase import units as u
import ase.md.velocitydistribution as vd

from hilde.structure.io import inform
from hilde.konstanten.einheiten import omega_to_THz
from hilde.helpers import talk


def create_samples(
    geometry,
    temperature,
    n_samples,
    force_constants,
    mc_rattle,
    quantum,
    deterministic,
    random_seed,
    format,
):
    """create samples for Monte Carlo sampling

    Parameters
    ----------
    geometry: str
        The input geometry file
    temperature: float
        The temperature in Kelvin
    n_samples: int
        The number of samples to create (default: 1)
    force_constants: str
        The filename of the file holding force constants for phonon rattle
    mc_rattle: bool
        If True use hiphive mc rattle
    quantum: bool
        If True use Bose-Einstein distribution instead of Maxwell-Boltzmann
    deterministic: bool
        If True create sample deterministically
    random_seed: int
        The seed the random number generator
    format: str
        The ASE file format for geometry files
    """

    atoms = read(geometry, format=format)
    inform(atoms, verbosity=0)

    seed = random_seed
    temp = temperature
    info_str = []

    if mc_rattle:
        try:
            from hiphive.structure_generation import mc_rattle
        except (ModuleNotFoundError, ImportError):
            exit("** hiphive needs to be installed to use mc_rattle.")
    else:
        if not temp:
            exit("** temperature needs to be given")

    if seed:
        talk(f"Random seed is {seed}")
        rng = np.random.RandomState(seed)
        info_str += [f"Random seed: {seed}"]
    else:
        rng = np.random

    if force_constants is not None:
        # if 3Nx3N shaped txt file:
        try:
            fc = np.loadtxt(force_constants)
        except ValueError:
            exit("other force constants not yet implemented")

        # Check dyn. matrix
        check_frequencies(atoms, fc)

        # collect arguments for PhononHarmonics
        phonon_harmonic_args = {
            "force_constants": fc,
            "quantum": quantum,
            "temp": temp * u.kB,
            "plus_minus": deterministic,
            "failfast": True,
            "rng": rng,
        }
        info_str += ["created from force constants", f"T = {temp} K"]
        talk(f"Use force constants from {force_constants} to prepare samples")

    else:
        mb_args = {"temp": temp * u.kB, "rng": rng}
        info_str += ["created from MB distrubtion", f"T = {temperature} K"]
        talk(f"Use Maxwell Boltzamnn to set up samples")

    for ii in range(n_samples):
        talk(f"Sample {ii:3d}:")
        sample = atoms.copy()

        if force_constants is not None:
            vd.PhononHarmonics(sample, **phonon_harmonic_args)

        elif mc_rattle:
            exit("mc_rattle not yet implemented")

        else:
            vd.MaxwellBoltzmannDistribution(sample, **mb_args)

        talk(f".. remove net momentum from sample and force temperature")
        talk(f".. temperature before cleaning: {sample.get_temperature():.3f}K")
        vd.force_temperature(sample, temp)
        vd.Stationary(sample)
        vd.ZeroRotation(sample)

        filename = f"{geometry}.{int(temp)}K"
        if n_samples > 1:
            filename += f".{ii+1:03d}"

        sample.write(filename, info_str=info_str, velocities=True, format=format)

        talk(f".. temperature in sample {ii}:     {sample.get_temperature():.3f}K")
        talk(f".. written to {filename}")


def main():
    """ main function """
    parser = argpars(description="Read geometry create supercell")
    parser.add_argument("geom", type=str, help="geometry input file")
    parser.add_argument("-T", "--temperature", type=int)
    parser.add_argument("-fc", "--force_constants")
    parser.add_argument("--mc_rattle", nargs="?", type=float, const=0.01, default=None)
    parser.add_argument("-n", "--n_samples", type=int, default=1, help="no. of samples")
    parser.add_argument("--quantum", action="store_true")
    parser.add_argument("--deterministic", action="store_true")
    parser.add_argument("--ignore_negative", action="store_false")
    parser.add_argument("--format", default="aims")
    parser.add_argument("--non_enforced_temp", action="store_true")
    parser.add_argument("--non_stationary", action="store_true")
    parser.add_argument("--random_seed", type=int, default=None)
    args = parser.parse_args()

    create_samples(
        args.geom,
        args.temperature,
        args.n_samples,
        args.force_constants,
        args.mc_rattle,
        args.quantum,
        args.deterministic,
        args.random_seed,
        args.format,
    )


if __name__ == "__main__":
    main()


def get_frequencies(atoms, force_constants):
    """ create dynamical matrix, return frequencies for sanity checks """
    masses = atoms.get_masses()
    # Build dynamical matrix
    rminv = (masses ** -0.5).repeat(3)
    dynamical_matrix = force_constants * rminv[:, None] * rminv[None, :]

    # Solve eigenvalue problem to compute phonon spectrum and eigenvectors
    w2_s, _ = np.linalg.eigh(dynamical_matrix)

    return w2_s * (omega_to_THz) ** 2


def check_frequencies(atoms, force_constants):
    """print lowest and highest frequencies obtained from force constants"""
    w2 = get_frequencies(atoms, force_constants)
    print("The first 6 frequencies:")
    for ii, freq in enumerate(w2[:6]):
        print(f" {ii + 1:4d}: {np.sign(freq) * np.sqrt(abs(freq))}")

    print("Highest 6 frequencies")
    for ii, freq in enumerate(w2[-6:]):
        print(f" {len(w2) - ii:4d}: {np.sign(freq) * np.sqrt(abs(freq))}")
