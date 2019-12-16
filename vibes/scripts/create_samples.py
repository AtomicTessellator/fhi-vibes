"""
Script to initialize positions and velocities with force constants.
Similar to canonical_sampling from TDEP.
"""
from argparse import ArgumentParser as argpars
import numpy as np

from ase.io import read
from ase import units as u
import vibes.ase.md.velocitydistribution as vd

from vibes.structure.io import inform
from vibes.helpers import talk
from vibes.harmonic_analysis.dynamical_matrix import get_frequencies


def generate_samples(
    atoms,
    temperature,
    n_samples,
    force_constants,
    rattle,
    quantum,
    deterministic,
    plus_minus,
    gauge_eigenvectors,
    ignore_negative,
    sobol,
    random_seed,
    propagate,
    format,
):
    """create samples for Monte Carlo sampling

    Args:
        geometry: input geometry file
        temperature: temperature in Kelvin
        n_samples: number of samples to create (default: 1)
        force_constants: filename of the file holding force constants for phonon rattle
        rattle: atoms.rattle
        quantum: use Bose-Einstein distribution instead of Maxwell-Boltzmann
        deterministic: create sample deterministically
        plus_minus: use +/-
        gauge_eigenvectors: make largest entry positive
        sobol: Use sobol numbers for the sampling
        random_seed: seed for random number generator
        propagate: propagate atoms according to velocities for this many fs
        format: The ASE file format for geometry files
    """

    inform(atoms, verbosity=0)

    seed = random_seed
    temp = temperature
    info_str = []

    if not rattle:
        if temp is None:
            exit("** temperature needs to be given")

    if not seed:
        seed = np.random.randint(2 ** 31)
        if sobol:
            seed = np.random.randint(2 ** 16)

    if sobol:
        from vibes.helpers.sobol import RandomState

        # create sobol generator with dimension 3N - 3
        # check that `nw` coincides with `nw` in `velocitydistribution.phonon_harmonics`
        nw = 3 * len(atoms) - 3
        rng = RandomState(dimension=nw, seed=seed, failsafe=False)
    else:
        rng = np.random.RandomState(seed)

    if force_constants is not None:
        # if 3Nx3N shaped force_constants:
        if np.any(np.array(force_constants.shape, dtype=int) != 3 * len(atoms)):
            exit("other force constants not yet implemented")

        # Check dyn. matrix
        check_frequencies(atoms, force_constants)

        # collect arguments for PhononHarmonics
        phonon_harmonic_args = {
            "force_constants": force_constants,
            "quantum": quantum,
            "temp": temp * u.kB,
            "failfast": True,
            "rng": rng,
            "deterministic": deterministic,
            "plus_minus": plus_minus,
            "gauge_eigenvectors": gauge_eigenvectors,
            "ignore_negative": ignore_negative,
        }

        info_str += ["created from force constants", f"T = {temp} K"]
        talk(f"Random seed: {seed}")
    else:
        mb_args = {"temp": temp * u.kB, "rng": rng}
        info_str += ["created from MB distrubtion", f"T = {temperature} K"]
        talk(f"Use Maxwell Boltzamnn to set up samples")

    info_str += [
        f"quantum:             {quantum}",
        f"deterministic:       {deterministic}",
        f"plus_minus:          {plus_minus}",
        f"gauge_eigenvectors:  {gauge_eigenvectors}",
        f"Sobol numbers:       {sobol}",
        f"Random seed:         {seed}",
    ]

    sample_list = []

    for ii in range(n_samples):
        talk(f"Sample {ii:3d}:")
        sample_info_str = info_str + [f"Sample number:       {ii + 1}"]
        sample = atoms.copy()

        if force_constants is not None:
            vd.PhononHarmonics(sample, **phonon_harmonic_args)

        elif rattle is not None:
            sample.rattle(rattle)
            sample_info_str = info_str + [f"Rattle with stdev:   {rattle}"]

        else:
            vd.MaxwellBoltzmannDistribution(sample, **mb_args)

        if force_constants is not None:
            d = np.ravel(sample.positions - atoms.positions)
            epot_ha = 0.5 * (force_constants @ d @ d)
            epot_ha_temp = epot_ha / u.kB / len(atoms) / 3 * 2

            ha_epot_str = f"{epot_ha:9.3f}eV ({epot_ha_temp:.2f}K)"

            sample_info_str += [f"Harmonic E_pot:    {ha_epot_str}"]
            talk(f".. harmonic potential energy:   {ha_epot_str})")

        talk(f".. temperature before cleaning: {sample.get_temperature():9.3f}K")
        talk(f".. remove net momentum from sample and force temperature")
        vd.force_temperature(sample, temp)
        vd.Stationary(sample)
        # vd.ZeroRotation(sample)

        if propagate:
            talk(f".. propagate positions for {propagate} fs")
            sample_info_str += [f"Propagated for:      {propagate} fs"]
            sample.positions += sample.get_velocities() * propagate * u.fs

        sample.info["info_str"] = sample_info_str
        sample_list.append(sample)

        talk(f".. temperature in sample {ii}:     {sample.get_temperature():9.3f}K")

    return sample_list


def create_samples(
    geometry,
    temperature,
    n_samples,
    force_constants,
    rattle,
    quantum,
    deterministic,
    plus_minus,
    gauge_eigenvectors,
    ignore_negative,
    sobol,
    random_seed,
    propagate,
    format,
):
    """create samples for Monte Carlo sampling

    Args:
        geometry: input geometry file
        temperature: temperature in Kelvin
        n_samples: number of samples to create (default: 1)
        force_constants: filename of the file holding force constants for phonon rattle
        mc_rattle: hiphive mc rattle
        quantum: use Bose-Einstein distribution instead of Maxwell-Boltzmann
        deterministic: create sample deterministically
        sobol: Use sobol numbers for the sampling
        random_seed: seed for random number generator
        format: The ASE file format for geometry files
        return_samples (bool): If True do not write the samples, but return them
    """

    atoms = read(geometry, format=format)
    inform(atoms, verbosity=0)

    fc = None
    if force_constants is not None:
        # if 3Nx3N shaped txt file:
        try:
            fc = np.loadtxt(force_constants)
        except ValueError:
            exit("other force constants not yet implemented")
        talk(f"\nUse force constants from {force_constants} to prepare samples")

    sample_list = generate_samples(
        atoms,
        temperature,
        n_samples,
        fc,
        rattle,
        quantum,
        deterministic,
        plus_minus,
        gauge_eigenvectors,
        ignore_negative,
        sobol,
        random_seed,
        propagate,
        format,
    )

    for ii, sample in enumerate(sample_list):
        talk(f"Sample {ii:3d}:")
        filename = f"{geometry}.{int(temperature):04d}K"
        if n_samples > 1:
            filename += f".{ii:03d}"

        sample.write(
            filename, info_str=sample.info.pop("info_str"), velocities=True, format=format
        )
        talk(f".. temperature in sample {ii}:     {sample.get_temperature():9.3f}K")
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
    parser.add_argument("--sobol", action="store_true")
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
        args.sobol,
        args.random_seed,
        args.format,
    )


if __name__ == "__main__":
    main()


def check_frequencies(atoms, force_constants):
    """print lowest and highest frequencies obtained from force constants"""
    w2 = get_frequencies(force_constants, masses=atoms.get_masses())

    print("The first 6 frequencies:")
    for ii, freq in enumerate(w2[:6]):
        print(f" {ii + 1:4d}: {freq}")

    print("Highest 6 frequencies")
    for ii, freq in enumerate(w2[:-7:-1]):
        print(f" {len(w2) - ii:4d}: {freq }")
