"""Provide Thermodynamics"""
import collections

import numpy as np
import pandas as pd
from ase import units as u

from vibes.helpers import talk


prefix = "td"


def _talk(msg, verbose: bool = True):
    return talk(msg, prefix=prefix, verbose=verbose)


print = _talk


def get_harmonic_properties(energies: list, temperature: list, classical: bool = False):
    """Compute harmonic energy, entropy, and free energy (U, S, F) at 1 temperature

    Args:
        energies: List of energies in eV (w/o translational modes!)
        temperature: temperature in K
        classical: enforce classical statistics

    Returns:
        namedtuple with colums "U, S, F" (energy, entropy, free energy)

    """
    T = temperature
    kT = T * u.kB

    N = len(energies)
    hw = energies
    weights = np.exp(-hw / kT)

    if classical:
        U = N * kT
        S = N * u.kB
        S -= u.kB * np.log(hw / kT).sum()
    else:
        # energy
        U = 0.5 * hw.sum()
        U += (hw / (weights ** -1 - 1)).sum()

        # entropy
        S = -u.kB * np.log(1 - weights).sum()
        S += 1 / T * (hw / (weights ** -1 - 1)).sum()

    # free energy
    F = U - T * S

    return collections.namedtuple("harmonic_properties", "U S F")(U, S, F)


def get_harmonic_properties_df(
    energies: list, temperatures: list, classical: bool = False, cutoff: float = 1e-6
) -> pd.DataFrame:
    """Compute harmonic energy, entropy, and free energy (U, S, F) at >= 1 temperatures

    Args:
        energies: List of energies in eV
        temperatures: list of temperatures in K
        classical: enforce classical statistics
        cutoff: ignore modes below cutoff (default: 0.001 meV)

    Returns:
        Dataframe with colums "T, U, S, F" (temperature, energy, entropy, free energy)

    """
    mask = energies > cutoff

    hw = energies[mask]

    print(f".. discard {len(energies)-len(hw)} modes (cutoff: {cutoff} eV)")

    rows = []
    for T in temperatures:
        row = get_harmonic_properties(hw, temperature=T, classical=classical)
        rows.append({"T": T, **row._asdict()})

    df = pd.DataFrame(rows).set_index("T")

    return df
