# encoding: utf-8
# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up velocity distributions such as Maxwell–Boltzmann.

Currently, only a few functions are defined, such as
MaxwellBoltzmannDistribution, which sets the momenta of a list of
atoms according to a Maxwell-Boltzmann distribution at a given
temperature.

COPYRIGHT notice: This file was originally distributed with ase3.18 under the LGPL
license
"""

import numpy as np
from ase import units
from ase.parallel import world

from vibes.helpers import warn


# define a ``zero'' temperature to avoid divisions by zero
eps_temp = 1e-12


class UnitError(Exception):
    """Exception raised when wrong units are specified"""


def force_temperature(atoms, temperature, unit="K"):
    """ force (nucl.) temperature to have a precise value

    Parameters:
    atoms: ase.Atoms
        the structure
    temperature: float
        nuclear temperature to set
    unit: str
        'K' or 'eV' as unit for the temperature
    """

    if unit == "K":
        E_temp = temperature * units.kB
    elif unit == "eV":
        E_temp = temperature
    else:
        raise UnitError("'{}' is not supported, use 'K' or 'eV'.".format(unit))

    if temperature > eps_temp:
        E_kin0 = atoms.get_kinetic_energy() / len(atoms) / 1.5
        gamma = E_temp / E_kin0
    else:
        gamma = 0.0
    atoms.set_momenta(atoms.get_momenta() * np.sqrt(gamma))


def _maxwellboltzmanndistribution(masses, temp, communicator=world, rng=np.random):
    # For parallel GPAW simulations, the random velocities should be
    # distributed.  Uses gpaw world communicator as default, but allow
    # option of specifying other communicator (for ensemble runs)
    xi = rng.standard_normal((len(masses), 3))
    communicator.broadcast(xi, 0)
    momenta = xi * np.sqrt(masses * temp)[:, np.newaxis]
    return momenta


def MaxwellBoltzmannDistribution(
    atoms, temp, communicator=world, force_temp=False, rng=np.random
):
    """Sets the momenta to a Maxwell-Boltzmann distribution.

    temp should be fed in energy units; i.e., for 300 K use
    temp=300.*units.kB. If force_temp is set to True, it scales the
    random momenta such that the temperature request is precise."""
    masses = atoms.get_masses()
    momenta = _maxwellboltzmanndistribution(masses, temp, communicator, rng)
    atoms.set_momenta(momenta)
    if force_temp:
        force_temperature(atoms, temperature=temp, unit="eV")


def Stationary(atoms, preserve_temperature=True):
    "Sets the center-of-mass momentum to zero."

    # Save initial temperature
    temp0 = atoms.get_temperature()

    p = atoms.get_momenta()
    p0 = np.sum(p, 0)
    # We should add a constant velocity, not momentum, to the atoms
    m = atoms.get_masses()
    mtot = np.sum(m)
    v0 = p0 / mtot
    p -= v0 * m[:, np.newaxis]
    atoms.set_momenta(p)

    if preserve_temperature:
        force_temperature(atoms, temp0)


def ZeroRotation(atoms, preserve_temperature=True):
    "Sets the total angular momentum to zero by counteracting rigid rotations."

    # Save initial temperature
    temp0 = atoms.get_temperature()

    # Find the principal moments of inertia and principal axes basis vectors
    Ip, basis = atoms.get_moments_of_inertia(vectors=True)
    # Calculate the total angular momentum and transform to principal basis
    Lp = np.dot(basis, atoms.get_angular_momentum())
    # Calculate the rotation velocity vector in the principal basis, avoiding
    # zero division, and transform it back to the cartesian coordinate system
    omega = np.dot(np.linalg.inv(basis), np.select([Ip > 0], [Lp / Ip]))
    # We subtract a rigid rotation corresponding to this rotation vector
    com = atoms.get_center_of_mass()
    positions = atoms.get_positions()
    positions -= com  # translate center of mass to origin
    velocities = atoms.get_velocities()
    atoms.set_velocities(velocities - np.cross(omega, positions))

    if preserve_temperature:
        force_temperature(atoms, temp0)


def n_BE(temp, omega):
    """Bose-Einstein distribution function.

    Args:
        temp: temperature converted to eV (*units.kB)
        omega: sequence of frequencies converted to eV

    Returns:
        Value of Bose-Einstein distribution function for each energy

    """

    omega = np.asarray(omega)

    # 0K limit
    if temp < eps_temp:
        n = np.zeros_like(omega)
    else:
        n = 1 / (np.exp(omega / (temp)) - 1)
    return n


def phonon_harmonics(
    force_constants,
    masses,
    temp,
    rng=np.random.rand,
    quantum=False,
    deterministic=False,
    plus_minus=False,
    gauge_eigenvectors=False,
    return_eigensolution=False,
    failfast=False,
    ignore_negative=True,
):
    r"""Return displacements and velocities that produce a given temperature.

    Parameters:

    force_constants: array of size 3N x 3N
        force constants (Hessian) of the system in eV/Å²
    masses: array of length N
        masses of the structure in amu
    temp: float
        Temperature converted to eV (T * units.kB)
    rng: function
        Random number generator function, e.g., np.random.rand
    quantum: bool
        True for Bose-Einstein distribution, False for Maxwell-Boltzmann
        (classical limit)
    deterministic: bool
        Displace atoms deterministically by fixing the phases
    plus_minus: bool
        Displace atoms with +/- the amplitude accoding to PRB 94, 075125
    gauge_eigenvectors: bool
        gauge eigenvectors such that largest entry is positive
    return_eigensolution: bool
        return eigenvalues and eigenvectors of the dynamical matrix
    failfast: bool
        If True, raise Error when phonon spectrum is not positive
    ignore_negative: bool
        If True freeze out the imaginary modes

    Returns:

        displacements, velocities generated from the eigenmodes,
        (optional: eigenvalues, eigenvectors of dynamical matrix)

    Purpose:

        Excite phonon modes to specified temperature.

    This excites all phonon modes randomly so that each contributes,
    on average, equally to the given temperature.  Both potential
    energy and kinetic energy will be consistent with the phononic
    vibrations characteristic of the specified temperature.

    In other words the system will be equilibrated for an MD run at
    that temperature.

    force_constants should be the matrix as force constants, e.g.,
    as computed by the ase.phonons module.

    Let X_ai be the phonon modes indexed by atom and mode, w_i the
    phonon frequencies, and let 0 < Q_i <= 1 and 0 <= R_i < 1 be
    uniformly random numbers.  Then

    .. code-block:: none


                    1/2
       _     / k T \     ---  1  _             1/2
       R  += | --- |      >  --- X   (-2 ln Q )    cos (2 pi R )
        a    \  m  /     ---  w   ai         i                i
                 a        i    i


                    1/2
       _     / k T \     --- _            1/2
       v   = | --- |      >  X  (-2 ln Q )    sin (2 pi R )
        a    \  m  /     ---  ai        i                i
                 a        i

    Reference: [West, Estreicher; PRL 96, 22 (2006)]
    """

    # Build dynamical matrix
    rminv = (masses ** -0.5).repeat(3)
    dynamical_matrix = force_constants * rminv[:, None] * rminv[None, :]

    # Solve eigenvalue problem to compute phonon spectrum and eigenvectors
    w2_s, X_is = np.linalg.eigh(dynamical_matrix)

    # First three modes are translational so ignore:
    last_ignore_mode = 3
    # if ignore_negative ignore all modes with frequency < 0.0 threshold (1e-3)
    if ignore_negative:
        last_ignore_mode = len(np.where(w2_s < 1e-3)[0])

    # Check for soft modes
    w2min = w2_s[last_ignore_mode:].min()
    if w2min < 0:
        msg = "Dynamical matrix has negative eigenvalues such as {}".format(w2min)
        if failfast:
            raise ValueError(msg)
        else:
            warn(msg, level=1)

    zeros = w2_s[last_ignore_mode - 3 : last_ignore_mode]
    worst_zero = np.abs(zeros).max()
    if worst_zero > 1e-3:
        msg = "Translational deviate from 0 significantly: {}".format(w2_s[:3])
        if failfast:
            raise ValueError(msg)
        else:
            warn(msg, level=1)

    nw = len(w2_s) - last_ignore_mode
    n_atoms = len(masses)

    w_s = np.sqrt(w2_s[last_ignore_mode:])
    X_acs = X_is[:, last_ignore_mode:].reshape(n_atoms, 3, nw)

    # Assign the amplitudes according to Bose-Einstein distribution
    # or high temperature (== classical) limit
    if quantum:
        hbar = units._hbar * units.J * units.s
        A_s = np.sqrt(hbar * (2 * n_BE(temp, hbar * w_s) + 1) / (2 * w_s))
    else:
        A_s = np.sqrt(temp) / w_s

    if ignore_negative:
        A_s[np.where(np.isnan(w_s))[0]] = 0.0
        w_s[np.where(np.isnan(w_s))[0]] = 0.0

    if plus_minus or gauge_eigenvectors:
        # gauge eigenvectors: largest value always positive
        for ii in range(X_acs.shape[-1]):
            vec = X_acs[:, :, ii]
            max_arg = np.argmax(abs(vec))
            X_acs[:, :, ii] *= np.sign(vec.flat[max_arg])

    if deterministic or plus_minus:
        # create samples by multiplying the amplitude with +/-
        # according to Eq. 5 in PRB 94, 075125

        if plus_minus:
            spread = (-1) ** np.arange(nw)
        else:
            spread = np.ones(nw)

        # Create velocities und displacements from the amplitudes and
        # eigenvectors
        A_s *= spread
        phi_s = 2.0 * np.pi * rng(nw)

        # Assign velocities, sqrt(2) compensates for missing sin(phi) in
        # amplitude for displacement
        v_ac = (w_s * A_s * np.sqrt(2) * np.cos(phi_s) * X_acs).sum(axis=2)
        v_ac /= np.sqrt(masses)[:, None]

        # Assign displacements
        d_ac = (A_s * X_acs).sum(axis=2)
        d_ac /= np.sqrt(masses)[:, None]

    else:
        # compute the gaussian distribution for the amplitudes
        # We need 0 < P <= 1.0 and not 0 0 <= P < 1.0 for the logarithm
        # to avoid (highly improbable) NaN.

        # Box Muller [en.wikipedia.org/wiki/Box–Muller_transform]:
        spread = np.sqrt(-2.0 * np.log(1.0 - rng(nw)))

        # assign amplitudes and phases
        A_s *= spread
        phi_s = 2.0 * np.pi * rng(nw)

        # Assign velocities and displacements
        v_ac = (w_s * A_s * np.cos(phi_s) * X_acs).sum(axis=2)
        v_ac /= np.sqrt(masses)[:, None]

        d_ac = (A_s * np.sin(phi_s) * X_acs).sum(axis=2)
        d_ac /= np.sqrt(masses)[:, None]

    if return_eigensolution:
        return d_ac, v_ac, w2_s, X_is
    # else
    return d_ac, v_ac


def PhononHarmonics(
    atoms,
    force_constants,
    temp,
    rng=np.random,
    quantum=False,
    deterministic=False,
    plus_minus=False,
    gauge_eigenvectors=False,
    failfast=True,
    ignore_negative=False,
):
    r"""Excite phonon modes to specified temperature.

    This will displace atomic positions and set the velocities so as
    to produce a random, phononically correct state with the requested
    temperature.

    Parameters:

    atoms: ase.atoms.Atoms() object
        Grumble
    force_constants: ndarray of size 3N x 3N
        Force constants for the the structure represented by atoms in eV/Å²
    temp: float
        Temperature in eV (T * units.kB)
    rng: Random number generator
        RandomState or other random number generator, e.g., np.random.rand
    quantum: bool
        True for Bose-Einstein distribution, False for Maxwell-Boltzmann
        (classical limit)
    deterministic: bool
        Displace atoms deterministically by fixing the phases
    plus_minus: bool
        Displace atoms with +/- the amplitude accoding to PRB 94, 075125
    gauge_eigenvectors: bool
        gauge eigenvectors such that largest entry is positive
    failfast: bool
        True for sanity checking the phonon spectrum for negative frequencies
        at Gamma.
    ignore_negative: bool
        If True freeze out the imaginary modes
    """

    # Receive displacements and velocities from phonon_harmonics()
    d_ac, v_ac = phonon_harmonics(
        force_constants=force_constants,
        masses=atoms.get_masses(),
        temp=temp,
        rng=rng.rand,
        deterministic=deterministic,
        plus_minus=plus_minus,
        quantum=quantum,
        gauge_eigenvectors=gauge_eigenvectors,
        failfast=failfast,
        ignore_negative=ignore_negative,
        return_eigensolution=False,
    )

    # Assign new positions (with displacements) and velocities
    atoms.positions += d_ac
    atoms.set_velocities(v_ac)
