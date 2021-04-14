""" high-level access to mode projection functionality """

import numpy as np
import scipy.linalg as la
from vibes.dynamical_matrix import fc2dynmat, get_frequencies
from vibes.helpers import Timer, lazy_property, warn


Timer.prefix = "mode_analysis"


class SimpleModeProjection:
    """provide tools to perform (simple) mode projection w/o using symmetry etc."""

    def __init__(self, atoms, force_constants, verbose=True):
        """instantiate

        Args:
            atoms (ase.Atoms): the reference structure
            force_constants (np.ndarray): the [3N, 3N] force constants
            vebose (bool): to request verbosity
        """
        self.atoms = atoms.copy()
        self.masses = np.asarray(self.atoms.get_masses())

        # verify shape of force_constants
        Na = len(self.masses)
        force_constants = np.squeeze(force_constants)
        is_shape = np.shape(force_constants)
        is_size = np.size(force_constants)
        fc_shape = (3 * Na, 3 * Na)
        Na = len(self.atoms)
        if is_shape == fc_shape:
            self.force_constants = np.asarray(force_constants)
        elif is_size == ((3 * Na) ** 2):
            self.force_constants = np.reshape(force_constants, fc_shape)
        else:
            print((3 * Na) ** 2)
            msg = f"FIXME Force constants shape {is_shape} not supported (yet)."
            warn(msg, level=2)

    @lazy_property
    def dynamical_matrix(self):
        return fc2dynmat(self.force_constants, self.masses)

    @lazy_property
    def eigenvectors(self):
        """ return eigenvectors """
        _, evecs = la.eigh(self.dynamical_matrix)
        return evecs

    @lazy_property
    def mode_projector(self):
        """eigenvectors for mode projection"""
        return self.eigenvectors.T

    @lazy_property
    def omegas(self):
        """ return angular frequencies """
        return get_frequencies(self.dynamical_matrix)

    def project(self, array, mass_weight=0.0, info=None):
        """perform mode projection on [Nt, Na, 3] shaped array

        mass_weight: Exponent for mass weighting when applying the projector
             0.0: No mass weighting
             0.5: for positions and velocities
            -0.5: for forces

            A_s = M ** (mass_weight) * P @ A_x

        Args:
            array (np.ndarray): positions, velocities or forces, ...
            mass_weight_prefactor (float): prefactor for mass weighting
            info (str): Info message for timer
        """
        timer = Timer(f"Project {info} onto modes")

        assert np.shape(array)[1:] == (len(self.atoms), 3)

        P = self.mode_projector

        _supported_mass_weights = (0, 0.5, -0.5)
        if mass_weight in _supported_mass_weights:
            M = self.masses[None, :, None] ** mass_weight
            array = M * array
        elif mass_weight is None:
            pass
        else:
            msg = f"`mass_weight` {mass_weight} not in {_supported_mass_weights}"
            warn(msg, level=2)

        result = np.array([P @ (f.flatten()) for f in np.asarray(array)])
        timer()

        return result
