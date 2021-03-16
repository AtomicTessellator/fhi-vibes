""" high-level access to mode projection functionality """

import numpy as np
import scipy.linalg as la
from ase import Atoms

from vibes.helpers import Timer, lazy_property, progressbar, warn
from vibes.helpers.displacements import get_dUdt, get_U
from vibes.helpers.force_constants import DynamicalMatrix

from .dynamical_matrix import fc2dynmat, get_frequencies


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


class HarmonicAnalysis:
    """ provide tools to perform harmonic analysis in periodic systems """

    def __init__(
        self,
        primitive: Atoms,
        supercell: Atoms,
        force_constants: np.ndarray,
        verbose: bool = False,
    ):
        """Initialize unit cell, supercell, force_constants and lattice points.

        Args:
            primitive: The unit cell structure
            supercell: The supercell structure
            force_constants: f.c. in phonopy format [Np, Ns, 3, 3]
            verbose: control verbosity level
        """

        self.dynamical_matrix_obj = DynamicalMatrix(
            force_constants=force_constants, primitive=primitive, supercell=supercell
        )

    @property
    def dmx(self):
        return self.dynamical_matrix_obj

    def get_Utsq(self, trajectory: list):
        """ Get the mode projected positions, weighted by mass.

        Args:
            trajectory: The trajectory to work over

        Returns:
            U_tsq, V_tsq: [Nt, Ns, Nq] arrays of projected displacements and velocities

        """
        print(f"Project trajectory onto modes:")
        proj = self.dmx.e_sqI
        shape = [len(trajectory), proj.shape[0], proj.shape[1]]
        Utsq = np.zeros(shape, dtype=complex)
        Vtsq = np.zeros(shape, dtype=complex)

        atoms0 = self.dmx.supercell

        for ii in progressbar(range(len(trajectory))):
            atoms = trajectory[ii]
            Utsq[ii] = proj @ get_U(atoms, atoms0=atoms0).flatten()
            Vtsq[ii] = proj @ get_dUdt(atoms).flatten()

        return Utsq, Vtsq

    def get_Ztsq(self, trajectory: list):
        """ Return the imaginary mode amplitude for [t, s, q]

        Args:
            trajectory: The trajectory to work over

        """
        omegas = self.dmx.w_sq
        U_tsq, V_tsq = self.get_Utsq(trajectory)

        Z_tsq = V_tsq - 1.0j * omegas[None, :, :] * U_tsq

        return Z_tsq

    def project(self, trajectory, times=None):
        """ perform mode projection for atoms objects in trajectory

        Args:
            trajectory: The trajectory to work over
            times: The times at each point in the trajectory

        Returns:
            A_tsq2: Amplitdues [Nt, Ns, Nq]
            phi_tsq: Angles [Nt, Ns, Nq]
            E_tsq: Energies [Nt, Ns, Nq]

        """

        timer = Timer("Perform mode analysis for trajectory")

        if isinstance(trajectory, Atoms):
            trajectory = [trajectory]

        w2_sq = self.dmx.w2_sq
        Z_tsq = self.get_Ztsq(trajectory)

        A_tsq2 = w2_sq[None, :, :] ** -1 * abs(Z_tsq) ** 2
        # phi_qst = get_phi_qst(U_tsq, V_tsq, self.omegas, in_times=times)

        E_tsq = 0.5 * w2_sq[None, :, :] * A_tsq2

        timer("project trajectory")

        return A_tsq2, E_tsq
