""" helper utilities:
    - FCCalculator for using force constants to compute forces
    - Logger for tracking custom MD """
from pathlib import Path

from ase.calculators.calculator import Calculator

from hilde import son
from hilde.harmonic_analysis.displacements import get_dR
from hilde.trajectory import input2dict


def get_F(dR, force_constants):
    """Compute force from force_constants @ displacement

    Parameters
    ----------
    dR: np.ndarray
        The displacement matrix
    force_constants: np.ndarray
        The Force constant Matrix

    Returns
    -------
    np.ndarray
        The harmonic forces
    """
    return -(force_constants @ dR.flatten()).reshape(dR.shape)


class FCCalculator(Calculator):
    """ Calculator that uses (2nd order) force constants to compute forces. """

    def __init__(self, ref_atoms, force_constants, **kwargs):
        """Initializor

        Parameters
        ----------
        ref_atoms: ase.atoms.Atoms
            Reference structure (where harmonic forces are zero)
        force_constant: np.ndarray
            The force constant matrix
        """
        super().__init__(**kwargs)
        self.implemented_properties = ["forces"]

        self.force_constants = force_constants
        self.atoms0 = ref_atoms

    def get_forces(self, atoms=None):
        """Get the harmonic forces

        Parameters
        ----------
        atoms: ase.atoms.Atoms
            displaced structure (only positions can be different w/rt ref_atoms)

        Returns
        -------
        np.ndarray
            The harmonic forces
        """
        dR = get_dR(atoms, self.atoms0)
        return get_F(dR, self.force_constants)


class MDLogger:
    """ MD logger class to write hilde trajectory files """

    def __init__(self, atoms, trajectory, metadata=None, overwrite=False):
        """initialize

        Parameters
        ----------
        atoms: ase.atoms.Atoms
            Atoms of the reference structure
        trajectory: str or Path
            path to the trajectory file
        metadata: dict
            metadata for the MD run
        overwrite: bool
            If true overwrite the trajectory file
        """

        if not metadata:
            metadata = {}

        self.trajectory = trajectory
        if Path(trajectory).exists() and overwrite:
            Path(trajectory).unlink()
            print(f"** {trajectory} deleted.")

        son.dump({**metadata, **input2dict(atoms)}, self.trajectory, is_metadata=True)

    def __call__(self, atoms, info=None):
        """Log the current step to the trajectory

        Parameters
        ----------
        atoms: ase.atoms.Atoms
            Atoms of the current step
        info: dict
            additional information to add to the update
        """
        if info is None:
            info = {}
        dct = {
            "atoms": {
                "info": info,
                "positions": atoms.positions,
                "velocities": atoms.get_velocities(),
            },
            "calculator": {
                "forces": atoms.get_forces(),
                "energy": atoms.get_kinetic_energy(),
            },
        }

        son.dump(dct, self.trajectory)
