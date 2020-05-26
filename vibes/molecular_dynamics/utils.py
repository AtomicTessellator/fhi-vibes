""" helper utilities:
    - FCCalculator for using force constants to compute forces
    - Logger for tracking custom MD """
from pathlib import Path

from ase.calculators.calculator import Calculator, PropertyNotImplementedError

from vibes import son
from vibes.helpers.converters import input2dict
from vibes.helpers.displacements import get_dR


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
    """ MD logger class to write vibes trajectory files """

    def __init__(self, atoms, trajectory_file, metadata=None, overwrite=False):
        """initialize

        Parameters
        ----------
        atoms: ase.atoms.Atoms
            Atoms of the reference structure
        trajectory_file: str or Path
            path to the trajectory file
        metadata: dict
            metadata for the MD run
        overwrite: bool
            If true overwrite the trajectory file
        """

        if not metadata:
            metadata = {}

        self.trajectory_file = trajectory_file
        if Path(trajectory_file).exists() and overwrite:
            Path(trajectory_file).unlink()
            print(f"** {trajectory_file} deleted.")

        son.dump(
            {**metadata, **input2dict(atoms)}, self.trajectory_file, is_metadata=True
        )

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

        try:
            stress = atoms.get_stress(voigt=False)
        except PropertyNotImplementedError:
            stress = None

        dct = {
            "atoms": {
                "info": info,
                "cell": atoms.cell[:],
                "positions": atoms.positions,
                "velocities": atoms.get_velocities(),
            },
            "calculator": {
                "forces": atoms.get_forces(),
                "energy": atoms.get_kinetic_energy(),
                "stress": stress,
            },
        }

        son.dump(dct, self.trajectory_file)
