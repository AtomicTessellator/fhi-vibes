"""Functions used to wrap around HiLDe Phonopy/Phono3py functions"""
from pathlib import Path

from vibes.helpers.dict import AttributeDict
from vibes.helpers.k_grid import update_k_grid
from vibes.relaxation.context import RelaxationContext
from vibes.settings import Settings


def run(atoms, calculator, kpt_density=None, relax_settings=None, fw_settings=None):
    """Creates a Settings object and passes it to the bootstrap function

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Atoms object of the primitive cell
    calculator: ase.calculators.calulator.Calculator
        Calculator for the force calculations
    kpt_density: float
        k-point density for the MP-Grid
    md_settings: dict
        kwargs for md setup
    fw_settings: dict
        FireWork specific settings

    Returns
    -------
    completed: bool
        True if the workflow completed
    """
    workdir = relax_settings.get("workdir", None)
    if workdir:
        workdir = Path(workdir)
        trajectory = workdir / "trajectory.son"
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        trajectory = "trajectory.son"
        workdir = Path(".")

    atoms.write(str(workdir / "geometry.in"), format="aims", scaled=True)

    if calculator.name.lower() == "aims":
        if kpt_density is not None:
            update_k_grid(atoms, calculator, kpt_density, even=True)
        calculator.parameters.pop("sc_init_iter", None)

    settings = Settings(settings_file=None)
    settings._settings_file = "relaxation.in"
    settings["relaxation"] = AttributeDict(relax_settings)

    if calculator.name.lower() == "aims":
        settings["control"] = AttributeDict(calculator.parameters.copy())
        settings["basissets"] = AttributeDict(
            {"default": calculator.parameters.pop("species_dir").split("/")[-1]}
        )

    settings["geometry"] = AttributeDict({"file": str(workdir / "geometry.in")})
    settings.write(f"{workdir}/relaxation.in")

    ctx = RelaxationContext(settings, workdir, trajectory)
    ctx.calculator = calculator

    return ctx.run()
