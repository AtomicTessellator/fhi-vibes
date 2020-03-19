"""Functions used to wrap around HiLDe Phonopy/Phono3py functions"""
from pathlib import Path

import numpy as np

from vibes.helpers.attribute_dict import AttributeDict
from vibes.helpers.converters import dict2atoms
from vibes.helpers.k_grid import k2d, update_k_grid
from vibes.helpers.numerics import get_3x3_matrix
from vibes.helpers.supercell import make_supercell
from vibes.molecular_dynamics.context import MDContext
from vibes.phonopy.postprocess import extract_results, postprocess
from vibes.phonopy.utils import remap_force_constants
from vibes.scripts.create_samples import generate_samples
from vibes.settings import Settings
from vibes.trajectory import reader


def run(atoms, calc, kpt_density=None, md_settings=None, fw_settings=None):
    """Creates a Settings object and passes it to the bootstrap function

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        Atoms object of the primitive cell
    calc: ase.calculators.calulator.Calculator
        Calculator for the force calculations
    kpt_density: float
        k-point density for the MP-Grid
    md_settings: dict
        kwargs for md setup
    fw_settings: dict
        FireWork specific settings

    Returns
    -------
    outputs: dict
        The output of vibes.phonopy.workflow.bootstrap for phonopy and phono3py
    """
    workdir = md_settings.get("workdir", None)
    if workdir:
        workdir = Path(workdir)
        trajectory = workdir / "trajectory.son"
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        trajectory = None
        workdir = Path(".")

    if "phonon_file" not in md_settings and not Path(workdir / "geometry.in").exists():
        sc_matrix = md_settings.pop("supercell_matrix", np.eye(3))
        supercell = make_supercell(atoms, sc_matrix)
        supercell.write(str(workdir / "supercell.in"), format="aims", scaled=True)
        atoms.write(str(workdir / "unitcell.in"), format="aims", scaled=True)
        supercell.write(str(workdir / "geometry.in"), format="aims", scaled=True)
    else:
        if Path(md_settings["phonon_file"]).parent.name == "converged":
            sc_list = Path(md_settings["phonon_file"]).parents[1].glob("sc_n*")
            sc_list = sorted(sc_list, key=lambda s: int(str(s).split("_")[-1]))
            md_settings["phonon_file"] = str(sc_list[-1] / "phonopy/trajectory.son")

        _, ph_metadata = reader(md_settings["phonon_file"], get_metadata=True)
        traj_dir = str(Path(md_settings["phonon_file"]).parent)
        traj_file = str(Path(md_settings["phonon_file"]).name)
        phonon = postprocess(traj_file, traj_dir, verbose=False)
        extract_results(phonon, output_dir=workdir / "phonopy_output")

        atoms_dict = ph_metadata["primitive"]["atoms"]
        ph_atoms = dict2atoms(atoms_dict, ph_metadata["calculator"], False)
        calc = ph_atoms.calc
        kpt_density = k2d(ph_atoms, calc.parameters["k_grid"])
        if not Path(workdir / "geometry.in").exists():
            sc = dict2atoms(ph_metadata["supercell"]["atoms"])
            sc_matrix = md_settings.pop(
                "supercell_matrix",
                get_3x3_matrix(ph_metadata["Phonopy"]["supercell_matrix"]),
            )
            supercell = make_supercell(ph_atoms, sc_matrix)

            fc = remap_force_constants(
                phonon.get_force_constants(), ph_atoms, sc, supercell, two_dim=True
            )

            atoms_md = generate_samples(
                supercell,
                md_settings["temperature"] * 1.05,
                n_samples=1,
                force_constants=fc,
                rattle=False,
                quantum=False,
                deterministic=False,
                plus_minus=False,
                gauge_eigenvectors=False,
                ignore_negative=True,
                sobol=False,
                random_seed=np.random.randint(2 ** 32),
                propagate=0,
                format="aims",
            )[0]

            supercell.write(str(workdir / "supercell.in"), format="aims", scaled=True)
            ph_atoms.write(str(workdir / "unitcell.in"), format="aims", scaled=True)
            info_str = atoms_md.info.pop("info_str")
            atoms_md.write(
                str(workdir / "geometry.in"),
                info_str=info_str,
                format="aims",
                scaled=True,
                velocities=True,
            )

    if kpt_density is not None:
        update_k_grid(atoms, calc, kpt_density, even=True)

    calc.parameters.pop("sc_init_iter", None)
    settings = Settings(settings_file=None)
    settings._settings_file = "md.in"
    settings["md"] = AttributeDict(md_settings)
    settings["basissets"] = AttributeDict(
        {"default": calc.parameters.pop("species_dir").split("/")[-1]}
    )
    settings["control"] = AttributeDict(calc.parameters.copy())
    settings["geometry"] = AttributeDict(
        {
            "file": str(workdir / "geometry.in"),
            "primitive": str(workdir / "unitcell.in"),
            "supercell": str(workdir / "supercell.in"),
        }
    )

    ctx = MDContext(settings, workdir, trajectory)

    return ctx.run()
