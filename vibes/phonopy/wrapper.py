"""
A leightweight wrapper for Phonopy()
"""

import json
from pathlib import Path

import numpy as np
from phonopy import Phonopy

from vibes import konstanten as const
from vibes.helpers import brillouinzone as bz
from vibes.helpers import talk, warn
from vibes.helpers.numerics import get_3x3_matrix
from vibes.materials_fp.material_fingerprint import (
    get_phonon_bs_fingerprint_phononpy,
    to_dict,
)
from vibes.phonopy.utils import get_supercells_with_displacements
from vibes.spglib.wrapper import map_unique_to_atoms
from vibes.structure.convert import to_Atoms, to_phonopy_atoms

from ._defaults import defaults


def prepare_phonopy(
    atoms,
    supercell_matrix,
    fc2=None,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    trigonal=defaults.is_trigonal,
    is_diagonal=defaults.is_diagonal,
    is_plusminus=defaults.is_plusminus,
    wrap=False,
):
    """Create a Phonopy object

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The primitive cell for the calculation
    supercell_matrix: np.ndarray
        The supercell matrix for generating the supercell
    fc2: np.ndarray
        The second order force constants
    displacement: float
        The magnitude of the phonon displacement
    symprec: float
        tolerance for determining the symmetry/spacegroup of the primitive cell
    trigonal: bool
        If True use trigonal displacements
    is_diagonal: bool
        If True use diagonal displacements
    is_plusminus: bool
        trigger to compute both + and - the displacement

    Returns
    -------
    phonon: phonopy.Phonopy
        The phonopy object corresponding to the parameters
    """
    ph_atoms = to_phonopy_atoms(atoms, wrap=wrap)
    supercell_matrix = get_3x3_matrix(supercell_matrix)

    phonon = Phonopy(
        ph_atoms,
        supercell_matrix=np.transpose(supercell_matrix),
        symprec=symprec,
        is_symmetry=True,
        factor=const.omega_to_THz,
    )

    phonon.generate_displacements(
        distance=displacement,
        is_plusminus=is_plusminus,
        # is_diagonal=False is chosen to be in line with phono3py, see
        # https://github.com/atztogo/phono3py/pull/15
        is_diagonal=is_diagonal,
        is_trigonal=trigonal,
    )

    if fc2 is not None:
        phonon.set_force_constants(fc2)

    return phonon


def preprocess(
    atoms,
    supercell_matrix,
    displacement=defaults.displacement,
    symprec=defaults.symprec,
    trigonal=defaults.is_trigonal,
    is_plusminus=defaults.is_plusminus,
    **kwargs,
):
    """Generate phonopy objects and return displacements as Atoms objects

    Parameters
    ----------
    atoms: ase.atoms.Atoms
        The primitive cell for the calculation
    supercell_matrix: np.ndarray
        The supercell matrix for generating the supercell
    displacement: float
        The magnitude of the phonon displacement
    symprec: float
        tolerance for determining the symmetry/spacegroup of the primitive cell
    trigonal: bool
        If True use trigonal displacements
    is_plusminus: bool
        use +/- displacements

    Returns
    -------
    phonon: phonopy.Phonopy
        The phonopy object with displacement_dataset, and displaced supercells
    supercell: ase.atoms.Atoms
        The undisplaced supercell
    supercells_with_disps: list of ase.atoms.Atoms
        All of the supercells with displacements
    """
    supercell_matrix = get_3x3_matrix(supercell_matrix)

    phonon = prepare_phonopy(
        atoms,
        supercell_matrix,
        displacement=displacement,
        symprec=symprec,
        trigonal=trigonal,
        is_plusminus=is_plusminus,
    )

    return get_supercells_with_displacements(phonon)


# TARP: This is depcricated and should not be used
# def get_force_constants(phonon, force_sets=None):
#     """ Take a Phonopy object, produce force constants from the given forces and
#         return in usable shape (3N, 3N) insated of (N, N, 3, 3)

#     """
#     n_atoms = phonon.get_supercell().get_number_of_atoms()

#     phonon.produce_force_constants(force_sets)

#     force_constants = phonon.get_force_constants()

#     if force_constants is not None:
#         # convert forces from (N, N, 3, 3) to (3*N, 3*N)
#         force_constants = (
#             phonon.get_force_constants().swapaxes(1, 2).reshape(2 * (3 * n_atoms,))
#         )
#         return force_constants
#     # else
#     raise ValueError("Force constants not yet created, specify force_sets.")


def get_dos(
    phonon,
    total=True,
    q_mesh=defaults.q_mesh,
    freq_min="auto",
    freq_max="auto",
    freq_pitch=None,
    tetrahedron_method=True,
    write=False,
    filename="total_dos.dat",
    force_sets=None,
    direction=None,
    xyz_projection=False,
):
    """Compute the DOS (and save to file)

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The phonopy object with calculated force constants if force_sets is None
    total: bool
        If True calculate the total density of states
    q_mesh: np.ndarray
        size of the interpolated q-point mesh
    freq_min: float
        minimum frequency to calculate the density of states on
    freq_max: float
        maximum frequency to calculate the density of states on
    freq_pitch: float
        spacing between frequency points
    tetrahedron_method: bool
        If True use the tetrahedron method to calculate the DOS
    write: bool
        If True save the DOS to a file
    filename: str
        Path to write to write the dos to
    force_sets: np.ndarray
        Force sets to calculate the force constants with
    direction: np.ndarray
        Specific projection direction. This is specified three values along basis
        vectors or the primitive cell. Default is None, i.e., no projection.
    xyz_projection: bool
        If True project along Cartesian directions

    Returns
    ------
    np.ndarray
        The total or projected phonon density of states
    """

    if force_sets is not None:
        phonon.produce_force_constants(force_sets)

    if total:
        phonon.run_mesh(q_mesh)

        if freq_max == "auto":
            freq_max = phonon.get_mesh()[2].max() * 1.05
        if freq_min == "auto":
            freq_min = phonon.get_mesh()[2].min()
            if freq_min < 0.0:
                freq_min *= 1.05
            else:
                freq_min = 0.0
        phonon.run_total_dos(
            freq_min=freq_min,
            freq_max=freq_max,
            freq_pitch=freq_pitch,
            use_tetrahedron_method=tetrahedron_method,
        )

        if write:
            phonon.write_total_dos()
            Path("total_dos.dat").rename(filename)

        return phonon.get_total_dos_dict()

    phonon.run_mesh(q_mesh, is_mesh_symmetry=False, with_eigenvectors=True)

    if freq_max == "auto":
        freq_max = phonon.get_mesh()[2].max() * 1.05

    phonon.run_projected_dos(
        freq_min=freq_min,
        freq_max=freq_max,
        freq_pitch=freq_pitch,
        use_tetrahedron_method=tetrahedron_method,
        direction=direction,
        xyz_projection=xyz_projection,
    )
    if write:
        phonon.write_projected_dos()
        Path("projected_dos.dat").rename(filename)
    return phonon.get_projected_dos_dict()


def set_bandstructure(phonon, paths=None):
    """Compute bandstructure for given path and attach to phonopy object

    Parameters
    ----------
    phonon: phonopy.Phonopy
        Phonopy object with calculated force constants if force_Sets is None
    paths: list of str
        List of high-symmetry point paths e.g. ['GXSYGZURTZ', 'YT', 'UX', 'SR']
    """
    if isinstance(paths, str):
        paths = [paths]

    bands, labels = bz.get_bands_and_labels(to_Atoms(phonon.primitive), paths)

    phonon.run_band_structure(bands, labels=labels)


def get_bandstructure(phonon, paths=None, force_sets=None):
    """Compute bandstructure for given path

    Parameters
    ----------
    phonon: phonopy.Phonopy
        Phonopy object with calculated force constants if force_Sets is None
    paths: list of str
        List of high-symmetry point paths e.g. ['GXSYGZURTZ', 'YT', 'UX', 'SR']
    force_sets: np.ndarray
        set of forces to calculate force constants with

    Returns
    -------
    band_structure_dict: dict
        dictionary of the band structure
    labels: list of str
        labels for the start and end of those force constants
    """
    if force_sets is not None:
        phonon.produce_force_constants(force_sets)

    _, labels = bz.get_bands_and_labels(to_Atoms(phonon.primitive), paths)
    set_bandstructure(phonon, paths)

    return (phonon.get_band_structure_dict(), labels)


def plot_thermal_properties(
    phonon, file="thermal_properties.pdf", t_step=20, t_max=1000, t_min=0
):
    """plot thermal properties to pdf file"""

    phonon.set_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)
    plt = phonon.plot_thermal_properties()

    try:
        plt.savefig(file)
    except (RuntimeError, FileNotFoundError):
        warn("saving the thermal properties not possible, latex probably missing?")


def plot_bandstructure(phonon, file="bandstructure.pdf", paths=None, force_sets=None):
    """Plot bandstructure for given path and save to file

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The phonopy object with calculated force constants if force_sets is None
    filename: str
        file name to store the pdf of the band structure
    paths: list of str
        List of high-symmetry point paths e.g. ['GXSYGZURTZ', 'YT', 'UX', 'SR']
    force_sets: np.ndarray
        Force sets to calculate the force constants with
    """

    set_bandstructure(phonon, paths)

    plt = phonon.plot_band_structure()

    try:
        plt.savefig(file)
    except (RuntimeError, FileNotFoundError):
        warn("saving the phonon dispersion not possible, latex probably missing?")


def plot_bandstructure_and_dos(
    phonon,
    q_mesh=defaults.q_mesh,
    partial=False,
    file="bandstructure_dos.pdf",
    run_mesh=True,
):
    """Plot bandstructure and PDOS

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The phonopy object with calculated force constants
    q_mesh: np.ndarray
        size of the interpolated q-point mesh
    partial: bool
        If True use projected density of states
    file: str
        Path to save the the plot to
    """

    set_bandstructure(phonon)

    if partial:
        if run_mesh:
            phonon.run_mesh(q_mesh, with_eigenvectors=True, is_mesh_symmetry=False)
            phonon.run_projected_dos(use_tetrahedron_method=True)
        pdos_indices = map_unique_to_atoms(phonon.get_primitive())
    else:
        if run_mesh:
            phonon.run_mesh(q_mesh)  # , with_eigenvectors=True,)
            phonon.run_total_dos(use_tetrahedron_method=True)
        pdos_indices = None

    plt = phonon.plot_band_structure_and_dos(pdos_indices=pdos_indices)

    try:
        plt.savefig(file)
    except (RuntimeError, FileNotFoundError):
        warn("saving the phonon DOS not possible, latex probably missing?")


def summarize_bandstructure(phonon, fp_file=None):
    """Print a concise symmary of the bandstructure fingerprint

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The phonopy object with calculated force constants if force_sets is None
    fp_file: str
        Path to save the fingerprint to

    Returns
    ------
    gamma_freq: np.ndarray
        The frequencies at the gamma point
    max_freq: float
        The maximum frequency at the gamma point
    """
    from vibes.konstanten.einheiten import THz_to_cm

    set_bandstructure(phonon)

    qpts = np.array(phonon.band_structure.qpoints).reshape(-1, 3)

    freq = np.array(phonon.band_structure.frequencies).reshape(qpts.shape[0], -1)

    gamma_freq = freq[np.where((qpts == np.zeros(3)).all(-1))[0][0]]
    max_freq = np.max(freq.flatten())

    if fp_file:
        talk(f"Saving the fingerprint to {fp_file}")
        fp = get_phonon_bs_fingerprint_phononpy(phonon, binning=False)
        fp_dict = to_dict(fp)
        for key, val in fp_dict.items():
            fp_dict[key] = val.tolist()
        with open(fp_file, "w") as outfile:
            json.dump(fp_dict, outfile, indent=4)

    mf = max_freq
    mf_cm = mf * THz_to_cm
    talk(f"The maximum frequency is: {mf:.3f} THz ({mf_cm:.3f} cm^-1)")
    talk(f"The frequencies at the gamma point are:")
    talk(f"              THz |        cm^-1")
    p = lambda ii, freq: print(f"{ii+1:3d}: {freq:-12.5f} | {freq*THz_to_cm:-12.5f}")
    for ii, freq in enumerate(gamma_freq[:6]):
        p(ii, freq)
    for _ in range(3):
        print("  .")
    for ii, freq in enumerate(gamma_freq[-3:]):
        p(len(gamma_freq) - 3 + ii, freq)
    return gamma_freq, max_freq


def get_animation(phonon, q_point, filename):
    """Gets the animation file at a q_point

    Parameters
    ----------
    phonon: phonopy.Phonopy
        The phonopy object with calculated force constants if force_sets is None
    q_point: np.ndarray
        q-point to write the animation file on
    filename: str
        Path to animation file output
    """
    return phonon.write_animation(q_point=q_point, filename=filename)


def get_debye_temperature(
    phonon=None, freq_pitch=5e-3, q_mesh=defaults.q_mesh, tetrahedron_method=True
):
    """Calculate the Debye Temperature from the Phonon Density of States

    Formulas taken from: J. Appl. Phys. 101, 093513 (2007)

    Args:
        phonon (phonopy.Phonopy): The phonon calculation
        freq_pitch (double): Energy spacing for calculating the total DOS
        q_mesh (np.ndarray): size of the interpolated q-point mesh
        tetrahedron_method (bool): If True use the tetrahedron method to calculate the DOS
    Returns:
        $\\Theta_P$ (float):
            Average phonon temperature
        $\\Theta_{D\\infty} (float):
            T -> $\\infty$ limiting magnitude of the Debye Temperature
        $\\Theta_D$ (float):
            The phonopy debye temperature as fitting the phonon DOS to Debye Model
    """
    dos = get_dos(
        phonon,
        q_mesh=q_mesh,
        freq_min=0.0,
        freq_pitch=freq_pitch,
        tetrahedron_method=tetrahedron_method,
    )
    ener = dos["frequency_points"] * const.THzToEv
    gp = dos["total_dos"]
    eps_p_1 = np.trapz(gp * ener, ener) / np.trapz(gp, ener)
    eps_p_2 = np.trapz(gp * ener ** 2.0, ener) / np.trapz(gp, ener)

    phonon.set_Debye_frequency()
    omgea_d = phonon.get_Debye_frequency() * 1e12 * np.pi * 2.0

    theta_p = eps_p_1 / const.kB
    theta_d_infty = np.sqrt(5.0 / 3.0 * eps_p_2) / const.kB
    theta_d = omgea_d * const.HBAR / (const.kB * const.EV)

    return theta_p, theta_d_infty, theta_d