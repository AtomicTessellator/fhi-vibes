"""Postprocessing for phonon clalculations"""
from pathlib import Path
from shutil import copyfile, rmtree

import numpy as np

from vibes.helpers.converters import atoms2dict
from vibes.helpers.k_grid import update_k_grid
from vibes.helpers.paths import cwd
from vibes.materials_fp.material_fingerprint import (
    get_phonon_dos_fingerprint_phononpy,
    fp_tup,
    scalar_product,
)
from vibes.phonon_db.row import phonon_to_dict
from vibes.phonopy.wrapper import preprocess as ph_preprocess
from vibes.settings import Settings
from vibes.structure.convert import to_Atoms_db
from vibes.trajectory import reader


def time2str(n_sec):
    """
    Converts a number of seconds into a time string
    Args:
        n_secs (int): A time presented as a number of seconds

    Returns
    -------
    time_str: str
        A string representing a specified time
    """
    secs = int(n_sec % 60)
    mins = int(n_sec / 60) % 60
    hrs = int(n_sec / 3600) % 24
    days = int(n_sec / 86400)
    return f"{days}-{hrs}:{mins}:{secs}"


def get_base_work_dir(wd):
    """
    Converts wd to be it's base (no task specific directories)
    Args:
        wd (str): Current working directory

    Returns:
        (str): The base working directory for the workflow
    """
    wd_list = wd.split("/")
    # remove analysis directories from path
    while "phonopy_analysis" in wd_list:
        wd_list.remove("phonopy_analysis")

    while "phono3py_analysis" in wd_list:
        wd_list.remove("phono3py_analysis")

    while "phonopy" in wd_list:
        wd_list.remove("phonopy")

    while "phono3py" in wd_list:
        wd_list.remove("phono3py")

    # Remove all "//" from the path
    while "" in wd_list:
        wd_list.remove("")

    # If starting from root add / to beginning of the path
    if wd[0] == "/":
        wd_list = [""] + wd_list

    # Remove "sc_natoms_???" to get back to the base directory
    if len(wd_list[-1]) > 10 and wd_list[-1][:10] == "sc_natoms_":
        return "/".join(wd_list[:-1])
    return "/".join(wd_list)


def get_memory_expectation(new_supercell, calc, k_pt_density, workdir):
    """Runs a dry_run of the new calculation and gets the estimated memory usage

    Parameters
    ----------
    new_supercell: ase.atoms.Atoms
        The structure to get the memory estimation for
    calc: ase.atoms.Calculator
        The ASE Calculator to be used
    k_pt_density: list of floats
        The k-point density in all directions
    workdir: str
        Path to working directory

    Returns
    -------
    float:
        The expected memory of the calculation, scaling based on empirical tests
    """
    settings = Settings()
    calc.parameters["dry_run"] = True
    calc.parameters.pop("use_local_index", None)
    calc.parameters.pop("load_balancing", None)
    calc.command = settings.machine.aims_command
    bs_base = settings.machine.basissetloc
    calc.parameters["species_dir"] = (
        bs_base + "/" + calc.parameters["species_dir"].split("/")[-1]
    )
    calc.parameters["compute_forces"] = False
    update_k_grid(new_supercell, calc, k_pt_density, even=True)
    new_supercell.set_calculator(calc)
    mem_expect_dir = workdir + "/.memory_expectation"
    with cwd(mem_expect_dir, mkdir=True):
        try:
            new_supercell.calc.calculate()
        except RuntimeError:
            calc.parameters["dry_run"] = False
        lines = open("aims.out").readlines()
    rmtree(mem_expect_dir)

    for line in lines:
        if "Maximum number of basis functions" in line:
            n_basis = int(line.split(":")[1])
            continue

        if "Number of Kohn-Sham states (occupied + empty)" in line:
            n_states = int(line.split(":")[1])
            continue

        if "Size of matrix packed + index [n_hamiltonian_matrix_size]" in line:
            n_hamiltonian_matrix_size = int(line.split(":")[1])
            break

        if "Number of k-point" in line:
            n_kpt = int(line.split(":")[1])
            continue
    n_spin = 1
    total_mem = n_basis * n_basis
    total_mem += n_hamiltonian_matrix_size
    total_mem += n_basis * n_states * n_spin * n_kpt
    total_mem *= 16 * 10
    return total_mem


def check_phonon_conv(dos_fp, prev_dos_fp, conv_crit):
    """
    Checks if the density of state finger prints are converged
    Args:
        dos_fp (MaterialsFingerprint): Current fingerprint
        prev_dos_fp (MaterialsFingerprint): Fingerprint of the previous step
        conv_crit (float): convergence criteria

    Returns:
        (bool): True if conv_criteria is met
    """
    if not isinstance(prev_dos_fp, fp_tup):
        for ll in range(4):
            prev_dos_fp[ll] = np.array(prev_dos_fp[ll])
        prev_dos_fp = fp_tup(
            prev_dos_fp[0], prev_dos_fp[1], prev_dos_fp[2], prev_dos_fp[3]
        )
    return (
        scalar_product(dos_fp, prev_dos_fp, col=1, pt=0, normalize=False, tanimoto=True)
        >= conv_crit
    )


def get_converge_phonon_update(
    workdir,
    trajectory,
    calc_times,
    ph,
    conv_crit=0.95,
    prev_dos_fp=None,
    init_workdir="./",
    **kwargs,
):
    """
    Check phonon convergence and set up future calculations after a phonon calculation
    Args:
        func (str): Path to the phonon analysis function
        func_fw_out (str): Path to this function
        args (list): list arguments passed to the phonon analysis
        fw_settings (dict): Dictionary for the FireWorks specific systems
        kwargs (dict): Dictionary of keyword arguments
            Mandatory Keys:
                outputs: The Phonopy object from post-processing
                serial (bool): If True use a serial calculation
                init_workdir (str): Path to the base phonon force calculations
                trajectory (str): trajectory file name

    Returns:
        (FWAction): Increases the supercell size or adds the phonon_dict to the spec
    """
    calc_time = np.sum(calc_times)

    traj = f"{workdir}/{trajectory}"

    _, metadata = reader(traj, True)
    calc_dict = metadata["calculator"]
    calc_dict["calculator"] = calc_dict["calculator"].lower()

    # Calculate the phonon DOS
    ph.set_mesh([45, 45, 45])
    if np.any(ph.get_frequencies([0.0, 0.0, 0.0]) < -1.0e-1):
        raise ValueError("Negative frequencies at Gamma, terminating workflow here.")
    if prev_dos_fp:
        de = prev_dos_fp[0][0][1] - prev_dos_fp[0][0][0]
        min_f = prev_dos_fp[0][0][0] - 0.5 * de
        max_f = prev_dos_fp[0][0][-1] + 0.5 * de
        ph.run_total_dos(freq_min=min_f, freq_max=max_f, freq_pitch=0.01)
    else:
        ph.run_total_dos(freq_pitch=0.01)

    # Get a phonon DOS Finger print to compare against the previous one
    n_bins = len(ph.get_total_dos_dict()["frequency_points"])
    dos_fp = get_phonon_dos_fingerprint_phononpy(ph, nbins=n_bins)

    # Get the base working directory
    init_workdir = get_base_work_dir(init_workdir)
    analysis_wd = get_base_work_dir(workdir)
    if prev_dos_fp:
        ph_conv = check_phonon_conv(dos_fp, prev_dos_fp, conv_crit)
    else:
        ph_conv = False

    # Check to see if phonons are converged
    if prev_dos_fp is not None and ph_conv:
        Path(f"{analysis_wd}/converged/").mkdir(exist_ok=True, parents=True)
        copyfile(traj, f"{analysis_wd}/converged/trajectory.son")
        if "prev_wd" in kwargs:
            Path(f"{analysis_wd}/converged_mn_1/").mkdir(exist_ok=True, parents=True)
            traj_prev = f"{kwargs['prev_wd']}/{trajectory}"
            copyfile(traj_prev, f"{analysis_wd}/converged_mn_1/trajectory.son")
        update_job = {
            "ph_dict": phonon_to_dict(ph),
            "ph_calculator": calc_dict,
            "ph_primitive": atoms2dict(to_Atoms_db(ph.get_primitive())),
            "ph_time": calc_time / len(ph.get_supercells_with_displacements()),
        }
        return True, update_job

    # Reset dos_fp to include full Energy Range for the material
    if prev_dos_fp:
        ph.set_total_DOS(tetrahedron_method=True)
        dos_fp = get_phonon_dos_fingerprint_phononpy(ph, nbins=201)

    # If Not Converged update phonons

    if "sc_matrix_original" not in kwargs:
        kwargs["sc_matrix_original"] = ph.get_supercell_matrix()

    ind = np.where(np.array(kwargs["sc_matrix_original"]).flatten() != 0)[0][0]
    if kwargs.get("sc_matrix_original", None) is not None:
        n_cur = int(
            round(
                ph.get_supercell_matrix().flatten()[ind]
                / np.array(kwargs["sc_matrix_original"]).flatten()[ind]
            )
        )
        sc_mat = (n_cur + 1) * np.array(kwargs["sc_matrix_original"]).reshape((3, 3))
    else:
        sc_mat = 2.0 * ph.get_supercell_matrix()

    displacement = ph._displacement_dataset["first_atoms"][0]["displacement"]
    disp_mag = np.linalg.norm(displacement)

    ratio = np.linalg.det(sc_mat) / np.linalg.det(ph.get_supercell_matrix())
    ph, _, _ = ph_preprocess(
        to_Atoms_db(ph.get_primitive()), sc_mat, displacement=displacement
    )

    if ph.get_supercell().get_number_of_atoms() > 500:
        time_scaling = 3.0 * ratio ** 3.0
    else:
        time_scaling = 3.0 * ratio

    expected_walltime = calc_time * time_scaling

    ntasks = int(np.ceil(ph.supercell.get_number_of_atoms() * 1.25))

    init_workdir += f"/sc_natoms_{ph.get_supercell().get_number_of_atoms()}"
    analysis_wd += f"/sc_natoms_{ph.get_supercell().get_number_of_atoms()}"

    displacement = ph._displacement_dataset["first_atoms"][0]["displacement"]
    disp_mag = np.linalg.norm(displacement)

    update_job = {
        "sc_matrix_original": kwargs["sc_matrix_original"],
        "supercell_matrix": sc_mat,
        "init_workdir": init_workdir,
        "analysis_wd": analysis_wd,
        "ntasks": ntasks,
        "expected_walltime": expected_walltime,
        # "expected_mem": expected_mem,
        "ph_calculator": calc_dict,
        "ph_primitive": atoms2dict(to_Atoms_db(ph.get_primitive())),
        "ph_supercell": atoms2dict(to_Atoms_db(ph.get_supercell())),
        "prev_dos_fp": dos_fp,
        "prev_wd": workdir,
        "displacement": disp_mag,
    }
    return False, update_job
