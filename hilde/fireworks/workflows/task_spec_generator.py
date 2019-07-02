"""Generates TaskSpec Objects"""
from hilde.fireworks.tasks.task_spec import TaskSpec


def gen_phonon_task_spec(func_kwargs, fw_settings=None):
    """Generate a parallel Phononpy or Phono3py calculation task

    Parameters
    ----------
    func_kwargs: dict
        The defined kwargs for func
    fw_settings: dict
        Settings used by fireworks to place objects in the right part of the MongoDB

    Returns
    -------
    TaskSpec
        The specification object of the task
    """
    if fw_settings is not None:
        fw_settings = fw_settings.copy()
    kwargs_init = {}
    kwargs_init_fw_out = {}
    preprocess_keys = {
        "ph_settings": ["supercell_matrix", "displacement"],
        "ph3_settings": ["supercell_matrix", "displacement", "cutoff_pair_distance"],
    }
    out_keys = ["walltime", "trajectory", "backup_folder", "serial"]
    for set_key in ["ph_settings", "ph3_settings"]:
        if set_key in func_kwargs:
            kwargs_init[set_key] = {}
            if "workdir" in func_kwargs[set_key]:
                wd = func_kwargs[set_key]["workdir"]
            else:
                wd = "."
            kwargs_init_fw_out[set_key] = {"workdir": wd}
            for key, val in func_kwargs[set_key].items():
                if key in preprocess_keys[set_key]:
                    kwargs_init[set_key][key] = val
                if key in out_keys:
                    kwargs_init_fw_out[set_key][key] = val
    if fw_settings and "kpoint_density_spec" in fw_settings:
        inputs = [fw_settings["kpoint_density_spec"]]
        args = []
        if "kpt_density" in func_kwargs:
            del func_kwargs["kpt_density"]
    elif "kpt_density" in func_kwargs:
        inputs = []
        args = [func_kwargs.pop("kpt_density")]
    else:
        inputs = []
        args = [None]

    return TaskSpec(
        "hilde.fireworks.tasks.phonopy_phono3py_functions.bootstrap_phonon",
        "hilde.fireworks.tasks.fw_out.phonons.post_init_mult_calcs",
        True,
        kwargs_init,
        inputs=inputs,
        args=args,
        func_fw_out_kwargs=kwargs_init_fw_out,
        make_abs_path=False,
    )


def gen_stat_samp_task_spec(func_kwargs, fw_settings=None):
    """Generate a Harmonic Analysis task

    Parameters
    ----------
    func_kwargs: dict
        The defined kwargs for func
    fw_settings: dict
        Settings used by fireworks to place objects in the right part of
                            the MongoDB

    Returns
    -------
    TaskSpec
        The specification object of the task
    """
    preprocess_keys = [
        "supercell_matrix",
        "phonon_file",
        "temperatures",
        "debye_temp_fact",
        "n_samples",
        "deterministic",
        "rng_seed",
    ]
    out_keys = ["walltime", "trajectory", "backup_folder", "serial"]
    kwargs_init = {}
    kwargs_init_fw_out = {}

    kwargs_init["stat_samp_settings"] = {}
    kwargs_init_fw_out["stat_samp_settings"] = {}

    if "workdir" in func_kwargs:
        wd = func_kwargs["workdir"]
    else:
        wd = "."

    kwargs_init_fw_out["stat_samp_settings"] = {"workdir": wd}
    kwargs_init["stat_samp_settings"] = {"workdir": wd}

    for key, val in func_kwargs.items():
        if key in preprocess_keys:
            kwargs_init["stat_samp_settings"][key] = val
        if key in out_keys:
            kwargs_init_fw_out["stat_samp_settings"][key] = val

    if fw_settings and "kpoint_density_spec" in fw_settings:
        inputs = [fw_settings["kpoint_density_spec"]]
        args = []
        if "kpt_density" in func_kwargs:
            del func_kwargs["kpt_density"]
    elif "kpt_density" in func_kwargs:
        inputs = []
        args = [func_kwargs.pop("kpt_density")]
    else:
        inputs = []
        args = [None]

    return TaskSpec(
        "hilde.fireworks.tasks.phonopy_phono3py_functions.bootstrap_stat_sample",
        "hilde.fireworks.tasks.fw_out.phonons.post_init_mult_calcs",
        True,
        kwargs_init,
        args=args,
        inputs=inputs,
        func_fw_out_kwargs=kwargs_init,
        make_abs_path=False,
    )


def gen_phonon_analysis_task_spec(
    func, func_kwargs, metakey, forcekey, timekey, make_abs_path=False
):
    """Generate a serial Phononpy or Phono3py calculation task

    Parameters
    ----------
    func: str
        The function path to the serial calculator
    func_kwargs: dict
        The defined kwargs for func
    metakey: str
        Key to find the phonon calculation's metadata to recreate the trajectory
    forcekey: str
        Key to find the phonon calculation's force data to recreate the trajectory
    timekey: str
        Key to find the time needed for the phonon supercell calculations
    make_abs_path: bool
        If True make the paths of directories absolute

    Returns
    -------
    TaskSpec
        The specification object of the task
    """
    if "workdir" in func_kwargs and "init_wd" not in func_kwargs:
        func_kwargs["init_wd"] = func_kwargs["workdir"]

    if "analysis_workdir" in func_kwargs:
        func_kwargs["workdir"] = func_kwargs["analysis_workdir"]
    elif "workdir" not in func_kwargs:
        func_kwargs["workdir"] = "."

    if "converge_phonons" in func_kwargs and func_kwargs["converge_phonons"]:
        func_out = "hilde.fireworks.tasks.fw_out.phonons.converge_phonons"
    else:
        func_out = "hilde.fireworks.tasks.fw_out.phonons.add_phonon_to_spec"

    if "trajectory" not in func_kwargs:
        func_kwargs["trajectory"] = "trajectory.son"
    task_spec_list = []
    task_spec_list.append(
        TaskSpec(
            "hilde.fireworks.tasks.phonopy_phono3py_functions.collect_to_trajectory",
            "hilde.fireworks.tasks.fw_out.general.fireworks_no_mods_gen_function",
            False,
            args=[func_kwargs["workdir"], func_kwargs["trajectory"]],
            inputs=[forcekey, metakey],
            make_abs_path=make_abs_path,
        )
    )
    task_spec_list.append(
        TaskSpec(
            "hilde.fireworks.tasks.phonopy_phono3py_functions.phonon_postprocess",
            func_out,
            False,
            args=[func],
            inputs=[timekey, "mem_use"],
            func_kwargs=func_kwargs,
            make_abs_path=make_abs_path,
        )
    )
    return task_spec_list


def gen_aims_task_spec(
    func_kwargs, func_fw_out_kwargs, make_abs_path=False, relax=True
):
    """Gets the task spec for an FHI-aims calculations

    Parameters
    ----------
    func_kwargs: dict
        The defined kwargs for func
    func_fw_outkwargs: dict
        The defined kwargs for fw_out
    make_abs_path: bool
        If True make the paths of directories absolute
    relax: bool
        If True it is a relaxation

    Returns
    -------
    TaskSpec
        The task_spec for the calculation
    """
    fw_out = "hilde.fireworks.tasks.fw_out.relax.check_aims_complete"
    if not relax:
        fw_out = "hilde.fireworks.tasks.fw_out.general.fireworks_no_mods"
    return TaskSpec(
        "hilde.fireworks.tasks.calculate_wrapper.wrap_calculate",
        fw_out,
        True,
        func_kwargs,
        func_fw_out_kwargs=func_fw_out_kwargs,
        make_abs_path=make_abs_path,
    )


def gen_kgrid_task_spec(func_kwargs, make_abs_path=False):
    """Gets the task spec for a k-grid optimization

    Parameters
    ----------
    func_kwargs: dict
        The defined kwargs for func
    make_abs_path: bool
        If True make the paths of directories absolute

    Returns
    -------
    TaskSpec
        The TaskSpec for the kgrid optimization
    """
    return TaskSpec(
        "hilde.k_grid.converge_kgrid.converge_kgrid",
        "hilde.fireworks.tasks.fw_out.optimizations.check_kgrid_opt_completion",
        True,
        func_kwargs,
        make_abs_path=make_abs_path,
    )
