"""Generate FWActions after post-processing statistical sampling calculations"""
from fireworks import FWAction

from vibes.fireworks.tasks.postprocess.statistical_sampling import get_sigma
from vibes.helpers.k_grid import k2d
from vibes.helpers.converters import dict2atoms
from vibes.trajectory.io import reader


def add_stat_samp_to_spec(func, func_fw_out, *args, fw_settings=None, **kwargs):
    """Add the phonon_dict to the spec

    Parameters
    ----------
    func: str
        Path to the phonon analysis function
    func_fw_out: str
        Path to this function
    fw_settings: dict
        Dictionary for the FireWorks specific systems
    kwargs: dict
        Dictionary of keyword arguments that must have the following objects
        workdir: str
            Working directory for the calculation
        trajectory: str
            filename for the trajectory

    Returns
    -------
    FWAction
        FWAction that adds the phonon_dict to the spec
    """
    traj = f"{kwargs['workdir']}/{kwargs['trajectory']}"

    sigma = get_sigma(traj)

    _, metadata = reader(traj, True)
    calc_dict = metadata["calculator"]
    calc_dict["calculator"] = calc_dict["calculator"].lower()
    if calc_dict["calculator"] == "aims":
        k_pt_density = k2d(
            dict2atoms(metadata["supercell"]),
            calc_dict["calculator_parameters"]["k_grid"],
        )
    else:
        k_pt_density = None
    qadapter = {}
    if fw_settings and "spec" in fw_settings:
        qadapter = fw_settings["spec"].get("_queueadapter", None)
    update_spec = {
        "stat_samp_calculator": calc_dict,
        "stat_samp_supercell": metadata["supercell"],
        "_queueadapter": qadapter,
        "sigma": sigma,
    }
    update_spec["kgrid"] = k_pt_density
    return FWAction(update_spec=update_spec)
