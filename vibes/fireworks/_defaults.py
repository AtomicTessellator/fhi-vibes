"""Default definitions for FireWorks"""
from vibes.settings import ConfigDict
from vibes.helpers.attribute_dict import AttributeDict as adict
from vibes._defaults import DEFAULT_FIREWORKS_FILE

SETTINGS = ConfigDict(config_files=[DEFAULT_FIREWORKS_FILE])

REMOTE_SETUP = SETTINGS.pop("remote_setup", {})
REMOTE_HOST_AUTH = SETTINGS.pop("remote_host_auth", {})
REMOTE_QUEUE_PARAM = SETTINGS.pop("remote_queue_param", {})
LAUNCH_PARAMS = SETTINGS.pop("launch_params", {})

FW_DEFAULTS = adict(
    {
        "launch_dir": REMOTE_SETUP.pop("launch_dir", "."),
        "remote_host": REMOTE_SETUP.pop("remote_host", None),
        "remote_config_dir": REMOTE_SETUP.pop("remote_config_dir", "~/.fireworks"),
        "remote_user": REMOTE_HOST_AUTH.pop("remote_user", None),
        "remote_password": REMOTE_HOST_AUTH.pop("remote_password", None),
        "njobs_queue": REMOTE_QUEUE_PARAM.pop("njobs_queue", 0),
        "njobs_block": REMOTE_QUEUE_PARAM.pop("njobs_block", 500),
        "nlaunches": LAUNCH_PARAMS.pop("nlaunches", 0),
        "sleep_time": LAUNCH_PARAMS.pop("sleep_time", None),
        "tasks2queue": LAUNCH_PARAMS.pop("tasks2queue", ""),
    }
)
