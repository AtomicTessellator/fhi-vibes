"""Default definitions for FireWorks"""
from hilde.settings import Settings
from hilde.helpers.attribute_dict import AttributeDict as adict

SETTINGS = Settings()
REMOTE_SETUP = SETTINGS.remote_setup if "remote_setup" in SETTINGS else {}
REMOTE_HOST_AUTH = SETTINGS.remote_host_auth if "remote_host_auth" in SETTINGS else {}
REMOTE_QUEUE_PARAM = (
    SETTINGS.remote_queue_param if "remote_queue_param" in SETTINGS else {}
)
LAUNCH_PARAMS = SETTINGS.launch_params if "launch_params" in SETTINGS else {}

FW_DEFAULTS = adict(
    {
        "launch_dir": (
            REMOTE_SETUP.launch_dir if "launch_dir" in REMOTE_SETUP else "."
        ),
        "remote_host": (
            REMOTE_SETUP.remote_host if "remote_host" in REMOTE_SETUP else None
        ),
        "remote_config_dir": (
            REMOTE_SETUP.remote_config_dir
            if "remote_config_dir" in REMOTE_SETUP
            else "~/.fireworks"
        ),
        "remote_user": (
            REMOTE_HOST_AUTH.remote_user if "remote_user" in REMOTE_HOST_AUTH else None
        ),
        "remote_password": (
            REMOTE_HOST_AUTH.remote_password
            if "remote_password" in REMOTE_HOST_AUTH
            else None
        ),
        "njobs_queue": (
            REMOTE_QUEUE_PARAM.njobs_queue if "njobs_queue" in REMOTE_QUEUE_PARAM else 0
        ),
        "njobs_block": (
            REMOTE_QUEUE_PARAM.njobs_block
            if "njobs_block" in REMOTE_QUEUE_PARAM
            else 500
        ),
        "nlaunches": (LAUNCH_PARAMS.nlaunches if "nlaunches" in LAUNCH_PARAMS else 0),
        "sleep_time": (
            LAUNCH_PARAMS.sleep_time if "sleep_time" in LAUNCH_PARAMS else None
        ),
        "tasks2queue": (
            LAUNCH_PARAMS.tasks2queue if "tasks2queue" in LAUNCH_PARAMS else ""
        ),
    }
)
