import os
import sys
import time

try:
    import fabric
    if int(fabric.__version__.split('.')[0]) < 2:
        raise ImportError
except ImportError:
    HAS_FABRIC = False
else:
    HAS_FABRIC = True

__authors__ = "Anubhav Jain, Shyue Ping Ong, Adapted by Thomas Purcell for use in HiLDe"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jan 14, 2013, Adaptation: October 31, 2018"

def convert_input_to_param(param_name, param, param_list):
    '''
    Converts a function input into a qlaunch parameter
    Args:
        param_name: str
            name of the parameter for qlaunch
        param:
            value of the parameter
    Returns: str
        command argument str for qlaunch
    '''
    if param is not None:
        param_list.append(f"--{param_name} {param}")
    return param_list
def qlaunch_remote(command,
                   maxjobs_queue=None,
                   maxjobs_block=None,
                   nlaunches=None,
                   sleep=None,
                   fw_ids=None,
                   wflow=None,
                   silencer=False,
                   reserve=False,
                   launcher_dir=None,
                   loglvl=None,
                   gss_auth=True,
                   remote_host='localhost',
                   remote_config_dir=['~/.fireworks'],
                   remote_user=None,
                   remote_password=None,
                   remote_shell='/bin/bash -l -c',
                   daemon=0):
    '''
    This function adapts the python definition of qlaunch in fireworks to a python function
    Args:
        command: str ("singleshot" or "rappidfire")
            Whether to do a single shot or rapidfire command
        maxjobs_queue: int
            maximum jobs to keep in queue for this user'
        maxjobs_block: int
            maximum jobs to put in a block
        nlaunches: int
            num_launches (int or "infinite"; default 0 is all jobs in DB)
        sleep: int
            sleep time between loops
        fw_ids: list of int
            specific fw_ids to run in reservation mode
        silencer: bool
            shortcut to mute log messages
        reserve: bool
            reserve a fw
        remote_host: str
            Remote host to exec qlaunch. Right now, only supports running from a config dir.
        remote_config_dir: list of str
            Remote config dir location(s). Defaults to ~/.fireworks. You can specify multiple
            locations if you have multiple configurations on the same cluster e.g., multiple
            queues or FireWorkers.
        remote_user: str
            Username to login to remote host.
        remote_password: str
            Password for remote host (if necessary). For best operation, it is recommended that
            you do passwordless ssh.
        remote_shell: str
            Shell command to use on remote host for running submission.
        daemon: int
            Daemon mode. Command is repeated every x seconds. Defaults to 0, which means non-daemon
            mode.
    '''
    assert remote_host is not "localhost"

    if not HAS_FABRIC:
        print("Remote options require the Fabric package v2+ to be installed!")
        sys.exit(-1)

    non_default = []
    convert_input_to_param("launch_dir", launcher_dir, non_default)
    convert_input_to_param("loglvl", loglvl, non_default)
    if command is "rapidfire":
        convert_input_to_param("maxjobs_queue", maxjobs_queue, non_default)
        convert_input_to_param("maxjobs_block", maxjobs_block, non_default)
        convert_input_to_param("nlaunches", nlaunches, non_default)
        convert_input_to_param("sleep", sleep, non_default)
        if fw_ids:
            non_default.append("--{} {}".format("firework_ids", fw_ids[0]))
            for fw in fw_ids[1:]:
                non_default[-1] += " {}".format(fw)
        if wflow:
            non_default.append("--{} {}".format("wflow", wflow.root_fw_ids[0]))
            for fw in wflow.root_fw_ids[1:]:
                non_default[-1] += " {}".format(fw)
    else:
        convert_input_to_param("fw_id", fw_id, non_default)
    non_default = " ".join(non_default)

    pre_non_default = []
    if silencer:
        pre_non_default.append("--silencer")
    if reserve:
        pre_non_default.append("--reserve")
    pre_non_default = " ".join(pre_non_default)
    print("qlaunch_hilde {} {} {}".format(pre_non_default, command, non_default))
    interval = daemon
    while True:
        for h in remote_host:
            with fabric.Connection(
                    host=h,
                    user=remote_user,
                    config=fabric.Config({'run': {'shell': remote_shell}}),
                    connect_kwargs={'password': remote_password,
                                    "gss_auth": gss_auth}) as conn:
                for r in remote_config_dir:
                    r = os.path.expanduser(r)
                    with conn.cd(r):
                        conn.run("qlaunch_hilde {} {} {}".format(
                            pre_non_default, command, non_default))
        if interval > 0:
            print("Next run in {} seconds... Press Ctrl-C to exit at any "
                  "time.".format(interval))
            time.sleep(daemon)
        else:
            break
