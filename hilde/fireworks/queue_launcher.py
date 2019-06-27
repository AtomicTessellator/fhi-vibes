"""A modified version of queue launcher to allow for a rapidfire over a single Workflow"""
import os
import glob
import time


from datetime import datetime

import numpy as np

from monty.os import cd, makedirs_p

from fireworks import FWorker
from fireworks.fw_config import (
    SUBMIT_SCRIPT_NAME,
    ALWAYS_CREATE_NEW_BLOCK,
    QUEUE_UPDATE_INTERVAL,
    QSTAT_FREQUENCY,
    RAPIDFIRE_SLEEP_SECS,
    QUEUE_JOBNAME_MAXLEN,
)
from fireworks.queue.queue_launcher import (
    # launch_rocket_to_queue,
    _njobs_in_dir,
    _get_number_of_jobs_in_queue,
    setup_offline_job,
)
from fireworks.utilities.fw_serializers import load_object
from fireworks.utilities.fw_utilities import (
    get_fw_logger,
    get_slug,
    log_exception,
    create_datestamp_dir,
)

from hilde.fireworks.combined_launcher import get_ordred_fw_ids
from hilde.fireworks._defaults import FW_DEFAULTS
from hilde.fireworks.workflows.firework_generator import get_time

__author__ = "Anubhav Jain, Michael Kocher, Modified by Thomas Purcell"
__copyright__ = "Copyright 2012, The Materials Project, Modified 2.11.2018"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Dec 12, 2012"


def launch_rocket_to_queue(
    launchpad,
    fworker,
    qadapter,
    launcher_dir=FW_DEFAULTS.launch_dir,
    reserve=False,
    strm_lvl="INFO",
    create_launcher_dir=False,
    fill_mode=False,
    fw_id=None,
):
    """ Submit a single job to the queue.

    Parameters
    ----------
    launchpad: LaunchPad
        LaunchPad for the launch
    fworker: FWorker
        FireWorker for the launch
    qadapter: QueueAdapterBase
        Queue Adapter for the resource
    launcher_dir: str
        The directory where to submit the job
    reserve: bool
        Whether to queue in reservation mode
    strm_lvl: str
        level at which to stream log messages
    create_launcher_dir: bool
        Whether to create a subfolder launcher+timestamp, if needed
    fill_mode: bool
        whether to submit jobs even when there is nothing to run (only in non-reservation mode)
    fw_id: int
        specific fw_id to reserve (reservation mode only)

    Raises
    ------
    RuntimeError
        If launch is not successful OR
        If queue script could not be submitted
    ValueError
        If the launch directory does not exist OR
        If in offline mode and not reservation mode OR
        If in Reservation Mode and not using a singleshot RocketLauncher OR
        If in Reservation Mode and Fill Mode is also requested OR
        If asking to launch a particular FireWork and not in Reservation Mode
    """
    fworker = fworker if fworker else FWorker()
    launcher_dir = os.path.abspath(launcher_dir)
    l_logger = get_fw_logger(
        "queue.launcher", l_dir=launchpad.logdir, stream_level=strm_lvl
    )

    l_logger.debug("getting queue adapter")
    qadapter = load_object(
        qadapter.to_dict()
    )  # make a defensive copy, mainly for reservation mode

    fw, launch_id = None, None  # only needed in reservation mode
    if not os.path.exists(launcher_dir):
        raise ValueError(
            "Desired launch directory {} does not exist!".format(launcher_dir)
        )

    if "--offline" in qadapter["rocket_launch"] and not reserve:
        raise ValueError(
            "Must use reservation mode (-r option) of qlaunch "
            "when using offline option of rlaunch!!"
        )

    if reserve and "singleshot" not in qadapter.get("rocket_launch", ""):
        raise ValueError(
            "Reservation mode of queue launcher only works for singleshot Rocket Launcher!"
        )

    if fill_mode and reserve:
        raise ValueError("Fill_mode cannot be used in conjunction with reserve mode!")

    if fw_id and not reserve:
        raise ValueError(
            "qlaunch for specific fireworks may only be used in reservation mode."
        )
    if fill_mode or launchpad.run_exists(fworker):
        launch_id = None
        try:
            if reserve:
                if fw_id:
                    l_logger.debug("finding a FW to reserve...")
                fw = launchpad._get_a_fw_to_run(
                    fworker.query, fw_id=fw_id, checkout=False
                )
                fw_id = fw.fw_id
                if not fw:
                    l_logger.info(
                        "No jobs exist in the LaunchPad for submission to queue!"
                    )
                    return False
                l_logger.info("reserved FW with fw_id: {}".format(fw.fw_id))
                # update qadapter job_name based on FW name
                job_name = get_slug(fw.name)[0:QUEUE_JOBNAME_MAXLEN]
                qadapter.update({"job_name": job_name})
                if "_queueadapter" in fw.spec:
                    l_logger.debug("updating queue params using Firework spec..")
                    if "queues" in qadapter:
                        nodes_needed = []
                        if "expected_mem" in fw.spec["_queueadapter"]:
                            expect_mem = fw.spec["_queueadapter"]["expected_mem"]
                        else:
                            expect_mem = 1e-10
                        for queue in qadapter["queues"]:
                            if "max_mem_per_rank" in queue:
                                accessible_mem = queue["max_mem_per_rank"]
                                if expect_mem > accessible_mem:
                                    nodes_needed.append(
                                        int(np.ceil(expect_mem / accessible_mem))
                                    )
                                else:
                                    nodes_needed.append(1)
                            else:
                                nodes_needed.append(0)
                        if "queue" in fw.spec["_queueadapter"]:
                            qu_ind = -1
                            for ii, queue in enumerate(qadapter["queues"]):
                                if queue["name"] == fw.spec["_queueadapter"]["queue"]:
                                    qu_ind = ii
                            if qu_ind < 0:
                                raise ValueError(
                                    "Queue Name not found for that resources"
                                )
                            if (
                                qadapter["queues"][qu_ind]["max_nodes"]
                                < nodes_needed[qu_ind]
                            ):
                                raise IOError(
                                    "Requested resource does not have enough memory to complete the job"
                                )
                            if "walltime" not in fw.spec["_queueadapter"]:
                                fw.spec["_queueadapter"]["walltime"] = qadapter[
                                    "queues"
                                ][qu_ind]["max_walltime"]
                        elif "walltime" in fw.spec["_queueadapter"]:
                            nodes = 1
                            if "nodes" in fw.spec["_queueadapter"]:
                                nodes = fw.spec["_queueadapter"]["nodes"]
                            sc_wt = nodes * get_time(
                                fw.spec["_queueadapter"]["walltime"]
                            )
                            min_acceptable_qu_ind = -1
                            min_acceptable_qu_wt = 1.0e15
                            min_acceptable_qu = None
                            for ii, queue in enumerate(qadapter["queues"]):
                                if nodes_needed[ii] > queue["max_nodes"]:
                                    continue
                                if (
                                    min_acceptable_qu_wt
                                    > get_time(queue["max_walltime"])
                                    * queue["max_nodes"]
                                    > sc_wt
                                ):
                                    min_acceptable_qu_wt = (
                                        get_time(queue["max_walltime"])
                                        * queue["max_nodes"]
                                    )
                                    min_acceptable_qu = queue["name"]
                                    min_acceptable_qu_ind = ii
                            if min_acceptable_qu is None:
                                raise ValueError("Job can't run on requested resource")
                            fw.spec["_queueadapter"]["queue"] = min_acceptable_qu
                            if sc_wt < get_time(
                                qadapter["queues"][min_acceptable_qu_ind][
                                    "max_walltime"
                                ]
                            ):
                                fw.spec["_queueadapter"]["nodes"] = 1
                            else:
                                fw.spec["_queueadapter"]["nodes"] = int(
                                    np.ceil(
                                        sc_wt
                                        / get_time(
                                            qadapter["queues"][min_acceptable_qu_ind][
                                                "max_walltime"
                                            ]
                                        )
                                    )
                                )
                                fw.spec["_queueadapter"]["walltime"] = qadapter[
                                    "queues"
                                ][min_acceptable_qu_ind]["max_walltime"]
                            if (
                                nodes_needed[min_acceptable_qu_ind]
                                > fw.spec["_queueadapter"]["nodes"]
                            ):
                                fw.spec["_queueadapter"]["nodes"] = nodes_needed[
                                    min_acceptable_qu_ind
                                ]
                        del qadapter["queues"]
                    qadapter.update(fw.spec["_queueadapter"])
                if "walltime" in qadapter:
                    for tt, task in enumerate(fw.spec["_tasks"]):
                        if "kwargs" in task and "walltime" in task["kwargs"]:
                            fw.spec["_tasks"][tt]["kwargs"]["walltime"] = (
                                get_time(qadapter["walltime"]) - 180.0
                            )
                        elif "calculate" in task["func"]:
                            fw.spec["_tasks"][tt]["args"][2]["walltime"] = (
                                get_time(qadapter["walltime"]) - 180.0
                            )
                            fw.spec["_tasks"][tt]["args"][3]["walltime"] = (
                                get_time(qadapter["walltime"]) - 180.0
                            )

                launchpad.update_spec([fw.fw_id], fw.spec)
                fw, launch_id = launchpad.reserve_fw(fworker, launcher_dir, fw_id=fw_id)

                # reservation mode includes --fw_id in rocket launch
                qadapter["rocket_launch"] += " --fw_id {}".format(fw.fw_id)

                # update launcher_dir if _launch_dir is selected in reserved fw
                if "_launch_dir" in fw.spec:
                    fw_launch_dir = os.path.expandvars(fw.spec["_launch_dir"])

                    if not os.path.isabs(fw_launch_dir):
                        fw_launch_dir = os.path.join(launcher_dir, fw_launch_dir)

                    launcher_dir = fw_launch_dir

                    makedirs_p(launcher_dir)

                    launchpad.change_launch_dir(launch_id, launcher_dir)
                elif create_launcher_dir:
                    # create launcher_dir
                    launcher_dir = create_datestamp_dir(
                        launcher_dir, l_logger, prefix="launcher_"
                    )
                    launchpad.change_launch_dir(launch_id, launcher_dir)

            elif create_launcher_dir:
                # create launcher_dir
                launcher_dir = create_datestamp_dir(
                    launcher_dir, l_logger, prefix="launcher_"
                )

            # move to the launch directory
            l_logger.info("moving to launch_dir {}".format(launcher_dir))

            with cd(launcher_dir):

                if "--offline" in qadapter["rocket_launch"]:
                    setup_offline_job(launchpad, fw, launch_id)

                l_logger.debug("writing queue script")
                with open(SUBMIT_SCRIPT_NAME, "w") as f:
                    queue_script = qadapter.get_script_str(launcher_dir)
                    f.write(queue_script)

                l_logger.info("submitting queue script")
                reservation_id = qadapter.submit_to_queue(SUBMIT_SCRIPT_NAME)
                if not reservation_id:
                    raise RuntimeError(
                        "queue script could not be submitted, check queue "
                        "script/queue adapter/queue server status!"
                    )
                if reserve:
                    launchpad.set_reservation_id(launch_id, reservation_id)
            return reservation_id
        except:
            log_exception(l_logger, "Error writing/submitting queue script!")
            if reserve and launch_id is not None:
                try:
                    l_logger.info(
                        "Un-reserving FW with fw_id, launch_id: {}, {}".format(
                            fw.fw_id, launch_id
                        )
                    )
                    launchpad.cancel_reservation(launch_id)
                    launchpad.forget_offline(launch_id)
                except:
                    log_exception(
                        l_logger, "Error unreserving FW with fw_id {}".format(fw.fw_id)
                    )

            return False

    else:
        l_logger.info("No jobs exist in the LaunchPad for submission to queue!")
        return (
            None
        )  # note: this is a hack (rather than False) to indicate a soft failure to rapidfire()


def rapidfire(
    launchpad,
    fworker,
    qadapter,
    launch_dir=FW_DEFAULTS.launch_dir,
    nlaunches=FW_DEFAULTS.nlaunches,
    njobs_queue=FW_DEFAULTS.njobs_queue,
    njobs_block=FW_DEFAULTS.njobs_block,
    sleep_time=FW_DEFAULTS.sleep_time,
    reserve=False,
    strm_lvl="CRITICAL",
    timeout=None,
    fill_mode=False,
    fw_ids=None,
    wflow_id=None,
):
    """ Submit many jobs to the queue.

    Parameters
    ----------
    launchpad: LaunchPad
        LaunchPad for the launch
    fworker: FWorker
        FireWorker for the launch
    qadapter: QueueAdapterBase
        Queue Adapter for the resource
    launch_dir: str
        directory where we want to write the blocks
    nlaunches: int
        total number of launches desired; "infinite" for loop, 0 for one round
    njobs_queue: int
        stops submitting jobs when njobs_queue jobs are in the queue, 0 for no limit
    njobs_block: int
        automatically write a new block when njobs_block jobs are in a single block
    sleep_time: int
        secs to sleep between rapidfire loop iterations
    reserve: bool
        Whether to queue in reservation mode
    strm_lvl: str
        level at which to stream log messages
    timeout: int
        # of seconds after which to stop the rapidfire process
    fill_mode: bool
        whether to submit jobs even when there is nothing to run (only in non-reservation mode)
    fw_ids: list of ints
        a list fw_ids to launch (len(fw_ids) == nlaunches)
    wflow_id: list of ints
        a list fw_ids that are a root of the workflow

    Raises
    ------
    ValueError
        If the luanch directory does not exist
    """
    if fw_ids and len(fw_ids) != nlaunches:
        print("WARNING: Setting nlaunches to the length of fw_ids.")
        nlaunches = len(fw_ids)
    sleep_time = sleep_time if sleep_time else RAPIDFIRE_SLEEP_SECS
    launch_dir = os.path.abspath(launch_dir)
    nlaunches = -1 if nlaunches == "infinite" else int(nlaunches)
    l_logger = get_fw_logger(
        "queue.launcher", l_dir=launchpad.logdir, stream_level=strm_lvl
    )

    # make sure launch_dir exists:
    if not os.path.exists(launch_dir):
        raise ValueError(
            "Desired launch directory {} does not exist!".format(launch_dir)
        )

    num_launched = 0
    start_time = datetime.now()
    try:
        l_logger.info("getting queue adapter")

        prev_blocks = sorted(
            glob.glob(os.path.join(launch_dir, "block_*")), reverse=True
        )
        if prev_blocks and not ALWAYS_CREATE_NEW_BLOCK:
            block_dir = os.path.abspath(os.path.join(launch_dir, prev_blocks[0]))
            l_logger.info("Found previous block, using {}".format(block_dir))
        else:
            block_dir = create_datestamp_dir(launch_dir, l_logger)
        while True:
            # get number of jobs in queue
            jobs_in_queue = _get_number_of_jobs_in_queue(
                qadapter, njobs_queue, l_logger
            )
            job_counter = 0  # this is for QSTAT_FREQUENCY option
            if wflow_id:
                wflow = launchpad.get_wf_by_fw_id(wflow_id[0])
                nlaunches = len(wflow.fws)
                fw_ids = get_ordred_fw_ids(wflow)

            while launchpad.run_exists(fworker, ids=fw_ids) or (
                fill_mode and not reserve
            ):
                if timeout and (datetime.now() - start_time).total_seconds() >= timeout:
                    l_logger.info("Timeout reached.")
                    break

                if njobs_queue and jobs_in_queue >= njobs_queue:
                    l_logger.info(
                        "Jobs in queue ({}) meets/exceeds "
                        "maximum allowed ({})".format(jobs_in_queue, njobs_queue)
                    )
                    break

                l_logger.info("Launching a rocket!")

                # switch to new block dir if it got too big
                if _njobs_in_dir(block_dir) >= njobs_block:
                    l_logger.info("Block got bigger than {} jobs.".format(njobs_block))
                    block_dir = create_datestamp_dir(launch_dir, l_logger)
                return_code = None
                # launch a single job
                if fw_ids or wflow_id:
                    return_code = launch_rocket_to_queue(
                        launchpad,
                        fworker,
                        qadapter,
                        block_dir,
                        reserve,
                        strm_lvl,
                        True,
                        fill_mode,
                        fw_ids[num_launched],
                    )
                else:
                    return_code = launch_rocket_to_queue(
                        launchpad,
                        fworker,
                        qadapter,
                        block_dir,
                        reserve,
                        strm_lvl,
                        True,
                        fill_mode,
                    )
                if return_code is None:
                    l_logger.info("No READY jobs detected...")
                    break
                elif not return_code:
                    raise RuntimeError("Launch unsuccessful!")

                if wflow_id:
                    wflow = launchpad.get_wf_by_fw_id(wflow_id[0])
                    nlaunches = len(wflow.fws)
                    fw_ids = get_ordred_fw_ids(wflow)
                num_launched += 1
                if nlaunches > 0 and num_launched == nlaunches:
                    l_logger.info(
                        "Launched allowed number of " "jobs: {}".format(num_launched)
                    )
                    break
                # wait for the queue system to update
                l_logger.info(
                    "Sleeping for {} seconds...zzz...".format(QUEUE_UPDATE_INTERVAL)
                )
                time.sleep(QUEUE_UPDATE_INTERVAL)
                jobs_in_queue += 1
                job_counter += 1
                if job_counter % QSTAT_FREQUENCY == 0:
                    job_counter = 0
                    jobs_in_queue = _get_number_of_jobs_in_queue(
                        qadapter, njobs_queue, l_logger
                    )

            if (
                (nlaunches > 0 and num_launched == nlaunches)
                or (
                    timeout and (datetime.now() - start_time).total_seconds() >= timeout
                )
                or (
                    nlaunches == 0
                    and not launchpad.future_run_exists(fworker, ids=fw_ids)
                )
            ):
                break

            l_logger.info(
                "Finished a round of launches, sleeping for {} secs".format(sleep_time)
            )
            time.sleep(sleep_time)
            l_logger.info("Checking for Rockets to run...")
    except:
        log_exception(l_logger, "Error with queue launcher rapid fire!")
