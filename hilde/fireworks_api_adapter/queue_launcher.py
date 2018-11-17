"""A modified version of queue launcher to allow for a rapidfire over a single Workflow"""
import os
import glob
import time
from datetime import datetime

from fireworks.fw_config import (
    ALWAYS_CREATE_NEW_BLOCK,
    QUEUE_UPDATE_INTERVAL,
    QSTAT_FREQUENCY,
    RAPIDFIRE_SLEEP_SECS,
)
from fireworks.queue.queue_launcher import (
    launch_rocket_to_queue,
    _njobs_in_dir,
    _get_number_of_jobs_in_queue,
)
from fireworks.utilities.fw_utilities import get_fw_logger, log_exception, create_datestamp_dir

from .combined_launcher import get_ordred_fw_ids

__author__ = "Anubhav Jain, Michael Kocher, Modified by Thomas Purcell"
__copyright__ = "Copyright 2012, The Materials Project, Modified 2.11.2018"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Dec 12, 2012"


def rapidfire(
    launchpad,
    fworker,
    qadapter,
    launch_dir=".",
    nlaunches=0,
    njobs_queue=0,
    njobs_block=500,
    sleep_time=None,
    reserve=False,
    strm_lvl="CRITICAL",
    timeout=None,
    fill_mode=False,
    fw_ids=None,
    wflow_id=None,
):
    """
    Submit many jobs to the queue.

    Args:
        launchpad (LaunchPad)
        fworker (FWorker)
        qadapter (QueueAdapterBase)
        launch_dir (str): directory where we want to write the blocks
        nlaunches (int): total number of launches desired; "infinite" for loop, 0 for one round
        njobs_queue (int): stops submitting jobs when njobs_queue jobs are in the queue, 0 for no
                           limit
        njobs_block (int): automatically write a new block when njobs_block jobs are in a single
                           block
        sleep_time (int): secs to sleep between rapidfire loop iterations
        reserve (bool): Whether to queue in reservation mode
        strm_lvl (str): level at which to stream log messages
        timeout (int): # of seconds after which to stop the rapidfire process
        fill_mode (bool): whether to submit jobs even when there is nothing to run (only in
                          non-reservation mode)
        fw_ids(list of ints): a list fw_ids to launch (len(fw_ids) == nlaunches)
        wflow_id(list of ints): a list fw_ids that are a root of the workflow
    """
    if fw_ids and len(fw_ids) != nlaunches:
        print("WARNING: Setting nlaunches to the length of fw_ids.")
        nlaunches = len(fw_ids)
    sleep_time = sleep_time if sleep_time else RAPIDFIRE_SLEEP_SECS
    launch_dir = os.path.abspath(launch_dir)
    nlaunches = -1 if nlaunches == "infinite" else int(nlaunches)
    l_logger = get_fw_logger("queue.launcher", l_dir=launchpad.logdir, stream_level=strm_lvl)

    # make sure launch_dir exists:
    if not os.path.exists(launch_dir):
        raise ValueError("Desired launch directory {} does not exist!".format(launch_dir))

    num_launched = 0
    start_time = datetime.now()
    try:
        l_logger.info("getting queue adapter")

        prev_blocks = sorted(glob.glob(os.path.join(launch_dir, "block_*")), reverse=True)
        if prev_blocks and not ALWAYS_CREATE_NEW_BLOCK:
            block_dir = os.path.abspath(os.path.join(launch_dir, prev_blocks[0]))
            l_logger.info("Found previous block, using {}".format(block_dir))
        else:
            block_dir = create_datestamp_dir(launch_dir, l_logger)
        while True:
            # get number of jobs in queue
            jobs_in_queue = _get_number_of_jobs_in_queue(qadapter, njobs_queue, l_logger)
            job_counter = 0  # this is for QSTAT_FREQUENCY option
            if wflow_id:
                wflow = launchpad.get_wf_by_fw_id(wflow_id[0])
                nlaunches = len(wflow.fws)
                fw_ids = get_ordred_fw_ids(wflow)
            while launchpad.run_exists(fworker, ids=fw_ids) or (fill_mode and not reserve):
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
                        launchpad, fworker, qadapter, block_dir, reserve, strm_lvl, True, fill_mode
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
                    l_logger.info("Launched allowed number of " "jobs: {}".format(num_launched))
                    break
                # wait for the queue system to update
                l_logger.info("Sleeping for {} seconds...zzz...".format(QUEUE_UPDATE_INTERVAL))
                time.sleep(QUEUE_UPDATE_INTERVAL)
                jobs_in_queue += 1
                job_counter += 1
                if job_counter % QSTAT_FREQUENCY == 0:
                    job_counter = 0
                    jobs_in_queue = _get_number_of_jobs_in_queue(qadapter, njobs_queue, l_logger)

            if (
                (nlaunches > 0 and num_launched == nlaunches)
                or (timeout and (datetime.now() - start_time).total_seconds() >= timeout)
                or (nlaunches == 0 and not launchpad.future_run_exists(fworker, ids=fw_ids))
            ):
                break

            l_logger.info("Finished a round of launches, sleeping for {} secs".format(sleep_time))
            time.sleep(sleep_time)
            l_logger.info("Checking for Rockets to run...")
    except:
        log_exception(l_logger, "Error with queue launcher rapid fire!")
