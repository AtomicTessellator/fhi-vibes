'''Used to combine a queue and rocket launcher based on which FireTasks are in a WorkFlow'''
import os
import glob
import time
from datetime import datetime

from monty.os import cd, makedirs_p

from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import launch_rocket
from fireworks.fw_config import SUBMIT_SCRIPT_NAME, ALWAYS_CREATE_NEW_BLOCK, QUEUE_RETRY_ATTEMPTS, \
    QUEUE_UPDATE_INTERVAL, QSTAT_FREQUENCY, RAPIDFIRE_SLEEP_SECS, QUEUE_JOBNAME_MAXLEN
from fireworks.queue.queue_launcher import launch_rocket_to_queue, _njobs_in_dir, \
    _get_number_of_jobs_in_queue, setup_offline_job
from fireworks.utilities.fw_serializers import load_object
from fireworks.utilities.fw_utilities import (get_fw_logger, log_exception, create_datestamp_dir,
                                              get_slug)

from hilde.fireworks_api_adapter.qlaunch_remote import qlaunch_remote
from hilde.tasks import fireworks as fw

def get_ordred_fw_ids(wflow):
    '''Gets an ordered (with respect to when jobs need to run) list of fws in a WorkFlow wflow'''
    fw_ids_ordered = wflow.leaf_fw_ids
    parent_links = wflow.links.parent_links
    for fw_id in fw_ids_ordered:
        parents = parent_links[fw_id] if fw_id in parent_links else []
        for parent in sorted(parents)[::-1]:
            if parent not in fw_ids_ordered:
                fw_ids_ordered.append(parent)
    return fw_ids_ordered[::-1]

def use_queue_launch(fire_work, tasks2queue):
    '''Determines if a particular FireWork should be ran on a cluster'''
    for task in fire_work.spec["_tasks"]:
        if task["func"] in tasks2queue:
            return True
    return False

def rapidfire(launchpad,
              fworker=None,
              qadapter=None,
              launch_dir='.',
              nlaunches=0,
              njobs_queue=0,
              njobs_block=500,
              sleep_time=None,
              reserve=False,
              strm_lvl='CRITICAL',
              timeout=None,
              fill_mode=False,
              fw_ids=None,
              wflow=None,
              tasks2queue=None,
              gss_auth=True,
              remote_host='localhost',
              remote_config_dir=None,
              remote_user=None,
              remote_password=None,
              remote_shell='/bin/bash -l -c', daemon=0):
    """
    Submit many jobs to the queue.

    Args:
        launchpad (LaunchPad)
        fworker (FWorker)
        qadapter (QueueAdapterBase)
        wflow (WorkFlow): the workflow this qlauncher is supposed to run
        launch_dir (str): directory where we want to write the blocks
        nlaunches (int): total number of launches desired; "infinite" for loop, 0 for one round
        njobs_queue (int): stops submitting jobs when njobs_queue jobs are in the queue, 0 for
                           no limit
        njobs_block (int): automatically write a new block when njobs_block jobs are in a
                           single block
        sleep_time (int): secs to sleep between rapidfire loop iterations
        reserve (bool): Whether to queue in reservation mode
        strm_lvl (str): level at which to stream log messages
        timeout (int): # of seconds after which to stop the rapidfire process
        fill_mode (bool): whether to submit jobs even when there is nothing to run (only in
            non-reservation mode)
        fw_ids(list of ints): a list fw_ids to launch (len(fw_ids) == nlaunches)
        wflow_id(list of ints): a list fw_ids that are a root of the workflow
    """
    if tasks2queue is None:
        tasks2queue = [fw.calculate.name]
    if remote_config_dir is None:
        remote_config_dir = ['~/.fireworks']
    r_args = [launchpad]
    r_kwargs = {"fworker": fworker, "strm_lvl": strm_lvl, "pdb_on_exception": False}
    if wflow:
        wflow_id = wflow.root_fw_ids
    if remote_host is 'localhost':
        qlaunch = launch_rocket_to_queue
        remote = False
        if not fworker or not qadapter:
            raise AttributeError("For a direct launch_rocket_to_queue fworker and qadapter " +
                                 "need to be specified")
        q_args = [launchpad, fworker, qadapter]
        q_kwargs = {"launcher_dir": launch_dir, "reserve": reserve, "strm_lvl": strm_lvl,
                    "create_launcher_dir": True, "fill_mode": fill_mode}
    else:
        qlaunch = qlaunch_remote
        remote = True
        q_args = ["rapidfire"]
        q_kwargs = {"reserve": reserve,
                    "gss_auth": gss_auth,
                    "remote_host": remote_host,
                    "remote_config_dir": remote_config_dir,
                    "remote_user": remote_user,
                    "remote_password": remote_password,
                    "remote_shell": remote_shell,
                    "daemon": daemon,
                    "maxjobs_queue": njobs_queue,
                    "nlaunches": 1,
                    "sleep": sleep_time}
        if launch_dir is not ".":
            q_kwargs["launcher_dir"] = launch_dir
        if strm_lvl is not "INFO":
            q_kwargs["loglvl"] = strm_lvl
        if njobs_block is not 500:
            q_kwargs["maxjobs_block"] = njobs_block
    if fw_ids and len(fw_ids) != nlaunches:
        print("WARNING: Setting nlaunches to the length of fw_ids.")
        nlaunches = len(fw_ids)
    sleep_time = sleep_time if sleep_time else RAPIDFIRE_SLEEP_SECS
    launch_dir = os.path.abspath(launch_dir)
    nlaunches = -1 if nlaunches == 'infinite' else int(nlaunches)
    l_logger = get_fw_logger('queue.launcher', l_dir=launchpad.logdir, stream_level=strm_lvl)

    # make sure launch_dir exists:
    if not os.path.exists(launch_dir):
        raise ValueError('Desired launch directory {} does not exist!'.format(launch_dir))

    num_launched = 0
    start_time = datetime.now()
    skip_check = False
    try:
        if remote_host is 'localhost':
            l_logger.info('getting queue adapter')

            prev_blocks = sorted(glob.glob(os.path.join(launch_dir, 'block_*')), reverse=True)
            if prev_blocks and not ALWAYS_CREATE_NEW_BLOCK:
                block_dir = os.path.abspath(os.path.join(launch_dir, prev_blocks[0]))
                l_logger.info('Found previous block, using {}'.format(block_dir))
            else:
                block_dir = create_datestamp_dir(launch_dir, l_logger)
            q_kwargs['launcher_dir'] = block_dir
        while True:
            if remote_host is 'localhost':
                # get number of jobs in queue
                jobs_in_queue = _get_number_of_jobs_in_queue(qadapter, njobs_queue, l_logger)
                job_counter = 0  # this is for QSTAT_FREQUENCY option
            if wflow_id:
                wflow = launchpad.get_wf_by_fw_id(wflow_id[0])
                nlaunches = len(wflow.fws)
                fw_ids = get_ordred_fw_ids(wflow)
            while (skip_check or launchpad.run_exists(fworker, ids=fw_ids) or
                   (fill_mode and not reserve)):

                if timeout and (datetime.now() - start_time).total_seconds() >= timeout:
                    l_logger.info("Timeout reached.")
                    break
                if remote_host is "localhost":
                    if njobs_queue and jobs_in_queue >= njobs_queue:
                        l_logger.info("Jobs in queue ({}) meets/exceeds "
                                      "maximum allowed ({})".format(jobs_in_queue,
                                                                    njobs_queue))
                        break
                    l_logger.info('Launching a rocket!')

                    # switch to new block dir if it got too big
                    if _njobs_in_dir(block_dir) >= njobs_block:
                        l_logger.info('Block got bigger than {} jobs.'.format(njobs_block))
                        block_dir = create_datestamp_dir(launch_dir, l_logger)
                return_code = None
                # launch a single job
                if fw_ids or wflow_id:
                    fw_id = fw_ids[num_launched]
                else:
                    fw_id = launchpad._get_a_fw_to_run(fworker.query, fw_id=None).fw_id
                use_queue = use_queue_launch(launchpad.get_fw_by_id(fw_id), tasks2queue)
                if use_queue:
                    rlaunch = qlaunch
                    args = q_args
                    kwargs = q_kwargs
                    if remote_host is "localhost":
                        kwargs['fw_id'] = fw_id
                    else:
                        kwargs['fw_ids'] = [fw_id]
                else:
                    rlaunch = launch_rocket
                    args = r_args
                    kwargs = r_kwargs
                    kwargs['fw_id'] = fw_id
                print(fw_id, fw_ids)
                return_code = rlaunch(*args, **kwargs)
                if wflow_id:
                    wflow = launchpad.get_wf_by_fw_id(wflow_id[0])
                    nlaunches = len(wflow.fws)
                    fw_ids = get_ordred_fw_ids(wflow)
                if use_queue:
                    num_launched += 1
                elif return_code is None:
                    l_logger.info('No READY jobs detected...')
                    break
                elif not return_code:
                    raise RuntimeError("Launch unsuccessful!")
                else:
                    num_launched += 1
                if nlaunches > 0 and num_launched == nlaunches:
                    l_logger.info('Launched allowed number of '
                                  'jobs: {}'.format(num_launched))
                    break
                # wait for the queue system to update
                if remote_host is "localhost":
                    l_logger.info('Sleeping for {} seconds...zzz...'.format(QUEUE_UPDATE_INTERVAL))
                    time.sleep(QUEUE_UPDATE_INTERVAL)
                if not use_queue and launchpad.run_exists(fworker, ids=fw_ids):
                    skip_check = True
                else:
                    skip_check = False
            if (nlaunches > 0 and num_launched == nlaunches) or \
                    (timeout and (datetime.now() - start_time).total_seconds()
                     >= timeout) or (nlaunches == 0 and not launchpad.future_run_exists(fworker)):
                break

            l_logger.info('Finished a round of launches, sleeping for {} secs'.format(sleep_time))
            time.sleep(sleep_time)
            l_logger.info('Checking for Rockets to run...')
    except:
        log_exception(l_logger, 'Error with queue launcher rapid fire!')
