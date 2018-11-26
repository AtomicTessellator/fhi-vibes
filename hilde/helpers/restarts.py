""" handling restarts of tasks and workflows """
import subprocess as sp


def restart(settings, verbose=True):
    """ restart a job according to the restart instructions in the settings """

    if "restart" in settings:
        sp.run(settings.restart.command.split())
        if verbose:
            print(f"Task restarted with {settings.restart.command}")
        return True
    else:
        if verbose:
            print("Task not completed, please inspect and rerun.")
        return False
