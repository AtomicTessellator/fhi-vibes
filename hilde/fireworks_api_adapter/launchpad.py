from fireworks.core.launchpad import LaunchPad
class LaunchPad_hilde(LaunchPad):
    '''
    The modified Launchpad that manges the FireWorks database
    '''
    def __init__(self, *args, **kwargs):
        super(LaunchPad_hilde, self).__init__(*args, **kwargs)
    def run_exists(self, fworker=None, ids=None):
        """
        Checks to see if the database contains any FireWorks of a given id set that are
        ready to run.

        Returns:
            bool: True if the database contains any FireWorks that are ready to run in a given set.
        """
        q = fworker.query if fworker else {}
        if ids:
            q["fw_id"] = {"$in": ids }
        return bool(self._get_a_fw_to_run(query=q, checkout=False))