import attr


@attr.s
class CliTracker:
    """click context object"""

    verbose = attr.ib(default=1)
