import attr

@attr.s
class CliTracker:
    verbose = attr.ib(default=1)
