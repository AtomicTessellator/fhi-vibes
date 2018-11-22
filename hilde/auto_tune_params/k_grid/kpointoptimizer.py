from ase.optimize.optimize import Dynamics
from ase.calculators.calculator import kptdensity2monkhorstpack
from copy import copy

from hilde.helpers.k_grid import d2k

class KPointOptimizer(Dynamics):
    def __init__(
        self,
        atoms,
        func=lambda x: x.calc.get_property('energy',x)/len(x),
        loss_func=lambda x: x,
        dfunc_min=1e-6,
        even=True,
        trajectory=None,
        logfile="-",
        loginterval=1,
        kpts_density_init=1.0,
    ):

        Dynamics.__init__(
            self, atoms, logfile=logfile, trajectory=trajectory, append_trajectory=True
        )

        self.kpts_density = kpts_density_init
        self.even = even
        kpts_initial = d2k(atoms, self.kpts_density, self.even)
        self.kpts = kpts_initial

        self.func = func
        self.loss_func = loss_func
        self.dfunc = 1e12
        self.ref = 0
        self.last = 0
        self.dfunc_min = dfunc_min
        self.n_atoms = len(atoms)

    @property
    def kpts(self):
        return self.atoms.calc.parameters.k_grid

    @kpts.setter
    def kpts(self, kp):
        self.atoms.calc.parameters.k_grid = kp

    def increase_kpts(self):
        # fkdev: wie vernuenftig loesen?!
        self.atoms.calc.results = {}
        while True:
            self.kpts_density += 1
            kpts = d2k(self.atoms, self.kpts_density, self.even)
            if kpts != self.kpts:
                self.kpts = kpts
                break

    def todict(self):
        return {"type": "kpoint-optimizer"}

    def irun(self, steps=100):
        self.ref = self.func(self.atoms)

        for _ in range(steps):
            val = self.func(self.atoms)
            self.dfunc = self.loss_func(self.last - val)
            self.call_observers()
            if self.dfunc < self.dfunc_min:
                yield True
                return True
            else:
                yield False
                self.increase_kpts()
                self.last = val

    def run(self, steps=100):
        for _ in self.irun(steps=steps):
            pass

