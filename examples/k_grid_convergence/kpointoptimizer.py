from ase.optimize.optimize import Dynamics
from ase.calculators.calculator import kptdensity2monkhorstpack
from copy import copy

class KPointOptimizer(Dynamics):
    def __init__(self, atoms, func, loss_func, dfunc_min=1e-4, even=True, trajectory=None,
                 logfile='-', loginterval=1):

        Dynamics.__init__(self, atoms, logfile=logfile, trajectory=trajectory,
                         append_trajectory=True)

        self.kpts_density = 1
        self.even = even
        kpts_initial = self.d2k()
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

    def d2k(self):
        """ Wrapper for kptdensity2monkhorstpack to return list """
        return list(kptdensity2monkhorstpack(self.atoms,
                                             kptdensity=self.kpts_density,
                                             even=self.even))

    def todict(self):
        return {'type': 'kpoint-optimizer'}

    def irun(self, nmax=100):
        self.ref = self.func(self.atoms)

        for _ in range(nmax):
            val = self.func(self.atoms)
            self.dfunc = self.loss_func(self.last - val)
            self.call_observers()
            yield
            if self.dfunc < self.dfunc_min:
                print('converged')
                return
            else:
                self.increase_kpts()
                self.last = val

    def run(self, nmax=100):
        for _ in self.irun(nmax=nmax):
            pass

    def increase_kpts(self):
        # fkdev: wie vernuenftig loesen?!
        self.atoms.calc.results = {}
        while True:
            self.kpts_density += 1
            kpts = self.d2k()
            if kpts != self.kpts:
                self.kpts = kpts
                break
