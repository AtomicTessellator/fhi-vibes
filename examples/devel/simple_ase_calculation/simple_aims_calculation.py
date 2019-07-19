from pebble import concurrent
from concurrent.futures import TimeoutError

from ase.io import read
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator
from hilde.helpers.paths import cwd
from hilde.tasks.calculate import calculate
from hilde.helpers import Timer

atoms = read("si.in", 0, "aims")

command = "orterun -n 4 aims.x"
species_dir = "/home/knoop/FHIaims/aimsfiles/species_defaults/light"

aims_settings = {
    "command": command,
    "species_dir": species_dir,
    "output_level": "MD_light",
    "relativistic": "atomic_zora scalar",
    "xc": "pw-lda",
    "k_grid": 3 * [2],
    "use_pimd_wrapper": ("localhost", 12345),
}

calc = Aims(**aims_settings, label="test")


# # option 1:
# atoms.calc = calc
# with cwd('tmp', mkdir=True):
#     atoms.calc.calculate()
# print(atoms.get_total_energy())
#
# # option 2:
# err = calculate(atoms, calc, 'tmp_hilde')
# if not err:
#     print(atoms.get_total_energy())


def launch_server(atoms, timeout=3):
    iocalc = atoms.calc

    @concurrent.process(timeout=timeout)
    def launch(iocalc):
        cmd = iocalc.calc.command.replace("PREFIX", iocalc.calc.prefix)
        iocalc.calc.write_input(atoms)
        iocalc.launch_server(cmd)

    return launch(iocalc).result()


def get_forces(atoms, timeout=3):
    @concurrent.process(timeout=timeout)
    def forces(atoms):
        return atoms.get_forces()

    return forces(atoms).result()


with SocketIOCalculator(calc, port=12345) as iocalc:
    atoms.calc = iocalc
    launch_server(atoms, timeout=3)

    try:
        f = get_forces(atoms)
    except TimeoutError as error:
        print("pebble")
        f = atoms.get_forces()
        print(f)
