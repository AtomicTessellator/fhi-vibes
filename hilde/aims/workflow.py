""" Provide a full highlevel phonopy workflow

    Input: geometry.in and settings.in
    Output: geometry.in.supercell and trajectory.son """

from ase.io import read

from hilde.tasks import calculate_socket
from hilde.helpers import talk
from hilde.helpers.converters import input2dict, atoms2dict
from hilde.helpers.restarts import restart

from .setup import setup_aims


def run_aims(ctx):
    """ high level function to run aims calculation

    Parameters:
       ctx (AimsContext): aims context

    """

    args = bootstrap(ctx)

    completed = calculate_socket(**args)

    if not completed:
        restart()
    else:
        talk("done.")


def bootstrap(ctx):
    """ load settings, prepare atoms and aims calculator

    Parameters:
        ctx (AimsContext): context object

    """

    # find geometries
    atoms_to_calculate = ctx.atoms_to_calculate

    if len(atoms_to_calculate) < 1:
        raise RuntimeError("no structures to compute.")

    calc = setup_aims(ctx)

    # save metadata
    metadata = input2dict(
        ctx.ref_atoms,
        calc=calc,
        settings=ctx.settings,
        primitive=ctx.primitive,
        supercell=ctx.supercell,
    )

    # save input files
    input_files = {}
    for file, atoms in zip(ctx.geometry_files, atoms_to_calculate):
        dct = {f"{file}": atoms2dict(atoms)}
        input_files.update(dct)

    metadata.update({"geometry_files": input_files})

    return {
        "atoms_to_calculate": atoms_to_calculate,
        "calculator": calc,
        "metadata": metadata,
        "workdir": ctx.workdir,
        "settings": ctx.settings,
        "backup_after_calculation": False,
    }
