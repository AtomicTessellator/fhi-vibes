"""hiphive utilities"""

import numpy as np
from ase import Atoms

from vibes.helpers import Timer as _Timer
from vibes.helpers import talk as _talk
from vibes.helpers import warn
from vibes.helpers.force_constants import ForceConstants as MyForceConstants
from vibes.helpers.geometry import inscribed_sphere_in_box


try:
    from hiphive import ClusterSpace, ForceConstantPotential, ForceConstants
    from hiphive import enforce_rotational_sum_rules as _enforce_sum_rules
    from hiphive.utilities import extract_parameters
except ModuleNotFoundError:
    warn(f"`hiphive` could not be loaded, enforcing sum rules not possible", level=1)


_prefix = "hiphive"


def talk(msg):
    return _talk(msg, prefix=_prefix)


def Timer(msg=None):
    return _Timer(message=msg, prefix=_prefix)


def enforce_rotational_sum_rules(
    force_constants: np.ndarray, primitive: Atoms, supercell: Atoms, decimals: int = 9,
):
    """enforce the rotational sum rules on phonopy object

    Args:
        force_constants: force constants in [Np, Ns, 3, 3] shape
        primitive: reference primitive structure
        supercell: reference supercell structure
        decimals: use to clean up positions to remove spurious 1s (hiphive bug)

    References:
        - hiphive.materialsmodeling.org/advanced_topics/rotational_sum_rules.html
        - hiphive.materialsmodeling.org/advanced_topics/force_constants_io.html
        - hiphive.materialsmodeling.org/moduleref/cluster_space.html
    """
    # clean positions
    primitive.positions += 10 ** (-decimals)
    primitive.wrap()
    primitive.positions -= 10 ** (-decimals)
    primitive.positions = primitive.positions.round(decimals=decimals)

    kw = {"primitive": primitive, "supercell": supercell}
    fc_obj = MyForceConstants(force_constants, **kw)

    cutoff = inscribed_sphere_in_box(supercell.cell) * 0.999

    # create ForceConstants
    fcs = ForceConstants.from_arrays(supercell, fc2_array=fc_obj.remapped_to_supercell)

    # create ClusterSpace
    cs = ClusterSpace(primitive, [cutoff])

    # extract parameters
    # kw = {"fit_method": "lasso", "alpha": .005}
    # talk(f"Use fit method `{kw}`")
    parameters = extract_parameters(fcs, cs)  # , **kw)

    # Apply sum rules
    sum_rules = ["Huang", "Born-Huang"]
    msg = f"Enforce {sum_rules} via `hiphive.enforce_rotational_sum_rules`"
    timer = Timer(msg)
    parameters_rot = _enforce_sum_rules(cs, parameters, sum_rules)
    timer()

    # create new ForceConstantPotential
    fcp_rot = ForceConstantPotential(cs, parameters_rot)

    # obtain new force constants
    fc_rot = fcp_rot.get_force_constants(supercell).get_fc_array(order=2)

    fc_rot_obj = MyForceConstants(fc_rot, **kw)

    # compare to old force constants
    dev = np.linalg.norm(fc_obj.array - fc_rot_obj.array)
    msg = ["Deviation to phonopy force_constants after enforcing sum rules:", f"{dev}"]
    talk(msg)

    # check Symmetry
    dev = np.linalg.norm(fc_rot_obj.remapped - fc_rot_obj.remapped.T)
    if dev > 1e-9:
        warn(f"Force constants not symmetric after enforcing sum rules by: {dev}")

    return fc_rot_obj.array
