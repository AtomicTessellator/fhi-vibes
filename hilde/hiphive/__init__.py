"""hiphive utilities"""
import numpy as np

from hiphive import ForceConstants, ForceConstantPotential, ClusterSpace
from hiphive import enforce_rotational_sum_rules as _enfore_sum_rules
from hiphive.utilities import extract_parameters

from hilde.structure.convert import to_Atoms
from hilde.helpers.geometry import inscribed_sphere_in_box
from hilde.helpers import warn, Timer as _Timer, talk as _talk
from hilde.phonopy.utils import remap_force_constants

_prefix = "hiphive"


def talk(msg):
    return _talk(msg, prefix=_prefix)


def Timer(msg=None):
    return _Timer(message=msg, prefix=_prefix)


def enforce_rotational_sum_rules(phonon):
    """enforce the rotational sum rules on phonopy object

    Reference:
        - [enfore rotational sum rules](https://hiphive.materialsmodeling.org/advanced_topics/rotational_sum_rules.html)
        - [ForceConstants from external sources](https://hiphive.materialsmodeling.org/advanced_topics/force_constants_io.html)
        - [create ClusterSpace](https://hiphive.materialsmodeling.org/moduleref/cluster_space.html)
    """
    primitive = to_Atoms(phonon.primitive)
    supercell = to_Atoms(phonon.supercell)
    fc = phonon.get_force_constants()

    Np = len(primitive)
    Na = len(supercell)
    p2s_map = phonon.primitive.get_primitive_to_supercell_map()
    cutoff = inscribed_sphere_in_box(supercell.cell) * 0.99

    reduced_fc = False
    if fc.shape[0] == Np:
        reduced_fc = True
        fc = remap_force_constants(fc, primitive, supercell)

    # create ForceConstants
    fcs = ForceConstants.from_arrays(supercell, fc2_array=fc)

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
    parameters_rot = _enfore_sum_rules(cs, parameters, sum_rules)
    timer()

    # create new ForceConstantPotential
    fcp_rot = ForceConstantPotential(cs, parameters_rot)

    # obtain new force constants
    fc_rot = fcp_rot.get_force_constants(supercell).get_fc_array(order=2)

    # compare to old force constants
    dev = np.linalg.norm(fc - fc_rot)
    msg = ["Deviation to phonopy force_constants after enforcing sum rules:", f"{dev}"]
    talk(msg)

    # check Symmetry
    H = fc_rot.swapaxes(1, 2).reshape(3 * Na, 3 * Na)
    dev = np.linalg.norm(H - H.T)
    if dev > 1e-9:
        warn(f"Force constants not symmetric after enforcing sum rules by: {dev}")

    if reduced_fc:
        fc_rot = fc_rot[p2s_map, :, :, :]

    phonon.set_force_constants(fc_rot)

    return True
