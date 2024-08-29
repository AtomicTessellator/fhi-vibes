from pathlib import Path
import numpy as np

from vibes.io import read, parse_force_constants
from vibes.dynamical_matrix import DynamicalMatrix

parent = Path(__file__).parent
materials = ["beta_ga2o3"]

for material in materials:
    print(f"Test {material}")

    primitive = read(parent / material / "geometry.in.primitive")
    supercell = read(parent / material / "geometry.in.supercell")
    fc_file = parent / material / "FORCE_CONSTANTS"
    fc = parse_force_constants(fc_file, two_dim=False)

    dmx_full = DynamicalMatrix(
        force_constants=fc,
        primitive=primitive,
        supercell=supercell,
        symmetry=False,
    )

    
    dmx_reduced = DynamicalMatrix(
        force_constants=fc,
        primitive=primitive,
        supercell=supercell,
        symmetry=True,
    )
    
    assert np.allclose(dmx_full.w_sq, dmx_reduced.w_sq)
    assert np.allclose(dmx_full.v_sqa_cartesian, dmx_reduced.v_sqa_cartesian)
    # Here are known mismatching need to check the rotation of eigenvectors
    # assert np.allclose(dmx_full._solution.D_qij, dmx_reduced._solution.D_qij)
    # assert np.allclose(dmx_full.e_isq, dmx_reduced.e_isq)
