#!/usr/bin/env python
"""Check the completion of an aims calculation"""

import balsam.launcher.dag as dag
from hilde.balsam.data_encoder import decode
from hilde.balsam.post.aims_calc_process import postprocess_aims

postprocess_aims(decode(dag.current_job.data))
