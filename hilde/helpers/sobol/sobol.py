""" Module to compute Sobol sequences as adapted from `SALib <github.com/SALib/SALib>`_
#==============================================================================
# The following code is based on the Sobol sequence generator by Frances
# Y. Kuo and Stephen Joe. The license terms are provided below.
#
# Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the names of the copyright holders nor the names of the
# University of New South Wales and the University of Waikato
# and its contributors may be used to endorse or promote products derived
# from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#==============================================================================
"""

import sys

import numpy as np
from .directions import directions
from hilde.helpers.utils import talk

if sys.version_info[0] > 2:
    long = int

n_samples = 500
lowest_startingpoint = 100
highest_startingpoint = 1000


def sample(N, D):
    """Generate (N x D) numpy array of Sobol sequence samples

    Parameters
    ----------
    N: int
        number of rows in the array
    D: int
        number of columns in the array

    Returns
    -------
    result: np.ndarray
        Sobol sequence sample array

    Raises
    -----
    ValueError
        If Not enough dimensions in directions OR
        If there are not enough bits for the sequence
    """
    scale = 31
    result = np.zeros([N, D])

    if D > len(directions) + 1:
        raise ValueError("Error in Sobol sequence: not enough dimensions")

    L = int(np.ceil(np.log(N) / np.log(2)))

    if L > scale:
        raise ValueError("Error in Sobol sequence: not enough bits")

    for i in range(D):
        V = np.zeros(L + 1, dtype=np.int64)

        if i == 0:
            for j in range(1, L + 1):
                V[j] = 1 << (scale - j)  # all m's = 1
        else:
            m = np.array(directions[i - 1], dtype=int)
            a = m[0]
            s = len(m) - 1

            # The following code discards the first row of the ``m`` array
            # Because it has floating point errors, e.g. values of 2.24e-314
            if L <= s:
                for j in range(1, L + 1):
                    V[j] = m[j] << (scale - j)
            else:
                for j in range(1, s + 1):
                    V[j] = m[j] << (scale - j)
                for j in range(s + 1, L + 1):
                    V[j] = V[j - s] ^ (V[j - s] >> s)
                    for k in range(1, s):
                        V[j] ^= ((a >> (s - 1 - k)) & 1) * V[j - k]

        X = long(0)
        for j in range(1, N):
            X ^= V[index_of_least_significant_zero_bit(j - 1)]
            result[j][i] = float(X / np.power(2, scale))

    return result


def index_of_least_significant_zero_bit(value):
    """Get the index of the least significant zero bit

    Parameters
    ----------
    value:
        value to check

    Returns
    -------
    index: int
        The index of the least significant zero bit
    """
    index = 1
    while (value & 1) != 0:
        value >>= 1
        index += 1

    return index


class RandomState:
    """ Similar to np.random.RandomState, but for Sobol sequences """

    def __init__(
        self,
        dimension=None,
        low=lowest_startingpoint,
        high=highest_startingpoint,
        startpoint=False,
        randomize=True,
        seed=None,
        failsafe=True,
    ):
        """ Initialize the QuasiRandomState for samples of specific dimension"""

        # print copyright
        self.copyright_notice()

        self.dimension = dimension
        self.low = low
        self.high = high
        self.startpoint = startpoint
        self.randomize = randomize
        self.seed = seed
        self.failsafe = failsafe

        # Choose the starting point of the Sobol sequence. Similar to a seed.
        if not self.seed:
            self.seed = np.random.randint(2 ** 32 - 1)

        talk(f"[sobol]: use random seed of {self.seed}")

        np.random.seed(seed)
        rng = np.random

        if not self.startpoint:
            if randomize:
                self.startpoint = rng.randint(low=low, high=high + low)
            else:
                self.startpoint = low

        talk(f"[sobol]: use startpoint of {self.startpoint}")

    def rand(self, nsamples, dimension=1):
        """create sobol numbers"""
        if dimension > 21200:
            raise ValueError(f"dimension must not exceed 21200 (is {dimension}).")

        startpoint = int(self.startpoint)
        self.startpoint += nsamples

        if self.dimension:
            if self.failsafe and dimension != self.dimension:
                msg = (
                    f"\nRandom numbers initialized with dimension: {self.dimension},"
                    + f" given: {dimension}."
                    + f"\nIf you know what you are doing, run with `failsafe=False`"
                )
                raise ValueError(msg)
            sequence = sample(startpoint + nsamples, self.dimension)
            # flatten, truncate and reshape
            sequence = sequence[startpoint:, :].flatten()[: nsamples * dimension]
            sequence.resize((nsamples, dimension))
        else:
            sequence = sample(startpoint + nsamples, dimension)[startpoint:, :]

        if self.randomize:
            # further randomize by adding 1 and taking modulo
            sequence = (sequence + np.random.rand()) % 1

        return sequence.squeeze()

    @staticmethod
    def copyright_notice():
        """print the copyright statement for SALib"""
        msg = (
            "\n",
            "# The Sobol numbers are created by a code that is based ",
            "# on the Sobol sequence generator by Frances Y. Kuo and Stephen Joe. ",
            "# Copyright (c) 2008, Frances Y. Kuo and Stephen Joe",
            "# All rights reserved.",
            "# More info: https://web.maths.unsw.edu.au/~fkuo/sobol/",
            "\n",
        )

        print("\n".join(msg), flush=True)
