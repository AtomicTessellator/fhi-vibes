""" Take a pickled phonopy object and print bandstructure """
from argparse import ArgumentParser

import json
import numpy as np

from hilde.helpers.pickletools import pread
from hilde.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy
from hilde.phonopy.wrapper import plot_bandstructure, get_bandstructure, get_dos
def main():
    """ main routine """
    parser = ArgumentParser(description="plot band structure from a phonopy object")
    parser.add_argument("file", help="pickled phonopy object")
    parser.add_argument("--fp_file", default=None)
    args = parser.parse_args()

    phonon = pread(args.file)
    get_bandstructure(phonon)

    qpts = np.array(phonon.band_structure.qpoints)
    qpts = qpts.reshape(-1, qpts.shape[-1])

    freq = np.array(phonon.band_structure.frequencies)
    freq = freq.reshape(-1, freq.shape[-1])

    gamma_freq = freq[np.where((qpts == np.zeros(3)).all(-1))[0][0]]
    max_freq = np.max(freq.flatten())

    if args.fp_file:
        print(f"Saving the fingerprint to {args.fp_file}")
        fp = get_phonon_bs_fingerprint_phononpy(phonon, binning=False)
        fp_dict = {}
        for freq, pt in zip(fp[0], fp[2]):
            fp_dict[pt] = freq.tolist()
        with open(args.fp_file, 'w') as outfile:
            json.dump(fp_dict, outfile, indent=4)
    print(f"The highest frequency is {max_freq}")
    print(f"The frequencies at the gamma point are:\n{gamma_freq}")

if __name__ == "__main__":
    main()
