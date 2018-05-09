#!/usr/bin/env python3

"""PCA of GRM computed on filtered data"""

import os
from configall import GRM_DIR, GCTA, GRM_CUTOFF, MYGCTA, NBPROC, USE_SBATCH
import preprocessing

def main():
    """Entry point if called as an executable"""
    in_grm_filtered = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF), 'all-' + str(GRM_CUTOFF))
    out_file_pca = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF), 'all-' + str(GRM_CUTOFF))
    nbpcs = 10

    if USE_SBATCH:
        rmode = "srun"
    else:
        rmode = "direct"

    preprocessing.run([MYGCTA, GCTA,
                       "--grm-bin", in_grm_filtered,
                       "--pca", str(nbpcs),
                       "--out", out_file_pca,
                       "--thread-num", str(NBPROC)],
                      mode=rmode,
                      slurm_par=["-c", str(NBPROC)])

if __name__ == '__main__':
    main()
