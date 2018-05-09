#!/usr/bin/env python3

'''Compute GRMs'''

import os
from configall import PRU_DIR, GRM_DIR, GRM_CUTOFF, GCTA, NBPROC, USE_SBATCH
import preprocessing

def main():
    """Entry point if called as an executable"""
    in_file_pruned = os.path.join(PRU_DIR, 'all')
    out_file_grm = os.path.join(GRM_DIR, 'grm-all', 'all')
    out_file_grm_filtered = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF),
                                         'all-' + str(GRM_CUTOFF))

    # create output directory if it doesn't exist
    if not os.path.exists(GRM_DIR):
        os.makedirs(GRM_DIR)

    ### Using All SNPs
    preprocessing.gcta_grm(in_file=in_file_pruned,
                           par_input='--bfile',
                           out_file=out_file_grm,
                           gcta=GCTA,
                           per_chr=True, #If True, splits by chromosome merges them
                           ncpus=NBPROC,
                           sbatch=USE_SBATCH)

    ### filtering related individuals
    preprocessing.gcta_grm_filter(in_file=out_file_grm,
                                  out_file=out_file_grm_filtered,
                                  grm_cutoff=GRM_CUTOFF,
                                  gcta=GCTA,
                                  per_chr=True, #If True, splits by chromosome merges them
                                  ncpus=NBPROC,
                                  srun=USE_SBATCH)


if __name__ == '__main__':
    main()
