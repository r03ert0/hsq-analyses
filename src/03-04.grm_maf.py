#!/usr/bin/env python3

"""GRMs using SNPs in 3 ranges of MAF"""

import os
import preprocessing
from configall import PRU_DIR, GRM_DIR, GCTA, NBPROC, KEEP_IND, USE_SBATCH

def main():
    """Entry point if called as an executable"""
    in_file_pruned = os.path.join(PRU_DIR, 'all')
    out_dir_grm_maf = os.path.join(GRM_DIR, 'grm-maf')
    maf_intervals = [(0.05, 0.20), (0.20, 0.35), (0.35, 0.50)]

    for maf_int in maf_intervals:
        print('Generating GRM for MAF interval:' + str(maf_int))
        maf_int_char = str(maf_int[0]) + '-' + str(maf_int[1])

        out_file_grm_mafint = os.path.join(out_dir_grm_maf, 'maf' + str(maf_int_char),
                                           'maf.' + str(maf_int_char))

        # Compute GRMs using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=out_file_grm_mafint,
                               gcta=GCTA,
                               per_chr=False,
                               ncpus=NBPROC,
                               other_gcta_par=["--maf", str(maf_int[0]),
                                               "--max-maf", str(maf_int[1]),
                                               "--keep", KEEP_IND],
                               sbatch=USE_SBATCH)

if __name__ == '__main__':
    main()
