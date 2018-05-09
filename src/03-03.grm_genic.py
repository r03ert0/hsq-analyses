#!/usr/bin/env python3

"""GRM using SNPs located in genic regions only"""

import os
import preprocessing
from configall import PRU_DIR, GRM_DIR, EXTERNAL_DIR, PLINK, GCTA, KEEP_IND, NBPROC, USE_SBATCH

def main():
    """Entry point if called as an executable"""
    #in_prefix = 'all'
    in_file_pruned = os.path.join(PRU_DIR, 'all')
    #in_dir_grm_filtered = os.path.join(GRM_DIR, 'all-0.025')
    in_genicregions = os.path.join(EXTERNAL_DIR, 'hg19genes.txt')
    in_geneids = os.path.join(EXTERNAL_DIR, 'genic.txt')
    out_dir_grm_genic = os.path.join(GRM_DIR, 'grm-genic')
    out_dir_grm_nongenic = os.path.join(GRM_DIR, 'grm-genic')

    for margin in [0, 20, 50]:
        ## Extract SNPs located in the genic regions
        preprocessing.plink_make_set(bfile=in_file_pruned,
                                     border=margin,
                                     subset=in_geneids,
                                     gene_set=in_genicregions,
                                     out_dir=out_dir_grm_genic,
                                     out_prefix='genic-margin'+str(margin),
                                     plink=PLINK
                                    )

        listsnp_margin = os.path.join(out_dir_grm_genic, 'genic-margin' + str(margin) + '.snplist')
        out_file_grm_genic_margin = os.path.join(out_dir_grm_genic, 'genic-margin' + str(margin),
                                                 'genic-margin' + str(margin))
        out_file_grm_nongenic_margin = os.path.join(out_dir_grm_nongenic, 'nongenic-margin'
                                                    + str(margin), 'nongenic-margin' + str(margin))


        ## Compute genic GRM using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=out_file_grm_genic_margin,
                               gcta=GCTA,
                               per_chr=False,
                               ncpus=NBPROC,
                               other_gcta_par=["--extract", str(listsnp_margin),
                                               "--keep", KEEP_IND],
                               sbatch=USE_SBATCH
                              )

        ## Compute non genic GRM using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=out_file_grm_nongenic_margin,
                               gcta=GCTA,
                               per_chr=False,
                               ncpus=NBPROC,
                               other_gcta_par=["--exclude", str(listsnp_margin),
                                               "--keep", KEEP_IND],
                               sbatch=USE_SBATCH
                              )

if __name__ == '__main__':
    main()
