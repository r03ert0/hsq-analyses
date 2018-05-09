#!/usr/bin/env python3

"""GRMs using lists of genes involved in various processes"""

import os
import preprocessing
from configall import (SELGENES_NEURODEV, SELGENES_CNSEXPRESSION, PRU_DIR, EXTERNAL_DIR, GCTA,
                       PLINK, GRM_DIR, NBPROC, KEEP_IND, USE_SBATCH)

def main():
    """Entry point if called as an executable"""
    in_margin = 50
    for namesel, filesel in [
            ('neurodev', SELGENES_NEURODEV),
            ('cnsexpression', SELGENES_CNSEXPRESSION)
        ]:
        print("Gene set name:", namesel)
        print("Gene set file:", filesel)

        in_file_pruned = os.path.join(PRU_DIR, 'all')
        in_genicregions = os.path.join(EXTERNAL_DIR, 'hg19genes.txt')
        out_dir_grm_genesel = os.path.join(GRM_DIR, 'grm-' + namesel)
        out_dir_grm_genic = os.path.join(out_dir_grm_genesel, namesel + '-margin' + str(in_margin))
        out_dir_grm_nongenic = os.path.join(out_dir_grm_genesel, 'non' + namesel + '-margin' +
                                            str(in_margin))

        if not os.path.exists(out_dir_grm_genesel):
            os.makedirs(out_dir_grm_genesel)

        ## extract SNPs located in the gene set
        preprocessing.plink_make_set(bfile=in_file_pruned,
                                     border=in_margin,
                                     subset=filesel,
                                     gene_set=in_genicregions,
                                     out_dir=out_dir_grm_genesel,
                                     out_prefix=namesel,
                                     plink=PLINK
                                    )

        listsnp_margin = os.path.join(out_dir_grm_genesel, namesel + '.snplist')

        ## Compute GRM using SNPs in those genes using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=os.path.join(out_dir_grm_genic, namesel),
                               gcta=GCTA,
                               per_chr=False,
                               ncpus=NBPROC,
                               other_gcta_par=["--extract", str(listsnp_margin),
                                               "--keep", KEEP_IND],
                               sbatch=USE_SBATCH)

        ## Compute GRM using SNPs not in those genes using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=os.path.join(out_dir_grm_nongenic, 'non'+namesel),
                               gcta=GCTA,
                               per_chr=False,
                               ncpus=NBPROC,
                               other_gcta_par=["--exclude", str(listsnp_margin),
                                               "--keep", KEEP_IND],
                               sbatch=USE_SBATCH)

if __name__ == '__main__':
    main()
