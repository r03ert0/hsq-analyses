#!/usr/bin/env python3
"""Compute heritablility partition per chromosome"""

import os
from configall import (GCTA, MYGCTA, HSQ_DIR, GRM_DIR, PHE_DIR, PHE_LIST, USE_SBATCH, GRM_CUTOFF,
                       NBPROC, QUANT_COVAR, QUAL_COVAR)
import preprocessing


#ipdb.set_trace()
#C-c C-z : open a python shell
#C-c C-c : run the content of the buffer in the opened python shell
#C-c C-r : run the selected region in the python shell
# Simple command
# subprocess.call(['ls', '-1'], shell=True)
# http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html

def main():
    """Entry point if called as an executable"""

    ## quantitative covariables
    qcovar_par = []
    for qcov in QUANT_COVAR:
        qcovar_par.append('--qcovar')
        qcovar_par.append(qcov)

    ## qualitative covariables
    covar_par = []
    for cov in QUAL_COVAR:
        covar_par.append('--covar')
        covar_par.append(cov)


    var_par = qcovar_par + covar_par


    ### ========= per chr =============

    in_dir_grm_perchr = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF), 'all-' +
                                     str(GRM_CUTOFF) + '-chr')
    out_hsq_perchr = os.path.join(HSQ_DIR, 'hsq-perchr')

    if not os.path.exists(out_hsq_perchr):
        os.makedirs(out_hsq_perchr)

    ## extract number of SNPs per chromosome
    in_filenb = open(os.path.join(out_hsq_perchr, 'perchr.nbSNPs.txt'), 'w')
    for chrom in list(range(1, 23)):
        nsnp = preprocessing.read_grm_bin_n(in_dir_grm_perchr + str(chrom))
        in_filenb.write('chr' + str(chrom) + ' ' + str(nsnp)  + '\n')
    in_filenb.close()

    for pheno in PHE_LIST:

        phenopath = os.path.join(PHE_DIR, pheno + '.txt')
        pars = var_par + ['--pheno', phenopath, '--reml']

        for chrom in list(range(1, 23)):

            in_file = in_dir_grm_perchr + str(chrom)
            out_file = os.path.join(out_hsq_perchr, 'chr' + str(chrom) + '.' + pheno)

            preprocessing.gcta_hsq(in_file=in_file,
                                   out_file=out_file,
                                   gcta=GCTA,
                                   mygcta=MYGCTA,
                                   other_gcta_par=pars,
                                   ncpus=NBPROC,
                                   sbatch=USE_SBATCH,
                                   sbatch_par_j="hsq-perchr")


    ## sum across chr
    in_file_allchr = open(os.path.join(out_hsq_perchr, 'perchr.test.txt'), 'w')
    for chrom in list(range(1, 23)):
        in_file_allchr.write(in_dir_grm_perchr + str(chrom) + '\n')
    in_file_allchr.close()

    for pheno in PHE_LIST:

        out_file = os.path.join(out_hsq_perchr, 'allchr' + '.' + pheno)

        phenopath = os.path.join(PHE_DIR, pheno + '.txt')
        pars = var_par + ['--pheno', phenopath, '--reml']

        preprocessing.gcta_hsq(in_file=in_file_allchr.name,
                               out_file=out_file,
                               gcta=GCTA,
                               mygcta=MYGCTA,
                               other_gcta_par=pars,
                               par_input='--mgrm-bin',
                               ncpus=NBPROC,
                               sbatch=USE_SBATCH,
                               sbatch_par_j="hsq-perchr")

if __name__ == '__main__':
    main()
