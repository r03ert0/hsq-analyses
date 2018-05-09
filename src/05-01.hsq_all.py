#!/usr/bin/env python3
"""Compute heritability for all SNPs"""

import os
from configall import (GCTA, MYGCTA, HSQ_DIR, GRM_DIR, PHE_DIR, PHE_LIST, NBPROC,
                       USE_SBATCH, QUANT_COVAR, QUAL_COVAR)
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
    #qcovar = [os.path.join(PHE_DIR, 'age.txt'),PCS]
    qcovar_par = []
    for qcov in QUANT_COVAR:
        qcovar_par.append('--qcovar')
        qcovar_par.append(qcov)

    ## qualitative covariables
    #covar = [os.path.join(PHE_DIR, 'sex.txt'),os.path.join(PHE_DIR, 'centre.txt')]
    covar_par = []
    for cov in QUAL_COVAR:
        covar_par.append('--covar')
        covar_par.append(cov)

    var_par = qcovar_par + covar_par


    ### ========= 1. All SNPS, 10 PCs =============

    in_file = os.path.join(GRM_DIR, 'grm-all', 'all')
    out_hsq = os.path.join(HSQ_DIR, 'hsq-all', 'all')


    for pheno in PHE_LIST:

        out_file = out_hsq +'.'+pheno
        print('running h^2 gcta estimation for phenotype: ' + pheno)
        pheno = os.path.join(PHE_DIR, pheno+'.txt')
        pars = var_par + ['--pheno', pheno, '--reml']

        preprocessing.gcta_hsq(in_file=in_file,
                               out_file=out_file,
                               gcta=GCTA,
                               other_gcta_par=pars,
                               ncpus=NBPROC,
                               mygcta=MYGCTA,
                               sbatch=USE_SBATCH,
                               sbatch_par_j="hsq-all")

    ## extract number of SNPs per chromosome
    with open(out_hsq + '.nbSNPs.txt', 'w') as in_filenb:
        nsnp = preprocessing.read_grm_bin_n(in_file)
        in_filenb.write('all' + ' ' + str(nsnp)  + '\n')



if __name__ == '__main__':
    main()
