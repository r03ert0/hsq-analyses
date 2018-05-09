#!/usr/bin/env python3
"""compute heritability without including PCs as covariates"""

#from future import print_function
import os
from configall import (GCTA, MYGCTA, HSQ_DIR, GRM_DIR, PHE_DIR, PHE_LIST, NBPROC, PCS,
                       USE_SBATCH, QUANT_COVAR, QUAL_COVAR)
import preprocessing
#import ipdb


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
        if qcov == PCS:
            continue
        qcovar_par.append('--qcovar')
        qcovar_par.append(qcov)

    ## qualitative covariables
    covar_par = []
    for cov in QUAL_COVAR:
        covar_par.append('--covar')
        covar_par.append(cov)


    var_par = qcovar_par + covar_par


    ### ========= 1bis. All SNPS, no PC =============

    in_file = os.path.join(GRM_DIR, 'grm-all/all')
    out_dir = os.path.join(HSQ_DIR, 'hsq-nopca', 'nopca')

    for pheno in PHE_LIST:
        out_file = out_dir+'.'+pheno
        print('running h^2 gcta estimation (no PCs) for phenotype: ' + pheno)
        pheno = os.path.join(PHE_DIR, pheno + '.txt')
        pars = var_par + ['--pheno', pheno, '--reml']

        preprocessing.gcta_hsq(in_file=in_file,
                               out_file=out_file,
                               gcta=GCTA,
                               mygcta=MYGCTA,
                               ncpus=NBPROC,
                               other_gcta_par=pars,
                               sbatch=USE_SBATCH,
                               sbatch_par_j="hsq-nopca")


if __name__ == '__main__':
    main()
