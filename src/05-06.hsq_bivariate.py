#!/usr/bin/env python3
"""bivariate heritability analysis"""

import os
from configall import (GCTA, MYGCTA, HSQ_DIR, GRM_DIR, PHE_DIR, PHE_LIST,
                       NBPROC, USE_SBATCH, QUANT_COVAR, QUAL_COVAR)
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
    ## Run if .hsq file was already computed ?
    # If TRUE, will compute again the hsq,
    # otherwise won't compute the hsq if the .hsq file is present
    overwrite = False

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



    ### ========= 2. All SNPS, 10 pcs, bivariate analysis =============

    out_biv = os.path.join(HSQ_DIR, 'hsq-biv')
    in_file = os.path.join(GRM_DIR, 'grm-all-0.025/all-0.025')
    treated = []

    for pheno1 in PHE_LIST:
        for pheno2 in PHE_LIST:

            treated.append(pheno1+'_'+pheno2)
            if (pheno1 != pheno2 and (not pheno1+'.'+pheno2 in treated) and
                    (not pheno2+'.'+pheno1 in treated)):

                pheno1path = os.path.join(PHE_DIR, pheno1 + '.txt')
                pheno2path = os.path.join(PHE_DIR, pheno2 + '.txt')
                phenopair = os.path.splitext(pheno1)[0] + '.' + os.path.splitext(pheno2)[0]

                out_biv_all = os.path.join(out_biv, phenopair)

                pars = var_par + ['--pheno', pheno1path,
                                  '--pheno', pheno2path,
                                  '--reml-maxit', str(200),
                                  '--reml-bivar', '--reml-bendV']

                if not os.path.isfile(out_biv_all+'.hsq') and not overwrite:
                    preprocessing.gcta_hsq(in_file=in_file,
                                           out_file=out_biv_all,
                                           gcta=GCTA,
                                           mygcta=MYGCTA,
                                           other_gcta_par=pars,
                                           ncpus=NBPROC,
                                           sbatch=USE_SBATCH,
                                           sbatch_par_j="hsq-biv",
                                           sbatch_par_p="common",
                                           sbatch_par_qos="normal")


                out_biv_rg0 = os.path.join(out_biv, 'all.rg=0.' + phenopair)
                pars_rg0 = pars + ['--reml-bivar-lrt-rg', str(0)]

                if not os.path.isfile(out_biv_rg0+'.hsq') and not overwrite:
                    preprocessing.gcta_hsq(in_file=in_file,
                                           out_file=out_biv_rg0,
                                           gcta=GCTA,
                                           mygcta=MYGCTA,
                                           ncpus=NBPROC,
                                           other_gcta_par=pars_rg0,
                                           sbatch=USE_SBATCH,
                                           sbatch_par_j="hsq-biv",
                                           sbatch_par_p="common",
                                           sbatch_par_qos="normal")


                out_biv_rg1 = os.path.join(out_biv, 'all.rg=1.' + phenopair)
                pars_rg1 = pars + ['--reml-bivar-lrt-rg', str(1)]

                if not os.path.isfile(out_biv_rg1+'.hsq') and not overwrite:
                    preprocessing.gcta_hsq(in_file=in_file,
                                           out_file=out_biv_rg1,
                                           gcta=GCTA,
                                           mygcta=MYGCTA,
                                           ncpus=NBPROC,
                                           other_gcta_par=pars_rg1,
                                           sbatch=USE_SBATCH,
                                           sbatch_par_j="hsq-biv",
                                           sbatch_par_p="common",
                                           sbatch_par_qos="normal")

if __name__ == '__main__':
    main()
