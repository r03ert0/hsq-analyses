#!/usr/bin/env python3
"""Compute heritability partition for genic / non-genic"""

import os
from configall import (GCTA, MYGCTA, GRM_DIR, PHE_DIR, PHE_LIST, HSQ_DIR,
                       NBPROC, USE_SBATCH, QUAL_COVAR, QUANT_COVAR)
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
        qcovar_par.append('--qcovar')
        qcovar_par.append(qcov)

    ## qualitative covariables
    covar_par = []
    for cov in QUAL_COVAR:
        covar_par.append('--covar')
        covar_par.append(cov)

    var_par = qcovar_par + covar_par


    ### ========= 3. Genic/non-genic  =============

    in_genic = os.path.join(GRM_DIR, 'grm-genic') # + -0.025 ?
    in_nongenic = os.path.join(GRM_DIR, 'grm-genic') # + -0.025 ?
    out_genic = os.path.join(HSQ_DIR, 'hsq-genic')

    if not os.path.exists(out_genic):
        os.makedirs(out_genic)

    in_filenb = open(os.path.join(out_genic, 'genic.nbSNPs.txt'), 'w')

    for margin  in [0, 20, 50]:

        out_genic_margin = os.path.join(out_genic, 'genic-margin' + str(margin))

        in_genic_margin = os.path.join(in_genic, 'genic-margin' + str(margin),
                                       'genic-margin' + str(margin))
        in_nongenic_margin = os.path.join(in_nongenic, 'nongenic-margin' + str(margin),
                                          'nongenic-margin' + str(margin))

        # input both genic and non-genic and genic grm for variance partitioning
        print(out_genic_margin + '.test.txt')
        in_file = open(out_genic_margin + '.test.txt', 'w')
        in_file.write(in_genic_margin + '\n')
        in_file.write(in_nongenic_margin)
        in_file.close()

        # save the number of SNPs associated with each subgroup
        lc_genic = preprocessing.read_grm_bin_n(in_genic_margin)
        lc_nongenic = preprocessing.read_grm_bin_n(in_nongenic_margin)
        in_filenb.write('genic-margin' + str(margin) + ' ' + str(lc_genic) + '\n')
        in_filenb.write('nongenic-margin' + str(margin) + ' ' + str(lc_nongenic) + '\n')

        for pheno in PHE_LIST:
            for lrt in [1, 2]:
                # --reml-lrt  1
                # Calculate the log likelihood of a reduce model with one or multiple genetic
                # variance components dropped from the full model and calculate the LRT and p-value.
                # By default, GCTA will always calculate and report the LRT for the first genetic
                # variance component, i.e. --reml-lrt 1, unless you re-specify this option,
                # e.g. --reml-lrt 2 assuming there are a least two genetic variance components
                # included in the analysis. You can also test multiple components simultaneously,
                # e.g. --reml-lrt 1 2 4. See FAQ #1 for more details.
                out_file = out_genic_margin + '.' + str(lrt) + '.' + pheno

                print(pheno)
                phenopath = os.path.join(PHE_DIR, pheno + '.txt')
                pars = var_par + ['--pheno', phenopath, '--reml-lrt', str(lrt), '--reml']

                preprocessing.gcta_hsq(in_file=in_file.name,
                                       out_file=out_file,
                                       gcta=GCTA,
                                       mygcta=MYGCTA,
                                       ncpus=NBPROC,
                                       other_gcta_par=pars,
                                       par_input='--mgrm-bin',
                                       sbatch=USE_SBATCH,
                                       sbatch_par_j="hsq-genic")


    in_filenb.close()


if __name__ == '__main__':
    main()
