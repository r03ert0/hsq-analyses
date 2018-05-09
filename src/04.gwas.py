#!/usr/bin/env python3

"""Compute GWASs
"""

import os
import pandas as pd
from configall import GWA_DIR, PLINK, PRU_DIR, PHE_DIR, PCS, MYPLINK, PHE_LIST, FIL_DIR, USE_SBATCH
import preprocessing

def main():
    """Entry point if called as an executable"""

    in_prefix = 'all'
    in_file_allsnps = os.path.join(FIL_DIR, in_prefix)
    out_dir_gwas_allsnps = os.path.join(GWA_DIR, 'gwas-all')
    log_dir_gwas_allsnps = os.path.join(out_dir_gwas_allsnps, 'log')
    in_file_prunedsnps = os.path.join(PRU_DIR, in_prefix)
    out_dir_gwas_prunedsnps = os.path.join(GWA_DIR, 'gwas-pruned')
    log_dir_gwas_prunedsnps = os.path.join(out_dir_gwas_prunedsnps, 'log')

    ## use sex, center and age as covariates

    ## filter individuals in centre.txt keeping only those in the genotype file
    fam_table = pd.read_table(in_file_allsnps + ".fam", delim_whitespace=True,
                              names=['FID', 'IID', 'PID', 'MID', 'Gender', 'Phenotype'])
    centre_table = pd.read_table(os.path.join(PHE_DIR, "centre.txt"), delim_whitespace=True,
                                 index_col=False)
    centre_table = centre_table[centre_table.IID.isin(fam_table.IID)]
    centre_table.to_csv(os.path.join(PHE_DIR, "centre.cov"), sep='\t', index=False)

    # slurm configuration
    if USE_SBATCH:
        smode = "sbatch"
    else:
        smode = "direct"

    ## All SNPs

    os.makedirs(log_dir_gwas_allsnps, exist_ok=True)

    for pheno in PHE_LIST:
        out_prefix = 'all.' + pheno
        preprocessing.run([MYPLINK, PLINK,
                           "--bfile", in_file_allsnps,
                           "--allow-no-sex", "--linear", "hide-covar",
                           "--pheno", os.path.join(PHE_DIR, pheno+".txt"),
                           "--qcovar", os.path.join(PHE_DIR, "age.txt"),
                           "--qcovar", PCS,
                           "--covar", os.path.join(PHE_DIR, "centre.cov"),
                           "--qcovar", os.path.join(PHE_DIR, "sex.txt"),
                           "--out", os.path.join(out_dir_gwas_allsnps, out_prefix)],
                          mode=smode,
                          slurm_par=["-J", "gwas",
                                     "-D", log_dir_gwas_allsnps])

    ## Pruned SNPs

    os.makedirs(log_dir_gwas_prunedsnps)

    for pheno in PHE_LIST:
        out_prefix = 'all.' + pheno
        preprocessing.run([MYPLINK, PLINK,
                           "--bfile", in_file_prunedsnps,
                           "--allow-no-sex", "--linear", "hide-covar",
                           "--pheno", os.path.join(PHE_DIR, pheno+".txt"),
                           "--qcovar", os.path.join(PHE_DIR, "age.txt"),
                           "--qcovar", PCS,
                           "--covar", os.path.join(PHE_DIR, "centre.cov"),
                           "--qcovar", os.path.join(PHE_DIR, "sex.txt"),
                           "--out", os.path.join(out_dir_gwas_prunedsnps, out_prefix)],
                          mode=smode,
                          slurm_par=["-J", "gwas",
                                     "-D", log_dir_gwas_prunedsnps])

if __name__ == '__main__':
    main()
