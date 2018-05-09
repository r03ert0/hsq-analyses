"""Make config
General configuration of paths and binaries
"""

import sys
import subprocess
import os
import preprocessing

# indicate where the code is running
CONFIG = 'gris'

if CONFIG == 'tars':
    ANNEX_DIR = '/pasteur/projets/policy01/cinq/rto/ukb-annex/'
    #subprocess.call("module load plink/1.90b2m", shell=True)
    #PLINK = '/local/gensoft2/exe/plink/1.90b2m/scripts/plink'
    PLINK = os.path.join(ANNEX_DIR, 'bin/plink_linux_x86_64/plink')
    GCTA = os.path.join(ANNEX_DIR, 'bin/gcta_1.91.3beta/gcta64')
    USE_SBATCH = False
elif CONFIG == 'gris':
    #ANNEX_DIR = '/Volumes/cinq/rto/ukb-annex/'
    ANNEX_DIR = '/Users/roberto/Documents/annex-ukb/'
    PLINK = '/Applications/_Geno/plink_mac-1.9-15aug2017/plink'
    GCTA = '/Applications/_Geno/gcta_1.91.4beta_mac/bin/gcta64'
elif CONFIG == 'bastet':
    ANNEX_DIR = '/Volumes/LaCie/ukb-annex'
    PLINK = '/Applications/plink_1.90b4.9/plink'
    GCTA = '/Applications/gcta_1.91.0beta_mac/bin/gcta64'
elif CONFIG == 'anne':
    ANNEX_DIR = '/Volumes/cinq/rto/ukb-annex/'
    PLINK = '/Users/abiton/tools/plink_mac/plink'
    GCTA = '/Users/abiton/tools/gcta_1.91.2beta_mac/bin/gcta64'
    USE_SBATCH = False
elif CONFIG == 'ra':
    pass

def getannexdir():
    """call git executable to find annex dir"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                          cwd=script_dir).stdout.decode().strip()

if not 'ANNEX_DIR' in globals():
    ANNEX_DIR = getannexdir()
if not 'PLINK' in globals():
    PLINK = 'plink'
if not 'GCTA' in globals():
    GCTA = 'gcta64'
if not 'USE_SBATCH' in globals():
    USE_SBATCH = False

# wrapper for plink
MYPLINK = os.path.join(ANNEX_DIR, "bin", "myplink", "myplink.sh")

# wrapper for gcta
MYGCTA = os.path.join(ANNEX_DIR, "bin", "mygcta", "mygcta.sh")


# indicate which dataset to use
#DATASET = 'ukb-9891'
#DATASET = 'test'
#DATASET = 'ukb-all'
#DATASET = 'IMAGEN123'
DATASET = 'ukb-20986'
#DATASET = 'ADNI12'

# print python version
print(sys.version)

# data subdirectories
RAW_DIR = os.path.join(ANNEX_DIR, 'data', 'raw')
DERIVED_DIR = os.path.join(ANNEX_DIR, 'data', 'derived')
EXTERNAL_DIR = os.path.join(ANNEX_DIR, 'data', 'external')

# derived data subdirectories
GEN_DIR = os.path.join(DERIVED_DIR, DATASET, '01.genotype')
FIL_DIR = os.path.join(DERIVED_DIR, DATASET, '02.filtered')
PRU_DIR = os.path.join(DERIVED_DIR, DATASET, '03.pruned')
GRM_DIR = os.path.join(DERIVED_DIR, DATASET, '04.grm')
GWA_DIR = os.path.join(DERIVED_DIR, DATASET, '05.gwas')
HSQ_DIR = os.path.join(DERIVED_DIR, DATASET, '06.hsq')

# GRM cutoff
GRM_CUTOFF = 0.025

# subjects to keep: those selected as unrelated (used for the --keep parameter of gcta)
KEEP_IND = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF), 'all-' + str(GRM_CUTOFF) + '.grm.id')
# Principal componants to be included as covariates in GWAS and heritability analyses
PCS = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF), 'all-' + str(GRM_CUTOFF) + '.eigenvec')

# number of processors
NBPROC = '3'

# gene lists
SELGENES_NEURODEV = os.path.join(EXTERNAL_DIR, 'neurodev.txt')
SELGENES_CNSEXPRESSION = os.path.join(EXTERNAL_DIR, 'cnsexpression.txt')

# directory containing all phenotypes
PHE_DIR = os.path.join(DERIVED_DIR, DATASET, '00.phenotype')

## ids of the phenotypes of interest, dir_pheno/id.txt will be used to find the phenotype file.
PHE_LIST = ["ICV", "brain", "accumbens", "amygdala", "caudate", "hippocampus",
            "pallidum", "putamen", "thalamus",  "height", "intelligence"]

QUANT_COVAR = [os.path.join(PHE_DIR, 'age.txt'), PCS]
QUAL_COVAR = [os.path.join(PHE_DIR, 'sex.txt'), os.path.join(PHE_DIR, 'centre.txt')]

# ===== Check that covariables and phenotypes are indeed available for this cohort  =====

# check that files exist
quant_sel = [covarfile for covarfile in QUANT_COVAR if os.path.exists(covarfile)]
QUANT_COVAR = quant_sel

qual_sel = [covarfile for covarfile in QUAL_COVAR if os.path.exists(covarfile)]
QUAL_COVAR = qual_sel

pheno_sel = [phenofile for phenofile in PHE_LIST if os.path.exists(os.path.join(PHE_DIR, phenofile + '.txt'))]
PHE_LIST= pheno_sel

# check that files do not contain only missing value,
# minimal number of non missing value accepted, arbitrary number for now
MIN_NONMISSING = 99

quant_sel = [covarfile for covarfile in QUANT_COVAR if preprocessing.pheno_count(covarfile)[0] > MIN_NONMISSING ]
QUANT_COVAR = quant_sel

qual_sel = [covarfile for covarfile in QUAL_COVAR if preprocessing.pheno_count(covarfile)[0] > MIN_NONMISSING ]
QUAL_COVAR = qual_sel

pheno_sel = []
for phenofile in PHE_LIST:
    if preprocessing.pheno_count(os.path.join(PHE_DIR, phenofile + '.txt'))[0] > MIN_NONMISSING:
        pheno_sel.append(phenofile)
PHE_LIST = pheno_sel


print("Phenotypes used for dataset " + DATASET + " are: " + ', '.join(PHE_LIST))
print("Qualitative covariables used for dataset " + DATASET + " are: " + ', '.join(QUAL_COVAR))
print("Quantitative covariables used for dataset " + DATASET + " are: " + ', '.join(QUANT_COVAR))

