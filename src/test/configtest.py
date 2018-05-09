"""Make config"""

import subprocess
import os

# indicate where the code is running
CONFIG = 'ra'

if CONFIG == 'tars':
    subprocess.call(["module load plink/1.90b2m"])
    ANNEX_DIR = '/pasteur/projets/policy01/cinq/rto/ukb-annex/'
    PLINK = '/local/gensoft2/exe/plink/1.90b2m/scripts/plink'
    GCTA = '/pasteur/projets/policy01/cinq/rto/ukb-annex/bin/gcta_1.90.1beta/gcta64'
elif CONFIG == 'gris':
    ANNEX_DIR = '/Volumes/cinq/rto/ukb-annex/'
    PLINK = '/Applications/_Geno/plink_mac-1.9-15aug2017/plink'
    GCTA = '/Applications/_Geno/gcta_1.90.1beta_mac/bin/gcta64'
elif CONFIG == 'bastet':
    ANNEX_DIR = '/Volumes/cinq/rto/ukb-annex/'
elif CONFIG == 'ra':
    ANNEX_DIR = os.path.expanduser('~/ukb-annex')
    PLINK = 'plink'
    GCTA = 'gcta64'

# indicate which dataset to use
DATASET = 'test'

# data subdirectories
RAW_DIR = os.path.join(ANNEX_DIR, 'data/raw')
DERIVED_DIR = os.path.join(ANNEX_DIR, 'data/derived')
EXTERNAL_DIR = os.path.join(ANNEX_DIR, 'data/external')

# directory containing all phenotypes
PHE_DIR = os.path.join(DERIVED_DIR, DATASET, '00.phenotype')

# derived data subdirectories
GEN_DIR = os.path.join(DERIVED_DIR, DATASET, '01.genotype')

# number of processors
NBPROC = 4

# plink input file
IN_FILE = os.path.join(RAW_DIR, "1000genomes",
                       "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes")

# number of extracted SNPs
N_SNP = 20000

# simulated phenotypes and their heritability
PHENOTYPES = {
    "height":.6,
    "intelligence":.3,
    "age":0,
    "ICV":.4,
    "brain":.5,
    "thalamus":.4,
    "caudate":.4,
    "putamen":.4,
    "pallidum":.4,
    "hippocampus":.4,
    "amygdala":.4,
    "accumbens":.4
}
