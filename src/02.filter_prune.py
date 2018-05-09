#!/usr/bin/env python3

"""Filtering/pruning of the SNP data"""

import os
from configall import PLINK, DERIVED_DIR, DATASET
import preprocessing

print(PLINK)

def main():
    '''Entry point if called as an executable'''
    in_dir = os.path.join(DERIVED_DIR, DATASET)
    in_prefix = 'all'

    preprocessing.filter_prune(
        in_dir=in_dir,
        in_prefix=in_prefix,
        plink=PLINK)

if __name__ == '__main__':
    main()
