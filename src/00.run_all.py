#!/usr/bin/env python3

""" Master script that runs in the appropriate order all the scripts necessary to
    reproduce the derived data.
    @Author Roberto Toro, 1 Octobre 2017
"""

# @todo Should exec using
# slurm --mem 10GB [file.py] |tee file.log.txt

## pylint: disable=C0103
# or else we have to put everything in a main function



import os
import subprocess
from configall import GRM_DIR, HSQ_DIR


subprocess.run("./01.download_external.py", check=True)

subprocess.run("./02.filter_prune.py", check=True)

subprocess.run("./03-01.grm_all.py", check=True)
subprocess.run("./03-02.grm_pca.py", check=True)
subprocess.run("./03-03.grm_genic.py", check=True)
subprocess.run("./03-04.grm_maf.py", check=True)
subprocess.run("./03-05.grm_gene_lists.py", check=True)

subprocess.run("./04.gwas.py", check=True)


subprocess.run("./05-01.hsq_all.py", check=True)
subprocess.run("./05-02.hsq_nopca.py", check=True)
subprocess.run("./05-03.hsq_genic.py", check=True)
subprocess.run("./05-04.hsq_maf.py", check=True)
subprocess.run("./05-05.hsq_gene_lists.py", check=True)
subprocess.run("./05-06.hsq_bivariate.py", check=True)
subprocess.run("./05-07.hsq_perchr.py", check=True)


# to be tested:
hsqsummarydir = os.path.join(HSQ_DIR, 'hsq-summary')
os.makedirs(hsqsummarydir, exist_ok=True)
subprocess.run(["bash", "06.summarize_hsq.sh", HSQ_DIR, GRM_DIR,
                os.path.join(hsqsummarydir, 'summary')], check=True)
