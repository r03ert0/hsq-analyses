#!/usr/bin/env python3

"""Download gene lists from the web"""

import os
from urllib.request import urlopen

from configall import EXTERNAL_DIR, SELGENES_CNSEXPRESSION, SELGENES_NEURODEV

import pymysql
#pymysql.install_as_MySQLdb()
#import MySQLdb

def main():
    """Entry point if called as an executable"""
    rs_file = os.path.join(EXTERNAL_DIR, 'refseq19.txt')
    genes_file = os.path.join(EXTERNAL_DIR, 'hg19genes.txt')
    genic_file = os.path.join(EXTERNAL_DIR, 'genic.txt')
    if not (os.path.exists(rs_file) and os.path.exists(genes_file) and os.path.exists(genic_file)):
        # get hg19 genes from UCSC server
        con = pymysql.connect(host='genome-mysql.cse.ucsc.edu', user='genome', db='hg19')
        cur = con.cursor()
        cur.execute("SELECT chrom,txStart,txEnd,strand,cdsStart,cdsEnd,name,name2 from refGene")
        tmp = cur.fetchall()
        con.close()

        # remove data in non straightforward chromosomes
        result = []
        for row in tmp:
            if row[0][0:4] != "chrU" and not '_' in row[0]:
                result.append(row)

        # save data in refseq19.txt
        file = open(rs_file, 'w')
        for row in result:
            file.write("%s\n" % ("\t".join(map(str, row))))
        file.close()

        # extract only chr, start bp, end bp and refseq gene names, save to hg19genes.txt
        file = open(genes_file, 'w')
        for row in result:
            file.write("%s\n" % ("\t".join(map(str, [row[0], row[1], row[2], row[7]]))))
        file.close()

        # save a list of all gene names to us in genic enrichment
        file = open(genic_file, 'w')
        sorted_result = []
        for row in result:
            sorted_result.append(row[7])
        sorted_result = sorted(sorted_result)
        for row in sorted_result:
            file.write("%s\n" % str(row))
        file.close()
    if not os.path.exists(SELGENES_CNSEXPRESSION):
        # get genes used for cns expression enrichment
        response = urlopen('https://github.com/r03ert0/ENIGMA-GCTA/raw/master/lists/'
                           'cnsexpression.txt')
        file = open(SELGENES_CNSEXPRESSION, 'w')
        file.write(response.read().decode("utf-8"))
        file.close()
    if not os.path.exists(SELGENES_NEURODEV):
        # get genes used for neurodev enrichment
        response = urlopen('https://github.com/r03ert0/ENIGMA-GCTA/raw/master/lists/neurodev.txt')
        file = open(SELGENES_NEURODEV, 'w')
        file.write(response.read().decode("utf-8"))
        file.close()

if __name__ == '__main__':
    main()
