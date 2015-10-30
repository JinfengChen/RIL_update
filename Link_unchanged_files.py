#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Link_unchanged_files.py

Link unchanged RILs from Sofia's results to Illumina_correct and genotypes_correct.
 
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    illumina = '/rhome/cjinfeng/Rice/RIL/Illumina_correct'
    genotype = '/rhome/cjinfeng/Rice/RIL/genotypes_correct/MSU_r7.corrected/'
    table_ucr     = '/rhome/cjinfeng/Rice/RIL/FC_raw/genotype_sample_illuminaID.txt'
    table_cornell = '/rhome/cjinfeng/Rice/RIL/FC_raw/cornell.rename.txt' 
    

if __name__ == '__main__':
    main()

