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

Check fastq and bam consistence in fastq_raw, illumina_correct and genotype_correct.

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#/rhome/cjinfeng/Rice/RIL/FC_raw/*_p1.fq
#/rhome/cjinfeng/Rice/RIL/genotypes_correct/MSU_r7.corrected/*.recal.bam
def extract_RIL_id(files, sufix):
    data = defaultdict(lambda : int())
    for fn in files:
        fn_base = os.path.basename(fn)
        ril_id  = re.sub(r'%s' %(sufix), r'', fn_base)
        data[ril_id] = 1
    return data


def extract_RIL_num(files, sufix):
    data = defaultdict(lambda : int())
    for fn in files:
        fn_base = os.path.basename(fn)
        ril_id  = re.sub(r'%s' %(sufix), r'', fn_base)
        unit    = re.split(r'_', ril_id)
        data[unit[0]] = 1
    return data



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

    fastq_raw_p1 = extract_RIL_id(glob.glob('/rhome/cjinfeng/Rice/RIL/FC_raw/fastq_raw/*_p1.fq'), '_p1.fq')
    fastq_raw_p2 = extract_RIL_id(glob.glob('/rhome/cjinfeng/Rice/RIL/FC_raw/fastq_raw/*_p2.fq'), '_p2.fq')
    illumina_p1  = extract_RIL_id(glob.glob('/rhome/cjinfeng/Rice/RIL/Illumina_correct/*/*_p1.fq'), '_p1.fq')
    illumina_p2  = extract_RIL_id(glob.glob('/rhome/cjinfeng/Rice/RIL/Illumina_correct/*/*_p2.fq'), '_p2.fq')
    bam          = extract_RIL_id(glob.glob('/rhome/cjinfeng/Rice/RIL/genotypes_correct/MSU_r7.corrected/*.recal.bam'), '.recal.bam')
    
    print 'fastq_raw: %s fq1, %s fq2' %(len(fastq_raw_p1.keys()), len(fastq_raw_p2.keys()))
    print 'illumina:  %s fq1, %s fq2' %(len(illumina_p1.keys()), len(illumina_p2.keys()))
    print 'bam:       %s bam' %(len(bam.keys()))

    #number of RILs
    bam_rils = extract_RIL_num(glob.glob('/rhome/cjinfeng/Rice/RIL/genotypes_correct/MSU_r7.corrected/*.recal.bam'), '.recal.bam')
    ofile = open('File_consistence.RILs.list', 'w')
    print 'Number of RILs: %s' %(len(bam_rils.keys()))
    for ril in sorted(bam_rils.keys()):
        print >> ofile, ril
    ofile.close() 

if __name__ == '__main__':
    main()

