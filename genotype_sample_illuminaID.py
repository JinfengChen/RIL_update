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
python genotype_sample_illuminaID.py

Go through all FC* directories, convert and print a table like genotype_sample_illuminaID.txt 
197	2	flowcell197_lane2_TGACCA	0	TGACCA	RIL116_0
197	2	flowcell197_lane2_AGTCAA	0	AGTCAA	RIL74_0
197	2	flowcell197_lane2_ATCACG	0	ATCACG	RIL112_0
197	2	flowcell197_lane2_TTAGGC	0	TTAGGC	RIL61_0
197	2	flowcell197_lane2_GCCAAT	0	GCCAAT	RIL119_0
197	2	flowcell197_lane2_ATGTCA	0	ATGTCA	RIL77_0
197	2	flowcell197_lane2_GGCTAC	0	GGCTAC	RIL67_0
197	2	flowcell197_lane2_CCGTCC	0	CCGTCC	RIL79_0

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#RIL39_1 8
def convert_id(infile, barcode):
    #print infile
    #RIL39_1 8
    barcode_ril = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                barcode_ril[barcode[unit[1]]] = unit[0]
        

    data = defaultdict(lambda : list())
    r = re.compile(r'flowcell(\d+)_lane(\d+)_pair1_(\w+)\.fastq')
    father_dir = os.path.dirname(infile)
    fastqs     = glob.glob('%s/flowcell*_lane*_pair*_*.fastq' %(father_dir))
    for fq in fastqs:
        if r.search(fq):
            m    = r.search(fq)
            fc   = m.groups(0)[0]
            lane = m.groups(0)[1]
            bc   = m.groups(0)[2]
            data[bc] = [fc, lane, 'flowcell%s_lane%s_%s' %(fc, lane, bc), '0', bc, barcode_ril[bc]]
            print '\t'.join(data[bc])

#1	ATCACG
def read_barcode(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[1]
    return data


def print_file(infile):
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            print line

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    r = re.compile(r'(FC\d+)_')  
    files = glob.glob('%s/FC*/genotype_sample_illuminaID.txt' %('/rhome/cjinfeng/Rice/RIL'))
    barcode = read_barcode('barcodes.txt')

    for f in sorted(files):
        if r.search(f):
            fc = r.search(f).groups(0)[0]
            print '#%s\n#%s' %(fc, f)
            if not fc in ['FC251', 'FC271', 'FC279']:    
                convert_id(f, barcode)
            else:
                print_file(f)
if __name__ == '__main__':
    main()

