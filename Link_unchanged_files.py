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

Link unchanged RILs to Illumina_correct and genotypes_correct.
Fastq from fastq_raw
Bam from Sofia's results 

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def createdir(outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir) 

#FC133
#/rhome/cjinfeng/Rice/RIL/FC133_RIL_39/genotype_sample_illuminaID.txt
#133     7       flowcell133_lane7_ACTTGA        0       ACTTGA  RIL39_1
#133     7       flowcell133_lane7_GATCAG        0       GATCAG  RIL39_2

#Fastq
#RIL103_0_GGTAGC_FC1213L5_p1.fq  RIL103_0_GGTAGC_FC1213L5_p2.fq
#Bam
#RIL152_0_ATTCCT_FC197L3.*
def update_RIL(infile, illumina, genotype, raw, bam_dir):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                #print line
                unit = re.split(r'\t',line)
                ril  = unit[5] if unit[5].startswith(r'RIL') else 'RIL%s' % (unit[5])
                ril_dir = '%s/%s' %(illumina, ril)
                createdir(ril_dir)
                #raw fastq
                fq1_raw = '%s/%s_R1.fastq' %(raw, unit[2].rstrip())
                fq2_raw = '%s/%s_R2.fastq' %(raw, unit[2].rstrip())
                if unit[2].startswith(r'flowcell'):
                    word = re.split(r'\_', unit[2]) 
                    fq1_raw = '%s/%s.fastq' %(raw, '_'.join([word[0], word[1], 'pair1', word[2]]))
                    fq2_raw = '%s/%s.fastq' %(raw, '_'.join([word[0], word[1], 'pair2', word[2]]))
                if not os.path.isfile(fq1_raw):
                    print 'file not found: %s' %(fq1_raw)
                #else:
                #    print 'file found: %s' %(fq1_raw)
                if not os.path.isfile(fq2_raw):
                    print 'file not found: %s' %(fq2_raw)
                #linked fastq with RIL id and barcode
                prefix = '%s_%s_FC%sL%s' %(ril, unit[4], unit[0], unit[1])
                fq1 = '%s/%s_p1.fq' %(raw, prefix)
                fq2 = '%s/%s_p2.fq' %(raw, prefix)
                if not os.path.isfile(fq1):
                    os.system('ln -s %s %s' %(fq1_raw, fq1))
                if not os.path.isfile(fq2):
                    os.system('ln -s %s %s' %(fq2_raw, fq2))
                #Sofia's bam
                bam = '%s/%s.recal.bai' %(bam_dir, prefix)
                bam_all = '%s/%s.*' %(bam_dir, prefix)
                #illumina/genotype_correct
                fq1_target = '%s/%s/%s_p1.fq' %(illumina, ril, prefix)
                fq2_target = '%s/%s/%s_p2.fq' %(illumina, ril, prefix)
                bam_target = '%s/%s.recal.bai' %(genotype, prefix)
                bam_all_target = '%s' %(genotype)
                #print '%s -> %s' %(fq1, fq1_target)
                #print '%s -> %s' %(fq2, fq2_target)
                #print '%s -> %s' %(bam, bam_target)
                if not os.path.isfile(fq1_target):
                    print fq1_target
                    os.system('ln -s %s %s' %(fq1, fq1_target))
                if not os.path.isfile(fq2_target):
                    print fq2_target
                    os.system('ln -s %s %s' %(fq2, fq2_target))
                if not os.path.isfile(bam_target):
                    print bam_target
                    #print 'ck: %s' %(line)
                    #os.system('ln -s %s %s' %(bam, bam_target))
                    os.system('ln -s %s %s' %(bam_all, bam_all_target))
    return data

#RIL168_0 switch with RIL187_0, 20151029
#0813    1       6204_N_RILLib83-230_ATCACG      1       ATCACG  187_0
#0813    1       6204_N_RILLib83-230_ATTCCT      27      ATTCCT  83_0
def update_cornell(infile):
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

    raw      = '/rhome/cjinfeng/Rice/RIL/FC_raw/fastq_raw'
    bam_sofia= '/rhome/cjinfeng/Rice/RIL/genotypes/MSU_r7.corrected'
    illumina = '/rhome/cjinfeng/Rice/RIL/Illumina_correct'
    genotype = '/rhome/cjinfeng/Rice/RIL/genotypes_correct/MSU_r7.corrected'
    table_ucr     = '/rhome/cjinfeng/Rice/RIL/FC_raw/genotype_sample_illuminaID.txt'
    table_cornell = '/rhome/cjinfeng/Rice/RIL/FC_raw/cornell.rename.txt' 
    
    update_RIL(table_ucr, illumina, genotype, raw, bam_sofia)
    update_RIL(table_cornell, illumina, genotype, raw, bam_sofia)

if __name__ == '__main__':
    main()

