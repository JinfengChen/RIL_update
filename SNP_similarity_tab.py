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
import subprocess
import multiprocessing as mp


def usage():
    test="name"
    message='''
Use .tab file stread of .SNP to calculate similarity

.tab:
#CHROM  POS     REF     RIL103_0_GGTAGC_FC1213L5
Chr1    31071   A       ./.
Chr1    31478   C       ./.
Chr1    31843   G       A/A


.SNP:
0100021547A     GN190   A
0100031071A     GN190   A
0100031478C     GN190   C



python SNP_similarity.py --cpu 32 > RILs_SNP.similarity
or
qsub -q highmem SNP_similarity.sh

#same ril
awk '$1~/RIL58_/ && $2~/RIL58_/' RILs_SNP.dupli_all.similarity
#all dupli
awk '$3>0.9' RILs_SNP.dupli_all.similarity |grep "NA" -v| less -S
awk '$3>0.9' RILs_SNP.dupli_all_number.similarity | grep "NA" -v > RILs_SNP.dupli_all_number.duplicate.table
#all
awk '$3>0.9' RILs_SNP.similarity | grep "NA" -v | awk '$1!=$2' | less -S
awk '$3>0.9' RILs_SNP.all_all_number.similarity | grep "NA" -v | awk '$1!=$2' > RILs_SNP.all_all_number.duplicate.table
#dupli vs dupli, no duplicate
awk '$3>0.9' RILs_SNP.dupli_dupli_number.similarity |grep "NA" -v | awk '$1!=$2' | less -S

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#0100031071A     GN278   G
def read_snp(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #pos  = int(unit[0][2:-1])
                pos   = unit[0]
                data[pos] = unit[2]
    return data

#Chr1    31843   G       A/A
def read_snp_tab(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                #print line
                unit = re.split(r'\t',line)
                #pos  = int(unit[0][2:-1])
                #pos   = unit[0]
                chrs = re.sub(r'Chr', r'', unit[0])
                pos  = '%02d' %(int(chrs)) + '%08d' %(int(unit[1])) + unit[2]
                bases= re.split(r'\/', unit[3])
                if bases[0] == bases[1] and bases[0] in ['A', 'T', 'C', 'G']:
                    data[pos] = bases[0]
                    #print pos, data[pos]
    return data




#return dict of lib -> SNP
#Maq.p1.map.pileup.SNP
def parse_bam_all(bam_list):
    data = defaultdict(lambda : str())
    for lib in sorted(re.split(r'\n', bam_list)):
        #print lib
        unit = re.split(r' ', lib)
        if len(unit) < 2:
            continue
        #print len(unit)
        #print '%s\t%s' %(unit[-3], unit[-1])
        bam = os.path.split(unit[-1])[1]
        bam = re.sub(r'.recal.bam', r'', bam)
        bam = re.sub(r'.bam', r'', bam)
        snp = unit[-3]
        #print snp
        snp = re.sub(r'\.bam', r'.Maq.p1.map.pileup.SNP', snp)
        #print snp, bam
        data[bam] = snp
    return data  

def snp_similarity(lib1, lib2, snp1, snp2):
    #print lib1, lib2
    #print snp1, snp2
    snp1_dict = read_snp_tab(snp1)
    snp2_dict = read_snp_tab(snp2)
    total = 0
    match = 0
    if len(snp1_dict.keys()) < 2000 or len(snp2_dict.keys()) < 2000:
        return [lib1, lib2, 'NA', 0, 0, len(snp1_dict.keys()), len(snp2_dict.keys())]
    for pos in snp1_dict.keys():
        if snp2_dict.has_key(pos):
            total += 1
            if snp1_dict[pos] == snp2_dict[pos]:
                match += 1
    sim = float(match)/float(total)
    return [lib1, lib2, sim, total, match, len(snp1_dict.keys()), len(snp2_dict.keys())]

def function_helper(args):
    return snp_similarity(*args)

##run function with parameters using multiprocess of #cpu
def mp_pool_function(function, parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(function, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list

#RIL30_0_ATCACG_FC153L7.genotype.tab
def parse_snp_tab(files):
    data = defaultdict(lambda : str())
    for fn in files:
        fn_base = os.path.basename(fn)
        ril_id  = re.sub(r'.genotype.tab', r'', fn_base)
        data[ril_id] = fn
    return data

 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-c', '--cpu')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
 
    if not args.cpu:
        args.cpu = 2

    #/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib
    #bam_all   = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib'), shell=True)
    #bam_dupli = subprocess.check_output('ls -all %s/*.bam' %('/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correction'), shell=True)

    #snp_dupli = parse_bam_all(bam_dupli)
    #snp_all   = parse_bam_all(bam_all)

    #testing function
    #read_snp_tab('../genotypes_correct/MSU_r7.corrected/RIL103_0_GGTAGC_FC1213L5.genotype.tab')

    snp_dupli = parse_snp_tab(glob.glob('/rhome/cjinfeng/Rice/RIL/genotypes_correct/MSU_r7.corrected/*.tab'))
    snp_all   = snp_dupli

    print 'Lib1\tLib2\tSimilarity\tTotal_Shared_SNP_Site\tTotal_Identical_SNP_Sites\tLib1_SNP\tLib2_SNP'
    parameters = []
    for lib_d in sorted(snp_dupli.keys()):
        snp_d = snp_dupli[lib_d]
        for lib_a in sorted(snp_all.keys()):
            snp_a   = snp_all[lib_a]
            #print snp_d, snp_a
            parameters.append([lib_d, lib_a, snp_d, snp_a])
            #snp_sim = snp_similarity(snp_d, snp_a)
            #print '%s\t%s' %(lib_d, lib_a)
    collect_list_list = mp_pool_function(function_helper, parameters, args.cpu)
    for result in collect_list_list:
        print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(result[0], result[1], result[2], result[3], result[4], result[5], result[6])

if __name__ == '__main__':
    main()

