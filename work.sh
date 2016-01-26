echo "prepare links of raw fastq from flowcell"
ln -s /bigdata/wesslerlab/shared/Rice/RIL/FC*/flowcell*.fastq ./
ln -s /bigdata/wesslerlab/shared/Rice/RIL/cornell_RILs_*/*.fastq ./
mkdir fastq_raw
mv *.fastq fastq_raw/

echo "merge barcode and ID information"
python genotype_sample_illuminaID.py > genotype_sample_illuminaID.raw.txt
cat /bigdata/wesslerlab/shared/Rice/RIL/cornell_RILs_*/cornell.rename.txt > cornell.rename.raw.txt

echo "edit genotype_sample_illuminaID.txt, cornell.rename.txt. copy lines need to reprocess to fastq_raw and run call_genotype"
echo "libraries need to achieve were list in achieve.list and comment out in these two txt file"
echo "we will redo the rename and genotype for these need to reprocess"
cd fasta_raw
bash call_genotypes.sh 
make all
echo "then we need to link these file that no need to change in, test if file already there, we will skip"
python Link_unchanged_files.py
python File_consistence.py

echo "ckeck RILs ID"
cut -f1 ~/BigData/00.RD/RILs/QTL_pipe/input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.ALL | grep "GN" | sed 's/GN-/RIL/' > File_consistence.RILs.275
python File_consistence.missing.py

echo "comfirm uniqueness"
qsub SNP_similarity_tab.sh
awk '$3>0.9' RILs_SNP.similarity |grep "NA" -v | awk '$1!=$2' | less -S
 
