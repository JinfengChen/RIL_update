#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l mem=60gb
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

python SNP_similarity_tab.py --cpu $PBS_NP > RILs_SNP.similarity

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

