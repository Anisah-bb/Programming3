#!/bin/bash

export paired_NGS_1=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq
export paired_NGS_2=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq


output=/students/2021-2022/master/Firdaws_DSLS/output/

mkdir -p ${output}
output_folder=${output}Firdaws_



# parallelise the task for seveal k-mer sizes

seq 29 2 31 | parallel -j16 "velveth ${output_folder}{} {} -longPaired -fastq -separate ${paired_NGS_1} ${paired_NGS_2} &&  velvetg ${output_folder}{} && cat ${output_folder}{}/contigs.fa | python3 assignment4.py -kmers {} >> /students/2021-2022/master/Firdaws_DSLS/output/kmers.csv"

output=/homes/fabadmus/programming3/Programming3/Assignment4/

mkdir -p ${output}

python3 best_kmer.py


#rm -r $output_folder*
