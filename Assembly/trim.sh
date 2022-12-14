#!/bin/bash

#SBATCH --job-name=trim
#SBATCH --account=co_minium
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --time=48:00:00

ls *_1.fq.gz | cut -d "_" -f -4 | sort -u | while read prefix; do
		/global/scratch/users/skyungyong/Software/TrimGalore-0.6.6/trim_galore \
		-q 20 \
		-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-a2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
		-j 10 --length 80 --paired $prefix\_1.fq.gz $prefix\_2.fq.gz ; 
done
