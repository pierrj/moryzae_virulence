#!/bin/bash

#SBATCH --job-name=velvet
#SBATCH --account=co_minium
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --qos=savio_lowprio
#SBATCH --time=72:00:00

#Test GUY11 assembly with velvet
#Hash seems to have an upper limit and won't inrease after ~31, 33

for k in $(seq 15 2 149); do 
  velveth hash_${k} ${k} \
  -shortPaired -fastq.gz \
  -separate /global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_1_val_1.fq.gz \
  /global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_2_val_2.fq.gz; 
done
