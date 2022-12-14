#!/bin/bash

#SBATCH --job-name=spades
#SBATCH --account=co_minium:
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --qos=savio_lowprio
#SBATCH --time=72:00:00

#Run GUY11 assemblies with spades with varying k-mer combinations
#Enable --isolate for high coverage data

cat k-mers.list | while read index kmer; do 
  spades.py -t 56 -k ${kmer} -o isolate_${index} --isolate \
  -1 /global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_1_val_1.fq.gz \
  -2 /global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_2_val_2.fq.gz; 
done

#Disable --isolate 

cat k-mers.list | while read index kmer; do 
  spades.py -t 56 -k ${kmer} -o noisolate_${index} \
  -1 /global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_1_val_1.fq.gz \
  -2 /global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_2_val_2.fq.gz; 
done
