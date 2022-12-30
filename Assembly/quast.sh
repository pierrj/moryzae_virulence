#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --account=ac_kvkallow
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00

#Quast version 5.2.0
#Directory: /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/orthogrouping/Assembly_renamed

files=$(ls *.fasta)
for file in $files
do
  prefix=$(echo $file | sed 's/.fasta//g')
  quast.py -t 20 -o ${prefix}\.Quast --fast ${file}
done
