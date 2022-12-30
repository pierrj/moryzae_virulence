#!/bin/bash
#SBATCH --job-name=bu_geno
#SBATCH --account=ac_kvkallow
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00

#BUSCO version 5.4.2
#Directory: /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/orthogrouping/Assembly_renamed
conda activate busco

files=$(ls *.fasta)
for file in $files
do
  fasta=$file
  prefix=$(echo $file | sed 's/.fasta//g')
  busco -c 20 -l sordariomycetes_odb10 -i ${fasta} -o ${prefix} -m genome --augustus --augustus_species magnaporthe_grisea
done
