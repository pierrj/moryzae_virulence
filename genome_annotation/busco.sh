#!/bin/bash
#SBATCH --job-name=bu_prot
#SBATCH --account=co_minium
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00

#BUSCO version 5.4.2
#Directory: /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/orthogrouping/all_proteomes_processed

conda activate busco

files=$(ls *.faa)
for file in $files
do
  fasta=$file
  prefix=$(echo $file | sed 's/.faa//g')
  busco -c 20 -l sordariomycetes_odb10 -i ${fasta} -o ${prefix} -m protein
done
