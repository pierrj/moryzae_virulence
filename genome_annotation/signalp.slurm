#!/bin/bash
#SBATCH --job-name=signalp
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

cd /global/scratch/users/pierrj/moryzae_virulence_project/signalp

ls /global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping/all_proteomes_processed/ > mapfile

if [ -d "jobqueue" ]; then
    rm jobqueue
fi

while read file; do
    echo "/global/scratch/users/pierrj/signalp-4.1/signalp -t euk -u 0.34 -U 0.34 -f short /global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping/all_proteomes_processed/$file > $file.signalp.out" >> jobqueue
done < mapfile

parallel -j ${SLURM_NTASKS} < jobqueue