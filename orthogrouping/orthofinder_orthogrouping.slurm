#!/bin/bash
#SBATCH --job-name=full_orthofinder_run
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

PROJECT_DIR=/global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping

cd $PROJECT_DIR

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module purge

orthofinder -oa -M msa -A mafft -T fasttree -t ${SLURM_NTASKS} -a ${SLURM_NTASKS} -n out -b orthofinder_out/Results_out/WorkingDirectory