#!/bin/bash
#SBATCH --job-name=fasttree_scos
#SBATCH --partition=savio3_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=48:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


PROJECT_DIR=/global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping

cd $PROJECT_DIR

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

module purge

SCO_DIR=orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_1/Single_Copy_Orthologue_Sequences
basename=$(basename ${SCO_DIR})
MSA=${basename}_msas_combined.trimmed.fasta

fasttree -gamma -out ${basename}_msas_combined.trimmed.fasttree ${basename}_msas_combined.trimmed.fasta