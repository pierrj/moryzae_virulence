#!/bin/bash
#SBATCH --job-name=generate_sco_msa
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
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

# ls ${SCO_DIR} > sco_mapfile

# if [ -d "${basename}_msas" ]; then
#     rm -r ${basename}_msas
# fi

# mkdir ${basename}_msas

# while read orthogroup; do
#     echo ${orthogroup}
#     mafft --maxiterate 1000 --globalpair --thread ${SLURM_NTASKS} ${SCO_DIR}/${orthogroup} > \
#         ${basename}_msas/${orthogroup}.msa
# done < sco_mapfile

/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/scratch/users/pierrj/moryzae_virulence/orthogrouping/concatenate_msas.py ${basename}_msas ${basename}_msas_combined.fasta

module load trimal

trimal -gt 0.8 -in ${basename}_msas_combined.fasta -out ${basename}_msas_combined.trimmed.fasta