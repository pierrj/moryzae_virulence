#!/bin/bash


PROJECT_DIR=/global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping

cd $PROJECT_DIR

module purge

module load python

conda activate /global/scratch/users/pierrj/conda_envs/mov_orthofinder

## copy over fungap outputs and process

if [ -d "all_gffs" ]; then
    rm -r all_gffs
fi
if [ -d "all_proteomes" ]; then
    rm -r all_proteomes
fi
mkdir all_gffs
mkdir all_proteomes

while read genome; do
    cp /global/scratch/users/pierrj/moryzae_virulence_project/genome_annotation/${genome}/fungap_out/fungap_out/fungap_out.gff3 all_gffs/${genome}_fungap_out.gff3
    cp /global/scratch/users/pierrj/moryzae_virulence_project/genome_annotation/${genome}/fungap_out/fungap_out/fungap_out_prot.faa all_proteomes/${genome}_fungap_out_prot.faa
done < /global/scratch/users/pierrj/moryzae_virulence_project/genome_annotation/strain_names

while read genome; do
    echo $genome
    awk '{print $1}' /global/scratch/users/pierrj/moryzae_virulence_project/genome_annotation/${genome}/fungap_out/fungap_out/fungap_out_prot.faa | grep t1
done < /global/scratch/users/pierrj/moryzae_virulence_project/genome_annotation/strain_names

orthofinder -op -S diamond_ultra_sens -f all_proteomes_corrected -o orthofinder | grep "diamond blastp" > jobqueue

N_NODES=50

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split -a 3 --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%03g" 1 ${N_NODES})
do
    echo $node
    sbatch -p savio2 --ntasks-per-node 24 --job-name=$node.blast --export=ALL,node=$node /global/home/users/pierrj/git/slurm/orthofinder_blast.slurm
done