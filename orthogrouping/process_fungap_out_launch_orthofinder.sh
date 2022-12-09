#!/bin/bash
PROJECT_DIR=/global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping

cd $PROJECT_DIR

conda activate /global/scratch/users/pierrj/conda_envs/orthofinder

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

cp -r /global/scratch/users/pierrj/moryzae_virulence_project/genome_annotation/Assembly .

/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python /global/scratch/users/pierrj/moryzae_virulence/orthogrouping/process_fungap_out_and_assemblies.py \
                                                                                                                                        Assembly \
                                                                                                                                        Assembly_renamed \
                                                                                                                                        all_gffs \
                                                                                                                                        all_gffs_processed \
                                                                                                                                        all_proteomes \
                                                                                                                                        all_proteomes_processed


## copy over mgrisae proteome
cp /global/scratch/users/pierrj/PAV_SV/PAV/re_gladieux_proteomes_fungap/all_proteomes_corrected/NI907_fungap_out_prot_filtered.faa all_proteomes_processed/

orthofinder -op -S diamond_ultra_sens -f all_proteomes_processed -n out -o orthofinder_out | grep "diamond blastp" > jobqueue

N_NODES=4

mv jobqueue jobqueue_old

shuf jobqueue_old > jobqueue

split -a 3 --number=l/${N_NODES} --numeric-suffixes=1 jobqueue jobqueue_

for node in $(seq -f "%03g" 1 ${N_NODES})
do
    echo $node
    sbatch -p savio4_htc -A co_minium --ntasks-per-node 56 --qos=minium_htc4_normal --job-name=$node.blast --export=ALL,node=$node /global/scratch/users/pierrj/moryzae_virulence/orthogrouping/orthofinder_blast.slurm
done

## once finished

sbatch -p savio4_htc -A co_minium --ntasks-per-node 56 --qos=minium_htc4_normal /global/scratch/users/pierrj/moryzae_virulence/orthogrouping/orthofinder_orthogrouping.slurm


# now generate MSAs and make the tree using fasttree
sbatch -p savio4_htc -A co_minium --ntasks-per-node 56 --qos=minium_htc4_normal /global/scratch/users/pierrj/moryzae_virulence/orthogrouping/generate_sco_msa.slurm
sbatch -p savio4_htc -A co_minium --ntasks-per-node 56 --qos=minium_htc4_normal /global/scratch/users/pierrj/moryzae_virulence/orthogrouping/fasttree_scos.slurm

## extract species tree here:
## and root in ITOL and output as nwk

