#!/bin/bash

cd /global/scratch/users/pierrj/moryzae_virulence_project/genome_annotation

cp -r /global/scratch/users/pierrj/fungap_runs/gladieux_all/template_run/ .

while read genome; do
    echo $genome
    if [ -d "${genome}" ]; then
        rm -r ${genome}
    fi
    mkdir ${genome}
    cd ${genome}
        cp -r ../template_run/fungap_out/ .
    cd ..
done < strain_names

while read genome; do
    sbatch --job-name=${genome}_run_fungap --export=genome=$genome /global/scratch/users/pierrj/moryzae_virulence/genome_annotation/run_fungap.slurm
done < strain_names

genome=KVK1
sbatch -p savio4_htc -A co_minium --qos=minium_htc4_normal --job-name=${genome}_run_fungap --export=genome=$genome /global/scratch/users/pierrj/moryzae_virulence/genome_annotation/run_fungap.slurm

## to relaunch failed jobs due to busco download error

squeue -u pierrj --format="%.100j" | tail -n +2 | awk '{print substr($1, 0,length($1)-11)}' > running_jobs

while read genome; do
    if grep -Fxq "$genome" running_jobs
    then
        echo "$genome is already running"
    else
    sbatch --job-name=${genome}_run_fungap --export=genome=$genome /global/scratch/users/pierrj/moryzae_virulence/genome_annotation/run_fungap.slurm
    fi
done < strain_names