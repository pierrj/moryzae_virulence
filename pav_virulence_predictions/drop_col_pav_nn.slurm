#!/bin/bash
#SBATCH --job-name=drop_col_pav_nn
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --account=co_minium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out


cd /global/scratch/users/pierrj/moryzae_virulence_project/pav_virulence_prediction

source activate /global/scratch/users/pierrj/conda_envs/pytorch/

jobqueue_file=jobqueue
OUTPUT_FILE=pav_drop_outs_nn_f1s.txt

if [ -f $OUTPUT_FILE ]; then
    rm $OUTPUT_FILE
fi

if [ -f $jobqueue_file ]; then
    rm $jobqueue_file
fi

N_CPUS=$SLURM_NTASKS

split --number=l/${N_CPUS} --numeric-suffixes=1 og_list og_list_

for i in $(seq -f "%02g" $N_CPUS); do
    echo "/global/scratch/users/pierrj/conda_envs/pytorch/bin/python /global/scratch/users/pierrj/moryzae_virulence_code/pav_virulence_predictions/drop_col_nn_virulence_predictions.py /global/scratch/users/pierrj/moryzae_virulence_project/pav_virulence_prediction og_list_${i}" >> $jobqueue_file
done

parallel -j ${SLURM_NTASKS} < ${jobqueue_file} >> ${OUTPUT_FILE}