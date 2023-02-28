#!/bin/bash
#SBATCH --job-name=pav_validation_parallel
#SBATCH --partition=savio
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out

PROJECT_DIR=/global/scratch/users/pierrj/moryzae_virulence_project/pav_validation
cd $PROJECT_DIR

CODE_DIRECTORY=/global/scratch/users/pierrj/moryzae_virulence_code/pav_validation/
LOST_GENOME_DIR=/global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping/Assembly_renamed # where the assemblies are
LOST_OG_DIR=/global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_1/Orthogroup_Sequences # where the OG_protein fastas are
ORTHOGROUPS_TSV=/global/scratch/users/pierrj/moryzae_virulence_project/orthogrouping/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_1/Orthogroups/Orthogroups.tsv # where the orthogroup tsv is
GENOME_SUFFIX=_final.scaffolds.fasta
E_VALUE=1e-10
PIDENT=55
QUERY_COV=55
HIT_COUNT=2
N_NODES=40
OUTPUT_FILE=${PROJECT_DIR}/pav_table
BLAST_DB=${PROJECT_DIR}/all_ogs_seqs.fasta
ABSENCES_FILE=absences_to_validate.txt
JOBQUEUE_FILE=jobqueue_pav_validation

source activate /global/scratch/users/pierrj/conda_envs/orthofinder

## generate orthogroups to validate
Rscript --vanilla $CODE_DIRECTORY/print_pavs_to_validate.r $ORTHOGROUPS_TSV $GENOME_SUFFIX absences_to_validate.txt

# concatenate all orthogroup sequences together into one file and make blast db
/global/scratch/users/pierrj/conda_envs/orthofinder/bin/python $CODE_DIRECTORY/make_single_file_from_og_dir.py ${LOST_OG_DIR} ${BLAST_DB}


makeblastdb -in ${BLAST_DB} -dbtype prot



if [ -d "${PROJECT_DIR}/pav_validation" ]; then
    rm -r ${PROJECT_DIR}/pav_validation
fi

mkdir ${PROJECT_DIR}/pav_validation

cd ${PROJECT_DIR}/pav_validation

if [ -f "jobqueue" ]; then
    rm jobqueue
fi

# validate missing ogs using tblastn and blastp
while read -r LOST_GENOME LOST_OG; do
    echo $CODE_DIRECTORY/tblastn_validation.sh -l ${LOST_OG_DIR}/${LOST_OG}.fa -g ${LOST_GENOME_DIR}/${LOST_GENOME} \
        -e ${E_VALUE} -p ${PIDENT} -q ${QUERY_COV} -c ${HIT_COUNT} -d ${BLAST_DB} -h ${CODE_DIRECTORY} >> jobqueue
done < ${ABSENCES_FILE}

if [ -f "${OUTPUT_FILE}" ]; then
    rm ${OUTPUT_FILE}
fi

if [ -f $JOBQUEUE_FILE ]; then
    rm $JOBQUEUE_FILE
fi

# parallelization stuff
split --number=l/${N_NODES} --numeric-suffixes=1 $JOBQUEUE_FILE ${JOBQUEUE_FILE}_

if [ -f slurm_ids ]; then
    rm slurm_ids
fi

for node in $(seq -f "%02g" 1 ${N_NODES})
do
    sbatch --parsable -p savio2 --ntasks-per-node 24 --job-name=$node.tblastn_validation --export=node=$node,OUTPUT_FILE=$OUTPUT_FILE,JOBQUEUE_FILE=$JOBQUEUE_FILE $CODE_DIRECTORY/gnu_parallel_multinode.slurm
done

if [ -f slurm_ids ]; then
    # wait for all jobs to finish
    sbatch -W --dependency=afterok:$(cat slurm_ids | head -c -1 | tr '\n' ':') $CODE_DIRECTORY/dummy_job.slurm
    wait
fi

## after everything is done, remove verbose outputs
awk -v OFS='\t' '{ if ($3 == "yes") {print $1, $2, $3} else {print $1, $2, "no"}}' ${OUTPUT_FILE} > ${OUTPUT_FILE}.simplified