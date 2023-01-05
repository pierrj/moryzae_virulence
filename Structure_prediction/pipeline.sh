#!/bin/bash

#Donwload proteomes from Uniprot
#70-15: https://www.uniprot.org/proteomes/UP000009058
#P131 : https://www.uniprot.org/proteomes/UP000011085
#Y34  : https://www.uniprot.org/proteomes/UP000011086

#Concatnate the proteins from Uniprot
cat uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2022.12.19-21.3* > uniprot.ref.fasta

#Concate all M. oryzae protein annotation sets from our study
cat ../orthogrouping/all_proteomes_processed/*.faa > Mo.fa

#Run Blast search against these two concatnated fasta files
module load blast #v2.7.1+
mkdir blastdb

makeblastdb -in Mo.fa -out blastdb/Mo -dbtype 'prot'
blastp -query uniprot.ref.fasta -db blastdb/Mo -max_target_seqs 5 -num_threads 52 -evalue 1e-10 \
       -max_hsps 1 -outfmt "6 std qlen slen" -out uniprot.ref.against.Mo.blast.out

#Based on the BLAST outputs, decide whether I need to predict the structures
#This will generate 'AF2.list' - there are already structures from the AF2 database
#and 'Predict.list' - I need to predict the structures
python choose_representative.py

#Get each sequence into a new folder to set up for AF2
#Rather use Biopython, since iterating for each sequence takes a while for the large fasta file....
cd ../Structure
less ../BLAST/Predict.list | awk '{print $2}' | sort -u | while read seq; do \
       mkdir ${seq} && awk -v seq=$seq -v RS=">" '$1 == seq {print RS $0; exit}' ../BLAST/Mo.fa > ${seq}\/${seq}\.fasta; \
done

#Predict the structures for these with alphafold
#Proteins > 800AA were predicted with GPUs with high memory (A40 for instance)
#Proteins > 1000AA were skipped, as they were too large
#I didn't search against the bfd database to save MSA computation time. 

conda activate alphafold
module load cuda/11.2
export PATH=$PATH:/global/scratch/users/skyungyong/Software/alphafold-msa/hh-suite-3.3.0/build/bin
export PATH=$PATH:/global/scratch/users/skyungyong/Software/alphafold-msa/hh-suite-3.3.0/build/scripts
prefix=$(less prefix.list | tr "\n" ",")

#These two scripts to collect MSAs. These are simply modified scripts of AF2
#These may or may not run independently of the AF2 packages (may need some more modificiation)
python /global/scratch/users/skyungyong/Software/alphafold-no-change_102522/alphafold/compute_msa._1_.py ${prefix}
python /global/scratch/users/skyungyong/Software/alphafold-no-change_102522/alphafold/compute_msa._2_.py ${prefix}

#Alphafold was run with the following commend
#for each {sequence}
#This script also had slight modification to search against pdb_seqres.txt instead of pdb70_databases
#To use recently available PDB structures as template
python run_alphafold.py --fasta_paths=${sequence} --model_preset=monomer --db_preset=full_dbs --output_dir=. \
                        --use_gpu_relax=True --use_precomputed_msas=True 

#Download pre-generated structures for M. oryzae from the AlphaFold database
#The TAXIDs are 242507 for 70-15, 1143189 for Y34, and 1143193 for P131
#This requires 'gsutil'

cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-242507-*.tar .
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-1143189-*.tar .
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-1143193-*.tar .

ls *.tar | while read tar; do tar xf $tar; done
gunzip *.cif.gz

