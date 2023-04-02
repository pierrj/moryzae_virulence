# Structure prediction

## AlphaFold
We will predict or download a structure for a representative sequence in each cluster. DeepMind's alphafold database hosts predicted structures for 70-15, P131 and Y34 strains. If our strains' proteins have good matches in this database, we can skip the prediction. Otherwise, we will predict the structures.  

First, donwload proteomes from Uniprot.  
70-15: https://www.uniprot.org/proteomes/UP000009058  
P131 : https://www.uniprot.org/proteomes/UP000011085  
Y34  : https://www.uniprot.org/proteomes/UP000011086  

We will search our proteomes against these three strains'  
`cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/BLAST`  

Concatnate the proteins from Uniprot  
`cat uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2022.12.19-21.3* > uniprot.ref.fasta`  

Concate all M. oryzae protein annotation sets from our study  
`cat ../orthogrouping/all_proteomes_processed/*.faa > Mo.fa`  

Run Blast search between these two concatnated fasta files  
<code>
module load blast #v2.7.1+   
mkdir blastdb  
makeblastdb -in Mo.fa -out blastdb/Mo -dbtype 'prot'  
blastp -query uniprot.ref.fasta -db blastdb/Mo -max_target_seqs 5 -num_threads 52 -evalue 1e-10  
       -max_hsps 1 -outfmt "6 std qlen slen" -out uniprot.ref.against.Mo.blast.out
</code>

Based on the BLAST outputs, decide whether we can download the existing structures or need to predict the structures.  
`python choose_representative.py`       

This will generate two files:  
[1] AF2.list: These 9446 sequences have predicted structures in the AF2 database with 98% or more sequence identity and 100% coverage. These will be downloaded. 
[2] Predict.list: These 5299 structures will be predicted with AlphaFold v2.2.2.

Later, we realized a few sequences in 'Predict.list' contained unknown sequence 'X'. This won't allow the relaxation step of Alphafold to run. These sequences had to be replaced as below:

BR0019_14_02091T0  -> replaced with gene_8576_NI907 (different length, no X)  
BR0019_1_00155T0   -> replaced with gene_5088_NI907 (different length, no X)  
BR0019_32_03871T0  -> replaced with gene_6814_NI907 (different length, no X)  
CD0073_61_05782T0  -> replaced with CH0452_89_07103T0 (same length, no X)  
CH0063_592_11157T0 -> replaced with CH0072_312_10461T0 (same length, no X)  
CH0333_1_00001T0   -> Removed a single 'X' at the very end of the sequences  
       
Get each sequence into a new folder to set up for AF2. For this iteration, using Biopython will be much faster instead of what is given here.
`cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures`  
<code>less ../BLAST/Predict.list | awk '{print $2}' | sort -u | while read seq; do \
&emsp;&emsp;mkdir ${seq} && awk -v seq=$seq -v RS=">" '$1 == seq {print RS $0; exit}' ../BLAST/Mo.fa > ${seq}\/${seq}\.fasta; \
done</code>  

Predict the structures for these with alphafold. Proteins > 800 AA were predicted with GPUs with high memory - e.g. A40. All proteins > 1000 AA were skipped, as they were too large. Moreover, to reduce computing time for the MSA generation, we didn't search against the bfd database.   

`conda activate alphafold`  
`module load cuda/11.2`  
`export PATH=$PATH:/global/scratch/users/skyungyong/Software/alphafold-msa/hh-suite-3.3.0/build/bin`   
`export PATH=$PATH:/global/scratch/users/skyungyong/Software/alphafold-msa/hh-suite-3.3.0/build/scripts`  

`ls -d * | grep -e "T0" -e "gene" | cut -d "_" -f 1 | sort -u > prefix.list`  
`prefix=$(less prefix.list | tr "\n" ",")`  

These two scripts to collect MSAs. These are simply modified scripts of AF2
These may or may not run independently of the AF2 packages (may need some more modificiation)
`python compute_msa._1_.py ${prefix}`  
`python compute_msa._2_.py ${prefix}`

Alphafold was run with the following commend for each {sequence}
This script also had slight modification to search against pdb_seqres.txt instead of pdb70_databases
To use recently available PDB structures as template
python run_alphafold.py --fasta_paths=${sequence} --model_preset=monomer --db_preset=full_dbs --output_dir=. \
                        --use_gpu_relax=True --use_precomputed_msas=True 

Download pre-generated structures for M. oryzae from the AlphaFold database
The TAXIDs are 242507 for 70-15, 1143189 for Y34, and 1143193 for P131
This requires 'gsutil'

<code>cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures-DB/AF2  
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-242507-*.tar .  
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-1143189-*.tar .  
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-1143193-*.tar .<\code>

ls *.tar | while read tar; do tar xf $tar; done
gunzip *.cif.gz

#I also need to collect MSAs for the sequences, the structures of which are available in the AF2 DB. 
#Generate the same folder architecture and collect MSAs
#Rather use Biopython, since iterating for each sequence takes a while for the large fasta file....
cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures-2
less ../BLAST/AF2.list | awk '{print $3}' | sort -u | while read seq; do \
       mkdir ${seq} && awk -v seq=$seq -v RS=">" '$1 == seq {print RS $0; exit}' ../BLAST/Mo.fa > ${seq}\/${seq}\.fasta; \
done

#Generate MSAs
conda activate alphafold
export PATH=$PATH:/global/scratch/users/skyungyong/Software/alphafold-msa/hh-suite-3.3.0/build/bin
export PATH=$PATH:/global/scratch/users/skyungyong/Software/alphafold-msa/hh-suite-3.3.0/build/scripts

ls -d * | grep -e "T0" -e "gene" | cut -d "_" -f 1 | sort -u > prefix.list
prefix=$(less prefix.list | tr "\n" ",")

python /global/scratch/users/skyungyong/Software/alphafold-no-change_102522/alphafold/compute_msa._1_.py ${prefix}
python /global/scratch/users/skyungyong/Software/alphafold-no-change_102522/alphafold/compute_msa._2_.py ${prefix}

##ESMfold 
Finally, for structures - regardless of they are predicted or downloaded - without plddt > 70,
#check if omegafold can predict better structures.
#for those sequences, omegafold was run as:
omegafold --model 1 {fasta} {output}
