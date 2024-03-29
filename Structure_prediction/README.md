# Structure prediction

## AlphaFold
We will predict or download a structure for a representative sequence in each cluster. DeepMind's alphafold database hosts predicted structures for 70-15, P131 and Y34 strains. If our strains' proteins have good matches in this database, we can skip the prediction. Otherwise, we will predict the structures.  

### Using existing structures  

First, donwload proteomes from Uniprot.  
70-15: https://www.uniprot.org/proteomes/UP000009058  
P131 : https://www.uniprot.org/proteomes/UP000011085  
Y34  : https://www.uniprot.org/proteomes/UP000011086  

We will search the annotation sets of our 70 strains against these three strains'  
```
cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/BLAST
```

Concatnate the proteins from Uniprot  
```
cat uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2022.12.19-21.3* > uniprot.ref.fasta
```

Concate all M. oryzae protein annotation sets from our study  
```
cat ../orthogrouping/all_proteomes_processed/*.faa > Mo.fa
```

Run Blast (v2.7.1+) search between these two concatnated fasta files  
```  
mkdir blastdb  
makeblastdb -in Mo.fa -out blastdb/Mo -dbtype 'prot'  
blastp -query uniprot.ref.fasta -db blastdb/Mo -max_target_seqs 5 -num_threads 52 -evalue 1e-10  
       -max_hsps 1 -outfmt "6 std qlen slen" -out uniprot.ref.against.Mo.blast.out
```

We have one representative sequence from 13387 OGs (one OG, N0.HOG0013128, is discarded as its two members are mostly composed of 'X'), 1256 unassigned genes, 133 phylogenetically misplaced genes, and 8 sequences that are somehow discarded by OrthoFinder. 14784 sequences' structures should be either downloaded or predicted. Based on the BLAST outputs, decide whether we can download the existing structures or need to predict the structures.  
```
python choose_representative.py
```  
  
  
This will generate two files:  
**AF2.list**: These 9447 sequences have predicted structures in the AF2 database with 98% or more sequence identity and 100% coverage. These will be downloaded.  

**Predict.list**: These 5337 structures will be predicted with AlphaFold v2.2.2.  

Download pre-generated structures for M. oryzae from the AlphaFold database. The TAXIDs are 242507 for 70-15, 1143189 for Y34, and 1143193 for P131. This requires 'gsutil'.

```
cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures-DB/AF2  
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-242507-*.tar .  
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-1143189-*.tar .  
gsutil -o "GSUtil:state_dir=/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures/AF2" \
       -m cp gs://public-datasets-deepmind-alphafold/proteomes/proteome-tax_id-1143193-*.tar .  
ls *.tar | while read tar; do tar xf $tar; done
gunzip *.cif.gz
```


We will generate folders to store any outputs associated with sequences in AF2.list. Biopython should be faster than the code below.  
```
cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures-2
less AF2.list | awk '{print $3}' | sort -u | while read seq; do \
     mkdir ${seq} && awk -v seq=$seq -v RS=">" '$1 == seq {print RS $0; exit}' ../BLAST/Mo.fa > ${seq}\/${seq}\.fasta; \
done
```

Convert the cif files into pdb files and store in each directory.  
```
python cif2pdb.py
```

Although we will not predict the structures for these sequences, we still need to generate their sequence profiles for sequence similarity searches at the clustering step. These sequence profiles are converted from the MSAs that come from AlphaFold, so let's run the MSA generation step. 

```
less AF2.list | awk '{print $3}' | cut -d "_" -f 1 | sort -u > prefix.list  
prefix=$(less prefix.list | tr "\n" ",")    
python compute_msa._1_.py ${prefix}    
python compute_msa._2_.py ${prefix}
```       
### Predicting protein structures  

We will predict 5337 sequences in Predict.list. Get each sequence into a new folder to set up for AF2. For this iteration, using Biopython will be much faster instead of what is given here. 
```
cd /global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures  
less Predict.list | awk '{print $2}' | sort -u | while read seq; do \
  mkdir ${seq} && awk -v seq=$seq -v RS=">" '$1 == seq {print RS $0; exit}' ../BLAST/Mo.fa > ${seq}\/${seq}\.fasta; \
done
```

Predict the structures for these with alphafold. Proteins > 800 AA were predicted with GPUs with high memory - e.g. A40. All proteins > 1000 AA were skipped, as they were too large. Moreover, to reduce computing time for the MSA generation, we didn't search against the bfd database. Collect the MSAs first.  

```
less Predict.list | awk '{print $2}' | cut -d "_" -f 1 | sort -u > prefix.list    
prefix=$(less prefix.list | tr "\n" ",")
python compute_msa._1_.py ${prefix}  
python compute_msa._2_.py ${prefix}
```

Alphafold was run with the following commend for each {sequence}. This script also had slight modification to search against pdb_seqres.txt instead of pdb70_databases to use recently available PDB structures as template.
```
python run_alphafold.py --fasta_paths=${sequence} --model_preset=monomer --db_preset=full_dbs --output_dir=. --use_gpu_relax=True --use_precomputed_msas=True
```



## ESMfold  
Finally, we will try to predict the structures with ESMfold which alphafold failed to model. If the plddt is not > 70, get the sequence and run ESMfold. ESMfold was run through Google Colab. Check out:
```
esmfold.ipynb
```
      
