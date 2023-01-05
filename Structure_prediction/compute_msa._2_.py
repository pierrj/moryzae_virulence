import os, sys
import subprocess
import multiprocessing
from joblib import Parallel, delayed

#This script is a modification of AF2 scripts
#Paths to excutables and databases
hhblits_binary_path = 'hhblits'
bfd_database_path = '/global/scratch2/groups/fc_kvkallow/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt'
uniclust30_database_path = '/global/scratch2/groups/fc_kvkallow/uniclust30/uniclust30_2018_08/uniclust30_2018_08'

def build_msa(item):

    os.chdir(item)
    
    #If MSA generation is not running in the current working directory
    if "alphafold.msa2.running" not in os.listdir(".") and "alphafold.msa2.done" not in os.listdir("."):

        #Mark running
        with open("alphafold.msa2.running", "w") as o: o.write(" ")
        
        #Get the sequence
        fa = os.getcwd() + "/" + item + ".fasta"

        #Generate the folders
        if not os.path.exists(item): os.mkdir(item)
        if not os.path.exists(f"{item}/msas"): os.mkdir(f"{item}/msas")
        os.chdir(f"{item}/msas")

        #Run HHblits
        cmd = [
          hhblits_binary_path,
          '-i', fa,
          '-cpu', str(4),
          '-oa3m', 'bfd_uniclust_hits.a3m',
          '-o', 'hhblits.out',
          '-n', str(5),
          '-e', str(0.001),
          '-maxseq', str(1000000),
          '-realign_max', str(100000),
          '-maxfilt', str(100000),
          '-min_prefilter_hits', str(1000)]

        #Run search 
        for db_path in [uniclust30_database_path]:  #Initially [bfd_database_path, uniclust30_database_path]: 
                                                    #Skip bfd to save computation time
           cmd.append('-d')
           cmd.append(db_path)

        print(' '.join(cmd))
        p = subprocess.Popen(cmd)
        p.wait()

        os.chdir("../..")

        #Mark done 
        with open("alphafold.msa2.done" , "w") as o: o.write("")
        os.remove("alphafold.msa2.running")

    os.chdir("..")

    prefix = sys.argv[1]

if "," in prefix:
  sequence  = []
  for p in prefix.split(","):
    sequence += [ s for s in os.listdir() if p in s ]
else:
  sequence = [ s for s in os.listdir(".") if prefix in s ]

Parallel(n_jobs = 8)(delayed(build_msa)(item) for item in sequence)
