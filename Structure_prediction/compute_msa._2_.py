import os, sys
import subprocess
import multiprocessing
from joblib import Parallel, delayed

hhblits_binary_path = 'hhblits'
bfd_database_path = '/global/scratch2/groups/fc_kvkallow/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt'
uniclust30_database_path = '/global/scratch2/groups/fc_kvkallow/uniclust30/uniclust30_2018_08/uniclust30_2018_08'

def build_msa(item):

    os.chdir(item)
    if "alphafold.msa2.running" not in os.listdir(".") and "alphafold.msa2.done" not in os.listdir("."):
        with open("alphafold.msa2.running", "w") as output_handle:
            output_handle.write(" ")

        fa = os.getcwd() + "/" + item + ".fasta"
        msa_output_dir = os.path.join(os.getcwd() + "/" + item, 'msas')
        if not os.path.exists(item):
            os.mkdir(item)
        os.chdir(item)
        if not os.path.exists("msas"):
            os.mkdir("msas")
        os.chdir("msas")

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

        for db_path in [uniclust30_database_path]:  #[bfd_database_path, uniclust30_database_path]: #bfd_database_path, uniclust30_database_path
           cmd.append('-d')
           cmd.append(db_path)

        print(' '.join(cmd))
        p = subprocess.Popen(cmd)
        p.wait()

        os.chdir("../..")
        os.remove("alphafold.msa2.running")
        with open("alphafold.msa2.done" , "w") as output_handle:
            output_handle.write("")

    os.chdir("..")

    prefix = sys.argv[1]

if "," in prefix:
  sequence  = []
  for p in prefix.split(","):
    sequence += [ s for s in os.listdir() if p in s ]
else:
  sequence = [ s for s in os.listdir(".") if prefix in s ]

Parallel(n_jobs = 8)(delayed(build_msa)(item) for item in sequence)
