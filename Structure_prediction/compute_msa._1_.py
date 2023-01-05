import os, sys
import subprocess
import multiprocessing
from joblib import Parallel, delayed

# Paths to excutables and databases
data_dir = "/global/scratch/users/skyungyong/Software/alphafold/Database"
uniref90_database_path = '/global/scratch/users/skyungyong/Software/alphafold/Database/uniref90/uniref90.fungalDB.added.fasta'
mgnify_database_path = '/global/scratch/users/skyungyong/Software/alphafold/Database/mgnify/mgy_clusters.fa'

# Takes a single fasta file as input. So run this function in the multi-threaded.

def compute_msa(input_fasta_path,
                msa_output_dir): # Run jackhmmer and hhsearch and return outputs

    with open(input_fasta_path) as f:
      input_fasta_str = f.read()

    input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)
    if len(input_seqs) != 1:
      raise ValueError(
          'More than one input sequence found in ' + input_fasta_path)
    input_sequence = input_seqs[0]
    input_description = input_descs[0]
    num_res = len(input_sequence)

    cmd_flags = [
          # Don't pollute stdout with Jackhmmer output.
          '-o', '/dev/null',
          '--noali',
          '--F1', str(0.0005),
          '--F2', str(0.00005),
          '--F3', str(0.0000005),
          '--incE', str(0.0001),
          # Report only sequences with E-values <= x in per-sequence output.
          '-E', str(0.0001),
          '--cpu', str(1),
          '-N', str(1),
      ]

    with utils.tmpdir_manager(base_dir='/tmp') as query_tmp_dir:
      sto_path = os.path.join(query_tmp_dir, 'output.uniref.sto')

# Run jackhmmer against uniref
      cmd = [jackhmmer_binary_path] + cmd_flags + ['-A', sto_path] + [input_fasta_path, uniref90_database_path]

      logging.info('Launching subprocess "%s"', ' '.join(cmd))
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      with utils.timing(
          f'Jackhmmer ({os.path.basename(uniref90_database_path)}) query'):
        _, stderr = process.communicate()
        retcode = process.wait()

      if retcode:
        raise RuntimeError(
            'Jackhmmer failed\nstderr:\n%s\n' % stderr.decode('utf-8'))

    # Get e-values for each target name
      tbl = ''
      with open(sto_path) as f:
       sto = f.read()

      jackhmmer_uniref90_result = dict(
        sto=sto,
        tbl=tbl,
        stderr=stderr,
        n_iter=1,
        e_value=0.0001)

      uniref90_out_path = os.path.join(msa_output_dir, 'uniref90_hits.sto')
      with open(uniref90_out_path, 'w') as f:
        f.write(jackhmmer_uniref90_result['sto'])

      uniref90_msa_as_a3m = parsers.convert_stockholm_to_a3m(
        jackhmmer_uniref90_result['sto'], max_sequences=10000)

# Run jackhmmer against mgnify

    with utils.tmpdir_manager(base_dir='/tmp') as query_tmp_dir:
      sto_path = os.path.join(query_tmp_dir, 'output.mgnify.sto')
      cmd = [jackhmmer_binary_path] + cmd_flags + ['-A', sto_path] + [input_fasta_path, mgnify_database_path]

      logging.info('Launching subprocess "%s"', ' '.join(cmd))
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      with utils.timing(
          f'Jackhmmer ({os.path.basename(mgnify_database_path)}) query'):
        _, stderr = process.communicate()
        retcode = process.wait()

      if retcode:
        raise RuntimeError(
            'Jackhmmer failed\nstderr:\n%s\n' % stderr.decode('utf-8'))


    # Get e-values for each target name
      tbl = ''
      with open(sto_path) as f:
        sto = f.read()

      jackhmmer_mgnify_result = dict(
        sto=sto,
        tbl=tbl,
        stderr=stderr,
        n_iter=1,
        e_value=0.0001)

      mgnify_out_path = os.path.join(msa_output_dir, 'mgnify_hits.sto')
      with open(mgnify_out_path, 'w') as f:
        f.write(jackhmmer_mgnify_result['sto'])

def msa_parallel(item):
    os.chdir(item)
    fa = item + ".fasta"
    outDir = os.getcwd()

    #If MSA generation is not running in the current working directory
    if "alphafold.msa.running" not in os.listdir(".") and "alphafold.msa.done" not in os.listdir("."):
        
        #Mark running to avoid redundant runs
        with open("alphafold.msa.running", "w") as o: o.write(" ")

        #The outputs will be stored in {sequence}/{sequence}/msas/ folder
        output_dir = os.path.join(outDir, item)

        #Generate folders
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        msa_output_dir = os.path.join(output_dir, 'msas')
        if not os.path.exists(msa_output_dir): os.makedirs(msa_output_dir)

        #Compute MSA
        compute_msa(fa, msa_output_dir)
        
        #Mark done
        with open("alphafold.msa.done" , "w") as o: o.write(" ")
        os.remove("alphafold.msa.running")
        
    os.chdir("..")


prefix = sys.argv[1]

if "," in prefix:
  sequences = []
  for p in prefix.split(","):
    sequences += [ s for s in os.listdir() if p in s ]
else:
  sequences = [ s for s in os.listdir(".") if prefix in s ]

Parallel(n_jobs = multiprocessing.cpu_count())(delayed(msa_parallel)(item) for item in sequences)
