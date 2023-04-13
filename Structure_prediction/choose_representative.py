from Bio import SeqIO
from collections import Counter
import logging
import os


def remove_single_x(sequence):
    """
    Remove a single 'X' from the start or end of a sequence, if present.
    """
    if sequence.startswith('X'):
        return sequence[1:]
    elif sequence.endswith('X'):
        return sequence[:-1]
    else:
        return sequence

def orthogroups():
  """Go through the OrthoFinder outputs and save the membership of each sequences
    Also, pick a representative sequence per OG to predict the structure for    """

  # File locations
  n0_orthogroup_path    = "../orthogrouping/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_6/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
  misplaced_genes_dir   = '../orthogrouping/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_6/Phylogenetically_Misplaced_Genes/'
  unassigned_genes_path = "../orthogrouping/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_1/Orthogroups/Orthogroups_UnassignedGenes.tsv"
  fasta_path            = "Mo.fa"

  # Concatnated fasta file
  fasta_handle = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

  memberships              = {}
  representative_seqs      = {}

  # Check orthogroup output to pick the representative sequences
  with open(n0_orthogroup_path, "r") as f:
    # Skip the first line
    next(f)

    for line in f:
      # Get all the members within that OG
      members = [ x.replace(",", "") for x in line.split()[3:]]

      # Move outgroup members (gene_*) to the end of the list
      ingroup_members  = [m for m in members if m.startswith('gene')]
      outgroup_members = [m for m in members if not m.startswith('gene')]
      members = ingroup_members + outgroup_members

      # Get length and orthogroup assignments
      lengths       = [ len(str(fasta_handle[member].seq)) for member in members ]
      orthogroup_id = line.split()[0] # orthogroups

      # store membership information for each sequence
      for member in members:
        memberships[member] = orthogroup_id

      # Because these are orthologous from closely related species
      # The gene lengths are supposedly the same in most cases
      # In case they differ, follow the majority rule
      cts = Counter(lengths)

      # Sort from most common occurence to the least
      # If on par, longer sequences are prioritized
      sorted_cts = sorted(cts.items(), key=lambda x: (-x[1], -x[0]))

      confirmed = False
      for item in sorted_cts:
        if not confirmed:
          target_length, occurrence = item
          target_members = [members[i] for i, length in enumerate(lengths) if length == target_length]

          for seq_id in target_members:
            sequence = str(fasta_handle[seq_id].seq)
            if 'X' in sequence:
              continue

            representative_seqs[orthogroup_id] = seq_id
            confirmed = True
            break

        if not confirmed:
          # None of the sequences without 'X' were suitable, so check for a single 'X'
          for seq_id in target_members:
            sequence = str(fasta_handle[seq_id].seq)
            if sequence.count('X') == 1:
                corrected_seq = remove_single_x(sequence)
                fasta_handle[seq_id].seq = corrected_seq
                representative_seqs[orthogroup_id] = seq_id
                confirmed = True
                break

      if not confirmed:
        # All sequences contain 'X', so discard this orthogroup
        logging.warning(f'In orthogroup {orthogroup_id}, all members contain X in their sequences. Ignore this orthogroup')
        for seq_id in members:
            logging.warning(f'{seq_id}  {fasta_handle[seq_id].seq}   {len(str(fasta_handle[seq_id].seq))}')


  # Misplaced genes
  with os.scandir(misplaced_genes_dir) as dir_entries:
    for entry in dir_entries:
        if entry.name.endswith('txt') and entry.is_file():
            with open(entry.path, 'r') as f:
                for line in f:
                    if 'T0' in line or 'gene' in line:
                        seqID = line.rstrip('\n')
                        memberships[seqID] = 's'
                        representative_seqs[seqID] = 's'

  # Singleton OGs
  with open(unassigned_genes_path, 'r') as f:
    next(f)  # skip header line
    for line in f:
        seqID = line.split()[1]
        memberships[seqID] = 's'
        representative_seqs[seqID] = 's'

  return memberships, representative_seqs, fasta_handle

def check_blast(memberships, representative_seqs, fasta_handle):
    # blast search output between our M. oryzae outputs and uniprot
    BLAST = "uniprot.ref.against.Mo.blast.out"

    AF2 = {}
    Predict = {}

    # Store sequences that are missing from the dictionary as singletons
    for line in open(BLAST, "r"):
        query = line.split()[1]
        try:
            OG = memberships[query]

        except KeyError:
            logging.warning(f'Key {query} is missing. Assigning as singleton')
            representative_seqs[query] = 's'
            OG = 's'

        if OG == 's':
            OG = query

        # Use the best match as representative
        identity = float(line.split()[3])
        query_len = int(line.split()[-2])
        hit_len = int(line.split()[-1])

        # Use stringent criteria to use the downloaded alphafold structures
        if identity >= 98 and query_len == hit_len:
            if OG not in AF2:
                AF2[OG] = [line.split()[0].split("|")[1], query, identity, query_len]
            else:
                if OG.startswith('N0'):
                    if query_len > AF2[OG][3]:
                        AF2[OG] = [line.split()[0].split("|")[1], query, identity, query_len]
                    elif query_len == AF2[OG][3] and identity > AF2[OG][2]:
                        AF2[OG] = [line.split()[0].split("|")[1], query, identity, query_len]
                elif identity > AF2[OG][2]:
                    AF2[OG] = [line.split()[0].split("|")[1], query, identity, query_len]

    # Store sequences that were not in AF2 as predictions
    for item in representative_seqs:
        if item not in AF2:
            Predict[item] = representative_seqs[item]

    # Write the sequences to files
    with open("AF2.list", "w") as o:
        for OG, info in AF2.items():
            rep_seq = info[0]
            seqID = info[1]
            sequence = fasta_handle[seqID].seq
            o.write(f"{OG}\t{rep_seq}\t{seqID}\t{sequence}\n")

    with open("Predict.list", "w") as o:
        for OG, seqID in Predict.items():
            rep_seq = 'NA'
            if OG.startswith('N0'):
                sequence = fasta_handle[seqID].seq
            else:
                sequence = fasta_handle[OG].seq
                seqID = OG
                OG    = 's'

            o.write(f"{OG}\t{rep_seq}\t{seqID}\t{sequence}\n")

memberships,refseq, fasta_handle = orthogroups()
check_blast(memberships, refseq, fasta_handle)
