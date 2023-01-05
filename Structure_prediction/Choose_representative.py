from Bio import SeqIO
from collections import Counter

def orthogroups():
  """Go through the OrthoFinder outputs and save the membership of each sequences
    Also, pick a representative sequence per OG to predict the structure for    """

  OG = "../orthogrouping/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_6/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
  UA = "../orthogrouping/orthofinder_out/Results_out/WorkingDirectory/OrthoFinder/Results_out_1/Orthogroups/Orthogroups_UnassignedGenes.tsv"
  fasta_handle = SeqIO.to_dict(SeqIO.parse("Mo.fa", "fasta"))

  memberships = {}
  refseq      = {}

  with open(OG, "r") as f:
    next(f)

    for line in f:
      members = [ x.replace(",", "") for x in line.split()[3:]]
      lengths = [ len(str(fasta_handle[m].seq)) for m in members ]
      o       = line.split()[0]

      #sequence ID = OG ID
      for m in members:
        memberships[m] = o

      #Because these are orthologous from closely related species
      #The gene lengths are supposedly the same
      #In case they differ, follow the majority rule
      cts = Counter(lengths)
      common = [ n for n,c in cts.most_common() if c == cts.most_common(1)[0][1] ]

      #if there is one winner, get the first sequence of that length
      if len(common) == 1:
        refseq[ o ] = [ members[x] for x in range(len(lengths)) if lengths[x] == cts.most_common(1)[0][0] ][0]

      #if on par, get the longest
      else:
        longest = [ members[x] for x in range(len(lengths)) if lengths[x] == max(common) ]

        if len(longest) == 1:
          refseq[ o ] = longest[0]
        else:
          refseq[ o ] = [ x for x in longest if "gene" not in x ][0] #remove outgroup if possible

  #Singleton OGs
  for line in open(UA, "r"):
    memberships[ line.split()[1] ] = line.split()[1]
    refseq[ line.split()[1] ]      = line.split()[1]

  return memberships, refseq

def check_blast(memberships, refseq):

  BLAST = "uniprot.ref.against.Mo.blast.out"

  AF2     = {}
  Predict = {}

  #If there are existing predicted structures from the AF2 DB,
  #with >= 98% sequence identity, and 100% coverage,
  #use that structure. Otherwise, predict one per OG

  for line in open(BLAST, "r"):
    #There are 19 sequences removed from phylogeny-based OGs by OrthoFinder
    #Igore them, since there are other informatic orthologs in the same OGs

    try:
      OG = memberships[ line.split()[1] ]
    except Exception:
      OG = "s"


    if not OG == "s":
      if float(line.split()[3]) >= 98: # sequence identity
        if line.split()[-1] == line.split()[-2]: # query/hit length

          if OG not in AF2:
            AF2[ OG ] = [line.split()[0].split("|")[1], line.split()[1]]

        else:
          if OG not in AF2:
            Predict[ OG ] = refseq[ OG ]

      else:
        if OG not in AF2:
          Predict[ OG ] = refseq[ OG ]

  with open("AF2.list", "w") as o:
     for a in AF2:
       o.write(f"{a}\t{AF2[a][0]}\t{AF2[a][1]}\n")

  with open("Predict.list", "w") as o:
     for a in Predict:
       o.write(f"{a}\t{Predict[a]}\n")

memberships,refseq = orthogroups()
check_blast(memberships, refseq)
