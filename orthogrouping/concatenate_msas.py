from Bio import SeqIO
import sys
import os
from os.path import join
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

msa_dir = sys.argv[1]
output = sys.argv[2]

genomes_dict = {}

strain_names = [
'BR0019_fungap_out_prot',
'CD0073_fungap_out_prot',
'CH0043_fungap_out_prot',
'CH0052_fungap_out_prot',
'CH0063_fungap_out_prot',
'CH0072_fungap_out_prot',
'CH0092_fungap_out_prot',
'CH0110_fungap_out_prot',
'CH0333_fungap_out_prot',
'CH0452_fungap_out_prot',
'CH0461_fungap_out_prot',
'CH0532_fungap_out_prot',
'CH0533_fungap_out_prot',
'CH0539_fungap_out_prot',
'CH0549_fungap_out_prot',
'CH0561_fungap_out_prot',
'CH0562_fungap_out_prot',
'CH0565_fungap_out_prot',
'CH0571_fungap_out_prot',
'CH0718_fungap_out_prot',
'CH1008_fungap_out_prot',
'CH1009_fungap_out_prot',
'CH1016_fungap_out_prot',
'CH1065_fungap_out_prot',
'CH1076_fungap_out_prot',
'CH1079_fungap_out_prot',
'CH1083_fungap_out_prot',
'CH1103_fungap_out_prot',
'CH1120_fungap_out_prot',
'CH1121_fungap_out_prot',
'CH1150_fungap_out_prot',
'CH1164_fungap_out_prot',
'CL0026_fungap_out_prot',
'FR0013_fungap_out_prot',
'GY0011_fungap_out_prot',
'HN0001_fungap_out_prot',
'IN0017_fungap_out_prot',
'IN0054_fungap_out_prot',
'IN0072_fungap_out_prot',
'IN0076_fungap_out_prot',
'IN0082_fungap_out_prot',
'IN0092_fungap_out_prot',
'IN0094_fungap_out_prot',
'IN0115_fungap_out_prot',
'IN0116_fungap_out_prot',
'LA0006_fungap_out_prot',
'LA0009_fungap_out_prot',
'LA0012_fungap_out_prot',
'LA0019_fungap_out_prot',
'LA0023_fungap_out_prot',
'MC0016_fungap_out_prot',
'MD0116_fungap_out_prot',
'MD0929_fungap_out_prot',
'NP0058_fungap_out_prot',
'NP0070_fungap_out_prot',
'PH0019_fungap_out_prot',
'PH0103_fungap_out_prot',
'PH0158_fungap_out_prot',
'PR0009_fungap_out_prot',
'SP0005_fungap_out_prot',
'SP0006_fungap_out_prot',
'TN0016_fungap_out_prot',
'TN0070_fungap_out_prot',
'US0032_fungap_out_prot',
'US0041_fungap_out_prot',
'US0090_fungap_out_prot',
'US0098_fungap_out_prot',
'VT0003_fungap_out_prot',
'VT0027_fungap_out_prot',
'VT0030_fungap_out_prot']

for strain in strain_names:
    if strain == "IN0092":
        continue
    genomes_dict[strain] = SeqRecord(Seq(""), id=strain)

genomes_dict['NI907'] = SeqRecord(Seq(""), id='NI907_fungap_out_prot_filtered')

msa_list = sorted(os.listdir(msa_dir))

for msa in msa_list:
    print(msa)
    msa_path = join(msa_dir, msa)
    for record in SeqIO.parse(msa_path, 'fasta'):
        if 'NI907' in record.id:
            genome = record.id.split("_")[2]
        else:
            genome = record.id.split("_")[0]
        genomes_dict[genome].seq += record.seq

with open(output, 'w') as handle:
    SeqIO.write(genomes_dict.values(), handle, 'fasta')