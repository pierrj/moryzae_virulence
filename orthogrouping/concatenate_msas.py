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
'BR0019',
'CD0073',
'CH0043',
'CH0052',
'CH0063',
'CH0072',
'CH0092',
'CH0110',
'CH0333',
'CH0452',
'CH0461',
'CH0532',
'CH0533',
'CH0539',
'CH0549',
'CH0561',
'CH0562',
'CH0565',
'CH0571',
'CH0718',
'CH1008',
'CH1009',
'CH1016',
'CH1065',
'CH1076',
'CH1079',
'CH1083',
'CH1103',
'CH1120',
'CH1121',
'CH1150',
'CH1164',
'CL0026',
'FR0013',
'GY0011',
'HN0001',
'IN0017',
'IN0054',
'IN0072',
'IN0076',
'IN0082',
'IN0092',
'IN0094',
'IN0115',
'IN0116',
'LA0006',
'LA0009',
'LA0012',
'LA0019',
'LA0023',
'MC0016',
'MD0116',
'MD0929',
'NP0058',
'NP0070',
'PH0019',
'PH0103',
'PH0158',
'PR0009',
'SP0005',
'SP0006',
'TN0016',
'TN0070',
'US0032',
'US0041',
'US0090',
'US0098',
'VT0003',
'VT0027',
'VT0030']

for strain in strain_names:
    genomes_dict[strain] = SeqRecord(Seq(""), id=strain)

genomes_dict['NI907'] = SeqRecord(Seq(""), id='NI907')

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