from BCBio import GFF
from os.path import join
from Bio import SeqIO
import sys
import os
import shutil

in_dir_assemblies = sys.argv[1]
out_dir_assemblies = sys.argv[2]
in_dir_gff = sys.argv[3]
out_dir_gff = sys.argv[4]
in_dir_prot = sys.argv[5]
out_dir_prot = sys.argv[6]

## create output dirs
if os.path.isdir(out_dir_assemblies):
    shutil.rmtree(out_dir_assemblies)
os.mkdir(out_dir_assemblies)
if os.path.isdir(out_dir_gff):
    shutil.rmtree(out_dir_gff)
os.mkdir(out_dir_gff)
if os.path.isdir(out_dir_prot):
    shutil.rmtree(out_dir_prot)
os.mkdir(out_dir_prot)

sample_ids = [
'KVK1',
'KVK2',
'KVK3',
'KVK4',
'KVK5',
'KVK6',
'KVK7',
'KVK8',
'KVK9',
'KVK10',
'KVK11',
'KVK12',
'KVK13',
'KVK14',
'KVK15',
'KVK16',
'KVK17',
'KVK18',
'KVK19',
'KVK20',
'KVK21',
'KVK22',
'KVK23',
'KVK24',
'KVK25',
'KVK26',
'KVK27',
'KVK28',
'KVK29',
'KVK30',
'KVK31',
'KVK32',
'KVK33',
'KVK34',
'KVK35',
'KVK36',
'KVK37',
'KVK38',
'KVK39',
'KVK40',
'KVK41',
'KVK42',
'KVK43',
'KVK44',
'KVK45',
'KVK46',
'KVK47',
'KVK48',
'KVK49',
'KVK50',
'KVK51',
'KVK52',
'KVK53',
'KVK54',
'KVK55',
'KVK56',
'KVK57',
'KVK58',
'KVK59',
'KVK60',
'KVK61',
'KVK62',
'KVK63',
'KVK64',
'KVK65',
'KVK66',
'KVK67',
'KVK68',
'KVK69',
'KVK70',
]

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

ids_to_name_dict = dict(zip(sample_ids, strain_names))
gene_name_translation_dict = {}
## get file lists
gff_list = os.listdir(in_dir_gff)
prot_list = os.listdir(in_dir_prot)

for id in sample_ids:
    print(id)
    if id == "KVK42": # contaminated samples
        continue
    strain_name = ids_to_name_dict[id]
    ## edit assembly file
    in_assembly = join(in_dir_assemblies, id+'_final.scaffolds.fasta')
    out_assembly = join(out_dir_assemblies, strain_name+'_final.scaffolds.fasta')
    record_list = list(SeqIO.parse(in_assembly, 'fasta'))
    with open(out_assembly, 'w') as corrected:
        for i in range(len(record_list)):
            record = record_list[i]
            record.id = strain_name + '_' + record.id.split('_')[1]
            record.description = ''
            SeqIO.write(record, corrected, 'fasta')
    ## gff now
    in_gff = join(in_dir_gff, id+'_fungap_out.gff3')
    in_handle = open(in_gff)
    out_gff = join(out_dir_gff, strain_name+'_fungap_out.gff3')
    # to rename protein sequences later...
    gene_name_translation_dict[strain_name] = {}
    gene_name_translation_dict_for_strain = gene_name_translation_dict[strain_name]
    with open(out_gff, "w") as out_handle:
        entries_dict = {}
        for rec in GFF.parse(in_handle):
            scaffold = rec.id.split('_')[1]
            rec.id = strain_name + '_' + scaffold
            ## prevents printing line that describes every single scaffold in the GFF
            rec.annotations = {}
            for feature in rec.features:
                if feature.type != "remark":
                    try:
                        original_feature_id = feature.id[:]
                        gene_number = str(original_feature_id.split('_')[1])
                        new_feature_id = strain_name + '_' + scaffold + '_' + gene_number
                        gene_name_translation_dict_for_strain[original_feature_id] = new_feature_id
                        ## yes you have to edit each subfeature individually...
                        feature.id = new_feature_id
                        feature.qualifiers['ID'] = [new_feature_id]
                        feature.qualifiers['Name'] = [new_feature_id]
                        feature.sub_features[0].id = [new_feature_id+'T0']
                        feature.sub_features[0].qualifiers['ID'] = [new_feature_id+'T0']
                        feature.sub_features[0].qualifiers['Parent'] = [new_feature_id]
                        for sub_feature in feature.sub_features[0].sub_features:
                            sub_feature.id = [new_feature_id+'T0']
                            sub_feature.qualifiers['Parent'] = [new_feature_id+'T0']
                            if len(sub_feature.qualifiers['ID']) > 1:
                                print('sub_feature too long')
                            sub_feature.qualifiers['ID'] = [new_feature_id+'T0'+'.'+sub_feature.qualifiers['ID'][0].split('.')[-1]]
                    except KeyError:
                        print(feature)
            entries_dict[int(scaffold)] = rec
        in_handle.close()
        ## this is so that the entries appear in scaffold order in the GFF
        for scaffold in sorted(entries_dict.keys()):
            GFF.write([entries_dict[scaffold]], out_handle)
    ## GFF writer likes to spawn a bunch of extra lines so just remove those and add a single one at the top
    with open(out_gff, "r") as input:
        with open("temp.txt", "w") as output:
            output.write('##gff-version 3\n')
            for line in input:
                if not line.strip("\n").startswith('##gff-version 3'):
                    output.write(line)
    # replace file with original name
    os.replace('temp.txt', out_gff)
    # proteome now
    in_prot = join(in_dir_prot, id+'_fungap_out_prot.faa')
    out_prot = join(out_dir_prot, strain_name+'_fungap_out_prot.faa')
    record_list = list(SeqIO.parse(in_prot, 'fasta'))
    with open(out_prot, 'w') as corrected:
        for i in range(len(record_list)):
            record = record_list[i]
            record.id = gene_name_translation_dict_for_strain[record.id.split('.t1')[0]]+'T0' ## rename records to have genome name in them
            record.description = ''
            if '*' in record.seq:
                if record.seq[-1] == '*': ## remove stop codon from end of sequences
                    record.seq = record.seq[:-1]
                    SeqIO.write(record, corrected, 'fasta')
            else:
                SeqIO.write(record, corrected, 'fasta')