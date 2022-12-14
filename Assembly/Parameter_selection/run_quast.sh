#!/bin/bash

#SBATCH --job-name=quast
#SBATCH --account=co_minium
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --qos=savio_lowprio
#SBATCH --time=72:00:00

#Run Quast to evaluate genome assembly quality against the reference genome (MagorGY11 from MycoCosm)
#Spades
ls ../../Spades/ | grep "isolate" | while read assembly; do 
  quast.py -t 56 -r MagorGY11_1_AssemblyScaffolds.fasta \
  -g mRNA:MagorGY11_1_GeneCatalog_20210817.gff3 --fast \
  -o ${assembly} ../../Spades/${assembly}/scaffolds.fasta; 
done

#Velvet
ls ../../Velvet/ | grep "hash" | while read assembly; do 
  quast.py -t 56 -r MagorGY11_1_AssemblyScaffolds.fasta \
  -g mRNA:MagorGY11_1_GeneCatalog_20210817.gff3 --fast \
  -o ${assembly} ../../Velvet/${assembly}/contigs.fa; 
done

#Abyss
ls ../../Abyss/ | grep "kc" | while read assembly; do 
  quast.py -t 56 -r MagorGY11_1_AssemblyScaffolds.fasta \
  -g mRNA:MagorGY11_1_GeneCatalog_20210817.gff3 --fast \
  -o ${assembly} ../../Abyss/${assembly}/MO-scaffolds.fa; 
done

cd ../NCBI

#Run Quast to evaluate genome assembly quality against the reference genome (GCA_002368485.1_ASM236848v1_genomic.fna from NCBI)
#Spades
ls ../../Spades/ | grep "isolate" | while read assembly; do 
  quast.py -t 56 -r GCA_002368485.1_ASM236848v1_genomic.fna \
  --fast \
  -o ${assembly} ../../Spades/${assembly}/scaffolds.fasta; 
done

#Velvet
ls ../../Velvet/ | grep "hash" | while read assembly; do 
  quast.py -t 56 -r GCA_002368485.1_ASM236848v1_genomic.fna \
  --fast \
  -o ${assembly} ../../Velvet/${assembly}/contigs.fa; 
done

#Abyss
ls ../../Abyss/ | grep "kc" | while read assembly; do 
  quast.py -t 56 -r GCA_002368485.1_ASM236848v1_genomic.fna \
  --fast \
  -o ${assembly} ../../Abyss/${assembly}/MO-scaffolds.fa; 
done
