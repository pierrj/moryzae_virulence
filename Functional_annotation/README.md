## InterProScan  

The representative per orthologous group was annotated with InterProScan v5.61-93.0.  
The databases scanned against are: PFAM v35.0, Gene3D v4.3.0, and SUPERFAMILY v1.75. 

<code>interproscan.sh -appl Pfam,SUPERFAMILY,Gene3D -cpu 40 -dp -i rep.fasta -o rep.tsv -f tsv</code>  

## SignalP  

SignalP v4.1 was used to predict the signal peptides.The parameters used for the run are given below.  
Also see <code>signalp.slurm</code>   

<code>signalp -t euk -u 0.34 -U 0.34 -f short > signalp.out </code>   
