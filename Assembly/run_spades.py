import os

for p in ["KVK" + str(x) for x in range(1,71)]:
  if "scaffolds.fasta" not in os.listdir(p):
    reads1 = sorted([ "-1 ../../Trim/" + r for r in os.listdir("../../Trim") if r.split("_")[0] == p and "1_val_1" in r])
    reads2 = sorted([ "-2 ../../Trim/" + r for r in os.listdir("../../Trim") if r.split("_")[0] == p and "2_val_2" in r])

    #Final parameter: k-mer = 127
    os.system(f"spades.py -o {p} -t 56 -k 127 " + " ".join(reads1) + "  " + " ".join(reads2))
