#!/bin/bash

#SBATCH --job-name=abyss
#SBATCH --account=co_minium
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --qos=savio_lowprio
#SBATCH --time=72:00:00

#Run GUY11 assemblies with abyss with varying kc and k

for kc in 2 3 4; do
        for k in $(seq 50 4 110); do
                mkdir k${k}-kc${kc}
                abyss-pe -C k${k}-kc${kc} \
                B=2G j=56 name=MO k=$k kc=$kc \
                in='/global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_1_val_1.fq.gz /global/scratch/users/skyungyong/CO_Pierre_MO/Trim/KVK35_CKDN220051608-1A_HFW35DSX5_L4_2_val_2.fq.gz'
        done
done

abyss-fac k*/MO-scaffolds.fa
