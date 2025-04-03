#!/bin/bash

#SBATCH -J CLSW_Bayescan2
#SBATCH -o CLSW_Bayescan2.out
#SBATCH -c 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vlhuerta@outlook.es
#SBATCH --time=0-30:00
#SBATCH --mem=32G

module load StdEnv/2020
module load bayescan/2.1


bayescan_2.1 CLSW_Bayescan2.txt -od ./ -o CLSW_outlier_fst -n 100000 -threads 8 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 700
