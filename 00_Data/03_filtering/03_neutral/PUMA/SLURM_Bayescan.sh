#!/bin/bash

#SBATCH -J PUMA_Bayescan
#SBATCH -o PUMA_Bayescan.out
#SBATCH -c 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vlhuerta@outlook.es
#SBATCH --time=0-40:00
#SBATCH --mem=40G

module load StdEnv/2020
module load bayescan/2.1


bayescan_2.1 PUMA_Bayescan.txt -od ./ -o PUMA_outlier_fst -n 100000 -threads 8 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 700
