#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --job-name=admix_bells
#SBATCH --output=_admix_out_%j
#SBATCH --error=_admix_error_%j
#SBATCH --partition=condo
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

module load miniconda3
source activate admixture-env

# By: Ethan Gyllenhaal
# Last updated: 1 Feb 2026
# Shell script for running admixture.
# Starts with non-12 formatted Plink ped and map files.
# I run this with a conda environment with plink and admixture installed.

# run Admixture
for K in 1 2 3 4 5; \
	do admixture -m EM --cv bevi_Z_focal_75_012.ped ${K} | tee Z_k_log${K}.out; done
# spit out K values
grep -h CV Z_k_log*.out > bevi_Z_kval

# run Admixture
for K in 1 2 3 4 5; \
	do admixture -m EM --cv bevi_noZ_focal_75_012.ped ${K} | tee noZ_k_log${K}.out; done
# spit out K values
grep -h CV noZ_k_log*.out > bevi_noZ_kval
