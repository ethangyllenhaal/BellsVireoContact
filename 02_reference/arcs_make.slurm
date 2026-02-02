#!/bin/bash

#SBATCH --job-name=arcs_make
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --partition=singleGPU
##SBATCH --mail-user=egyllenhaal@unm.edu

# Load conda environment with arcs in it
module load miniconda3/4.8.2-pilj
source activate assembly-env

cd $SLURM_SUBMIT_DIR

# run arcs-make with 16 threads
arcs-make arks draft=bevi_ref_pseudohap reads=bevi_longranger/outs/barcoded t=16

# set base name
f=bevi_ref_pseudohap_

# run links
LINKS -f $f -s empty.fof -k 20 -b $f.c5_e30000_r0.05 -l 5 -t 2 -a 0.3

