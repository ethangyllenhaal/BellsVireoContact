#!/bin/bash

#SBATCH --job-name=ragtag
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --partition=singleGPU
#SBATCH --output=_ragtag.out
##SBATCH --mail-user=egyllenhaal@unm.edu

# load conda env with ragtag
module load miniconda3/4.8.2-pilj
source activate assembly-env

cd $SLURM_SUBMIT_DIR

# command for keeping unmapped contigs and fragments
#ragtag.py scaffold crow_ref/hocr_ref.fna arks/bevi_ref_pseudohap_c5_m50-10000_k30_r0.05_e30000_z500_l5_a0.3.scaffolds.fa \
#    -t 16 -u

# commands for making a fake chromosome 0 of unmatched bits
ragtag.py scaffold crow_ref/hocr_ref.fna arks/bevi_ref_pseudohap_c5_m50-10000_k30_r0.05_e30000_z500_l5_a0.3.scaffolds.fa \
    -t 16 -u -C -o ./ragtag_chr0
