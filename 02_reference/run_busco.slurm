#!/bin/bash

#SBATCH --job-name=busco_test
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --partition=singleGPU
##SBATCH --mail-user=egyllenhaal@unm.edu

# load environment with BUSCO and others
# note that BUSCO can be hard to install, there is a yml in this tutorial I made:
## https://github.com/UNM-CARC/QuickBytes/blob/master/genome_evaluation.md
module load miniconda3/4.8.2-pilj
source activate assembly-env

cd $SLURM_SUBMIT_DIR

# busco commands for a variety of references
# Supernova output
busco -i bevi_ref_pseudohap.fasta -o busco_original -f --cpu 16 --mode genome -l passeriformes_odb10
# arks output
busco -i arks/bevi_ref_pseudohap_c5_m50-10000_k30_r0.05_e30000_z500_l5_a0.3.scaffolds.fa \
    -o busco_arks -f --cpu 16 --mode genome -l passeriformes_odb10
# full ragtag
busco -i ragtag_output/ragtag.scaffold.fasta -o busco_ragtag -f --cpu 16 --mode genome -l passeriformes_odb10
# chromosome 0 ragtag
busco -i ragtag_chr0/ragtag.scaffold.fasta -o busco_chr0 -f --cpu 16 --mode genome -l passeriformes_odb10
# ragtag throwing out chromosome 0
busco -i ragtag_chr0/ragtag.scaffold_nochr0.fasta -o busco_nochr0 -f --cpu 16 --mode genome -l passeriformes_odb10
