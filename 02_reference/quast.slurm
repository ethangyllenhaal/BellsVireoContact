#!/bin/bash

#SBATCH --job-name=quast
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --partition=singleGPU
#SBATCH --output=_quast.out
##SBATCH --mail-user=egyllenhaal@unm.edu

# environment with only quast 
module load miniconda3/4.8.2-pilj
source activate quast-env

# run a few batches of quast

# supernova output
quast $SLURM_SUBMIT_DIR/bevi_ref_pseudohap.fasta.gz -o $SLURM_SUBMIT_DIR/quast_original \
    --large -t 16
# arks output
quast $SLURM_SUBMIT_DIR/arks/bevi_ref_pseudohap_c5_m50-10000_k30_r0.05_e30000_z500_l5_a0.3.scaffolds.fa \
    -o $SLURM_SUBMIT_DIR/quast_arks --glimmer --large -t 16
# full ragtag
quast $SLURM_SUBMIT_DIR/ragtag_output/ragtag.scaffold.fasta -o $SLURM_SUBMIT_DIR/quast_ragtag --large -t 16
# chromosome 0 ragtag
quast $SLURM_SUBMIT_DIR/ragtag_chr0/ragtag.scaffold.fasta -o $SLURM_SUBMIT_DIR/quast_ragtag_chr0 --large -t 16
# ragtag without chr0
quast /projects/shared/tutorials/quickbytes/GATK/sagegrouse_reference.fa -o $SLURM_SUBMIT_DIR/quast_output --large --threads 8
