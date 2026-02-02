#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --job-name=gatk_bells
#SBATCH --output=_pixy_out_%j
#SBATCH --error=_pixy_error_%j
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

# load env with vcftools, samtools, and pixy
module load miniconda3
source activate pixy-env

cd $SLURM_SUBMIT_DIR
src=$SLURM_SUBMIT_DIR

## Preparing vcf for pixy input, renaming a messed up sample name
zcat $src/invar/invar_vcf.vcf.gz | sed 's/{}/Vireo_bellii_UTEP3350/g' | gzip > $src/invar/fixed_invar_vcf.vcf.gz

# process in a way that retains invariants
vcftools --gzvcf $src/invar/fixed_invar_vcf.vcf.gz --remove-indels --remove-filtered-all \
    --min-meanDP 4 --max-meanDP 30 --minDP 2 --recode --stdout | bgzip -c > $src/analysis_vcfs/pixy_input.vcf.gz
tabix $src/analysis_vcfs/pixy_input.vcf.gz

# run pixy
pixy --stats pi fst dxy \
     --vcf analysis_vcfs/pixy_input.vcf.gz \
     --window_size 50000 \
     --populations popmap_pixy \
     --bypass_invariant_check 'no' \
     --output_folder analysis_files/pixy \
     --output_prefix pixy_output \
     --n_cores 32
