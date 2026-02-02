#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=cover
#SBATCH --output=_cover_out_%j
#SBATCH --error=_cover_error_%j
#SBATCH --partition=debug
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

module load miniconda3/4.10.3-an4v
source activate gatk-env

cd $SLURM_SUBMIT_DIR

vcftools --vcf analysis_vcfs/filtered_snps.vcf --remove-indels --remove-filtered-all --depth \
         --out ./analysis_files/bevi_coverage

vcftools --vcf analysis_vcfs/filtered_snps.vcf --remove-indels --remove-filtered-all --depth \
	 --out ./analysis_files/bevi_coverage
