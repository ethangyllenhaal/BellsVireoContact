#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --job-name=filter
#SBATCH --output=filter_out_%j
#SBATCH --error=filter_error_%j
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

# load environment with VCFtools
module load miniconda3
source activate gatk-env

cd $SLURM_SUBMIT_DIR

# Make and then rename a lot of VCFs! Summary above each one before command.

# main thinned 75% complete VCF without Z
vcftools --vcf analysis_vcfs/filtered_snps.vcf --keep sample_list_focal \
	 --remove-indels --remove-filtered-all --not-chr CM022210.2_RagTag \
	 --min-alleles 2 --max-alleles 2 --mac 3 --max-meanDP 20 --min-meanDP 4 \
	 --max-missing .75 --thin 10000 --minDP 2 --recode
mv out.recode.vcf analysis_vcfs/bevi_noZ_focal_75.vcf

# main thinned 75% complete VCF with only Z, but without a few odd individuals
vcftools --vcf analysis_vcfs/filtered_snps.vcf --keep sample_list_focal \
	 --remove-indels --remove-filtered-all --chr CM022210.2_RagTag \
	 --remove-indv Vireo_bellii_SDSU2911 --remove-indv Vireo_bellii_LSUB51626 --remove-indv Vireo_bellii_SDSU2915 \
	 --min-alleles 2 --max-alleles 2 --mac 3 --max-meanDP 20 --min-meanDP 4 \
	 --max-missing .75 --thin 10000 --minDP 2 --recode
mv out.recode.vcf analysis_vcfs/bevi_Z_focal_75.vcf
