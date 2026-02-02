#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --job-name=filter
#SBATCH --output=filter_out_%j
#SBATCH --error=filter_error_%j
#SBATCH --partition=condo
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

# load environment with vcftools
module load miniconda3
source activate gatk-env

cd $SLURM_SUBMIT_DIR
dir=$SLURM_SUBMIT_DIR

# filter and make new VCF with vcftools first without then with only Z
vcftools --vcf analysis_vcfs/filtered_snps.vcf --keep sample_list_eems \
         --remove-indels --remove-filtered-all --not-chr CM022210.2_RagTag \
         --min-alleles 2 --max-alleles 2 --mac 3 --max-meanDP 20 --min-meanDP 4 \
         --max-missing .75 --thin 10000 --minDP 2 --recode
mv out.recode.vcf analysis_vcfs/bevi_noZ_eems_75.vcf

vcftools --vcf analysis_vcfs/filtered_snps.vcf --keep sample_list_eems \
         --remove-indels --remove-filtered-all --chr CM022210.2_RagTag \
         --remove-indv Vireo_bellii_SDSU2911 --remove-indv Vireo_bellii_LSUB51626 --remove-indv Vireo_bellii_SDSU2915 \
         --min-alleles 2 --max-alleles 2 --mac 3 --max-meanDP 20 --min-meanDP 4 \
         --max-missing .75 --thin 10000 --minDP 2 --recode
mv out.recode.vcf analysis_vcfs/bevi_Z_eems_75.vcf

# load environment with plink2
conda deactivate
source activate plink2-env

## autosomes

# Set up a re-named VCF to make it compatible with PLINK/bed2diffs
# Just replaced chromosome names with numbers
cp analysis_vcfs/bevi_noZ_eems_75.vcf analysis_vcfs/bevi_noZ_rename_eems_75.vcf
while IFS=',' read -ra rename; do
    old="${rename[0]}"
    new="${rename[1]}"
    sed -i "s/$old/$new/g" analysis_vcfs/bevi_noZ_rename_eems_75.vcf
done < intervals_rename

# Make input with PLINK2
plink2 --chr-set 43 --vcf analysis_vcfs/bevi_noZ_rename_eems_75.vcf \
       --make-bed --out $dir/analysis_files/eems/bevi_noZ_eems_75
# make eems input with bed2diffs script
/users/egyllenh/eems/bed2diffs/src/bed2diffs_v1 --bfile analysis_files/eems/bevi_noZ_eems_75 --nthreads 1

## Z chromosome

# Set up a re-named VCF to make it compatible with PLINK/bed2diffs
# Just replaced chromosome names with numbers
cp analysis_vcfs/bevi_Z_eems_75.vcf analysis_vcfs/bevi_Z_rename_eems_75.vcf
while IFS=',' read -ra rename; do
    old="${rename[0]}"
    new="${rename[1]}"
    sed -i "s/$old/$new/g" analysis_vcfs/bevi_Z_rename_eems_75.vcf
done < intervals_rename

# Make input with PLINK2
plink2 --chr-set 42 --vcf analysis_vcfs/bevi_Z_rename_eems_75.vcf \
       --make-bed --out analysis_files/eems/bevi_Z_eems_75
# make eems input with bed2diffs script
/users/egyllenh/eems/bed2diffs/src/bed2diffs_v1 --bfile analysis_files/eems/bevi_Z_eems_75 --nthreads 1
