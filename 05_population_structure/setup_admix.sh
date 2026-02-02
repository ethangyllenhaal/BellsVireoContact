#!/bin/bash

# By: Ethan Gyllenhaal
# Last updated: 1 Feb 2026
# script for making admixture input with 5 steps per VCF (with only or without Z)
## 1) run VCFtools to make plink format output
## 2) Make Plink 012 output
## 3) Replace holder integers with individual names
## 4) Re-arrange in desired order
## 5) Make a copy of plink's map file for new sorted output

vcftools --vcf analysis_vcfs/bevi_noZ_focal_75.vcf --plink --out bevi_noZ_focal_75
plink --file bevi_noZ_focal_75 --recode12 --out bevi_noZ_focal_75_012
sh sed_admixture.sh bevi_noZ_focal_75
sh greps_admix.sh bevi_noZ_focal_75
cp bevi_noZ_focal_75_012.map bevi_noZ_focal_75_012_sort.map

vcftools --vcf analysis_vcfs/bevi_Z_focal_75.vcf --plink --out bevi_Z_focal_75
plink --file bevi_Z_focal_75 --recode12 --out bevi_Z_focal_75_012
sh sed_admixture.sh bevi_Z_focal_75
sh greps_admix.sh bevi_Z_focal_75
cp bevi_Z_focal_75_012.map bevi_Z_focal_75_012_sort.map
