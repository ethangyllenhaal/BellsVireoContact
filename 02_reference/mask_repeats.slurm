#!/bin/bash

#SBATCH --job-name=repeat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --partition=condo
#SBATCH --output=_repeat_out_%j
#SBATCH --error=_repeat_error_%j
#SBATCH --mail-user=egyllenhaal@unm.edu

# load repeat masker env
module load miniconda3
source activate repeat-env

cd $SLURM_SUBMIT_DIR
dir=$SLURM_SUBMIT_DIR

# build repeat database
BuildDatabase -name $dir/repeatmask_hopper/Vireo-bellii $dir/ragtag_chr0.fa

# run RepeatModeler
RepeatModeler -threads 32 -database $dir/repeatmask_hopper/Vireo-bellii 

## RepeatModeler took a while, so ran this script in chunks, below after RepeatModeler finished

# RepeatModeler didn't get past round 4, so limited subset for useage
RepeatClassifier -consensi $dir/rmod_best_30Mar24_round4/consensi.fa \
                 -stockholm $dir/rmod_best_30Mar24_round4/families.stk

cd rmod_best_30Mar24_round4

# run RepeatMasker with denovo repeats
RepeatMasker -pa 32 -xsmall -lib consensi.fa.classified \
	     -dir $dir/repeatmask/Vireo-bellii_denovo \
	     $dir/ragtag_chr0.fa
# run it again on the denovo output with Zebra Finch
RepeatMasker -pa 32 -xsmall -species "zebra finch" \
    -dir $dir/repeatmask/Vireo-bellii_combined \
    $dir/repeatmask/Vireo-bellii_denovo/ragtag_chr0.fa.masked

# make a copy
cp $dir/repeatmask/Vireo-bellii_combined/ $dir/bells_ragtag_masked.fa
