#!/bin/bash

#SBATCH --ntasks=38
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=psmc_bells
#SBATCH --output=_psmc_out_%j
#SBATCH --error=_psmc_error_%j
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

module load miniconda3
module load parallel

source $(which env_parallel.bash)
source activate psmc-env
export LANG=C

cd $SLURM_SUBMIT_DIR
src=$SLURM_SUBMIT_DIR
scontrol show hostname > node_list_${SLURM_JOB_ID}
reference=$src/reference/ragtag_chr0

# standard pipeline for making PSMC input
cat sample_list_ref | env_parallel -j 1 --sshloginfile node_list_${SLURM_JOB_ID} \
    'mkdir $src/psmc/{} $src/psmc/{}/boot
     bcftools mpileup -Q 30 -q 30 -Ovu -f ${reference}.fa \
         $src/alignments/{}_dedup.bam --threads 8 | \
     bcftools call -c --threads 8 | \
     vcfutils.pl vcf2fq -d 10 -D 50 -Q 30 > $src/psmc/{}/{}_consensus.fq
     ~/psmc/utils/fq2psmcfa $src/psmc/{}/{}_consensus.fq > \
         $src/psmc/{}/{}.psmcfa
     ~/psmc/psmc -N25 -t10 -r5 -p "32*2+4+6" \
         -o $src/psmc/{}/{}.psmc \
         $src/psmc/{}/{}.psmcfa'

# run PSMC bootstraps in parallel
env_parallel -j 8 --sshloginfile node_list_${SLURM_JOB_ID} \
    '~/psmc/psmc -N25 -t10 -r5 -b -p "32*2+4+6" \
         -o $src/psmc/{1}/boot/{1}_r{2}.psmc $src/psmc/{1}/{1}.psmcfa' \
    ::: $(cat sample_list_psmc) ::: $(seq 50)

# combine bootstaps
env_parallel -j 1 --sshloginfile ./node_list_${SLURM_JOB_ID} \
    'cat $src/psmc/{}/{}.psmc $src/psmc/{}/boot/{}_r*.psmc > $src/psmc/{}/{}_combined.psmc' \
    ::: $(cat sample_list_psmc)
