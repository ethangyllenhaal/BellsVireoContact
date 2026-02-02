#!/bin/bash

#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --job-name=gatk_bells
#SBATCH --output=_gatk_out_%j
#SBATCH --error=_gatk_error_%j
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

# load gatk environment and env_parallel
module load miniconda3
module load parallel
source activate gatk-env
source $(which env_parallel.bash)

# make shortcut for working dir and reference
cd $SLURM_SUBMIT_DIR
src=$SLURM_SUBMIT_DIR
reference=${src}/reference/ragtag_chr0

# make node list
scontrol show hostname > $src/node_list_${SLURM_JOB_ID}

# make genomics databases
cat $src/intervals_redo.list | env_parallel -j 1 --sshloginfile $src/node_list_${SLURM_JOB_ID} \
    'mkdir $src/gendb_temp/{}
    interval_list=""
    # loop to generate list of sample-specific intervals to combine
    while read sample; do
          interval_list="${interval_list}-V ${src}/gvcfs/${sample}/${sample}_{}_raw.g.vcf.gz "
          done < $src/sample_list
    gatk --java-options "-Xmx44g" GenomicsDBImport \
         ${interval_list} \
         --genomicsdb-workspace-path $src/genomics_databases/{} \
         --tmp-dir $src/gendb_temp/{} \
         -L {}'


# Run GenotypeGVCFs on each interval GVCF, keeping all sites
cat $src/intervals_redo.list | env_parallel -j 1 --sshloginfile $src/node_list_${SLURM_JOB_ID} \
    'gatk --java-options "-Xmx44g" GenotypeGVCFs \
    	  -R ${reference}.fa \
	  -V gendb://$src/genomics_databases/{} \
	  -O $src/invar/intervals/{}_genotyped.vcf.gz \
          -all-sites \
          -L {}'

# This point onward doesn't use parallel, best practices would be to stop here and restart with a single node

# Make a file with a list of paths for GatherVcfs to use
> $src/invar/gather_list
while read interval; do
      echo ${src}/invar/intervals/${interval}_genotyped.vcf.gz >> \
      $src/invar/gather_list
done < $src/intervals_nochr0.list

# Run GatherVcfs
gatk GatherVcfs \
     -I $src/invar/gather_list \
     -O $src/invar/invar_vcf.vcf.gz

# Index the gathered VCF
gatk IndexFeatureFile \
     -I $src/invar/invar_vcf.vcf.gz

