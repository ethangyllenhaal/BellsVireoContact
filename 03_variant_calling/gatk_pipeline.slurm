#!/bin/bash

#SBATCH --ntasks=5
#SBATCH --cpus-per-task=8
#SBATCH --time=30:00:00
#SBATCH --job-name=gatk_bells
#SBATCH --output=_gatk_out_%j
#SBATCH --error=_gatk_error_%j
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

# load GNU parallel module (specifically with env parallel)
# also load an environment with gatk4, samtools, bwa, and picard
module load miniconda3
module load parallel
source activate gatk-env
source $(which env_parallel.bash)

# change to directory and makes a shortcut variable for dir and reference
cd $SLURM_SUBMIT_DIR
src=$SLURM_SUBMIT_DIR
reference=${src}/reference/ragtag_chr0

# silences some warning messages
export LANG=C

# make node list for parallelization
scontrol show hostname > $src/node_list_${SLURM_JOB_ID}

# indexing reference
bwa index -p $reference ${reference}.fa
samtools faidx ${reference}.fa -o ${reference}.fa.fai
picard CreateSequenceDictionary \
       R=${reference}.fa \
       O=${reference}.dict

# Trimming section, only need to do once
# assumes you have the trimmomatic fasta in your ref folder
adapters=$src/reference/TruSeq3-PE.fa
cat $src/sample_list | env_parallel -j 8 --sshloginfile $src/node_list_${SLURM_JOB_ID} \
    'read1=$src/raw_reads/{}_R1_001.fastq.gz
    read2=$src/raw_reads/{}_R2_001.fastq.gz
    paired_r1=$src/clean_reads/{}_paired_R1.fastq.gz
    paired_r2=$src/clean_reads/{}_paired_R2.fastq.gz
    unpaired_r1=$src/clean_reads/{}_unpaired_R1.fastq.gz
    unpaired_r2=$src/clean_reads/{}_unpaired_R2.fastq.gz
    # the minimum read length accepted, we do the liberal 30bp here
    min_length=30
    trimmomatic PE -threads 1 \
	$read1 $read2 $paired_r1 $unpaired_r1 $paired_r2 $unpaired_r2 \
	ILLUMINACLIP:${adapters}:2:30:10:2:keepBothReads \
	LEADING:3 TRAILING:3 MINLEN:${min_length}'

# Section for alignment and marking duplicates. 
# Note we parallelize such that BWA uses exactly one node.
# Then, we have a number of jobs equal to the number of nodes requested.
# MarkDuplicates can take a while depending on the reference and number of reads
## you may need to break it into steps for certain files

cat $src/sample_list_trim | env_parallel -j 1 --sshloginfile $src/node_list_${SLURM_JOB_ID} \
    'bwa mem \
    	 -t 8 -M \
	 -R "@RG\tID:{}\tPL:ILLUMINA\tLB:{}\tSM:{}" \
	 $reference \
	 $src/clean_reads/{}_paired_R1.fastq.gz \
	 $src/clean_reads/{}_paired_R2.fastq.gz \
	 | samtools view -b -F 4 - > $src/alignments/{}_unsort.bam
     mkdir $src/alignments/dedup_temp/{}
     mkdir $src/alignments/dedup_temp/spark/{}
     gatk MarkDuplicatesSpark \
	 -I $src/alignments/{}_unsort.bam \
	 --tmp-dir $src/alignments/dedup_temp \
         --conf "spark.local.dir=$src/alignments/dedup_temp/spark/{}" \
         --spark-master local[*] \
	 -O $src/alignments/{}_dedup.bam'

# make chromosome interval list, manually removed chr0 from it later
cut -f 1 ${reference}.fa.fai > $src/intervals.list

# Scatter-gather HaploType Caller, probably the most likely to need checkpoints.
# This can take a lot of different forms, this one is best for large files.
# Note that if the job terminates early, you have to redo the whole sample with this set up
# That shouldn't be too big of a problem, tho

env_parallel -j 32 --sshloginfile $src/node_list_${SLURM_JOB_ID} \
      'mkdir ${src}/gvcfs/{1}
      gatk --java-options "-Xmx3g" HaplotypeCaller \
           -R ${reference}.fa \
	   -I $src/alignments/{1}_dedup.bam \
	   -O $src/gvcfs/{1}/{1}_{2}_raw.g.vcf.gz \
	   -L {2} \
	   -ERC GVCF' ::: $(cat $src/sample_list_missdedup) ::: $(cat $src/intervals_nochr0.list)

# run genomicsdbimport for all samples
# This runs it scatter-gather
## the path to each interval (genomics_databases/{interval}) must be empty
cat $src/intervals_nochr0.list | env_parallel -j 4 --sshloginfile $src/node_list_${SLURM_JOB_ID} \
    'mkdir $src/gendb_temp/{}
    interval_list=""
    # loop to generate list of sample-specific intervals to combine
    while read sample; do
          interval_list="${interval_list}-V ${src}/gvcfs/${sample}/${sample}_{}_raw.g.vcf.gz "
          done < $src/sample_list
    gatk --java-options "-Xmx12g" GenomicsDBImport \
         ${interval_list} \
         --genomicsdb-workspace-path $src/genomics_databases/{} \
         --tmp-dir $src/gendb_temp/{} \
         -L {}'

# Run GenotypeGVCFs on each interval GVCF
cat $src/intervals_nochr0.list | env_parallel --sshloginfile $src/node_list_${SLURM_JOB_ID} \
    'gatk --java-options "-Xmx6g" GenotypeGVCFs \
    	  -R ${reference}.fa \
	  -V gendb://$src/genomics_databases/{} \
	  -O $src/combined_vcfs/intervals/{}_genotyped.vcf.gz'

# This point onward doesn't use parallel, best practices would be to stop here and restart with a single node

# Make a file with a list of paths for GatherVcfs to use
> $src/combined_vcfs/gather_list
while read interval; do
      echo ${src}/combined_vcfs/intervals/${interval}_genotyped.vcf.gz >> \
      $src/combined_vcfs/gather_list
done < $src/intervals_nochr0.list

# Run GatherVcfs
gatk GatherVcfs \
     -I $src/combined_vcfs/gather_list \
     -O $src/combined_vcfs/combined_vcf.vcf.gz

# Index the gathered VCF
gatk IndexFeatureFile \
     -I $src/combined_vcfs/combined_vcf.vcf.gz

# Select and filter variants
gatk SelectVariants \
     -R ${reference}.fa \
     -V $src/combined_vcfs/combined_vcf.vcf.gz \
     -select-type SNP \
     -O $src/combined_vcfs/raw_snps.vcf.gz

gatk SelectVariants \
     -R ${reference}.fa \
     -V $src/combined_vcfs/combined_vcf.vcf.gz \
     -select-type INDEL \
     -O $src/combined_vcfs/raw_indel.vcf.gz

gatk VariantFiltration \
     -R ${reference}.fa \
     -V $src/combined_vcfs/raw_snps.vcf.gz \
     -O $src/analysis_vcfs/filtered_snps.vcf  \
     -filter "DP < 4" --filter-name "DP_filter" \
     -filter "QUAL < 30.0" --filter-name "Q_filter" \
     -filter "QD < 2.0" --filter-name "QD_filter" \
     -filter "FS > 60.0" --filter-name "FS_filter" \
     -filter "MQ < 40.0" --filter-name "MQ_filter"

gatk VariantFiltration \
     -R ${reference}.fa \
     -V $src/combined_vcfs/raw_snps.vcf.gz \
     -O $src/analysis_vcfs/filtered_snps_gatkDP2.vcf  \
     -filter "DP < 2" --filter-name "DP_filter" \
     -filter "QUAL < 30.0" --filter-name "Q_filter" \
     -filter "QD < 2.0" --filter-name "QD_filter" \
     -filter "FS > 60.0" --filter-name "FS_filter" \
     -filter "MQ < 40.0" --filter-name "MQ_filter"
# At this point you run vcftools to filter again, I'll include some examples!
