#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --job-name=mtdna_bells
#SBATCH --output=_ND2_out_%j
#SBATCH --error=_ND2_error_%j
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenhaal@unm.edu

# load conda env with bwa/picard/bcftools and env_parallel
module load miniconda3
module load parallel
source activate bcftools-env
source $(which env_parallel.bash)
src=$SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

# reference base name
reference=${src}/bells_ref_SDSU2873

# rest combined fasta
> bells_combined_ND2.fasta

# index reference
bwa index -p $reference ${reference}.fasta

# loop for running across samples
# pipeline order is allign, sort, call variants, index, generating consensus, and combining
while read sample; do
    bwa mem \
         -t 24 \
         $reference \
         /carc/scratch/projects/andersen2016005/EFG_stuff/bells_genomes/clean_reads/${sample}_paired_R1.fastq.gz \
         /carc/scratch/projects/andersen2016005/EFG_stuff/bells_genomes/clean_reads/${sample}_paired_R2.fastq.gz \
         | samtools view -b -F 4 - > $src/bams/${sample}_unsort_ND2.bam
     picard SortSam I=$src/bams/${sample}_unsort_ND2.bam \
            O=$src/bams/${sample}_sort_ND2.bam \
            SORT_ORDER=coordinate
     bcftools mpileup -Q 30 -q 30 -Ou -f $reference.fasta $src/bams/${sample}_sort_ND2.bam --threads 8 | \
              bcftools call -mv -Oz --threads 24 -o $src/vcfs/${sample}_ND2.vcf.gz
     bcftools index $src/vcfs/${sample}_ND2.vcf.gz
     cat ${reference}.fasta | bcftools consensus $src/vcfs/${sample}_ND2.vcf.gz \
         > $src/output/${sample}_consensus_ND2.fasta
     echo \>${sample}_ND2 >> $src/bells_combined_ND2.fasta
     tail -n +2 $src/output/${sample}_consensus_ND2.fasta >> $src/bells_combined_ND2.fasta
done < ../sample_list
