# Asymmetric gene flow across a desert contact zone in a riparian songbird

Last updated 2 February 2026

Repository for scripts and input files used in a study of a contact zone of Bell's Vireos.

This README describes the scripts and data files used for the afforementioned project. Each section corresponds to a zipped directory in the Dryad repository and folder on GitHub. \*\*\* denotes files in the Dryad only (mostly VCF or alignment files).

## Sampling table and map (01_sampling)

Raw reads are housed at NCBI's Sequencing Read Archive (PRJNAXXX).

bells_sampling_coords - Sampling coordinates, with columns of: sample name, latitude, longitude, mtDNA haplotype, and geographic cluster.

map_BEVI.R - R Script for making sampling maps used in the paper.

## Reference assembly and annotation (02_reference)

supernova.pbs - Torque script for running SuperNova to assemble the reference (dates the project a bit!).

arcs_make.slurm - Slurm script to run the arcs-make and LINKS on the fragmented reference.

ragtag.slurm - Slurm script to make pseudoreference genome from Arks output with RagTag.

run_busco.slurm - Script for running BUSCO with a variety of references. Note that BUSCO has had issues on conda in the past, there is a YML that worked for me in this tutorial I put together while working on this: https://github.com/UNM-CARC/QuickBytes/blob/master/genome_evaluation.md.

quast.slurm - Script for running quast quality checks.

\*\*\*ragtag_chr0.fa - Reference with unmatched scaffolds combined into "chromosome 0".

mask_repeats.slurm - Script for annotating and masking repeats.

gemoma.slurm - Run reference-based annotation.

final_annotation.gff - Genome feature file output by GeMoMa.

## Primary variant calling pipeline (03_variant_calling)

gatk_pipeline.slurm - The big pipeline! Run in chunks as parts finished, using GNU parallel for parallelization.

filter_vcf.slurm - Script for making main two VCF subsets.

coverage.slurm - Script for estimating coverage based on lightly filtered called SNPs.

## Mitochondrial DNA genotyping (04_mitochondrial)

mtdna_bells.slurm - Slurm script for making mtDNA fasta.

\*\*\*bells_combined_ND2.fasta - FASTA file for all samples' ND2.

## Population structure (05_population_structure)

There are four different analyses we include, in four subsections outlined below. Ones used across categories are immediately below.

\*\*\*bevi_noZ_focal_75.vcf - VCF used for PCA, admixture, and pairwise FST estimates, excludes sex chromosomes.

\*\*\*bevi_Z_focal_75.vcf - As above, but with only the Z chromosome included and three individuals removed that were flagged as unusual for the Z in PCA analyses.


### PCA

BEVI_PCA.R - R script for making genomic PCA plots, including specific regions. For the specific regions, density plots of PC scores are also included.

popmap.csv/popmap_Z.csv - Population maps for relevant sets of individuals, used for labeling clusters in the PCA.

\*\*\*bevi_\[Z13M/Chr1-16M/Chr12-8M\]_focal_75.vcf - VCFs for three specific regions investigated in depth in Figure 4, otherwise like noZ/Z filtering schemes described at start of section.


### Admixture

setup_admix.sh - Script for making admixture input, starting from a VCF then modifying it into a sorted file. Relies on two helper scripts:

* sed_admixture.sh - Replaces numeric placeholders with sample names in 012 format file
* greps_admix.sh - Prints lines in desired order for plotting

\*\*\*bevi_\*\_focal_75_012_sort.ped - Plink style input for admixture for Z or autosomal datasets.

\*\*\*bevi_\*\_focal_75_012_sort.map - Plink style region map for admixture for Z or autosomal datasets.

run_admixture_bells.slurm - Slurm script for running admixture across K values.

bevi_admix_plot.R - R script for plotting admixture output.


### FST

drive_fst.sh - Script with commands used to run pairwise FST analyses.

run_fst.sh/run_fst_Z.sh - Scripts for running pairwise FST analyses (a bit outdated, but get the job done).

pops.zip/pops_Z.zip - Zipped directories of simple text files of samples from given populations.


### EEMS

make_eems_input.slurm - Slurm script for filtering main VCFs then making eems input.

make_chains.sh - Simple slurm script for copying chain1's parameters to 19 other chains.

run_eems_\*.slurm - Slurm script for running EEMS chains without and with only Z as a Slurm array job.

The next are a set of input files for EEMS, each one has a version of only the Z and everything but the Z

\*\*\*bevi_\*_eems_75.bed - BED file with genetic data.

\*\*\*bevi_\*_eems_75.bim - BIM file outlining variants

bevi_\*_eems_75.coord - File of sample coordinates.

\*\*\*bevi_\*_eems_75.count - List of allele counts.

\*\*\*bevi_\*_eems_75.diffs - Genetic distances for eems.

\*\*\*bevi_\*_eems_75.fam - FAM metadata file.

bevi_\*_eems_75.order - Sample name orders.

bevi_\*_eems_75.outer - EEMS outer bounds.

params_\*_chain1 - Parameter file that is then cloned per chain.

## Effective population size over time (06_demography)

run_psmc.slurm - Script for making PSMC input (standard pipeline from package's github) and then running PSMC, including its bootstraps.

## Divergence and diversity across the genome (07_pixy)

invar_gatk.slurm - Pipeline used for making invariant site VCFs

pixy_bevi.slurm - Slurm script for running pixy.

popmap_pixy - Population map used for pixy analyses.

pixy_output_\*.txt - Output from pixy for Dxy, Fst, and pi.

div_ratio_plots.R - Script for exploring primarily diversity ratio, but also other pixy-related outputs (e.g., pi vs Dxy, genome-wide weighted mean pi)

## Population genomics simualtions (08_simulations)

fst_div_ratio.slim - SLiM 5 script for running simple WF simulations of divergence followed by asymmetric gene flow.

make_params.R - R script for making parameters for job array to loop over.

run_ratio.slurm - Slurm script for running replicates as a job array.

\*\*\*simulation_output.zip - Zipped directory of output.

sim_calc_div_fst.R - R script for plotting and exploring output. Includes a section comparing autosomal and Z chromosome values.
