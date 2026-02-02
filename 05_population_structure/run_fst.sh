#!/bin/bash

# Made by Ethan F. Gyllenhaal (egyllenhaal@unm.edu)
# Last updated: 7Jan2026
#
# This is a shell script used as a wrapper for generating a pairwise Fst table for a given set of populations using vcftools.
# Takes one argument (output file), with two other inputs I manually set here (see below).
# One is a VCF file with all relevant samples included.
# The other input is a directory of text files. Each text files has a population name and names of all samples (names as defined in VCF).
### here "pops"
#
# I did this the "lazy" way, and didn't make a command line argument for anything but the output file.
# As a matter of principle, I always go with the exact code used to run original analysis as opposed to adding features post-analysis.
# Future implementations would also take the population directory as input. 
#
# Run like: "sh pairwise_Fst_vcftools.sh /path/to/input.vcf /path/to/output.txt "

# Just to fill space, really only needs to be a tab 
printf "%s\t" "Pops" > $2
> $2.list
# Prints population file names across the top row
for pop in `ls pops`; do 
	printf "%s\t" $pop >> $2 
done
# Adds a newline
echo >> $2 

# Main for loop for making table, nested loop doing every pairwise comparison of populations defined by user.
# Essentially: "For each population, print its name then calculate and print pairwise for each other population."
# Currently relies on internal consistency in order. I don't think this ever ever not the case, but could be remedied by adding a population list input file.
for pop1 in `ls pops`; do
	# prints input file name
	printf "%s\t" $pop1 >> $2
	# internal loop for performing pairwise calculation
	for pop2 in `ls pops`; do
		# Only calculates pairwise Fst if populations aren't the same.
		if [ $pop1 != $pop2 ]
		then
			# Uses VCFtools to calculate pairwise Fst. However, the output isn't simply a number, and needs to be processed (this is why the output is used).
			vcftools --vcf $1 --weir-fst-pop pops/$pop1 --weir-fst-pop pops/$pop2 --out fst_output
			# Locates and outputs the weighted Fst value for the comparison.
			# Note that this doesn't purge this temporary file afterwards.
			output=$(grep 'weighted' ./fst_output.log | awk '{printf ("%s\t", $7)}')
			printf "%s\t" $output >> $2
			echo "$pop1 $pop2\t$output" >> $2.list
		else
			# Prints a "-" as a placeholder for comparisons of the same population.
			# Note that this could be some other value (e.g. nucleotide diversity, Fis) in an analysis.
			printf "%s\t" "-" >> $2
		fi
	done
	# adds a newline after each "pop1" iteration.
	echo >> $2
done

# Removed the .txt after file names. This can/should be changed if you named files differently (e.g. PopName.pop).
sed -i 's/.txt//g' $2
