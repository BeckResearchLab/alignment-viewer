#!/bin/bash

# This script is meant to be used in conjunction with GNU parallel to run
# many align_coverage_calls scripts at once for different read-reference pairs.
# The inputs are a list of .fastq.gz raw read files, a list of .fsa fastq
# reference files, and a root output directory.
# The script first checks to make sure that all the fasta files have been
# indexed, and indexes them if it does not find a .bwt file.
# Next the script generates a tree of directories in the root output directory.
# The output directories are organized first by reference genome, then by
# fastq raw reads.
# Finally the script generates a task list, echoing the ./align_coverage_calls
# script for each of the fastq-reference pairs.
# This output can be re-directed to a text file, and GNU parallel can then be
# used to execute each of the calls in parallel.

# "usage: $0 <fastq_directory> <reference_fasta_directory> <root_output_directory>"
# "example: ./generate_align_tasklist.sh ../scratchspace/fastq_data/ ../scratchspace/joes_reference_genomes/fasta/ ../scratchspace/test_output/"

# check arguments and assign variables
if [ $# -eq 3 ]; then
		fastq_dir=$1
		ref_dir=$2
		out_dir=$3
else
		echo "usage: $0 <fastq_directory> <reference_fasta_directory> <root_output_directory>"
		echo "example: ./generate_align_tasklist.sh ../scratchspace/fastq_data/ ../scratchspace/joes_reference_genomes/fasta/ ../scratchspace/test_output/"
		exit 1
fi

# check that all reference sequences are indexed
# if there is no .bwt file, use bwa to index the reference file
for reference in $ref_dir*.fsa
do
	if [ ! -e "$reference.bwt" ]; then
		bwa index $reference
	fi
done

# make output directory for each fastq-reference pair, and echo script call
for reference in $ref_dir*.fsa
do
	org_dir=$out_dir$(basename $reference .fsa)/
	mkdir $org_dir
	for fastq in $fastq_dir*Metagenome*/Raw_Data/*.fastq.gz
	do
		output_dir=$org_dir$(basename $fastq .fastq.gz)/
		mkdir $output_dir
		echo "$fastq $reference $output_dir"
	done
done
