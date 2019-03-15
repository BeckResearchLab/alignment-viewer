#!/bin/bash

# This script takes a single .fastq.gz file and a single .fsa fasta file,
# as well as an output directory as inputs.
# It then aligns the fastq reads to the fasta referernce with bwa mem,
# converts the output to .bam, sorts, and indexes with samtools,
# calculates genome coverage with bedtools, and finally makes variant calls
# based on the reference fasta using bcftools

# "usage: $0 <fastq_file> <reference_fasta_file> <output_directory>"
# "example: ./align_coverage_calls.sh ../scratchspace/fastq_data/test_Metagenome_1/Raw_Data/sample_1.fastq.gz ../scratchspace/joes_reference_genomes/fasta/NC_012968.fsa ../scratchspace/test_output/"

if [ $# -eq 3 ]; then
		fastq=$1
		reference=$2
		dir_out=$3
		output_handle=$dir_out$(basename $fastq .fastq.gz)_$(basename $reference .fsa)
else
		echo "usage: $0 <fastq_file> <reference_fasta_file> <output_directory>"
		echo "example: ./align_coverage_calls.sh ../scratchspace/fastq_data/test_Metagenome_1/Raw_Data/sample_1.fastq.gz ../scratchspace/joes_reference_genomes/fasta/NC_012968.fsa ../scratchspace/test_output/"
		exit 1
fi

echo "aligning $fastq to reference file $reference"
bwa mem -M -t 28 $reference $fastq > $output_handle.sam

echo "converting sam to bam"
samtools view -bT $reference -o $output_handle.unsorted.bam -@28 $output_handle.sam

echo "sorting bam"
samtools sort -o $output_handle.sorted.bam -@28 $output_handle.unsorted.bam

echo "indexing bam"
samtools index $output_handle.sorted.bam

echo "calculating genome coverage"
bedtools genomecov -bga -ibam $output_handle.sorted.bam > $output_handle.coverage.bed

echo "calling genome variants"
bcftools mpileup -f $reference $output_handle.sorted.bam | bcftools call -cv --ploidy 1 -o $output_handle.vcf
