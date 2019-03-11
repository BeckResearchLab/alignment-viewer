#!/bin/bash
# for Dave - what is this? ^^

# this script aligns all fastq files in an input directory tree against
# all reference fasta genome files in an input directory tree
# and outputs all alignments to a third input directory location

# usage: $0 <fastq_root_directory> <genome_fasta_root_directory> <output_directory>
# example: ./all_by_all_align.sh ../scratchspace/fastq_data/ ../scratchspace/reference_fastas/ ../scratchspace/test_output/

if [ "$1" != "" ]; then
		fastq_dir="$1"
		echo "using '$fastq_dir' as the fastq root directory"
else
		echo "usage: $0 <fastq_root_directory> <genome_fasta_root_directory> <output_directory>"
		echo "example: ./all_by_all_align.sh ../scratchspace/fastq_data/ ../scratchspace/reference_fastas/ ../scratchspace/test_output/"
		exit 1
fi
if [ "$2" != "" ]; then
		fasta_dir="$2"
		echo "using '$fasta_dir' as the reference genome fasta directory"
else
		echo "usage: $0 <fastq_root_directory> <genome_fasta_root_directory> <output_directory>"
		echo "example: ./all_by_all_align.sh ../scratchspace/fastq_data/ ../scratchspace/reference_fastas/ ../scratchspace/test_output/"
		exit 1
fi
if [ "$3" != "" ]; then
		output_dir="$3"
		echo "using '$output_dir' as the bam output directory"
else
		echo "usage: $0 <fastq_root_directory> <genome_fasta_root_directory> <output_directory>"
		echo "example: ./all_by_all_align.sh ../scratchspace/fastq_data/ ../scratchspace/reference_fastas/ ../scratchspace/test_output/"
		exit 1
fi

# check that all reference sequences are indexed
# question: how to do full tree search?

for file in $fasta_dir*.fsa
do
	if [ ! -e "$file.bwt" ]; then
		# echo "indexing reference genome fasta files"
		bwa index $file
	fi
done


# run all alignments
for fastq in $fastq_dir*genome*/Raw_Data/*.fastq.gz
# have to finagle here b/c of weird file organization of data on waffle
do
	for reference in $fasta_dir*.fsa
	do
		output_handle=$output_dir$(basename $fastq .fastq.gz)_$(basename $reference .fsa)
		if [ ! -e "$output_handle.sorted.bam" ]; then
			# echo "aligning $(basename $fastq) to reference $(basename $reference)"
			echo "bwa mem -M -t 28 $reference $fastq > $output_handle.sam"
			# echo "converting sam to bam"
			echo "samtools view -bT $reference -o $output_handle.unsorted.bam $output_handle.sam"
			# the -T flag here is for a fasta reference file - still need to figure out what this is doing
			# rm $output_handle.sam
			# echo "sorting bam"
			echo "samtools sort -o $output_handle.sorted.bam $output_handle.unsorted.bam"
			# rm $output_handle.unsorted.bam
			# echo "indexing bam"
			echo "samtools index $output_handle.sorted.bam"
		fi
	done
done
