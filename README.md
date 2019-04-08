# alignment-viewer

## purpose
The purpose of this tool is to provide a simplified way to examine and explore coverage and mutations in sequence alignments


## implementation
The alignment-viewer tool was built around the metagenome dataset of Lake Washington sediment microbes used in Janet Matsen's thesis dissertation that can be found [here](https://digital.lib.washington.edu/researchworks/handle/1773/39975).
The tool is built as a jupyter notebook to allow visual eploration of different genes and regions. 


## inputs
The tool looks for one directory per metagenome/reference pair. 

Each data directory (`base_directory/referenceID/sampleID/`) must contain at least:
1. a `.coverage.bed` file containing coverage depth data along the genome
2. a `.vcf` file containing sequence variant data along the genome

Additionally, the base directory must contain two `.csv` reference files:
1. `aligned_isolate_genomes.csv` - contains the human-readable reference genomes and the names of the corresponding directories
2. `fastq_sample_lookup.csv` - contains metadata about the metagenomes included in the analysis

_see the [miscellaneous](https://github.com/blasks/alignment-viewer/tree/master/miscellaneous) directory for examples of the `.csv` reference files_

Details on generating these files and the directory organization can be found below, and in the documentation for the bash scripts in the [scripts](https://github.com/blasks/alignment-viewer/tree/master/scripts) directory.


## usage
**dependencies:**

The notebook depends on [pybedtools](https://daler.github.io/pybedtools/), among other packages. I found some dependency conflicts on my machine between samtools and pybedtools. A `.yml` file with the specifications for the environment that worked for me can be found in the [miscellaneous](https://github.com/blasks/alignment-viewer/tree/master/miscellaneous) directory.

**data organization:**

Organized data directories were generated using the [generate_align_tasklist.sh](https://github.com/blasks/alignment-viewer/blob/master/scripts/generate_align_tasklist.sh) script that can be found in the [scripts](https://github.com/blasks/alignment-viewer/tree/master/scripts) directory.
The directories should be organized in a nested fashion as follows:

```
base_directory
|   aligned_isolate_genomes.csv
|   fastq_sample_lookup.csv
|
└───referenceID (directory named with reference genome ID)
    |
    └───sampleID (directory named with metagenome sample ID)
            sampleID_referenceID.coverage.bed
            sampleID_referenceID.vcf
            sampleID_referenceID.extension (other data files, such as .bam)
```

**file generation:**

Alignment output files (`.coverage.bed` & `.vcf`) were generated using the [align_coverage_call.sh](https://github.com/blasks/alignment-viewer/blob/master/scripts/align_coverage_call.sh) script that can be found in the [scripts](https://github.com/blasks/alignment-viewer/tree/master/scripts) directory. 
This workflow runs the following steps:

1. uses bwa to run a burrows-wheeler alignment, outputting a `.sam` file
2. uses samtools to convert the `.sam` to a `.bam` file
3. uses samtools to sort the `.bam` file, outputting a `.sorted.bam` file
4. uses samtools to index the `.sorted.bam` file
5. uses bedtools to calculate coverage across the alignment, outputting a `.coverage.bed` file
6. uses bcftools to calculate sequence variants, outputting a `.vcf` file

**running the notebook:**

1. Import packages
2. Point the notebook to the correct base directory for the data by editing the `base_directory` global variable
3. Run the code cell to define the core functionality
4. Render the user input widgets and use them to identify the region of interest
5. Check that the .coverage.bed and .vcf files exist in the expected location for the data of interest
6. Use the notebook to render the visualization and explore the data of interest


