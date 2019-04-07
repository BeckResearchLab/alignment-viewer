# alignment-viewer

### purpose
The purpose of this tool is to provide a simplified way to examine and explore coverage and mutations in sequence alignments

### implementation
The alignment-viewer tool was built around the metagenome dataset of Lake Washington sediment microbes used in Janet Matsen's thesis dissertation that can be found [here](https://digital.lib.washington.edu/researchworks/handle/1773/39975).
The tool is built as a jupyter notebook to allow visual eploration of different genes and regions. 

### inputs
The tool looks for one directory per metagenome/reference pair, containing at least:
1) a `.coverage.bed` file containing coverage depth data along the genome
2) a `.vcf` file containing sequence variant data along the genome

Details on generating these files and the directory organization can be found below, and in the documentation for the bash scripts in the [scripts](https://github.com/blasks/alignment-viewer/tree/master/scripts) directory

### usage
**dependencies:**
The notebook depends on [pybedtools](https://daler.github.io/pybedtools/), among other packages. I found some dependency conflicts on my machine between samtools and pybedtools. A `.yml` file with the specifications for the environment that worked for me can be found in the [miscellaneous](https://github.com/blasks/alignment-viewer/tree/master/miscellaneous) directory.

**file generation:**
Alignment output files (`.coverage.bed` & `.vcf`) were generated using the [align_coverage_call.sh](https://github.com/blasks/alignment-viewer/blob/master/scripts/align_coverage_call.sh) script that can be found in the [scripts](https://github.com/blasks/alignment-viewer/tree/master/scripts) directory. 
This workflow runs the following steps:
1) uses bwa to run a burrows-wheeler alignment, outputting a `.sam` file
2) uses samtools to convert the `.sam` to a `.bam` file
3) uses samtools to sort the `.bam` file, outputting a `.sorted.bam` file
4) uses samtools to index the `.sorted.bam` file
5) uses bedtools to calculate coverage across the alignment, outputting a `.coverage.bed` file
6) uses bcftools to calculate sequence variants, outputting a `.vcf` file

**running the notebook:**
- Point the notebook to the correct base directory for the data by changing the `base_directory` variable in the second cell
- Run the core functionality and user input cells
- Use the user input widgets to specify the data of interest
- Run the final cell to render the visualization to explore the data of interest

