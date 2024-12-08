## RNA-Seq Data Processing Makefile
This README provides a simple guide on how to use the Makefile to automate RNA-Seq data processing, alignment, and count matrix generation.

#### Prerequisites
To use this pipeline, you need the following tools:

**fasterq-dump:** For downloading RNA-Seq reads from SRA.

**hisat2:** For aligning RNA-Seq reads to the reference genome.

**samtools:** For handling alignment files.

**featureCounts:** For counting reads aligned to genes.

**micromamba or conda:** For managing environments and running tools.

**R:** To generate the final count matrix.

**curl:** For downloading the reference genome and annotations.

#### Running the Workflow Prepare Files:

Place the Makefile and design.csv (which contains sample accessions) in the same directory.

Ensure the reference genome URL is correctly set in the Makefile.

#### Install Dependencies:

Make sure all tools are installed or use micromamba to manage dependencies.

#### Run the Pipeline: 

To run the entire RNA-Seq pipeline, use:

	make all
This will:

Download RNA-Seq reads.

Download the reference genome and annotations.

Align the reads.

Generate the count matrix.

#### Makefile Targets
**usage:** Displays usage instructions.

**all:** Runs the entire RNA-Seq pipeline (download reads, align, count).

**download_reads:** Downloads RNA-Seq reads from SRA using fasterq-dump.

**download_reference:** Downloads the reference genome and annotations from NCBI.

**align:** Aligns RNA-Seq reads to the reference genome using HISAT2.

**generate_counts:** Converts aligned files (SAM to BAM) and counts reads using featureCounts.

**count_matrix:** Generates a final count matrix from the count files using an R script.

**clean:** Removes intermediate files to clean up the directory.

#### Notes
The pipeline handles multiple samples defined in the design.csv file.
The reference genome and annotations are downloaded from NCBI using curl.
Adjust the generate_count_matrix.R script to customize how count data is processed.

#### Troubleshooting

Ensure that all required tools are installed and available.

Verify the design.csv file and SRA accessions are correct.

Check if the reference genome files are downloaded properly.