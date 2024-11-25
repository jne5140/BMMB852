## How to run HW12.md (Variant calling workflow Makefile)

This repository provides a Makefile for a complete workflow to download a reference genome, process sequencing data, and perform variant calling and annotation.

### Prerequisites

**Required Software**

Ensure the following software is installed and accessible in your menv environment:

BWA

FastQC

MultiQC

samtools

fastp

SRA Toolkit

snpEff

bcftools

**Required Environment**

menv

Activate environment before running the Makefile:

	micromamba activate menv  

 
### File/Directory	Description

**Makefile**	
- Contains all commands for running the workflow

**genome.fna**	
- Reference genome in FASTA format

**reads/**	
- Directory for storing raw and processed reads

**bam/**	
- Directory for storing BAM files and statistics

**reports/**	
- Directory for storing QC and MultiQC reports


### Step-by-Step Workflow

Run the workflow by using the following commands.

1. Download the Reference Genome

		make genome  

	This downloads and extracts the reference genome, saving it as genome.fna.

2. Simulate Reads

		make simulate  
	Generates simulated reads using wgsim for testing the workflow.

3. Download Reads from SRA

		make download  
	Downloads a subset of reads from the specified SRA accession (SRR).

4. Trim Reads

		make trim  
	Trims adapters and low-quality bases from raw reads and generates QC reports.

5. Index the Reference Genome

		make index  
	Creates a BWA index for the reference genome.

6. Align Reads

		make align  
	Aligns trimmed reads to the reference genome and outputs a sorted BAM file.

7. Generate BAM Index

		make index_bam  
	Creates BAM index files required for downstream analysis.

8. Generate Alignment Statistics

		make stats  
	Generates alignment statistics for the BAM file.

9. Filter BAM Files

		make filter  
	Filters the BAM file for high-quality, paired-end reads.

10. Call Variants

		make variants  
	Calls variants from the BAM file and saves them in VCF format.

11. Annotate Variants

		make annotate  
	Annotates variants using snpEff.

12. Run All Steps

		make all  
	Executes the entire workflow, from downloading the genome to annotating variants.

### Outputs
**Genome File:** genome.fna

**Trimmed Reads:** Stored in reads/

**BAM Files:** Stored in bam/

**Variant File:** bam/variants.vcf

**Annotated Variant File:** bam/
variants_snpeff.vcf











