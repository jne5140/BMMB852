## BMMB 852 HW#7

**Jessica Eckard**


This assignment requires writing a Makefile and a markdown report. Use the recommended practices in the ⭐ AI-powered Makefiles document to guide your work.

In the previous assignments you were asked to write a scripts to simulate reads. Then in the second assignment you were asked to obtain and trim reads for a realistic dataset.

Merge both scripts into a single Makefile. The file should have the following targets:

usage - Print the help

genome - Download the genome

simulate - Simulate reads for the genome

download - Download reads from SRA

trim - Trim reads

clean - Remove all files

Add additional targets as needed.

Have your makefile generate fastqc reports upon downloading and trimming data.

The report should explain how the makefile works and what each target does.

Use variables to set up the target of the reads. Here are some recommendations:

R1 - The first read file 

R2 - The second read file 

GENOME - The genome file 

N - The number of reads to simulate or download

_________________________

Makefile is a way to create an automated list of tasks where each individual task is sectioned and therefore can be individually tested and modified. 
### **Variables**

The variables section defines filenames that can be customizeable. 

- **GENOME_URL:** URL for downloading the genome


- **GENOME_FILE:** Name of the output genome file


- **READS_DIR:** Directory for storing read files


- **REPORTS_DIR:** Directory for storing reports


- **SRR:** The accession number for the SRA dataset


- **N:** Number of reads to simulate


- **R1 and R2:** Filenames for simulated reads


- **TRIMMED_R1 and TRIMMED_R2:** Filenames for trimmed reads


- **ADAPTER:** The adapter sequence used for the trimming process

- **ENV_NAME:** Name of micromamba environment used
### **usage:**

A quick reference guide on the sections in the makefile. 

- **genome:** downloads the genome data and unzips it into a directory. Moves the FASTA file to the current directory and cleans up temporary files.

- **simulate:** Creates the directory for storing reads. Uses wgsim to simulate reads from the downloaded genome file. Generates two FASTQ files.

- **download:** Creates the reads directory. Uses fastq-dump to download reads from the SRA using the accession number. Saves reads in the READS_DIR.

- **trim:** Runs FastQC on the raw reads to check their quality. Saves the reports in the REPORTS_DIR. Uses fastp to trim low-quality bases and remove adapter sequences. Runs FastQC again on the trimmed reads to assess their quality and generate reports. Runs multiqc to compile the FastQC reports into a single summary report.

- **qc:** Creates the reports directory if it doesn't exist. Runs FastQC on both raw and trimmed reads to generate quality reports.

- **clean:** Removes the genome file, reads directory, and reports directory.

- **all:** Executes all the primary steps above (genome download, read simulation, downloading reads, trimming, and quality control) in order.
