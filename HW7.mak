# Variables
GENOME_URL="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000013425.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
GENOME_FILE=genome.fna
READS_DIR=reads
REPORTS_DIR=reports
SRR=SRR519926
N=10000
R1=$(READS_DIR)/$(SRR)_1.fastq
R2=$(READS_DIR)/$(SRR)_2.fastq
TRIMMED_R1=$(READS_DIR)/$(SRR)_1.trimmed.fastq
TRIMMED_R2=$(READS_DIR)/$(SRR)_2.trimmed.fastq
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
ENV_NAME=menv

# Makefile directory
usage:
	@echo "Available commands:"
	@echo "  make genome      - Download a simple genome file"
	@echo "  make simulate    - Simulate reads from the genome"
	@echo "  make download    - Download some reads from SRA"
	@echo "  make trim        - Trim the downloaded reads (removing low-quality bases)"
	@echo "  make qc          - Check the quality of the reads"
	@echo "  make clean       - Remove all generated files"

# Download the genome
genome:
	curl -L -o genome.zip $(GENOME_URL)
	unzip genome.zip -d genome_dir
	mv genome_dir/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna $(GENOME_FILE)
	rm -rf genome_dir genome.zip
	@echo "Downloaded and extracted genome to $(GENOME_FILE)"

# Simulate reads for the genome
simulate:
	mkdir -p $(READS_DIR)
	wgsim -N $(N) -1 150 -2 150 -e 0.01 -d 100 -s 50 $(GENOME_FILE) $(R1) $(R2)
	@echo "Simulated reads saved to $(R1) and $(R2)"

# Download reads
download:
	mkdir -p $(READS_DIR)
	fastq-dump -X $(N) --split-files -O $(READS_DIR) $(SRR)
	@echo "Downloaded reads to $(READS_DIR)"


# Trim reads with fastp
trim:
	@echo "Setting up environment..."
	micromamba run -n $(ENV_NAME) bash -c
	mkdir -p $(READS_DIR) $(REPORTS_DIR)

	@echo "Running FastQC on raw reads..."
	micromamba run -n $(ENV_NAME) fastqc -q -o $(REPORTS_DIR) $(R1) $(R2)

	@echo "Trimming reads with fastp..."
	micromamba run -n $(ENV_NAME) fastp --adapter_sequence=$(ADAPTER) --cut_tail \
	      -i $(R1) -I $(R2) -o $(TRIMMED_R1) -O $(TRIMMED_R2)

	@echo "Running FastQC on trimmed reads..."
	micromamba run -n $(ENV_NAME) fastqc -q -o $(REPORTS_DIR) $(TRIMMED_R1) $(TRIMMED_R2)

	@echo "Running MultiQC on reports directory..."
	micromamba run -n $(ENV_NAME) multiqc -o $(REPORTS_DIR) $(REPORTS_DIR)

	@echo "Trimmed reads saved to $(TRIMMED_R1) and $(TRIMMED_R2)"

# Run quality control using FastQC
qc:
	mkdir -p $(REPORTS_DIR)
	fastqc $(R1) $(R2) $(TRIMMED_R1) $(TRIMMED_R2) -o $(REPORTS_DIR)
	@echo "FastQC reports are in $(REPORTS_DIR)"

# Clean up all files
clean:
	rm -f $(GENOME_FILE)
	rm -rf $(READS_DIR)
	rm -rf $(REPORTS_DIR)
	@echo "Cleaned up all files"

# Run everything
all: genome simulate download trim qc
	@echo "All steps are done!"