# Variables
READS_DIR=reads
ALIGN_DIR=alignment
COUNTS_DIR=counts
GENOME_DIR=genome_index
DESIGN_FILE=design.csv
REFERENCE_GENOME_URL=https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_002853715.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED
REFERENCE_GENOME=GCF_002853715.1_genomic.fna
REFERENCE_GFF=GCF_002853715.1_genomic.gff
ENV_NAME=rnaseq_env
COUNT_MATRIX=counts_matrix.csv

# Phony targets
.PHONY: all usage download_reads align generate_counts count_matrix clean

# Usage instructions
usage:
	@echo "Usage:"
	@echo "  make all            - Run the full RNA-Seq pipeline"
	@echo "  make download_reads - Download RNA-Seq reads from SRA"
	@echo "  make count_matrix   - Generate final count matrix"
	@echo "  make clean          - Remove intermediate files"

# Run the entire workflow
all: download_reads download_reference align generate_counts count_matrix
	@echo "All steps completed successfully!"

# Download RNA-Seq reads using design.csv
download_reads:
	mkdir -p $(READS_DIR)
	while IFS=',' read -r sample accession; do \
		[ $$sample = "Sample_Name" ] || fasterq-dump --split-files -O $(READS_DIR) $$accession; \
	done < $(DESIGN_FILE)
	@echo "Downloaded RNA-Seq reads to $(READS_DIR)"

# Download the reference genome and annotation files
download_reference:
	mkdir -p $(GENOME_DIR)
	curl -L $(REFERENCE_GENOME_URL) -o $(GENOME_DIR)/genome.tar.gz
	tar -xzvf $(GENOME_DIR)/genome.tar.gz -C $(GENOME_DIR)
	mv $(GENOME_DIR)/$(REFERENCE_GENOME) $(GENOME_DIR)/genome.fasta
	mv $(GENOME_DIR)/$(REFERENCE_GFF) $(GENOME_DIR)/genome.gff
	@echo "Downloaded and extracted reference genome and annotations"

# Align the reads to the reference genome using HISAT2
align:
	mkdir -p $(ALIGN_DIR)
	hisat2-build $(GENOME_DIR)/genome.fasta $(GENOME_DIR)/genome_index
	while IFS=',' read -r sample accession; do \
		[ $$sample = "Sample_Name" ] || micromamba run -n $(ENV_NAME) hisat2 -x $(GENOME_DIR)/genome_index -U $(READS_DIR)/$$accession_1.fastq -S $(ALIGN_DIR)/$$sample.sam; \
	done < $(DESIGN_FILE)
	@echo "Alignment completed"

# Generate counts from aligned BAM files using featureCounts
generate_counts:
	mkdir -p $(COUNTS_DIR)
	while IFS=',' read -r sample accession; do \
		[ $$sample = "Sample_Name" ] || samtools view -bS $(ALIGN_DIR)/$$sample.sam > $(ALIGN_DIR)/$$sample.bam; \
		samtools sort $(ALIGN_DIR)/$$sample.bam -o $(ALIGN_DIR)/$$sample_sorted.bam; \
		featureCounts -a $(GENOME_DIR)/genome.gff -o $(COUNTS_DIR)/$$sample_counts.txt $(ALIGN_DIR)/$$sample_sorted.bam; \
	done < $(DESIGN_FILE)
	@echo "Counts generated for all samples"

# Generate the final count matrix from the count files
count_matrix: generate_counts
	Rscript generate_count_matrix.R $(COUNTS_DIR) $(COUNT_MATRIX)
	@echo "Final count matrix saved to $(COUNT_MATRIX)"

# Clean intermediate files
clean:
	rm -rf $(READS_DIR) $(ALIGN_DIR) $(COUNTS_DIR) $(GENOME_DIR) *.tar.gz
	@echo "Removed intermediate files"
