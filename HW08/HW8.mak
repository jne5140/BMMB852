# Variables
GENOME_URL="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000005845.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
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
BAM_DIR=bam
SIM_BAM=$(BAM_DIR)/simulated_reads.bam
SRA_BAM=$(BAM_DIR)/sra_reads.bam
GENOME_INDEX=genome_index
SIM_BAM=$(BAM_DIR)/simulated_reads.bam
SRA_BAM=$(BAM_DIR)/sra_reads.bam
SIM_STATS=$(BAM_DIR)/simulated_reads.stats.txt
SRA_STATS=$(BAM_DIR)/sra_reads.stats.txt

# Directory
usage:
	@echo "Available commands:"
	@echo "  make genome      - Download the genome file"
	@echo "  make simulate    - Simulate reads for the genome"
	@echo "  make download    - Download reads from SRA"
	@echo "  make trim        - Trim reads with fastp"
	@echo "  make index       - Create an index for the reference genome"
	@echo "  make align       - Align the reads to the reference genome"
	@echo "  make index_bam   - Create BAM index files"
	@echo "  make stats       - Generate alignment statistics"
	@echo "  make clean       - Remove all generated files"

# Download the genome
genome:
	curl -L -o genome.zip $(GENOME_URL)
	unzip genome.zip -d genome_dir
	mv genome_dir/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna $(GENOME_FILE)
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
	mkdir -p $(READS_DIR) $(REPORTS_DIR)

	@echo "Running FastQC on raw reads..."
	micromamba run -n $(ENV_NAME) fastqc -q -o $(REPORTS_DIR) $(R1) $(R2)

	@echo "Trimming reads with fastp..."
	micromamba run fastp --adapter_sequence=$(ADAPTER) --cut_tail \
	      -i $(R1) -I $(R2) -o $(TRIMMED_R1) -O $(TRIMMED_R2)

	@echo "Running FastQC on trimmed reads..."
	micromamba run fastqc -q -o $(REPORTS_DIR) $(TRIMMED_R1) $(TRIMMED_R2)

	@echo "Running MultiQC on reports directory..."
	micromamba run multiqc -o $(REPORTS_DIR) $(REPORTS_DIR)

	@echo "Trimmed reads saved to $(TRIMMED_R1) and $(TRIMMED_R2)"

# Create an index for the reference genome
index:
	mkdir -p $(GENOME_INDEX)
	micromamba run bwa index -p $(GENOME_INDEX)/genome $(GENOME_FILE)
	@echo "Index created in $(GENOME_INDEX)"

# Align reads to the reference genome
align: 
	mkdir -p $(BAM_DIR)
	micromamba run bwa mem $(GENOME_INDEX)/genome $(TRIMMED_R1) $(TRIMMED_R2) | samtools view -Sb - | samtools sort -o $(SRA_BAM)
	micromamba run bwa mem $(GENOME_INDEX)/genome $(R1) $(R2) | samtools view -Sb - | samtools sort -o $(SIM_BAM)
	@echo "Aligned reads saved to $(SRA_BAM) and $(SIM_BAM)"

# Create BAM index files
index_bam:
	micromamba run samtools index $(SRA_BAM)
	micromamba run samtools index $(SIM_BAM)
	@echo "Index files created for BAM files"

# Clean up all files
clean:
	rm -f $(GENOME_FILE)
	rm -rf $(READS_DIR)
	rm -rf $(REPORTS_DIR)
	rm -rf $(BAM_DIR)
	rm -rf $(GENOME_INDEX)*
	@echo "Cleaned up all files"

# Generate alignment statistics
stats:
	micromamba run samtools flagstat $(SIM_BAM) > $(SIM_STATS)
	@echo "Alignment statistics for simulated reads saved to $(SIM_STATS)"
	micromamba run samtools flagstat $(SRA_BAM) > $(SRA_STATS)
	@echo "Alignment statistics for SRA reads saved to $(SRA_STATS)"

# Run everything
all: genome simulate download trim index align index_bam stats
	@echo "All steps are done!"
