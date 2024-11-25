# Variables
GENOME_URLS=$(shell cut -d ',' -f2 design.csv)
GENOME_FILES=$(GENOME_URLS:%.fna=%_genome.fna)
READS_DIR=reads
REPORTS_DIR=reports
BAM_DIR=bam
SIM_BAM=$(BAM_DIR)/simulated_reads.bam
SRA_BAM=$(BAM_DIR)/sra_reads.bam
FILTERED_BAM=$(BAM_DIR)/filtered_reads.bam
FILTERED_STATS=$(BAM_DIR)/filtered_reads.stats.txt
VARIANTS_VCF=$(BAM_DIR)/variants.vcf
SNPEFF_DB=Escherichia_coli_str_k_12_substr_mg1655
ENV_NAME=menv

# Commands to print usage
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
	@echo "  make filter      - Filter BAM file and generate statistics for filtered reads"
	@echo "  make variants    - Call variants on aligned reads"
	@echo "  make annotate    - Annotate variants using snpEff"
	@echo "  make all         - Run all steps"

# Download the genome for each sample
genome: 
	$(GENOME_FILES): 
	curl -L -o $@ $(word 2,$(shell grep $(basename $@) design.csv))
	@echo "Downloaded and saved genome to $@"

# Simulate reads for each genome
simulate: 
	simulate_%: $(GENOME_FILES)
	mkdir -p $(READS_DIR)
	wgsim -N 10000 -1 150 -2 150 -e 0.01 -d 100 -s 50 $^ $(READS_DIR)/$(basename $@)_1.fastq $(READS_DIR)/$(basename $@)_2.fastq
	@echo "Simulated reads for $@ saved."

# Download SRA reads for each sample
download: 
	mkdir -p $(READS_DIR)
	fastq-dump -X 10000 --split-files -O $(READS_DIR) $(SRR)
	@echo "Downloaded reads to $(READS_DIR)"

# Trim reads for each genome
trim:
	mkdir -p $(REPORTS_DIR)
	micromamba run -n $(ENV_NAME) fastqc -q -o $(REPORTS_DIR) $(READS_DIR)/*.fastq
	micromamba run -n $(ENV_NAME) fastp --adapter_sequence=$(ADAPTER) --cut_tail \
	      -i $(R1) -I $(R2) -o $(TRIMMED_R1) -O $(TRIMMED_R2)
	micromamba run -n $(ENV_NAME) multiqc -o $(REPORTS_DIR) $(REPORTS_DIR)
	@echo "Trimmed reads saved to $(TRIMMED_R1) and $(TRIMMED_R2)"

# Create index for the reference genome for each sample
index: 
	mkdir -p $(GENOME_INDEX)
	micromamba run -n $(ENV_NAME) bwa index -p $(GENOME_INDEX)/$(basename $@) $@
	@echo "Index created for genome $@"

# Align reads for each genome
align: 
	mkdir -p $(BAM_DIR)
	micromamba run -n $(ENV_NAME) bwa mem $(GENOME_INDEX)/$(basename $@) $(TRIMMED_R1) $(TRIMMED_R2) | \
	samtools view -Sb - | samtools sort -o $(SRA_BAM)
	@echo "Aligned reads saved to $(SRA_BAM)"

# Generate BAM index
index_bam:
	micromamba run -n $(ENV_NAME) samtools index $(SRA_BAM)
	@echo "Index files created for BAM files"

# Generate alignment statistics for each sample
stats:
	micromamba run -n $(ENV_NAME) samtools flagstat $(SRA_BAM) > $(SRA_STATS)
	@echo "Alignment statistics saved to $(SRA_STATS)"

# Filter BAM and generate filtered statistics
filter:
	micromamba run -n $(ENV_NAME) samtools view -b -h -q 10 -f 0x2 -F 0x100 $(SRA_BAM) -o $(FILTERED_BAM)
	micromamba run -n $(ENV_NAME) samtools flagstat $(FILTERED_BAM) > $(FILTERED_STATS)
	@echo "Filtered BAM saved to $(FILTERED_BAM) and statistics saved."

# Call variants for each sample
variants:
	micromamba run -n $(ENV_NAME) bcftools mpileup -f $(GENOME_FILE) $(SRA_BAM) | \
	micromamba run -n $(ENV_NAME) bcftools call -mv -Ov -o $(VARIANTS_VCF)
	@echo "Variants called for $(SRA_BAM) and saved to $(VARIANTS_VCF)"

# Annotate variants using snpEff
annotate:
	micromamba run -n $(ENV_NAME) snpEff -Xmx4g $(SNPEFF_DB) $(VARIANTS_VCF) > $(SNPEFF_OUTPUT)
	@echo "Variants annotated with snpEff. Output saved."

# Run everything for each sample
all: 
	$(GENOME_FILES) simulate download trim index align index_bam stats filter variants annotate
	@echo "All steps complete for each genome."
