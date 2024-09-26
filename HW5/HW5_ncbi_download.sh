#!/bin/bash

# Define the NCBI download URL
NCBI_URL="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

# Define output directories
OUTPUT_DIR="ncbi_genome_download"
FASTQ_OUTPUT_DIR="$OUTPUT_DIR/fastq_output"
mkdir -p $OUTPUT_DIR $FASTQ_OUTPUT_DIR

# do not change anything below this line
#------------------------------------------------------------------
# Download and unzip dataset
echo "Downloading genome files..."
curl -L -o $OUTPUT_DIR/genome_dataset.zip "$NCBI_URL"
echo "Unzipping dataset..."
unzip -o $OUTPUT_DIR/genome_dataset.zip -d $OUTPUT_DIR

# Locate the FASTA file
FASTA_FILE=$(find $OUTPUT_DIR -type f -name "*.fna" | head -n 1)

if [ -f "$FASTA_FILE" ]; then
    GENOME_SIZE=$(grep -v "^>" "$FASTA_FILE" | wc -c)
    CHROMOSOME_NAME=$(grep "^>" "$FASTA_FILE" | head -n 1 | sed 's/^>//')
    CHROMOSOME_LENGTH=$(awk '/^>/ {if (seq) print length(seq); seq=""; print; next} {seq=seq $0} END {print length(seq)}' "$FASTA_FILE" | tail -n 1)

    echo "FASTA file located: $FASTA_FILE"
    echo "Genome size: $GENOME_SIZE bp"
    echo "Chromosome: $CHROMOSOME_NAME, Length: $CHROMOSOME_LENGTH bp"

    # Generate simulated FASTQ files
    COVERAGE=10
    READ_LENGTH=150
    echo "Generating simulated FASTQ files with $COVERAGEx coverage..."
    wgsim -N $((COVERAGE * GENOME_SIZE / READ_LENGTH)) -1 $READ_LENGTH -2 $READ_LENGTH -e 0.01 -d 100 -s 50 "$FASTA_FILE" "$FASTQ_OUTPUT_DIR/simulated_R1.fastq" "$FASTQ_OUTPUT_DIR/simulated_R2.fastq"

    for FASTQ_FILE in "$FASTQ_OUTPUT_DIR"/simulated_R{1,2}.fastq; do
        if [ -f "$FASTQ_FILE" ]; then
            AVG_READ_LENGTH=$(awk 'NR % 4 == 2 {total += length($0); count++} END {if (count > 0) print total / count; else print 0}' "$FASTQ_FILE")
            FASTQ_SIZE=$(du -h "$FASTQ_FILE" | cut -f1)

            gzip "$FASTQ_FILE"
            COMPRESSED_FILE="$FASTQ_FILE.gz"
            COMPRESSED_SIZE=$(du -h "$COMPRESSED_FILE" | cut -f1)

            echo "FASTQ file: $COMPRESSED_FILE"
            echo "Avg read length: $AVG_READ_LENGTH bp"
            echo "Original size: $FASTQ_SIZE"
            echo "Compressed size: $COMPRESSED_SIZE"
        else
            echo "FASTQ file not found!"
        fi
    done

    # Check the total number of reads generated
    TOTAL_READS_R1=$(grep -c "^@" "$FASTQ_OUTPUT_DIR/simulated_R1.fastq.gz")
    TOTAL_READS_R2=$(grep -c "^@" "$FASTQ_OUTPUT_DIR/simulated_R2.fastq.gz")
    TOTAL_READS=$((TOTAL_READS_R1 + TOTAL_READS_R2))

    echo "Total reads generated: $TOTAL_READS"
    echo "Total original size: $(du -h "$FASTQ_OUTPUT_DIR" | cut -f1)"
    echo "Total compressed size: $(du -h "$FASTQ_OUTPUT_DIR" | cut -f1)"
else
    echo "FASTA file not found!"
fi

echo "All tasks completed."