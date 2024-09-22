#!/bin/bash

# To make the script executable input: chmod +x ncbi_genome_download.sh
# To run input: ./ncbi_genome_download.sh

# Define the NCBI download URL (click download on the data, check all of the boxes, go to URL header on specific dataset, copy URL here)
NCBI_URL="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000146045.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

# Define output directory
OUTPUT_DIR="ncbi_genome_download"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Download the dataset zip file
echo "Downloading genome and annotation files from NCBI..."
curl -L -o $OUTPUT_DIR/genome_dataset.zip "$NCBI_URL"

# Unzip the downloaded dataset
echo "Unzipping the dataset..."
unzip -o $OUTPUT_DIR/genome_dataset.zip -d $OUTPUT_DIR

# Find the GFF3 file
GFF_FILE=$(find $OUTPUT_DIR -type f -name "*.gff*" | head -n 1)

# Check if the GFF file exists
if [[ -z "$GFF_FILE" ]]; then
  echo "Error: No GFF file found in the downloaded dataset."
  exit 1
fi

# Create separate GFF files for 'gene' and 'cds' labels
# Only extract lines where the 3rd column exactly equals "gene"
echo "Extracting gene annotations..."
awk '$3 == "gene"' $GFF_FILE > $OUTPUT_DIR/gene.gff

# Only extract lines where the 3rd column exactly equals "CDS"
echo "Extracting CDS annotations..."
awk '$3 == "CDS"' $GFF_FILE > $OUTPUT_DIR/cds.gff

# Check if the gene.gff and cds.gff files were successfully created and populated
if [[ -s $OUTPUT_DIR/gene.gff ]]; then
  echo "gene.gff created successfully with gene information."
else
  echo "Warning: gene.gff is empty or could not find 'gene' entries."
fi

if [[ -s $OUTPUT_DIR/cds.gff ]]; then
  echo "cds.gff created successfully with CDS information."
else
  echo "Warning: cds.gff is empty or could not find 'CDS' entries."
fi

echo "All tasks completed. Files saved in $OUTPUT_DIR."

