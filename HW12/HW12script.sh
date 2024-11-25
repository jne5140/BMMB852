#!/bin/bash

# Path to your CSV file with Sample Names and Genome URLs
CSV_FILE="design.csv"

# Read the CSV file and loop through each line (skipping the header)
tail -n +2 $CSV_FILE | while IFS=, read -r SAMPLE GENOME_URL; do
    echo "Processing $SAMPLE with URL: $GENOME_URL"
    
    # Pass the current sample's genome URL to make
    MAKE_CMD="make -f HW12.mak GENOME_URL=\"$GENOME_URL\""

    # Execute the "genome" step to download the genome file for each sample
    echo "Running make genome for $SAMPLE"
    $MAKE_CMD genome

    # Optionally, run other steps (e.g., simulate, download, etc.)
    # For example, to run "simulate" after downloading the genome:
    echo "Running make simulate for $SAMPLE"
    $MAKE_CMD simulate

    echo "Processing for $SAMPLE complete!"
done

echo "All samples processed."
