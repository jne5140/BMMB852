Create directory for SRA data
mkdir -p sra_data
echo "Downloading SRA data from ENA for S. aureus (Accession: SRR2584863)..."

Download the FASTQ file for S. aureus data from ENA (SRA via FTP)
wget -P sra_data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz

Check if the download was successful
if [ $? -ne 0 ]; then
    echo "Error downloading the FASTQ file."
    exit 1
fi

echo "FASTQ file downloaded successfully."
echo "Downloaded files:"
ls -lh sra_data/

Perform FastQC quality control analysis
echo "Running FastQC on the downloaded data..."
fastqc sra_data/SRR2584863_2.fastq.gz -o ./ 

Check if FastQC was successful
if [ $? -ne 0 ]; then
    echo "Error running FastQC."
    exit 1
fi

List FastQC results
echo "FastQC analysis completed successfully."
echo "FastQC results:"
ls -lh *.html *.zip

Open FastQC for visualization
echo "Opening FastQC ..."
fastqc

Open FastQC report manually
