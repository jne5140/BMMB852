# BMMB 852 HW # 5

#### Jessica Eckard

This assignment requires writing a script and a markdown report.
Both files should be committed and published in your repository.

Tips

Some tools may produce lots of unnecessary output. Only include the relevant bits.
Select a smaller genome (< 100Mb) or a subset (a single chromosome) of a large genome to make the analysis more manageable. For example, if you want to select the human genome, use only chromosome 22.

You may reuse code written for previous assignments.

Your script should be runnable and produce the output required to answer each question. 

You don't have to make the script report the answers eloquently. The script job is to generate the type of data that you can use to answer the questions.

1.	Select a genome, then download the corresponding FASTA file.

I used an S. aureus genome from NCBI Datasets.

> NCBI RefSeq: GCF_000013425.1

Report:

•	The size of the file

•	The total size of the genome

•	The number of chromosomes in the genome

•	The name(id) and length of each chromosome in the genome.

> FASTA file located: ncbi_genome_download/ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna
FASTA file size: 2.8M
Genome size (bp): 2856629
Number of chromosomes: 1
Chromosome name and length (first):
Chromosome: NC_007795.1 Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete genome, Length: 2821361 bp`


2. Generate a simulated FASTQ output for a sequencing instrument of your choice. Set the parameters so that your target coverage is 10x.

Report:

•	How many reads have you generated?

•	What is the average read length?

•	How big are the FASTQ files?

•	Compress the files and report how much space it saves.

> FASTQ file: ncbi_genome_download/fastq_output/simulated_R2.fastq.gz
Avg read length: 150 bp
Original size: 64M
Compressed size: 12M
Total reads generated: 723`

Discuss whether you could get the same coverage with different parameter settings (read length vs. read number).

Coverage = (number of reads  x  length of each read in base pairs)  /  size of the genome in base pairs

Because of this equation, it is possible to get the same coverage with different parameter settings. There could be a larger number of reads, but then the length of each read would need to be smaller and vice versa. 

3. How much data would be generated when covering the yeast,  the drosophila and the human genomes?

 Assumptions:

Coverage = 10x

Read length = 150 bp


Yeast

> NCBI RefSeq: GCF_000146045.2

Saccharomyces cerevisiae S288C genome assembly R64 - NCBI - NLM (nih.gov)
FASTQ file: ncbi_genome_download/fastq_output/simulated_R2.fastq.gz
Avg read length: 150 bp
Original size: 208M
Compressed size: 21M
Total reads generated: 683`

Drosophila

> NCBI RefSeq: GCF_000001215.4

Drosophila melanogaster genome assembly release 6 plus ISO1 MT - NCBI - NLM (nih.gov)
FASTQ file: ncbi_genome_download/fastq_output/simulated_R2.fastq.gz
Avg read length: 150 bp
Original size: 1.5G
Compressed size: 146M
Total reads generated: 8308`

Human

> NCBI RefSeq: GCF_000001405.40

Homo sapiens genome assembly GRCh38.p14 - NCBI - NLM (nih.gov)
FASTQ file: ncbi_genome_download/fastq_output/simulated_R2.fastq.gz
Avg read length: 150 bp
Original size: 7.0G
Compressed size: 671M
Total reads generated: 40406
`

Now imagine that instead of your genome, each instrument generated reads that cover the yeast, the drosophila and the human genome at 30x coverage. You don't have to run the tool!

Using the information for the genome you've studied in the previous points, for each of the organisms estimate the size of the FASTA file that holds the genome, the number of FASTQ reads needed for 30x, and the size of the FASTQ files before and after compression.

Assumptions

Read length: 150 bp

Coverage: 30x

Compression ratio: 10x

FASTA file size is roughly proportional to genome size, assume these are equal

100 bytes per read

Yeast

Genome size: 12 Mbp

FASTA size: 12 MB

(30 x 12,000,000) / 150 = 2,400,000 reads

Before compression: 2,400,000 x 100 = 240 MB

After compression: 240 / 4 = 60 MB

Drosophila

Genome size: 180 Mbp

FASTA size: 180 MB

(30 x 180,000,000) / 150 = 36,000,000 reads

Before compression: 36,000,000 x 100 = 3.6 GB

After compression: 3,600 / 4 = 900 MB

Human 

Genome size: 3,200 Mbp

FASTA size: 3.2 GB

(30 x 3,200,000,000) / 150 = 640,000,000 reads

Before compression: 640,000,000 x 100 = 64 GB

After compression: 64,000 / 4 = 16 GB





Penalty Warning

Do not commit the FASTA or FASTQ files to the repository!

The penalty may sound a bit too harsh, but I want to send a message and emphasize the importance of NOT COMMITTING large data files to the repository.

