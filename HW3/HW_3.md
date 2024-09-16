# **BMMB 852 HW #3**

## Jessica Eckard

**1. Reformat your previous assignment. Check the discussion thread: List of student repositories and verify that your repository passes the requirements. Rewrite your previous homework as a Markdown file. Commit the file next to the current file.**

Done, file is located in [Github repository](https://github.com/jne5140/BMMB852)


**2. Visualize the GFF file of your choice.**

**Using a resource of your choice, download the genome and annotation files for an organism of your choice.
(We recommend a smaller genome to make things go faster and to look at a simpler GFF file)**


`datasets download genome accession GCF_000005845.2`




`datasets download genome accession GCF_000005845.2 --include gff3,cds,protein,rna,genome`

**Use IGV to visualize the annotations relative to the genome.**

See attached image in canvas

**Separate intervals of type "gene" into a different file. If you don't have genes, pick another feature.**



`cat ncbi_dataset/data/GCF_000005845.2/genomic.gff | head`

`cat ncbi_dataset/data/GCF_000005845.2/genomic.gff | awk ' $3=="gene" { print $0 }'`

`cat ncbi_dataset/data/GCF_000005845.2/genomic.gff | awk ' $3=="gene" {print $0 }' > ncbi_dataset/data/GCF_000005845.2/genomic.gff`

`cat ncbi_dataset/data/GCF_000005845.2/genomic.gff | awk ' $3=="CDS" {print $0 }' > ncbi_dataset/data/GCF_000005845.2/cds.gff`

 **Load up the simplified GFF file into IGV as a separate track.**

See attached image in canvas

**Zoom in to see the sequences, line up CDS (coding sequences) with the translation table in IGV**


See attached image in canvas

**Verify that coding sequences start with a start codon and end with a stop codon.**

See attached images on canvas


 **Manually generate a simple GFF file that can be visualized in your genome.**

`cat GCF_000005845.2_ASM584v2_genomic.fna | head -1>NC_000913.3`

`NC_000913.3	.	.	1000	2000	.	.	.	.`

**Load your manually generated GFF file as a separate track in IGV.**

See attached image

Report your findings in text, and provide a the code for the commands you ran in your markdown report.

Report your findings in text, and provide a few relevant screenshots.
You can attach the screenshots to CANVAS the assignment or even put them into the repository and link them to your markdown report.
Your submission needs to list the link to the markdown file in your repository so that reviews can access it.
