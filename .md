# **BMMB 852 HW #3**

## Jessica Eckard

**1. Reformat your previous assignment. Check the discussion thread: List of student repositories and verify that your repository passes the requirements. Rewrite your previous homework as a Markdown file. Commit the file next to the current file.**

Done, file is located in [Github repository](https://github.com/jne5140/BMMB852)


**2. Visualize the GFF file of your choice.**

**Using a resource of your choice, download the genome and annotation files for an organism of your choice.
(We recommend a smaller genome to make things go faster and to look at a simpler GFF file)**


`datasets download genome accession GCF_000005845.2
New version of client (16.28.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    1.38MB valid zip structure -- files not checked
Validating package [================================================] 100% 5/5
(bioinfo)`




`datasets download genome accession GCF_000005845.2 --include gff3,cds,protein,rna,genome
New version of client (16.28.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    4.2MB valid zip structure -- files not checked
Validating package [================================================] 100% 8/8
(bioinfo)`

**Use IGV to visualize the annotations relative to the genome.**

See attached image in canvas

**Separate intervals of type "gene" into a different file. If you don't have genes, pick another feature.**



`cat ncbi_dataset/data/GCF_000005845.2/genomic.gff | head
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
#!genome-build ASM584v2
#!genome-build-accession NCBI_Assembly:GCF_000005845.2
##sequence-region NC_000913.3 1 4641652
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=511145
NC_000913.3     RefSeq  region  1       4641652 .       +       .ID=NC_000913.3:1..4641652;Dbxref=taxon:511145;Is_circular=true;Name=ANONYMOUS;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=K-12;substrain=MG1655
NC_000913.3     RefSeq  gene    190     255     .       +       .ID=gene-b0001;Dbxref=ASAP:ABE-0000006,ECOCYC:EG11277,GeneID:944742;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;gene_synonym=ECK0001;locus_tag=b0001
NC_000913.3     RefSeq  CDS     190     255     .       +       0ID=cds-NP_414542.1;Parent=gene-b0001;Dbxref=UniProtKB/Swiss-Prot:P0AD86,GenBank:NP_414542.1,ASAP:ABE-0000006,ECOCYC:EG11277,GeneID:944742;Name=NP_414542.1;gbkey=CDS;gene=thrL;locus_tag=b0001;orig_transcript_id=gnl|b0001|mrna.NP_414542;product=thr operon leader peptide;protein_id=NP_414542.1;transl_table=11
(bioinfo)`

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

`code demo.gff
(bioinfo)
cat GCF_000005845.2_ASM584v2_genomic.fna | head -1>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
(bioinfo)`

`NC_000913.3	.	.	1000	2000	.	.	.	.`

**Load your manually generated GFF file as a separate track in IGV.**

See attached image

Report your findings in text, and provide a the code for the commands you ran in your markdown report.

Report your findings in text, and provide a few relevant screenshots.
You can attach the screenshots to CANVAS the assignment or even put them into the repository and link them to your markdown report.
Your submission needs to list the link to the markdown file in your repository so that reviews can access it.
