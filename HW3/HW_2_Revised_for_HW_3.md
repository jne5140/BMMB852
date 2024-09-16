# BMMB 852 Applied Bioinformatics HW #2

**Jessica Eckard**

Instructions:
The Ensemble Public FTP site provides biological information for many species in various formats. For this assignment, we will investigate files in the GFF3 format found at:
http://ftp.ensembl.org/pub/current_gff3/

GFF3 is a tab-delimited file format described in (or use ChatGPT for a short summary):
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

Select an organism, download its corresponding GFF file then investigate this file with Unix command line tools.




```python
cd bmmb852l3
```





```python
wget https://ftp.ensembl.org/pub/current_gff3/salmo_trutta/Salmo_trutta.fSalTru1.1.112.gff3.gz
```
Find answers to the following questions:

**1.	Tell us a bit about the organism.**

Salmo trutta is the scientific name for brown trout. It is a species of salmonid ray-finned fish that is endemic to Europe and North Africa and has become invasive. 

**2.	How many features does the file contain?**


```python
cat Salmo_trutta.fSalTru1.1.112.gff3 | grep -v '#' | wc -l
3537901
```
3,537,901 features

```python
cat Salmo_trutta.fSalTru1.1.112.gff3 | grep -v '#' > trutta.gff3

cat trutta.gff3 | cut -f 3 | sort | uniq -c | sort -rn
1574331 exon
1521266 CDS
 116381 mRNA
 100180 biological_region
  97513 five_prime_UTR
  71832 three_prime_UTR
  43935 gene
   4441 ncRNA_gene
   2437 lnc_RNA
   1441 region
    854 snRNA
    757 rRNA
    702 snoRNA
    581 pseudogenic_transcript
    581 pseudogene
    400 miRNA
    132 V_gene_segment
     49 transcript
     39 D_gene_segment
     22 tRNA
     20 scRNA
      5 J_gene_segment
2	Y_RNA
```

**3.	How many sequence regions (chromosomes) does the file contain?**

```python
cat Salmo_trutta.fSalTru1.1.112.gff3 | head
##gff-version 3
##sequence-region   1 1 81542925
##sequence-region   10 1 46595448
##sequence-region   11 1 22956941
##sequence-region   12 1 97529106
##sequence-region   13 1 91488822
##sequence-region   14 1 86253525
##sequence-region   15 1 66900148
##sequence-region   16 1 61349433
##sequence-region   17 1 59764774

cat Salmo_trutta.fSalTru1.1.112.gff3 | grep '##sequence-region' | wc -l
1441
```

**4.	How many genes are listed for this organism?**

```python
cat trutta.gff3 | cut -f 3 | sort | uniq -c | sort -rn
1574331 exon
1521266 CDS
 116381 mRNA
 100180 biological_region
  97513 five_prime_UTR
  71832 three_prime_UTR
  43935 gene
   4441 ncRNA_gene
   2437 lnc_RNA
   1441 region
    854 snRNA
    757 rRNA
    702 snoRNA
    581 pseudogenic_transcript
    581 pseudogene
    400 miRNA
    132 V_gene_segment
     49 transcript
     39 D_gene_segment
     22 tRNA
     20 scRNA
      5 J_gene_segment
      2 Y_RNA
```
43,935 genes listed

**5.	What are the top-ten most annotated feature types (column 3) across the genome?**

```python
cat trutta.gff3 |cut -f 3 | sort-uniq-count-rank | head
1574331 exon
1521266 CDS
116381  mRNA
100180  biological_region
97513   five_prime_UTR
71832   three_prime_UTR
43935   gene
4441    ncRNA_gene
2437    lnc_RNA
1441    region
```

**6. Having analyzed this GFF file, does it seem like a complete and well-annotated organism?**

This organism seems complete and well-annotated, judging by the number of sequencing regions and genes listed. Of course, this determination is relative, and more sequencing could always be done to have an even more complete genome. 

**7. Share any other insights you might note.**

It’s interesting that there were only two small noncoding RNAs ‘Y_RNA’ found so far. I would think that that number should be higher. 

**8. Create a text file that shows both how you downloaded the data and how you generated each of your results.**

**9. Commit the file to your GitHub repository that you created for this course.**

[(link to repository)](https://github.com/jne5140/BMMB852)

Note that future assignments may ask someone else to repeat your findings. Build your report with repeatability/reproducibility in mind.




