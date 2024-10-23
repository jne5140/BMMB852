# BMMB 852 HW#9

#### Jessica Eckard


This assignment requires writing a Makefile and a markdown report. 
Use the Makefile developed for your previous assignment, which contains targets for creating a BAM file.
Use the BAM file generated using reads from the SRA.
In your report, show the commands and the answers for each of these questions:

1.	How many reads did not align with the reference genome?

2.	How many primary, secondary, and supplementary alignments are in the BAM file?

3.	How many properly-paired alignments on the reverse strand are formed by reads contained in the first pair?

4.	Make a new BAM file that contains only the properly paired primary alignments with a mapping quality of over 10. 

5.	Compare the flagstats for your original and your filtered BAM file.

__________________________________
	micromamba run samtools flagstat sra_reads.bam > sra_reads.stats.txt

 ![image](https://github.com/user-attachments/assets/2df7632b-2a9e-438e-8a4c-19b4111a5621)


1.	Primary reads: 18,118

	Primary mapped reads (main alignments): 18,000

	18,118 – 18,000 = **118 reads**

	118 reads did not align with the reference genome

2.	**Primary reads: 18,118**
	
	**Secondary reads: 0**

	**Supplementary reads: 12**

3.	Number of properly paired read1 alignments on the reverse strand

		samtools view -f 0x2 -f 0x10 -f 0x40 sra_reads.bam | wc -l
	**4,443 reads**

4.	New BAM file

		samtools view -b -h -q 10 -f 0x2 -F 0x100 sra_reads.bam -o filtered_reads.bam

5.	Comparison

		micromamba run samtools flagstat filtered_reads.bam > filtered_reads.stats.txt
	
	**Original SRA reads stats:**

	![image](https://github.com/user-attachments/assets/d37cc4e9-6559-48fc-9241-39ae60867cbe)


	**Filtered SRA reads stats:**

	(Insert Image)

 	**Total Reads**

	18,130 original – 17,712 filtered = 418
	
	After filtering, there are 418 fewer reads in total. There were 418 reads with a mapping quality of 10 or lower.

	**Primary Alignments:**

	18,118 original – 17,706 filtered = 412

	There are 412 fewer primary reads after filtering.

	**Supplementary Alignments:**

	12 original – 6 filtered = 6

	There are 6 fewer supplementary alignments after filtering. 

	**Mapped Reads:**

	18,012 original (99.35%)
	
	17,712 filtered (100%)

	This indicates that unaligned or low-quality aligned reads have been removed.

	**Properly Paired Reads:**

	17,862 original (98.59%)
	
	17,706 filtered (100%)

	In the filtered set, all reads are properly paired (100%), while the original set had a slightly lower percentage (98.59%).

	**Singletons:**

	108 original (0.60%)
	
	0 filtered (0%)

	The filtering process removed all singletons, leaving only properly paired reads.

	**Reads Mapped to a Different Chromosome:**

	0 original
	
	0 filtered
	
	Both the original and filtered datasets have no reads mapped to a different chromosome.

	**Paired Reads:**

	Read1: 9,059 original – 8,855 filtered = 204

	Read2: 9,059 original – 8,851 filtered = 208

	There are slightly fewer read1 and read2 pairs in the filtered set, indicating that some low-quality pairs were removed.

	
	Overall, the filtering improved the dataset by removing low-quality, unpaired, and unmapped reads. 

