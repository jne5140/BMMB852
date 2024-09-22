# BMMB 852 HW #4
## Jessica Eckard


Your submission should be a link to a folder that contains a script and a report written in Markdown, all committed to your repository.

**Part 1: Write a script**

In the previous lecture, you wrote a Markdown report. Rewrite all the code from that report as a Bash script. Your script should be:
1.	Reusable: separate variable definitions from the code.

2.	Documented: include comments that explain what each part of the script does.

Run your script on your original data and verify that it works.

	My script works on my data.
![image](https://github.com/user-attachments/assets/3d337c60-4ea3-494d-817c-9be07cd23a76)




 
You were also assigned to review someone else's report. Now, run your script on their data. If the script is reusable, you can replace your variables with theirs and run the script with minimal edits.

	I ran my script on Tram Ha’s data on S. cerevisiae. 

	Species: Saccharomyces cerevisiae

	RefSeq: GCF_000146045.2

	Resouce: NCBI

Add more functions to the script that also print some of their results. Were you able to reproduce their results? Make a note in the report.
 
	These are my results from Tram Ha’s data visualized in IGV. It should match their data, however I can’t tell as their images/results are on CANVAS and not on their GitHub repository. 

![image](https://github.com/user-attachments/assets/2169af9b-206a-43ca-9fa2-014937a71d02)





Commit the script to your repository.
	
	Done!
Part 2: Make use of ontologies

	Species: Escherichia coli

	RefSeq: GCF_000005845.2

	Resouce: NCBI

[	Link to dataset used](https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000005845.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED)



Sequence Ontology
1.	Choose a feature type from the GFF file and look up its definition in the sequence ontology.

`./ncbi_genome_download.sh
cat genomic.gff | grep -v '#' > genomic.gff3`

`cat genomic.gff3 | cut -f 3 | sort | uniq -c | sort -rn`

	I chose the ‘exon’ feature type from the GFF file. 
`URL=https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so-edit.obo`

`wget $URL`

`cat so-edit.obo | grep 'name: exon$' -B 1 -A 6`

	Definition of exon: "A region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing." [SO:ke]
2.	Find both the parent terms and children nodes of the term.

...

	Parent Terms:

`grep -A 10 "id: SO:0000147" sequence_ontology.obo | grep "is_a"`

	Transcript region is a parent term to exon.
	
...

	Children Nodes:
`grep -B 2 'is_a: SO:0000147' sequence_ontology.obo`

	Children nodes of exon are coding exon, noncoding exon, single exon, and interior exon. 

3.	Provide a short discussion of what you found.

...

	According to the Sequence Ontology (SO) for exon (SO:0000147), its definition is "A region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing." This shows that exons play a key role in translation, protein synthesis, and gene expression. The parent term for exon is transcript region (SO:0000833) which indicates that exons make up the transcript region in RNA transcription. The children nodes for exon are coding exon (SO:0000196), noncoding exon (SO:0001792), single exon (SO:0000842), and interior exon (SO:0000198). These terms show that there are different types of exons depending on their functions and locations within a transcript. 

Gene Ontology
1.	Identify a CC, MF, and BP term in the gene ontology relevant to your organism with a method of your choice.

[Gene Ontology Resource](https://www.geneontology.org/)

	CC: Cytoplasmic Region GO:0099568

	MF: ATP Binding GO:0005524
	
	BP: Glycolytic Process GO:0006096
2.	Explain the term and show its parent terms and children nodes.

...

	(Parent terms are listed above, children nodes are listed below)

...

	Cytoplasmic Region: “Any (proper) part of the cytoplasm of a single cell of sufficient size to still be considered cytoplasm.” Source: GOC:dos
![image](https://github.com/user-attachments/assets/d12defc4-4a7c-49be-bf79-fcf076ed057b)



	ATP Binding: “Binding to ATP, adenosine 5'-triphosphate, a universally important coenzyme and enzyme regulator.” Source: ISBN:0198506732
	 
![image](https://github.com/user-attachments/assets/effd320f-f294-4d64-b94b-64d622adb307)

 
	Glycolytic Process: “The chemical reactions and pathways resulting in the breakdown of a carbohydrate into pyruvate, with the concomitant production of a small amount of ATP and the reduction of NAD(P) to NAD(P)H. Glycolysis begins with the metabolism of a carbohydrate to generate products that can enter the pathway and ends with the production of pyruvate. Pyruvate may be converted to acetyl-coenzyme A, ethanol, lactate, or other small molecules.” Source: GOC:dph, GOC:bf, Wikipedia:Glycolysis, ISBN:0879010479, ISBN:0716720094, ISBN:0201090910
	 
![image](https://github.com/user-attachments/assets/e2177e65-38c8-4059-880c-7b276cf39b2c)

 
3.	Find genes that are annotated with the term. List the genes.

...

	Cytoplasmic Region (main genes): rpoA (RNA polymerase subunit alpha), gapA (Glyceraldehyde 3-phosphate dehydrogenase), groEL (Chaperonin GroEL)

...

	ATP Binding (main genes): dnaK (Chaperone DnaK), gyrA (DNA gyrase subunit A), recA (Recombination protein RecA)

...

	Glycolytic Process (main genes): gapA (Glyceraldehyde 3-phosphate dehydrogenase), pfkA (6-phosphofructokinase), fbaA (Fructose-bisphosphate aldolase)

4.	Discuss how well the genome seems to be annotated and whether the terms you found are broad or narrowly specific.

...

	The genome of Escherichia coli appears to be well annotated since essential genes listed above (dnaK, etc.) are known and well-defined. The terms I chose provide a good overview of cellular components/processes in E. coli. 

	CC: Cytoplasmic Region (GO:0099568)
	-	This is a broad term, as it includes numerous cytoplasmic components. This can be related to many proteins.
	
	MF: ATP Binding (GO:0005524)
	-	This is a relatively broad term but it focuses on a specific function, ATP Binding. This process is important to numerous cellular functions.

	BP: Glycolytic Process (GO:0006096)
	-	This is a more specific term, but is still a common process in E. coli. 
Submit the link to the folder containing your script and the report.

	Done!

