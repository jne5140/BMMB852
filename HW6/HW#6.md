# BMMB 852 HW #6

#### Jessica Eckard

This assignment requires writing a script and a markdown report.

Write a script that downloads data from the SRA or ENA database and performs a quality control analysis on the data.

When searching for data use the genome you selected for the read simulation assignment as the target organism.

Identify a bad sequencing dataset. You may need to evaluate multiple SRR numbers to find one with poor quality.

**Data downloaded:**

**Database:** SRA

**Organism:** *S. aureus* (same as the read simulation assignment)

**Accession #:** SRR2584863

### 1.	Write a script to download data from the SRA database.

Done, script (HW#6_genome_analysis.sh) is located in the HW#6 folder in the GitHub repository. 

(Run HW#6_trim_data.sh right after, then use the finalized data to answer the questions below)

Files used for analysis:

![image](https://github.com/user-attachments/assets/d4d47426-e86a-478c-a7e8-d44b4dc74ff4)



### 2.	Evaluate the quality of the downloaded data.
 
![image](https://github.com/user-attachments/assets/c96758b8-ecd2-4edd-a2ed-fe2692dc6994)


Overall, per base sequence quality, along with per tile sequence quality and per base sequence content all are of poor quality. For the purpose of this assignment, I will focus on improving the per base sequence quality, which is shown below:

![image](https://github.com/user-attachments/assets/a9c04a65-7235-44ea-97b1-7055df6657cc)

 
### 3.	Improve the quality of the reads in the dataset.
Done, script (HW#6_trim_data.sh) is located in the HW#6 folder in the GitHub repository.

### 4.	Evaluate the quality again and document the improvements.
 
![image](https://github.com/user-attachments/assets/7337e6c5-9b48-4aea-94d0-eb1f7d53cbbe)

![image](https://github.com/user-attachments/assets/7e3f83cf-53b9-4151-abe3-6b32b161fcdf)



The quality after trimming has improved as there are no errors/fails, just warnings. As you can see, the trimming cut all of the worst-quality base sequences in red. 

--------------

Your script should include all necessary code. The markdown report should explain the data and publication it corresponds to and present the results of your analysis.

Commit the two files to your repository and submit the link to the folder containing the files.

Tip: Older datasets are often of lower quality the sequencing instruments were less reliable in the past.
