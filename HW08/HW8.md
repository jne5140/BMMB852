# BMMB 852 HW8
### Jessica Eckard

This assignment requires writing a Makefile and a markdown report.

Use the Makefile developed for your previous assignment that contains rules for simulating reads and obtaining reads from SRA.

Make a new version of the Makefile that includes rules for aligning reads to a reference genome.

Add a new target called index and align to the Makefile. The index target should create an index for the reference genome. The align target should align the reads to the reference genome.

Visualize the resulting BAM files for your simulated reads and for the reads downloaded from SRA.

Generate alignment statistics for the reads from both sources, simulated and SRA.
Briefly describe the differences between the two datasets.

___________________________________________________________________________________________________________________________________

Done, file is located in Github repository.

![image](https://github.com/user-attachments/assets/0c18fa9f-6439-4941-a586-74507d82309e)

### Simulated Reads:
![image](https://github.com/user-attachments/assets/bb30d0c5-d5b9-41ec-adea-fe6c1db2defe)


### SRA Reads:
![image](https://github.com/user-attachments/assets/23ea81ba-100d-47cc-b389-d893f330d062)

![image](https://github.com/user-attachments/assets/023d5dec-2e2d-4e2c-b2ad-e47d62369d44)


The simulated reads have a slightly lower mapping rate (97.70%) compared to the SRA reads (99.35%), which indicates the SRA data better aligns to the reference genome. The SRA dataset has fewer singletons (0.60%) than the simulated reads (2.29%), suggesting it captures real sequencing variations better. Properly paired reads are also higher in the SRA dataset (98.59%) versus the simulated reads (95.20%). Overall, while the simulated reads are more of a controlled representation, the SRA reads are more refletive of actual biological data.
