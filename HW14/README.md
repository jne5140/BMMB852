# README for HW14 Makefile: Differential Expression Analysis

#### Jessica Eckard


### Files and Directories
- **`differential_expression.R`**: The R script that performs differential expression analysis, PCA plotting, and heatmap generation.

- **`simulated_count_matrix.csv`**: Input file containing the count matrix.

- **`design.csv`**: Input file describing the experimental design.

- **`output/`**: Directory where results (e.g., plots, differential expression results) are saved.

### Prerequisites

1. **R**: Ensure R is installed and configured to use the custom library path: `~/Rlibs`.
3. **Input Files**: Place `simulated_count_matrix.csv` and `design.csv` in the working directory if they are not already there.

## How to Use the Makefile

### Targets

#### 1. **`all`**
Runs the entire workflow, including differential expression analysis, PCA plot generation, and heatmap creation.
```bash
make all
```

#### 2. **`differential_expression_results`**
Runs the differential expression analysis using the input files.
```bash
make differential_expression_results
```


#### 3. **`pca`**
Generates a PCA plot based on the normalized count data.
```bash
make pca
```


#### 4. **`heatmap`**
Generates a heatmap of the top differentially expressed genes.
```bash
make heatmap
```


#### 5. **`clean`**
Removes all files from the `output/` directory.
```bash
make clean
```

#### 6. **`usage`**
Displays usage information for the Makefile.
```bash
make usage
```

### Makefile Variables
- **`R_SCRIPT`**: Specifies the R script used for the analysis (`differential_expression.R`).

- **`RESULTS_DIR`**: Directory where the output files are saved (`output`).

- **`SIMULATED_COUNT_MATRIX`**: Input count matrix file (`simulated_count_matrix.csv`).

- **`DESIGN_FILE`**: Input experimental design file (`design.csv`).

- **`R_LIBS_USER`**: Specifies the path to the R library directory (`~/Rlibs`).


### Notes
- Ensure all input files are correctly formatted and located in the working directory.
- If any errors occur, check the R script and input data for issues.
- The PCA plot and heatmap rely on successful completion of the differential expression analysis.

-------

Your submission should include a readme, a makefile, and a design file.

You may re-use code and genomes from your previous submission.

Perform a differential expression analysis of a count matrix.

Take a count matrix, this count matrix may be one generated in a previous assignment or one that is simulated.

Use method to identify genes/transcripts that show differential expression.

Draw a PCA plot and a heatmap for these genes.

Discuss the results. How many genes have you found. What kind of expression levels can you observe. How reliable does your data seem?

------

I simulated a count matrix using the `simulate_count_matrix.R` script found in this folder. I ran

      Rscript simulate_count_matrix.R

on my terminal, which gave me the `simulated_count_matrix.csv` and `design.csv` files in my working directory.

I used the limma (Linear Models for Microarray Data) package for identifying genes/transcripts that show differential expression.

(PCA plot image here)

(heatmap image here)

This command in R shows that none of the simulated gene values in the matrix were statistically significant (differentially expressed), as the output was 0.

      results <- read.csv("output/differential_expression_results.csv", row.names = 1)
      significant_genes <- subset(results, adj.P.Val < 0.05)
      cat("Number of significant genes:", nrow(significant_genes), "\n")

I then decided to examine the most significantly changed genes:

      top_genes <- head(results[order(results$adj.P.Val), ], 10)
      print(top_genes)
      write.csv(top_genes, "output/top_genes_summary.csv")

which gave me this output:

(insert image here)

There were a total of 1000 simulated genes analyzed. The expression levels show a range of values, with some genes showing substantial changes in expression, such as Gene793 (upregulated with a logFC of 892) and Gene238 (downregulated with a logFC of -497). Most genes have very low p-values, suggesting that the differential expression results are statistically significant, although a few genes, like Gene352 and Gene422, have adjusted p-values slightly above the typical significance threshold of 0.05 (0.106 and 0.107, respectively). This indicates that the results for these genes may be less reliable due to the large amount f genes analyzed. Overall, the data seems reliable-ish, but further validation (such as biological replication) is recommended to confirm the findings.

The PCA plot conditions appear to be properly separated, further indicating that the data is valid.

The heatmap also reflects that the data is reliable as there is clear clustering. 


**Note:** Thank you Dr. Albert for being such a wonderful professor! Although I struggled a bit with this class as I have no coding background whatsoever and thought this would be a fun elective, I learned so much and can honestly say this class is one of my favorites and will be one of the ones I'll remember when looking back on my college years. Thank you so much and have a wonderful break!
