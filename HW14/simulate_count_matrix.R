#!/usr/bin/env Rscript

# Set parameters for simulation
set.seed(42)
n_genes <- 1000  # Number of genes
n_samples <- 6   # Number of samples
conditions <- rep(c("Control", "Treatment"), each = n_samples / 2)

# Simulate counts
simulate_counts <- function(base_mean, dispersion, n_samples) {
  matrix(rnbinom(n_samples * length(base_mean), mu = base_mean, size = 1/dispersion),
         nrow = length(base_mean))
}

base_mean <- runif(n_genes, 10, 1000)  # Base mean for counts
dispersion <- 0.1

control_counts <- simulate_counts(base_mean, dispersion, n_samples / 2)
treatment_counts <- simulate_counts(base_mean * runif(n_genes, 0.5, 1.5), dispersion, n_samples / 2)

# Combine counts and add labels
count_matrix <- cbind(control_counts, treatment_counts)
rownames(count_matrix) <- paste0("Gene", 1:n_genes)
colnames(count_matrix) <- paste0("Sample", 1:n_samples)

# Save to CSV
write.csv(count_matrix, "simulated_count_matrix.csv", row.names = TRUE)

# Create a design file
design <- data.frame(SampleName = colnames(count_matrix), Condition = conditions)
write.csv(design, "design.csv", row.names = FALSE)

cat("Simulated count matrix and design file generated successfully.\n")

