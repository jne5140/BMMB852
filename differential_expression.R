# Set library path to custom location
.libPaths("~/Rlibs")

# Load necessary libraries
library(limma)
library(pheatmap)

# Check if the output folder exists, and create it if not
if (!dir.exists("output")) {
  dir.create("output")
}

# Read input data
simulated_count_matrix <- read.csv("simulated_count_matrix.csv", row.names = 1)
design_matrix <- read.csv("design.csv", row.names = 1)

# Create design matrix
design <- model.matrix(~ Condition, data = design_matrix)

# Linear modeling for differential expression
fit <- lmFit(simulated_count_matrix, design)
fit <- eBayes(fit)
results <- topTable(fit, adjust = "BH", number = nrow(simulated_count_matrix))

# Save results to a file
write.csv(results, "output/differential_expression_results.csv")

# Perform PCA
normalized_counts <- voom(simulated_count_matrix, design)
print(dim(normalized_counts$E))  # Debug: Check dimensions of normalized data
print(head(normalized_counts$E))  # Debug: Preview normalized data

# Filter out zero-variance rows
filtered_data <- normalized_counts$E[rowSums(normalized_counts$E) > 0, ]
print(dim(filtered_data))  # Debug: Check dimensions after filtering

# Perform PCA
pca <- prcomp(t(filtered_data))  # Perform PCA on filtered data

# Assign colors based on conditions
conditions <- as.factor(design_matrix$Condition)
print(levels(conditions))  # Debug: Check condition levels
colors <- rainbow(length(levels(conditions)))[conditions]

# Plot PCA
png("output/pca_plot.png")
plot(pca$x[, 1:2],
     col = colors,
     pch = 16,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of Differentially Expressed Genes")
legend("topright", legend = levels(conditions), col = rainbow(length(levels(conditions))), pch = 16)
dev.off()

# Generate Heatmap for top 50 differentially expressed genes
top_genes <- rownames(results)[1:50]  # Use row names if `results$ID` is empty
heatmap_data <- simulated_count_matrix[top_genes, , drop = FALSE]
if (nrow(heatmap_data) > 0) {
    png("output/heatmap.png")
    pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
    dev.off()
    message("Heatmap saved to output/heatmap.png")
} else {
    message("Heatmap data is empty, skipping heatmap generation.")
}


# Check if all genes exist in the matrix
if (all(top_genes %in% rownames(simulated_count_matrix))) {
  heatmap_data <- simulated_count_matrix[top_genes, ]
  print(heatmap_data)  # Debugging: Print heatmap data
} else {
  stop("Some top genes are not present in the count matrix.")
}

# Check if heatmap_data is not empty before plotting
if (nrow(heatmap_data) > 0) {
  png("output/heatmap.png")
  pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
  dev.off()
} else {
  message("Heatmap data is empty, skipping heatmap generation.")
}
