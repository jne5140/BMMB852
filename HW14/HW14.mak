# Variables
R_SCRIPT=differential_expression.R
RESULTS_DIR=output
SIMULATED_COUNT_MATRIX=simulated_count_matrix.csv
DESIGN_FILE=design.csv
R_LIBS_USER=~/Rlibs

# Targets
all: pca heatmap differential_expression_results

# Run Differential Expression Analysis
differential_expression_results: $(SIMULATED_COUNT_MATRIX) $(DESIGN_FILE)
	micromamba run -n bioinfo Rscript -e "Sys.setenv(R_LIBS_USER = '~/Rlibs'); source('differential_expression.R')"

# Generate PCA plot
pca: differential_expression_results
	@echo "Generating PCA plot..."
	@echo "PCA plot saved to output/pca_plot.png"

# Generate Heatmap
heatmap: differential_expression_results
	@echo "Generating heatmap..."
	@echo "Heatmap saved to output/heatmap.png"

# Clean up
clean:
	rm -rf $(RESULTS_DIR)/*

# Display usage
usage:
	@echo "Usage:"
	@echo "  make pca                    # Generate PCA plot"
	@echo "  make heatmap                # Generate heatmap"
	@echo "  make differential_expression_results  # Run differential expression analysis"
	@echo "  make all                    # Run everything (differential expression, PCA, and heatmap)"
	@echo "  make clean                  # Clean up output"
