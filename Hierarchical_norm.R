library(ggplot2)
library(preprocessCore)

# Prompt user for the gene data file path
gene_data_file <- readline(prompt = "Enter the path to the gene data file (e.g., '~/desktop/df_no0_final'): ")
gene_data <- read.delim(gene_data_file, header = TRUE, sep = "\t")
gene_column <- gene_data$GENE
gene_data <- gene_data[, -1]

# Standardize the data
standardized_data <- t(apply(gene_data, 1, scale))
gene_data <- cbind(GENE = gene_column, as.data.frame(standardized_data))
colnames(gene_data) <- c("GENE", "X0hr", "X1hr", "X3hr", "X24hr")

# Set row names
rownames_gene_data <- gene_data[,1]
gene_data <- gene_data[,-1]
rownames(gene_data) <- rownames_gene_data

# Compute WCSS for the elbow method
wcss <- vector()  # Initialize vector to store WCSS values

# Loop to calculate WCSS for each number of clusters
for (i in 1:30) {
  kmeans_model <- kmeans(gene_data, centers = i)
  wcss[i] <- kmeans_model$tot.withinss  # Store WCSS for each number of clusters
}

# Plot the elbow curve
plot(1:length(wcss), wcss, type = "b", pch = 19, xlab = "Number of Clusters", 
     ylab = "Within-Cluster Sum of Squares (WCSS)", 
     main = "Elbow Method for Determining Number of Clusters")

# Hierarchical clustering
dist_matrix <- dist(gene_data)
hc <- hclust(dist_matrix)

# Cut tree to get clusters
clusters <- cutree(hc, k = 10)

# Add clusters to the gene data frame
gene_data_df <- data.frame(Gene = rownames(gene_data), gene_data)
gene_data_df$Cluster <- as.factor(clusters)

# Reshape the data for plotting
library(tidyr)
library(dplyr)

gene_data_long <- gene_data_df %>%
  pivot_longer(cols = -c(Gene, Cluster), names_to = "Sample", values_to = "Expression")

# Define the order of time points
time_order <- c("X0hr", "X1hr", "X3hr", "X24hr")
gene_data_long$Sample <- factor(gene_data_long$Sample, levels = time_order)

# Create the plot
combined_plot <- gene_data_long %>%
  filter(Cluster %in% unique(clusters)) %>%
  ggplot(aes(x = Sample, y = Expression, group = Gene, color = Gene)) +  # Mapping Gene to color
  geom_line(show.legend = FALSE) +  # Suppress legend
  facet_wrap(~Cluster, scales = "free") +
  ggtitle("MOA Profiles for Each Cluster") +
  xlab("Samples") +
  ylab("Gene Expression") +
  theme_minimal()

# Show combined plot
print(combined_plot)

# Prompt user for the output file path to save the plot
output_file <- readline(prompt = "Enter the path to save the combined plot (e.g., '~/desktop/combined_plot.svg'): ")

# Save the plot
ggsave(output_file, plot = combined_plot, device = "svg")
