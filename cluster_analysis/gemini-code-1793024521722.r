library(Biostrings)
library(tidyverse)
library(umap)
library(cluster)
library(RColorBrewer)

# ==========================================
# 1. LOAD SEQUENCES
# ==========================================
# Replace 'centromere_monomers.fasta' with your actual file path
fasta_file <- "centromere_monomers.fasta"
seqs <- readDNAStringSet(fasta_file)

cat("Loaded", length(seqs), "sequences.\n")

# ==========================================
# 2. FEATURE EXTRACTION
# ==========================================
# Feature 1: Sequence Length
seq_lengths <- width(seqs)

# Feature 2: Dinucleotide Frequencies (normalized proportions)
# This captures local sequence variations and motifs characteristic of different groups
di_freq <- dinucleotideFrequency(seqs, as.prob = TRUE)

# Combine length and sequence features into a single data frame
feature_df <- data.frame(
  Length = seq_lengths,
  as.data.frame(di_freq)
)

# ==========================================
# 3. FEATURE SCALING & DIMENSIONALITY REDUCTION
# ==========================================
# Scale features so that length (in base pairs) doesn't overly dominate 
# the proportional dinucleotide frequencies
scaled_features <- scale(feature_df)

# Run UMAP to project high-dimensional features into 2D space
set.seed(42)
umap_res <- umap(scaled_features)

# Create a plotting data frame
plot_df <- data.frame(
  UMAP1 = umap_res$layout[, 1],
  UMAP2 = umap_res$layout[, 2],
  Length = seq_lengths
)

# ==========================================
# 4. CLUSTERING (K-means, k = 2)
# ==========================================
# Since you noted they form two groups, we specify centers = 2
set.seed(42)
kmeans_res <- kmeans(scaled_features, centers = 2, nstart = 25)
plot_df$Cluster <- as.factor(kmeans_res$cluster)

# Print cluster summary sizes
print(table(plot_df$Cluster))

# ==========================================
# 5. VISUALIZATION
# ==========================================

# Plot 1: UMAP Clustering Scatter Plot
p_umap <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Centromeric Monomer Clustering",
    subtitle = "Based on Length and Dinucleotide Composition",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Group"
  ) +
  scale_color_brewer(palette = "Set1")

print(p_umap)

# Plot 2: Length Distribution Comparison between the two groups
p_length <- ggplot(plot_df, aes(x = Cluster, y = Length, fill = Cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.3, color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Length Distribution per Group",
    subtitle = "Validating length separation between clusters",
    x = "Cluster Group",
    y = "Monomer Length (bp)",
    fill = "Group"
  ) +
  scale_fill_brewer(palette = "Set1")

print(p_length)