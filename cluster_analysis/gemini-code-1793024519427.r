if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
install.packages(c("tidyverse", "umap", "cluster", "RColorBrewer"))