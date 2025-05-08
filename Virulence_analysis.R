# Virulence analysis

## Set the working space

```R
# Libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(pheatmap)
library(tidyr)

# Set directory and read files
rute <- ".../directory/"
files <- list.files(path = rute, pattern = "VFDB.csv", full.names = T)
data <- do.call(rbind, lapply(files, function(file) {
  sample_id <- sub("_VFDB\\.csv$", "", basename(file))
  df <- read.csv(file, header = F, sep = "\t")
  df$sample <- sample_id
  return(df)
}))

colnames(files) <- c("qseqid", "sseqid", "stitle", "pident", "evalue", 
                     "staxids", "sample")
data$sample <- str_remove(data$sample, "eICh24_") # remove part of the name
content <- str_extract_all(data$stitle, "\\[([^]]+)\\]")
data$title1 <- sapply(content, `[`, 1)
data$title2 <- sapply(content, `[`, 2)


# percentage of identity
ggplot(data, aes(x = sample, y = pident))+
  geom_boxplot(fill = "blue", color = "black")+
  labs(title = "Identity",
       x = "Sample",
       y = "Percentage of identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# e-values
data$log_evalue <- -log10(data$evalue)
ggplot(data, aes(x = sample, y = log_evalue)) +
  geom_boxplot(fill = "orange", color = "black") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red")+ # P-value: 0.00001
  labs(title = "e-values",
       x = "Sample",
       y = "-log10(e-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```
## Frequency for each virulence factor
```R
# Calculate total frequency for each virulence factor
top_fact <- data %>%
  filter(evalue <= 0.0000001) %>%
  count(title1, sort = TRUE) %>%  # Count total occurrences
  top_n(100, n)  # Select top 100 most frequent

# Count occurrences of each factor per sample
heatmap_data <- data %>%
  filter(evalue <= 0.001, title1 %in% top_fact$title1) %>%  # Filter
  group_by(title1, sample) %>%  # Group by virulence factor and sample
  summarise(frecuency = n(), .groups = "drop") %>%  # Count occurrences
  pivot_wider(names_from = sample, values_from = frecuency, values_fill = 0) # Convert to wide format

# Convert to matrix (excluding the first column, which is the factor name)
matrix_heatmap <- as.matrix(heatmap_data[,c(-1)])

# Assign row names
rownames(matrix_heatmap) <- heatmap_data$title1

# write the table in csv
write.csv(matrix_heatmap, "virulence_factors.csv")

# Plot heatmap
pheatmap(matrix_heatmap,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Colors from low to high
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize = 25,
         fontsize_col = 30,
         fontsize_row = 25)
grid.text("Z-scores of frequency", 
          x = 0.92, y = 0.95, gp = gpar(fontsize = 20, fontface = "bold"))

# Save at 40x30
```
