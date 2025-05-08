# Assessment of taxonomic diversity

## Set the working space
```R
# CRAN libraries

library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(tibble)
library(tidyr)
library(zoo)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(permute)
library(lattice)
library(vegan)
library(boot)

# Biocmanager libraries
library(rentrez)
library(taxize)
library(DESeq2)


# Working directory
setwd("your/working/directory/")

# Directory containing files
dir_path <- "Directory/with/kreports/"


# Set ENTREZ API Key globally
ent_ap_key = "your_key"
options(entrez_api_key = ent_ap_key)
getOption("entrez_api_key")

# directory containing files
dir_path <- ".../Centrifuge_kreport"

# list with file paths
file_list <- list.files(path = dir_path, full.names = TRUE, pattern = "kreport")

# function to load kreport files
read_kreport <- function(file_path){
  read.delim(file = file_path, 
             header = FALSE, sep = "\t", 
             stringsAsFactors = FALSE,
             col.names = c("percentage", 
                           "reads", 
                           "children_reads", 
                           "rank", 
                           "taxid", 
                           "name"))
}

# function to extract the filename without extension and add it as a new column
process_files <- function(file_processed){
  d <- read_kreport(file_processed)
  e <- basename(file_processed)
  e <- gsub("kreport_filtered_eICh24_|_bacteria","", e)
  d$sample <- tools::file_path_sans_ext(e)
  return(d)
}

# apply process_files function to all files in file_list
all_data <- lapply(file_list, process_files)

# combine all into a single table
combined_data <- do.call(rbind, all_data)

# Filter for orders only (rank == 'O') and keep those with reads > 0
filtered_data_1 <- combined_data %>%
  filter(rank == "O" & reads > 0 )

# trim and format the table
filtered_data_1$name <-  trimws(filtered_data_1$name)

```

## Convert names to taxIDs, and retrieve desired orders
```R

filtered_data_2 <- unique(filtered_data_1$name) # obtain unique values
filtered_data_2[102] <- "Bacteroidetes Incertae sedis" # the name give problems, so it has to be correctly formatted

# get IDs
taxids_ori <- get_uid(filtered_data_2)
sum(is.na(taxids_ori)) # if there are NAs, check filtered_data_2 and correct names

# Get classification for each taxID
classification_data <- classification(
  taxids_ori,
  db = "ncbi",
  apikey = ent_ap_key
)

sum(is.na(names(classification_data))) # if any NA, review from filtered_data_2

# Retrieve desired filter (e.g., viruses, bacteria, families)
long <- length(classification_data)

superkingdoms_taxid <- lapply(1:long, function(i) {
  # Check if element i is a data frame
  if (is.data.frame(classification_data[[i]])) {
    # Add a column with the key (name of the list) to the dataframe
    df <- classification_data[[i]]
    df$taxid <- names(classification_data)[i]  # Add identifier (taxid)
    return(df)  # Return dataframe with 'taxid' column
  } else {
    return(NULL)  # If not a data frame, return NULL
  }
}) %>%
  # Combine all tables into a single data frame
  bind_rows() %>%
  # Select columns 'name', 'id', 'rank' and 'taxid'
  select(name, id, rank, taxid)

# number of bacteria orders in the samples
sum(superkingdoms_taxid$name == "Bacteria")

# Filter by bacteria superkingdom the IDs
bacteria_id <- superkingdoms_taxid %>%
  filter(name == "Bacteria")

```

## Filter all the data by the ID of the domain previously obtained and add metadata for comparisions
```R
# Filter orders (genus, class, family, order) by domain (virus, bacteria)
filtered_data <- filtered_data_1 %>%
  filter(taxid %in% bacteria_id$taxid)

# Trimm sample names
filtered_data <- filtered_data %>%
  mutate(sample = str_remove(sample, "filtered_" ))
filtered_data$name <- trimws(filtered_data$name)

# Add dates
m <- paste0("eICh24_", sprintf("%02d", seq(1, 48)))
f <- c(
  "2023_2", "2023_2", "2023_2", "2023_2", "2023_3", "2023_3", "2023_2", "2023_1",
  "2023_1", "2023_1", "2023_10", "2023_2", "2023_10", rep("2024_4", 35))
dates <- data.frame(date = f, sample = m)
filtered_data <- merge(filtered_data, dates, by = "sample", all.x = TRUE)
Sys.setlocale("LC_TIME", "C")
filtered_data$date <- as.yearmon(filtered_data$date, format = "%Y_%m")

# Calculate absolute abundance matrix
abundance_matrix <- filtered_data %>%
  group_by(sample, name) %>%
  summarise(total_percentage = sum(reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = total_percentage)
abundance_matrix$name <- trimws(abundance_matrix$name)

# Get phenodata
feno <- filtered_data %>%
  select(sample, date) %>%
  distinct()
feno <- as.data.frame(feno) # convert to data.frame
row.names(feno) <- feno$sample # add sample name as row names
feno <-  feno %>%
  mutate(date = str_remove(date, " ")) # trimm spaces
feno$date <- as.factor(feno$date) # make sure the date column is a factor

# Add seasson as a factor
vector1 <-c(rep("Hibernation", 13), rep("Breeding season", 35)) 
feno$date2 <- as.factor(vector1)
feno <- feno[,-1]

```

## Differential abundance analysis

```R
# Change the matrix format to data.frame
abundance_matrix <- as.data.frame(abundance_matrix)
row.names(abundance_matrix) <- abundance_matrix$name 
abundance_matrix <- abundance_matrix[,-1]

# Fit DESeq model
dds <- DESeqDataSetFromMatrix(countData = abundance_matrix,
                              colData = feno,
                              design = ~ date2)  # ~1 means no factor model

# Normalize absolute abundance matrix using DESeq
dds <- DESeq(dds)

```
# Taxa that change among seasons
```R
resultsNames(dds)
res <- results(dds, name = "date2_Breeding_season_vs_Hibernation_season", alpha = 0.05)
dif_spec <- res[!is.na(res$pvalue) & res$pvalue < 0.05,]
dif_spec

# Change reference level to compare between months, e.g. 
dds$date <- relevel(dds$date, ref = "Oct2023")  # Change reference to "Oct2023"
dds <- DESeq(dds, fitType = "local")  # Refit the model

# Now compare "Apr2024" vs "Oct2023"
res <- results(dds, name = "date_Apr2024_vs_Oct2023")
dif_spec <- res[!is.na(res$padj) & res$padj < 0.05,]
write.csv(dif_spec, "res.csv")

```
## Principal Component Analysis

```R
# Variance Stabilizing Transformation (VST)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean")
rlog_counts <- assay(vsd)

# PCA
pca_data <- prcomp(t(rlog_counts))

# Create DataFrame with PCA coordinates and seasons
pca_df <- data.frame(pca_data$x)
pca_df$season <- feno$date2

# Rename rows
row.names(pca_data$x) <- c(sprintf("%02d", seq(1, 48)))

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = season)) +
  geom_point(size = 4) +
  labs(title = "", x = "Principal component 1", 
       y = "Principal component 2", ) +
  geom_label(aes(label = row.names(pca_data$x)), vjust = -1, hjust = 0.5, size = 4) + 
  theme_minimal()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         plot.background = element_blank(),
         axis.line = element_line(color = "black"),
         legend.position = "none",
         axis.title.x = element_text(face = "bold", size = 20),
         axis.title.y = element_text(face = "bold", size = 20))

# Save plot at 16x10
```

## Diversity Analyses

```R

# Change the format of abundance matrix
abundance_matrix_2 <- filtered_data %>%
  group_by(sample, name) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0) %>%
  column_to_rownames("name")
dat <- t(abundance_matrix_2)
```

### Rarefaction

``` R
# Function that generates rarefaction points for each sample
rarefaction_df <- function(mat, step = 100, max_depth = NULL) {
  results <- list()
  for (i in 1:nrow(mat)) {
    sample_name <- rownames(mat)[i]
    counts <- mat[i, ]
    counts <- counts[which(counts > 0)]
    total_reads <- sum(counts)
    max_reads <- if (is.null(max_depth)) total_reads else min(max_depth, total_reads)
    depths <- seq(1, max_reads, by = step)
    richness <- sapply(depths, function(d) {
      rarefy(counts, sample = d)
    })
    df <- data.frame(
      Sample = sample_name,
      Depth = depths,
      Richness = richness
    )
    results[[i]] <- df
  }
  do.call(rbind, results)
}

# Minimum depth
min_depth <- 6000

# Generate ggplot dataframe
rarefaction_data <- rarefaction_df(dat, step = 100, max_depth = min_depth)

ggplot(rarefaction_data, aes(x = Depth, y = Richness, color = Sample)) +
  geom_line(size = 1) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(color = "black"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title = element_text(size = 20, face = "bold")) +
  labs(title = "",
       x = "Number of reads",
       y = "Estimated richness (number of orders)") 

# Save plot at 16x10
```

### Alpha diversity

```R
# Apply Total Sum Scaling (TSS) - Convert to relative abundance
read_tss <- sweep(dat, 1, rowSums(dat), FUN = "/")

# Verify row sums are equal to 1
rowSums(read_tss)  

# Shannon index
shannon_index <- diversity(t(read_tss), index = "shannon")

# Pielou's evenness index (J)
richness <- specnumber(t(read_tss))
pielou_evenness <- shannon_index / log(richness)

# Simpson index
simpson_index <- 1-(diversity(t(read_tss), index = "simpson"))

# Inverse Simpson index
inv_simpson_index <- diversity(t(read_tss), index = "invsimpson")

# Combine results into a dataframe
diversity_indices <- data.frame(Shannon = shannon_index,
                                Simpson = simpson_index,
                                InvSimpson = inv_simpson_index,
                                Pielou = pielou_evenness)

# Show results
print(diversity_indices)

# Change row names
row.names(diversity_indices) <- c(sprintf("%02d", seq(1,48)))
write.csv(diversity_indices,"diversity_indices.csv")

# Plotting
shanon_g <- ggplot(diversity_indices, aes(x = feno$date2, y = shannon_index)) +
  geom_violin(trim = TRUE, fill = "#31688EFF", linewidth = 0.1) + 
  geom_boxplot(width = 0.1, color = "grey", fill = NA) +  
  theme_minimal() + 
  labs(x = "", y = "Shanon index") + 
  theme(axis.text.x = element_text(angle = 0, face = "bold",
                                   hjust = 0.5, size = 12),  
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank()) +
  stat_compare_means(method = "wilcox.test", label.x = 1.3)

# Estimate differences in abundances
wilcox.test(diversity_indices$Shannon[1:13], diversity_indices$Shannon[14:48],
            exact = FALSE, conf.int = TRUE)

simpson_g <- ggplot(diversity_indices, aes(x = feno$date2, y = simpson_index)) +
  geom_violin(trim = TRUE, fill = "#35B779FF", linewidth = 0.1) + 
  geom_boxplot(width = 0.1, color = "grey", fill = NA) +  
  theme_minimal() + 
  labs(x = "", y = "Simpson index") + 
  theme(axis.text.x = element_text(angle = 0, face = "bold",
                                   hjust = 0.5, size = 12),  
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank()) +
  stat_compare_means(method = "wilcox.test", label.x = 1.3)

wilcox.test(diversity_indices$Simpson[1:13], diversity_indices$Simpson[14:48],
            exact = FALSE, conf.int = TRUE)

InvSimpson_g <- ggplot(diversity_indices, aes(x = feno$date2, 
                                              y = inv_simpson_index)) +
  geom_violin(trim = TRUE, fill = "#FDE725FF", linewidth = 0.1) + 
  geom_boxplot(width = 0.1, color = "grey", fill = NA) +  
  theme_minimal() + 
  labs(x = "", y = "Simpson inverse index") + 
  theme(axis.text.x = element_text(angle = 0, face = "bold",
                                   hjust = 0.5, size = 12),  
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank()) +
  stat_compare_means(method = "wilcox.test", label.x = 1.3)

wilcox.test(diversity_indices$InvSimpson[1:13], diversity_indices$InvSimpson[14:48],
            exact = FALSE, conf.int = TRUE)

pielou_evenness_g <- ggplot(diversity_indices, 
                            aes(x = feno$date2, y = pielou_evenness)) +
  geom_violin(trim = TRUE, fill = "#440154ff", linewidth = 0.1) + 
  geom_boxplot(width = 0.1, color = "grey", fill = NA) +
  theme_minimal() + 
  labs(x = "", y = "Pielou evenness (J)") + 
  theme(axis.text.x = element_text(angle = 0, face = "bold",
                                   hjust = 0.5, size = 12),  
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank()) +
  stat_compare_means(method = "wilcox.test", label.x = 1.3)

wilcox.test(diversity_indices$Pielou[1:13], diversity_indices$Pielou[14:48],
            exact = FALSE, conf.int = TRUE)

## Layout: combine plots
layout <- (shanon_g / simpson_g) | (InvSimpson_g / pielou_evenness_g)
final_plot <- layout + plot_layout(widths = c(1,1), heights = c(1,1)) +
  plot_annotation(tag_levels = "A")
print(final_plot)

# Save at 16 X 10
```


### Beta diversity

```R
# TSS normalization - Total Sum Scaling - and hellinger transformation
tss_matrix <- sweep(abundance_matrix, 1, rowSums(abundance_matrix), "/")
hellinger_data <- decostand(tss_matrix, method = "hellinger")

# Beta diversity analysis
bray_dist <- vegdist(t(hellinger_data), method = "bray") # abundances
jaccard_dist <- vegdist(t(hellinger_data), method = "jaccard") # species presence
euclidean_dist <- dist(t(hellinger_data), method = "euclidean") # normalized abundances, very sensitive to extreme differences

bray_matrix <- as.matrix(bray_dist)
jaccard_matrix <- as.matrix(jaccard_dist)
euclidean_matrix <- as.matrix(euclidean_dist)

write.csv(bray_matrix, "bray_matrix.csv", row.names = T)
write.csv(jaccard_matrix, "jaccard_matrix.csv", row.names = T)
write.csv(euclidean_matrix, "euclidean_matrix.csv", row.names = T)

# Visualization PCoA (Principal Coordinates Analysis)
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)

# Calculate explained variance
var_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)

# Change row names
row.names(pcoa_res$points) <- c(sprintf("%02d", seq(1,48)))

# Plot
plot(pcoa_res$points, col = feno$date2, pch = 19,
     xlab = paste0("PCoA 1 (", var_explained[1], "%)"),
     ylab = paste0("PCoA 2 (", var_explained[2], "%)"),
     main = "")
text(pcoa_res$points[,1], pcoa_res$points[,2], 
     labels = rownames(pcoa_res$points), pos = 2, cex = 0.8)

# Ellipses by date group
ordiellipse(pcoa_res$points, feno$date2, kind = "sd", draw = "polygon",
            col = c("grey","red"), alpha = 50, label = FALSE,
            border = c("grey","red"))

# Plot at 16 X 10

# Statistical analysis: PERMANOVA (Permutational Multivariate Analysis of Variance)
adonis2(bray_dist ~ date2, data = feno, permutations = 999)
```
