# Differential Abundance of Taxa, Metabolic Pathways, and Virulence Factors Across Two Sample Groups Identified by Multivariate Clustering Methods"
# Including PCA-GMM of functional abundances, Bray-Curtis PoCA, and PCA of VST-Transformed Abundance Data 

```R
# CRAN libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(tibble)
library(tidyr)
library(zoo)
library(ggpubr)
library(ggExtra)
library(compositions)
library(edgeR)

# Biocmanager
library(vsn)
library(Maaslin2)
```

### Start from the results of the Bacteria_taxa.R analysis.

```R

### TAXA ANALYSIS ###

# Filter gender by domain (bacteria)
filtered_data <- filtered_data_1 %>%
  filter(taxid %in% bacteria_id$taxid)

# Modify Data Format
filtered_data <- filtered_data %>%
  mutate(sample = str_remove(sample, "filtered_" ))
filtered_data$name <- trimws(filtered_data$name)

## Define the samples of group G1
G1_samples <- c(12, 24, 26, 27, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 42, 44, 45, 46, 48)
all_samples <- 1:48
group_vector <- ifelse(all_samples %in% G1_samples, "G1", "G2")
group_vector
m <- paste0("", sprintf("%02d", seq(1, 48)))
grp <- data.frame(sample = m, group = group_vector)
filtered_data <- merge(filtered_data, grp, by = "sample", all.x = TRUE)

# Create the absolute abundance matrix
abundance_matrix <- filtered_data %>%
  group_by(sample, name) %>%
  summarise(total_percentage = sum(reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = total_percentage)
abundance_matrix$name <- trimws(abundance_matrix$name)

# Get phenodata
feno <- filtered_data %>%
  select(sample, group) %>%
  distinct()
feno <- as.data.frame(feno)
row.names(feno) <- feno$sample
feno$group <- as.factor(feno$group)
abundance_matrix <- as.data.frame(abundance_matrix)
row.names(abundance_matrix) <- abundance_matrix$name 
abundance_matrix <- abundance_matrix[,-1]

######### Differential abundance analysis ##########
# Normalization by TSS
tss_normalized_data <- sweep(abundance_matrix, 2, 
                             colSums(abundance_matrix), FUN = "/")
colSums(tss_normalized_data)

# Transformation  by log1p
joined_data_3 <- log1p(tss_normalized_data)


# Select samples, define directory, run Maaslin2, etc.
# Total number of samples (columns) in the dataset
n_muestras <- ncol(joined_data_3)
results_list <- vector("list", n_muestras)

set.seed(123)

for(i in 1:n_muestras){
  # Define jackknife indices: exclude the i-th sample
  jackknife_indices <- setdiff(1:n_muestras, i)
  
  # Subset of data without the i-th sample
  subset_data <- joined_data_3[, jackknife_indices]
  
  # Adjust metadata to include only the samples used
  subset_metadata <- feno[rownames(feno) %in% colnames(subset_data), ]

  # Define an output directory for each iteration
  out_dir <- paste0("maslin/maaslin2_jackknife_", i)
  
  tryCatch({
    # Run Maaslin2 with the subset
    fit <- Maaslin2(
      input_data = subset_data,
      input_metadata = subset_metadata,
      output = out_dir,
      fixed_effects = c("group"),
      normalization = "NONE",
      transform = "NONE",
      cores = 6,
      plot_heatmap = FALSE,
      heatmap_first_n = FALSE,
      plot_scatter = FALSE,
      analysis_method = "CPLM" 
    )
    res <- read.table(file.path(out_dir, "all_results.tsv"), header = TRUE, sep = "\t")
    if(nrow(res) > 0){
      res$iteration <- i
      write.table(res, file = file.path(out_dir, "filtered_results.tsv"), sep = "\t", row.names = FALSE)
    }
    
    # Free memory
    rm(fit, res)
    gc()
    
  }, error = function(e){
    warning("Error in jackknife iteration ", i, ": ", e$message)
  })
}

# Combine all results into a single data frame
# List all 'filtered_results.tsv' files in output directories
result_files <- list.files("maslin/", pattern = "filtered_results.tsv", 
                                     full.names = TRUE, recursive = TRUE)

# Read all files and combine them
results_list <- lapply(result_files, read.delim)
combined_results <- do.call(rbind, results_list)

# ACAT function to combine p-values
acat <- function(pvals, weights = NULL) {
  # Prevent issues with p-values exactly 0 or 1
  pvals[pvals < 1e-32] <- 1e-32
  pvals[pvals > 1 - 1e-32] <- 1 - 1e-32
  
  if(is.null(weights)){
    weights <- rep(1/length(pvals), length(pvals))
  }
  
  # Calculate ACAT statistic
  t_stat <- sum(weights * tan((0.5 - pvals) * pi))
  
  # Compute combined p-value using inverse tangent transformation
  p_combined <- 0.5 - atan(t_stat) / pi
  
  # Ensure p-value is between 0 and 1
  p_combined <- min(max(p_combined, 0), 1)
  
  return(p_combined)
}

# Get list of features present in each iteration
features_por_iteracion <- lapply(results_list, function(res) {
  if(!is.null(res)) unique(res$feature) else NULL
})

# Intersection of all features
features_comunes <- Reduce(intersect, features_por_iteracion)

# Filter combined_results to keep only common features
combined_results_filtrados <- combined_results %>% 
  filter(feature %in% features_comunes)

# Now proceed with summary and combination of p-values
combined_pvalues_acat <- combined_results_filtrados %>%
  group_by(feature) %>%
  summarise(
    k = n(),
    acat_p = acat(pval)
  ) %>%
  mutate(acat_p_adj = p.adjust(acat_p, method = "BH"))

coef_summary <- combined_results_filtrados %>%
  group_by(feature) %>%
  summarise(
    mean_coef = mean(coef),
    sd_coef = sd(coef),
    q_2.5_coef = quantile(coef, probs = 0.025),
    q_97.5_coef = quantile(coef, probs = 0.975)
  ) %>%
  left_join(combined_pvalues_acat %>% 
              select(feature, acat_p_adj), by = "feature")

print(head(coef_summary))
write.csv(coef_summary ,"Maaslin2_results.csv")

# Plot
# Assuming results have coefficient columns
ggplot(coef_summary, aes(x = mean_coef, y = -log10(acat_p_adj))) +
  geom_point(aes(color = acat_p_adj < 0.05), size = 2) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue")+ 
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"))+
  labs(x = "Log₁₀ fold change of log1p-transformed TSS-normalized abundances", 
       y = "-log10(adjusted p-value) (FDR via ACAT)", 
       color = "Significance") +
  scale_color_manual(values = c("red", "black"))
```


### VIRULENCE FACTOR ANALYSIS ###
```R
# Set directory and read files from VFDB
files_rute <- "../Virulencia_results"

files <- list.files(path = files_rute, pattern = "VFDB.csv", full.names = T)
data <- do.call(rbind, lapply(files, function(file) {
  sample_id <- sub("_VFDB\\.csv$", "", basename(file))
  df <- read.csv(file, header = F, sep = "\t")
  df$sample <- sample_id
  return(df)
}))

# Change format
colnames(data) <- c("qseqid", "sseqid", "stitle", "pident", "evalue", 
                     "staxids", "sample")
data$sample <- str_remove(data$sample, "eICh24_") # remove part of the name
contenidos_corchetes <- str_extract_all(data$stitle, "\\[([^]]+)\\]")
data$title1 <- sapply(contenidos_corchetes, `[`, 1)
data$title2 <- sapply(contenidos_corchetes, `[`, 2)
data$virulence_id <- gsub("\\s+", "_", gsub("[^[:alnum:]_]", "", data$title1))
head(data)

# Count occurrences of each factor per sample
datos_abundance <- data %>%
  filter(evalue <= 0.001) %>%  # Filter
  group_by(virulence_id, sample) %>%  # Group by virulence factor and sample
  summarise(frecuencia = n(), .groups = "drop") %>%  # Count occurrences
  pivot_wider(names_from = sample, values_from = frecuencia, values_fill = 0) # Convert to wide format

datos_abundance <- as.data.frame(datos_abundance)
row.names(datos_abundance) <- datos_abundance$virulence_id
datos_abundance <- datos_abundance[,-1]

# Filter by minimum prevalence (20% of samples)
prevalencia_min <- 0.2
n_samples <- ncol(datos_abundance)
prevalencia <- rowSums(datos_abundance > 0) / n_samples
datos_filtrados <- datos_abundance[prevalencia >= prevalencia_min, ]

# Minimum mean abundance
abundancia_min <- 10
abundancia_media <- rowMeans(datos_filtrados)
datos_filtrados <- datos_filtrados[abundancia_media >= abundancia_min, ]

# Filter variances
varianzas <- apply(datos_filtrados, 1, var)
datos_filtrados <- datos_filtrados[varianzas < quantile(varianzas, 0.99), ]

######### Differential abundance analysis ##########
# Pseudocount to avoid log(0)
pseudo <- 0.5 * calcNormFactors(datos_filtrados, method="TMM") / colSums(datos_filtrados)
datos_pseudo <- sweep(datos_filtrados, 2, pseudo, "+")

# CLR. Transpose because compositions::clr expects rows = samples
clr_data <- clr(t(datos_pseudo))

# Return to MaAsLin2 expected format (features x samples)
joined_data_2 <- as.data.frame(t(clr_data))

is.data.frame(joined_data_2)  # Must be TRUE
colnames(joined_data_2) <- as.character(colnames(tss_normalized_data_v))
rownames(feno) <- as.character(rownames(feno))
all(colnames(joined_data_2) %in% rownames(feno))  # Should return TRUE

# Select samples, define output directory, run Maaslin2, etc.
# Total number of samples (columns) in the dataset
n_muestras <- ncol(joined_data_2)
results_list_virulence <- vector("list", n_muestras)

set.seed(123)  # For reproducibility

for(i in 1:n_muestras){
  # Define jackknife indices: exclude the i-th sample
  jackknife_indices <- setdiff(1:n_muestras, i)
  
  # Subset of data without the i-th sample
  subset_data <- joined_data_2[, jackknife_indices]
  
  # Adjust metadata to include only the used samples
  subset_metadata <- feno[rownames(feno) %in% colnames(subset_data), ]
  
  # Define an output directory for each iteration
  out_dir <- paste0("maslin_virulence/maaslin2_jackknife_", i)
  
  tryCatch({
    # Run Maaslin2 with the subset
    fit <- Maaslin2(
      input_data = subset_data,
      input_metadata = subset_metadata,
      output = out_dir,
      fixed_effects = c("group"),
      normalization = "NONE",
      transform = "NONE",
      cores = 6,
      plot_heatmap = FALSE,
      heatmap_first_n = FALSE,
      plot_scatter = FALSE,
      analysis_method = "LM"
    )
    
    res <- read.table(file.path(out_dir, "all_results.tsv"), header = TRUE, sep = "\t")
    if(nrow(res) > 0){
      res$iteration <- i
      write.table(res, file = file.path(out_dir, "filtered_results.tsv"), sep = "\t", row.names = FALSE)
    }
    
    # Free memory
    rm(fit, res)
    gc()
    
  }, error = function(e){
    warning("Error in jackknife iteration ", i, ": ", e$message)
  })
}

# Combine all results into a single data frame
# List all 'filtered_results.tsv' files in output directories
result_files_virulence <- list.files("maslin_virulence/", pattern = "filtered_results.tsv", 
                           full.names = TRUE, recursive = TRUE)

# Read all files and combine them
results_list_virulence <- lapply(result_files_virulence, read.delim)
combined_results_virulence <- do.call(rbind, results_list_virulence)


# Get list of features present in each iteration
features_por_iteracion_virulence <- lapply(results_list_virulence, function(res) {
  if(!is.null(res)) unique(res$feature) else NULL
})

# Intersection of all features
features_comunes_virulence <- Reduce(intersect, features_por_iteracion_virulence)

# Filter combined_results to keep only common features
combined_results_filtrados_virulence <- combined_results_virulence %>% 
  filter(feature %in% features_comunes_virulence)

# Now proceed with summary and p-value combination
combined_pvalues_acat_virulence <- combined_results_filtrados_virulence %>%
  group_by(feature) %>%
  summarise(
    k = n(),
    acat_p = acat(pval)
  ) %>%
  mutate(acat_p_adj = p.adjust(acat_p, method = "BH"))

coef_summary_virulence <- combined_results_filtrados_virulence %>%
  group_by(feature) %>%
  summarise(
    mean_coef = mean(coef),
    sd_coef = sd(coef),
    q_2.5_coef = quantile(coef, probs = 0.025),
    q_97.5_coef = quantile(coef, probs = 0.975)
  ) %>%
  left_join(combined_pvalues_acat_virulence %>% 
              select(feature, acat_p_adj), by = "feature")

# Clean data
coef_summary_virulence <- coef_summary_virulence %>%
  mutate(clean_feature = gsub("^X\\.", "", feature),
         clean_feature = gsub("\\.+", ".", clean_feature),
         clean_feature = gsub("\\.$", "", clean_feature)) %>%
  filter(clean_feature %in% rownames(joined_data_2))


print(head(coef_summary_virulence))
write.csv(coef_summary_virulence ,"Maaslin2_results_virulence.csv")

# Plot
# Assuming results have coefficient columns
ggplot(coef_summary_virulence, aes(x = mean_coef, y = -log10(acat_p_adj))) +
  geom_point(aes(color = acat_p_adj < 0.05), size = 2) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue")+ 
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"))+
  labs(x = "Group-associated change from CLR-transformed data", 
       y = "-log10(adjusted p-value) (FDR via ACAT)", color = "Significance") +
  scale_color_manual(values = c("red", "black"))
```

### METABOLIC PATHWAYS ANALYSIS ###
```R
# Assume the file already contains normalized data and pathabundance_cpm.tsv in the same directory
joined_data_4 <- read.delim("pathabundance_cpm.tsv", header = TRUE, sep = "\t", 
                          row.names = 1)
colnames(joined_data_4) <- str_remove(colnames(joined_data_4), "_combined_Abundance")
colnames(joined_data_4) <-  str_remove(colnames(joined_data_4), "eICh24_")
filtro <- grep("UNINTEGRATED", rownames(joined_data_4)) # Remove Unintegrated
joined_data_4 <- joined_data_4[-c(1,filtro),] # Remove Unmapped

# # Filter by minimum prevalence (20% of samples)
prevalencia_min <- 0.2
n_samples_m <- ncol(joined_data_4)
prevalencia_m <- rowSums(joined_data_4 > 0) / n_samples
datos_filtrados_m <- joined_data_4[prevalencia_m >= prevalencia_min, ]

# Minimum mean abundance
abundancia_min <- 10
abundancia_media_m <- rowMeans(joined_data_4)
datos_filtrados_m <- datos_filtrados_m[abundancia_media >= abundancia_min, ]

# Filter variances
varianzas <- apply(datos_filtrados_m, 1, var)
datos_filtrados_m <- datos_filtrados_m[varianzas < quantile(varianzas, 0.99), ]

######### Differential abundance analysis ##########

# TSS normalization
tss_normalized_data_m <- sweep(datos_filtrados_m, 2, 
                             colSums(datos_filtrados_m), FUN = "/")
colSums(tss_normalized_data_m)
joined_data <- tss_normalized_data_m


is.data.frame(joined_data)  # Must be TRUE
colnames(joined_data) <- as.character(colnames(joined_data))
rownames(feno) <- as.character(rownames(feno))
all(colnames(joined_data) %in% rownames(feno))  # Should return TRUE


# Select samples, define output directory, run Maaslin2, etc.
# Total number of samples (columns) in the dataset
n_muestras <- ncol(joined_data)
results_list_metab <- vector("list", n_muestras)

set.seed(123)  # For reproducibility

for(i in 1:n_muestras){
  # Define jackknife indices: exclude the i-th sample
  jackknife_indices <- setdiff(1:n_muestras, i)
  
  # Subset of data without the i-th sample
  subset_data <- joined_data[, jackknife_indices]
  
  # Adjust metadata to include only the used samples
  subset_metadata <- feno[rownames(feno) %in% colnames(subset_data), ]
  
  # Define an output directory for each iteration
  out_dir <- paste0("maslin_metab/maaslin2_jackknife_", i)
  
  tryCatch({
    # Run Maaslin2 with the subset
    fit <- Maaslin2(
      input_data = subset_data,
      input_metadata = subset_metadata,
      output = out_dir,
      fixed_effects = c("group"),
      normalization = "NONE",
      transform = "NONE",
      cores = 6,
      plot_heatmap = FALSE,
      heatmap_first_n = FALSE,
      plot_scatter = FALSE,
      analysis_method = "LM"
    )
    
    res <- read.table(file.path(out_dir, "all_results.tsv"), header = TRUE, sep = "\t")
    if(nrow(res) > 0){
      res$iteration <- i
      write.table(res, file = file.path(out_dir, "filtered_results.tsv"), sep = "\t", row.names = FALSE)
    }
    
    # Free memory
    rm(fit, res)
    gc()
    
  }, error = function(e){
    warning("Error in jackknife iteration ", i, ": ", e$message)
  })
}

# Combine all results into a single data frame
# List all 'filtered_results.tsv' files in output directories
result_files_metab <- list.files("maslin_metab/", pattern = "filtered_results.tsv", 
                                     full.names = TRUE, recursive = TRUE)

# Read all files and combine them
results_list_metab <- lapply(result_files_metab, read.delim)
combined_results_metab <- do.call(rbind, results_list_metab)

# Get list of features present in each iteration
features_por_iteracion_metab <- lapply(results_list_metab, function(res) {
  if(!is.null(res)) unique(res$feature) else NULL
})

# Intersection of all features
features_comunes_metab <- Reduce(intersect, features_por_iteracion_metab)

# Filter combined_results to keep only common features
combined_results_filtrados_metab <- combined_results_metab %>% 
  filter(feature %in% features_comunes_metab)

# Now proceed with summary and p-value combination
combined_pvalues_acat_metab <- combined_results_filtrados_metab %>%
  group_by(feature) %>%
  summarise(
    k = n(),
    acat_p = acat(pval)
  ) %>%
  mutate(acat_p_adj = p.adjust(acat_p, method = "BH"))

coef_summary_metab <- combined_results_filtrados_metab %>%
  group_by(feature) %>%
  summarise(
    mean_coef = mean(coef),
    sd_coef = sd(coef),
    q_2.5_coef = quantile(coef, probs = 0.025),
    q_97.5_coef = quantile(coef, probs = 0.975)
  ) %>%
  left_join(combined_pvalues_acat_metab %>% 
              select(feature, acat_p_adj), by = "feature")

print(head(coef_summary_metab))
write.csv(coef_summary_metab ,"Maaslin2_results_metab.csv")

# Plot
# Assuming results have coefficient columns
ggplot(coef_summary_metab, aes(x = mean_coef, y = -log10(acat_p_adj))) +
  geom_point(aes(color = acat_p_adj < 0.05), size = 2) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue")+ 
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"))+
  labs(x = "Difference in relative abundance (TSS-normalized, no transformación)", 
       y = "-log10(adjusted p-value) (FDR via ACAT)", color = "Significance") +
  scale_color_manual(values = c("red", "black"))
```

# Save filtered results 
```R
filtered_metab <- coef_summary_metab %>%
  filter(
    abs(mean_coef) >= 0.003,                # Large effect
    acat_p_adj < 0.001,                     # FDR significance
    sd_coef < 1,                            # Low variability
    q_2.5_coef * q_97.5_coef > 0            # Confidence interval does not cross 0 (both ends same sign)
  )

filtered_virulence <- coef_summary_virulence %>%
  filter(
    abs(mean_coef) >= 3,                    # Large effect
    acat_p_adj < 0.001,                     # FDR significance
    sd_coef < 1,                            # Low variability
    q_2.5_coef * q_97.5_coef > 0            # Confidence interval does not cross 0 (both ends same sign)
  )

filtered_taxon <- coef_summary %>%
  filter(
    abs(mean_coef) >= 0.3,                  # Large effect
    acat_p_adj < 0.001,                     # FDR significance
    sd_coef < 1,                            # Low variability
    q_2.5_coef * q_97.5_coef > 0            # Confidence interval does not cross 0 (both ends same sign)
  )

write.csv(filtered_taxon, "taxa.csv")
write.csv(filtered_virulence ,"virulence.csv")
write.csv(filtered_metab ,"metab.csv")
```
