# Metabolic pathways analyses

## Set the working space

```R
# Libraries
library(tidyverse)
library(pheatmap)
library(vegan)
library(Maaslin2)
library(zoo)
library(stringr)
library(grid)
library(patchwork)
library(tibble)
library(viridis)
library(RColorBrewer)
library(forcats)
library(stats)
library(mclust)
library(ggplotify)
library(patchwork)
library(readr)
library(cowplot)
library(Maaslin2)

### Load HUMAnN data (combined file)
# Assume the file already contains normalized data
joined_data <- read.delim("pathabundance_cpm.tsv", header = TRUE, sep = "\t", 
                          row.names = 1)
colnames(joined_data) <- str_remove(colnames(joined_data), "_combined_Abundance")
colnames(joined_data) <-  str_remove(colnames(joined_data), "eICh24_")
filtro <- grep("UNINTEGRATED", rownames(joined_data)) # Remove Unintegrated
joined_data <- joined_data[-c(1,filtro),] # Remove Unmapped
joined_data_2 <- joined_data # copy for species-level abundance
joined_data_3 <- joined_data # copy for maaslin2
```
## Set metadata

```R
# Filter by total abundances and exclude species-specific abundances
filtro2 <- grep("\\|", row.names(joined_data))
joined_data <- joined_data[-filtro2,]

# Metadata
m <- sprintf("%02d", seq(1, 48))
f <- c(
  "2023_2", "2023_2", "2023_2", "2023_2", "2023_3", "2023_3", "2023_2", "2023_1",
  "2023_1", "2023_1", "2023_10", "2023_2", "2023_10", rep("2024_4", 35))
dates <- data.frame(date = f, sample = m) # data.frame with dates
dates$date <- as.yearmon(dates$date, format = "%Y_%m") # change format of dates
dates <-  dates %>%
  mutate(date = str_remove(date, " ")) # remove spaces
dates$date <- as.factor(dates$date) # convert to factor
row.names(dates) <- dates$sample # row names as sample names

# Add seasons
vector1 <-c(rep("Hibernation_season", 13), rep("Breeding_season", 35))
dates$date2 <- as.factor(vector1)
metadata <- dates[,-2] # change names

### Data transformation (log transform for visualization)
joined_data_norm <- joined_data %>% 
  mutate_all(function(x) log10(x + 1))

joined_data_norm2 <- joined_data_norm # copy for plots
```
             
## Identify the 100 most abundant metabolic pathways across all samples

```R
## change the order of pathways by decreasing order
top_pathways <- rowMeans(joined_data_norm) %>% 
  sort(decreasing = TRUE) %>% 
  head(100)

# split metadata            
hiber_sea <- metadata %>% 
  filter(date2 == "Hibernation_season")

bree_sea <- metadata %>% 
  filter(date2 == "Breeding_season")

# select by hibernation season
selected_data <- joined_data_norm[rownames(joined_data_norm) %in% 
                                    names(top_pathways), 
                                  colnames(joined_data_norm) %in%
                                   rownames(hiber_sea) ]

# select by breeding season
selected_data_2 <- joined_data_norm2[rownames(joined_data_norm2) %in% 
                                    names(top_pathways), 
                                  colnames(joined_data_norm2) %in%
                                    rownames(bree_sea) ]

# write csv with results
write.csv(selected_data, "rutes.csv")
```

#### Here we add the OG manually to rutes.csv. Make sure to define a column named id with the route values ​​(e.g. 1CMET2-PWY). There should be 100 routes in total.
#### Delete the route definition and change the internal periods to hyphens to match the routes.csv table: e.g., ANAEROFRUCAT.PWY..homolactic.fermentation -> ANAEROFRUCAT-PWY

             
## Maaslin2 analysis

```R
### Subsampling to compare seasons
# Filter by total abundances and exclude species-specific abundances
filtro3 <- grep("\\|", row.names(joined_data_3))
joined_data_3 <- joined_data_3[-filtro3,]


# Total number of samples (columns) in the dataset
n_samples <- ncol(joined_data_3)
results_list <- vector("list", n_samples)

set.seed(123)  # For reproducibility
print(dim(fit$results)) # see dimensions

# Function that select samples, define directory and run Maaslin2
for(i in 1:n_samples){
  # Define jackknife indices: exclude the i-th sample
  jackknife_indices <- setdiff(1:n_samples, i)
  # Subset of data without the i-th sample
  subset_data <- joined_data_3[, jackknife_indices]
  # Adjust metadata to include only used samples
  subset_metadata <- metadata[rownames(metadata) %in% colnames(subset_data), ]
  # Define output directory for each iteration
  out_dir <- paste0("maslin/maaslin2_jackknife_", i)
  tryCatch({
    # Run Maaslin2 with the subset
    fit <- Maaslin2(
      input_data = subset_data,
      input_metadata = subset_metadata,
      output = out_dir,
      fixed_effects = c("date2"),
      normalization = "TSS",
      transform = "LOG",
      cores = 5,
      plot_heatmap = FALSE,
      heatmap_first_n = FALSE,
      plot_scatter = FALSE,
      analysis_method = "CPLM"
    )
    # Read results
    res <- read.table(file.path(out_dir, "all_results.tsv"), header = TRUE, sep = "\t")
    if(nrow(res) > 0){
      res$iteration <- i
      results_list[[i]] <- res
    } else {
      warning("Jackknife iteration ", i, ": results file is empty.")
      results_list[[i]] <- NULL
    }
  }, error = function(e){
    warning("Error in jackknife iteration ", i, ": ", e$message)
    results_list[[i]] <- NULL
  })
}


# Combine all results into one data frame
combined_results <- do.call(rbind, results_list)

# ACAT function to combine p-values
acat <- function(pvals, weights = NULL) {
  # Avoid issues with exact 0 or 1 p-values
  pvals[pvals < 1e-32] <- 1e-32
  pvals[pvals > 1 - 1e-32] <- 1 - 1e-32
  if(is.null(weights)){
    weights <- rep(1/length(pvals), length(pvals))
  }

  # Calculate ACAT statistic
  t_stat <- sum(weights * tan((0.5 - pvals) * pi))
  # Calculate combined p-value from inverse tangent transformation
  p_combined <- 0.5 - atan(t_stat) / pi
  # Ensure p-value is between 0 and 1
  p_combined <- min(max(p_combined, 0), 1)
  return(p_combined)
}

# Get feature list per iteration
features_per_iteration <- lapply(results_list, function(res) {
  if(!is.null(res)) unique(res$feature) else NULL
})

# Intersection of all features
features_common <- Reduce(intersect, features_per_iteration)

# Filter combined_results to keep only common features
combined_results_filtered <- combined_results %>% 
  filter(feature %in% features_common)

# Proceed with summary and combining p-values
combined_pvalues_acat <- combined_results_filtered %>%
  group_by(feature) %>%
  summarise(
    k = n(),
    acat_p = acat(pval)
  ) %>%
  mutate(acat_p_adj = p.adjust(acat_p, method = "BH"))

# Combine results
coef_summary <- combined_results_filtered %>%
  group_by(feature) %>%
  summarise(
    mean_coef = mean(coef),
    sd_coef = sd(coef),
    q_2.5_coef = quantile(coef, probs = 0.025),
    q_97.5_coef = quantile(coef, probs = 0.975)
  ) %>%
  left_join(combined_pvalues_acat %>% 
              select(feature, acat_p_adj), by = "feature")

# inspect results
print(head(coef_summary))

# write results
write.csv(coef_summary ,"Maaslin2_results.csv")

# Volcanoplot
ggplot(coef_summary, aes(x = mean_coef, y = -log10(acat_p_adj))) +
  geom_point(aes(color = acat_p_adj < 0.05), size = 2) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue")+ 
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"))+
  labs(x = "Regression coefficient", y = "-log10(ACAT)", color = "Significance") +
  scale_color_manual(values = c("red", "black"))
# Save 10 X 16
```
#### Make sure to define a column named id with the route values ​​(e.g., 1CMET2-PWY). There should be 100 routes in total. 
#### Delete the route definition and change the internal periods to hyphens to match the routes.csv table: e.g., ANAEROFRUCAT.PWY..homolactic.fermentation -> ANAEROFRUCAT-PWY
#### Route PWY-6527 has values ​​of 0 according to the fitted model for mean_coef, sd_coef, q_2.5_coef, q_97.5_coef, and acat_p_adj, so maaslin2 automatically removes them from the results table. 

## Heat map plot

```R

# change order of data (100 metabolic pathways) from hibernation season by custom order
custom_order <- c("01", "02","03","04","05","06", "07", "08", "09",
                  "10", "11","12","13") 
selected_data <- selected_data[, match(custom_order, colnames(selected_data))]

# Read data with OG anotation and Maaslin2 results
maslin_res <- read.csv("Maaslin2_results.csv", header = T, sep = ";")
rutes_ <- read.csv("rutes.csv", header = TRUE, sep = ";")

# join both results by id
rutes <- inner_join(rutes_, maslin_res, by = "id")

# Annotation for fist plot
annotation <- data.frame( Molecule = trimws(as.factor(rutes$Molecule)),
                          GO2 = as.factor(rutes$lvl2),
                          coef = rutes$mean_coef,
                          sig = rutes$acat_p_adj) # to data.frame
row.names(annotation) <- row.names(selected_data) # add row names as the rute

# Sort rows
order <- order(annotation$Molecule)
selected_data <- selected_data[order, ]
annotation <- annotation[order, ]

# Merge tables
# Hibernation season
annotation_t <- tibble::rownames_to_column(annotation, var = "rowname")
selected_data_t <- tibble::rownames_to_column(selected_data, var = "rowname")
merged_data <- merge(annotation_t, selected_data_t, by = "rowname")
merged_data <- merged_data[order, ]
# write csv
write.csv(merged_data, "hibernation_order.csv")


# Annotation for second plot
annotation_2 <- data.frame( Molecule = trimws(as.factor(rutes$Molecule)),
                           GO2 = as.factor(rutes$lvl2),
                           coef = rutes$mean_coef,
                           sig = rutes$acat_p_adj)
row.names(annotation_2) <- row.names(selected_data) # add row names as the rute

# Sort rows
order_2 <- order(annotation_2$Molecule)
selected_data_2 <- selected_data_2[order_2, ]
annotation_2 <- annotation_2[order_2, ]

# Merge tables
# Breeding season
annotation_t_2 <- tibble::rownames_to_column(annotation_2, var = "rowname")
selected_data_2_t <- tibble::rownames_to_column(selected_data_2, var = "rowname")
merged_data_2 <- merge(annotation_t_2, selected_data_2_t, by = "rowname")
merged_data_2 <- merged_data_2[order_2, ]
# write csv
write.csv(merged_data_2, "breeding_order.csv")


### Visualization: Heatmap of the most abundant pathways
# hibernation season
fig1 <- pheatmap(selected_data,
                 color = colorRampPalette(c( "#005a32","white","#fc4e2a"))(9),
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 fontsize = 25,
                 fontsize_col = 30,
                 fontsize_row = 25,
                 show_rownames = FALSE,
                 cellwidth = 32,
                 cellheight = 16,
                 border_color = "white",
                 legend = FALSE,
                 annotation_legend = FALSE,
                 annotation_row = subset(annotation, select = c(-coef, -sig)),
                 annotation_colors = list(
                   Molecule = c("Amino acids" = "#0c2c84", "Aromatic acids" = "#225ea8", 
                                "Cell wall components" = "#1d91c0", "Cofactors" = "#41b6c4", 
                                "Isoprenoids" = "#7fcdbb", "Lipids" = "#c7e9b4",
                                "Nitrogen compounds" = "#fed976", "Nucleotides" = "#feb24c", 
                                "Organic acids" = "#fd8d3c", "Sugars" = "#fc4e2a", 
                                "Sulfate" = "#e31a1c", "Vitamins" = "#b10026"),
                   GO2 = c("Biosynthesis" = "#005a32", "Central_pathways" = "#238b45", 
                           "Degradation" = "#41ab5d", "Fermentation" = "#74c476", 
                           "Photosynthesis" = "#a1d99b", "Recycling" = "#c7e9c0"))
# breeding season
fig2 <- pheatmap(selected_data_2,
                 color = colorRampPalette(c( "#005a32","white","#fc4e2a"))(9),
                 cluster_rows = FALSE,
                 cluster_cols = FALSE, 
                 fontsize = 25,
                 fontsize_col = 30,
                 fontsize_row = 25,
                 show_rownames = FALSE,
                 cellwidth = 32,
                 cellheight = 16,
                 border_color = "white",
                 annotation_row = subset(annotation_2, select = c(-coef,-sig)),
                 annotation_colors = list(
                   Molecule = c("Amino acids" = "#0c2c84", "Aromatic acids" = "#225ea8", 
                                "Cell wall components" = "#1d91c0", "Cofactors" = "#41b6c4", 
                                "Isoprenoids" = "#7fcdbb", "Lipids" = "#c7e9b4",
                                "Nitrogen compounds" = "#fed976", "Nucleotides" = "#feb24c", 
                                "Organic acids" = "#fd8d3c", "Sugars" = "#fc4e2a", 
                                "Sulfate" = "#e31a1c", "Vitamins" = "#b10026"),
                   GO2 = c("Biosynthesis" = "#005a32", "Central_pathways" = "#238b45", 
                           "Degradation" = "#41ab5d", "Fermentation" = "#74c476", 
                           "Photosynthesis" = "#a1d99b", "Recycling" = "#c7e9c0")))
# title text
grid.text("log10(x + 1) of the abundance", 
          x = 0.98, y = 0.99, gp = gpar(fontsize = 20, fontface = "bold"))


# Combine plots using patchwork
fig1_gg <- as.ggplot(fig1)
fig2_gg <- as.ggplot(fig2)

```

## Functional diversity assessment using Shannon index

```R

# diversity index
diversity_shannon <- diversity(t(joined_data), index = "shannon")
diversity_simpson <- diversity(t(joined_data), index = "simpson")

# Convert to data frame for plotting
diversity_df <- data.frame(Sample = colnames(joined_data), 
                           Shannon = diversity_shannon, Simpson = diversity_simpson)

# Calculate Pielou's J
diversity_df$richness <- colSums(joined_data)
diversity_df$Pielou <- diversity_df$Shannon / log(diversity_df$richness)

# Bar plot of functional diversity for hibernation season
fig3_data <- diversity_df[colnames(selected_data),]

# plot
fig3 <- ggplot(fig3_data, aes(x = Sample, y = Shannon, fill = Shannon)) +
  geom_bar(stat = "identity", width = 0.88, position = position_dodge(width = 0)) +
  theme_minimal() +
  labs(y = "Shannon Index", x ="") + 
  theme(
    axis.text.x = element_text(size = 20), 
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
    aspect.ratio = 1/2.7,
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )


# Bar plot of functional diversity for breeding season
fig4_data <- diversity_df[colnames(selected_data_2),]

# plot
fig4 <- ggplot(fig4_data, aes(x = Sample, y = Shannon, fill = Shannon)) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") +
  theme_minimal() +
  labs(y = "Shannon Index", x = "") + 
  theme(
    axis.text.x = element_text(size = 20), 
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
    aspect.ratio = 1/7.1,
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    plot.margin = margin(0, 0, 0, 0)
  )


pielou_ind <- ggplot(fig4_data, aes(x = Sample, y = Pielou, fill = Pielou)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Sample", y = "Pielou (J´)", fill = "Pielou J´")

layout <- fig3 / fig4
# Save at 16 X10
```
## Maaslin2 log-fold change
```R
# Coefficients from maaslin2 model in the same order as the breeding season figure
rutes_ord <- rutes[order_2, ]

# plot
fig5_gg <- ggplot(rutes_ord, aes(x = mean_coef, y = seq_along(mean_coef))) +
  geom_col(aes(fill = acat_p_adj < 0.01), orientation = "y") +  # Horizontal bars
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey"),
                    name = "Significant") +  # Colors
  theme_minimal() +
  theme(axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)) +
  labs(x = "log-fold change", y = NULL)  # Labels


# Combine figures in one layout
space <- plot_spacer()
combined <- fig3 + fig4 + space + fig1_gg + fig2_gg  + fig5_gg +
  plot_layout(ncol = 3, widths = c(2, 4, 1), nrow = 2, heights = c(1,13)) +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 20))

# Show combined plot
print(combined) 
# Save at 40 X 30

```

## Species contributing to the most abundant pathways

```R
# Filter (those containing "|")
df <- joined_data_2[filtro2,]
df <- rownames_to_column(df, var = "Pathway")

# Split into two columns: Pathway and Species
df_species <- df %>%
  separate(Pathway, into = c("Ruta", "Especie"), sep = "\\|", remove = TRUE)

# Transform dataframe to long format
df_long <- df_species %>%
  pivot_longer(cols = -c(Ruta, Especie), names_to = "Muestra", values_to = "Abundancia")

# Select a specific pathway to plot
row_means <- rowMeans(selected_data_2) %>%
  sort(decreasing = T)
head(row_means, 28)

# Most abundant pathway
ruta_elegida <- "PWY-7229: superpathway of adenosine nucleotides de novo biosynthesis I" 
df_filtrado <- df_long %>%
  filter(Ruta == ruta_elegida)

# Select only taxa with contribution >= 2
df_filtrado <- df_filtrado %>%
  filter(Abundancia >= 2) %>%
  mutate(Abun_norm = log10(Abundancia + 1))

# change format of the taxa names
df_filtrado <- df_filtrado %>%
  separate(Especie, into = c("Genero", "Sp"), sep = "s\\__", remove = TRUE)

# Colors
my_colors <- c('#fdae61','#276419','#b8e186','#7fbc41',
               '#543005','#8c510a','#313695','#4575b4',
               '#40004b','#762a83','#9970ab','#c2a5cf',
               '#bf812d','#dfc27d',"#737373","#000000",
               '#f6e8c3','#b10026','#c7eae5','#80cdc1',
               '#35978f','#01665e','#003c30','#8e0152',
               '#f46d43','#4d9221','#fee090','#ffffbf',
               '#c51b7d','#de77ae','#f1b6da','#fde0ef',
               '#e31a1c','#e7298a','#f7f7f7','#e6f5d0',
               '#e0f3f8','#abd9e9',"#d9d9d9", '#74add1',
               "yellow", "red", "blue", "pink", "green", "purple")

# Create bar chart
g1 <- ggplot(df_filtrado, aes(x = fct_reorder(Muestra, Abun_norm, .fun = sum, .desc = TRUE), 
                        y = Abun_norm, fill = Sp)) +
  stat_summary(geom = "area", fun = sum, aes(group = 1), 
               fill = "#034e7b", alpha = 0.4) +
  geom_bar(stat = "identity", position = "stack", 
           width = 0.7) +
  scale_fill_manual(values = my_colors)+
  theme_minimal() +
  labs(title = paste("Species contribution (>= 2) to", ruta_elegida),
       x = "Sample",
       y = "log10(CPM abundance + 1)",
       fill = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black"))
# Save at 16 X10



# Pathway 2, this can be used also for:
# "ARGSYNBSUB-PWY: L-arginine biosynthesis II (acetyl cycle)"
# "CALVIN-PWY: Calvin-Benson-Bassham cycle"

ruta_elegida <- "PWY-3841: folate transformations II (plants)"
 
df_filtrado <- df_long %>%
  filter(Ruta == ruta_elegida)

# Select only taxa with contribution >= 2
df_filtrado <- df_filtrado %>%
  filter(Abundancia >= 2) %>%
  mutate(Abun_norm = log10(Abundancia + 1))

# Change format of taxa name
df_filtrado <- df_filtrado %>%
  separate(Especie, into = c("Genero", "Sp"), sep = "s\\__", remove = TRUE)

# Colors
my_colors <- c("#737373","#000000", '#fee090','#4575b4',
               '#40004b','#762a83','#9970ab','#c2a5cf',
               '#543005','#8c510a','#bf812d','#dfc27d',
               '#f6e8c3','#b10026','#c7eae5','#4d9221',
               '#35978f','#01665e','#003c30','#8e0152',
               '#9970ab','#fdae61','#313695','#e31a1c',
               '#c51b7d','#de77ae','#f1b6da','#fde0ef',
               '#ffffbf','#e7298a','#f7f7f7','#e6f5d0',
               '#80cdc1','#276419','#b8e186','#7fbc41',
               'yellow','#abd9e9',"#d9d9d9", '#74add1',
               '#e0f3f8','#CC0000','#f46d43','#FFC711', "grey")

# Plot
g2 <- ggplot(df_filtrado, aes(x = fct_reorder(Muestra, Abun_norm, .fun = sum, .desc = TRUE), 
                              y = Abun_norm, fill = Sp)) +
  stat_summary(geom = "area", fun = sum, aes(group = 1), 
               fill = "#034e7b", alpha = 0.4) +
  geom_bar(stat = "identity", position = "stack", 
           width = 0.7) +
  scale_fill_manual(values = my_colors)+
  theme_minimal() +
  labs(title = paste("Species contribution (>= 2) to", ruta_elegida),
       x = "Sample",
       y = "log10(CPM abundance + 1)",
       fill = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black"))
# Save at 16x10

# Combined
combinado2 <- g1 + g2 + g3 + g4 +
  plot_layout(ncol = 2, nrow = 2 ) +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size = 20))

# Show combined plot
print(combinado2) 
# Save at 30 X 22
```


### Clustering analysis (PCA based on functional abundances)
```R

# PCA
pca_result <- prcomp(t(joined_data_norm), scale. = TRUE)
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Show explained variance
print(explained_variance)
sum(explained_variance)

# GMM clustering (Gaussian Mixture Models)
# Description: Models data as a mixture of Gaussian distributions.
# Unlike k-means, GMM allows clusters of different sizes, orientations, and densities.
set.seed(123)
pca_data <- as.data.frame(pca_result$x)

# Cluster analysis
gmm_result <- Mclust(pca_data[, c("PC1", "PC2")])
pca_data$group <- as.factor(gmm_result$classification)
gmm_result$BIC

# PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) + 
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2, linewidth = 1) +
  geom_text(aes(label = rownames(pca_data)), color = "black", hjust = 1.25) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(title = paste("Log-likelihood: ", round(gmm_result$loglik,2),
                "
Bayesian Information Criterion: ", round(gmm_result$bic),
                "
 (Next top BIC: VEV -646.05     VVV -649.85     VVE -651.43"), 
       x = "Principal Component 1 (34%)", y = "Principal Component 2 (11.4%)", color = "Group")
# Export to 16X 10

```
