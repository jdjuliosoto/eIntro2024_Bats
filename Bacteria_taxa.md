# Obtain bacteria orders plots

## Prepare the data

```R
# BiocManager libraries
library(treeio)
library(taxize)
library(rentrez)
library(ggtree)
library(tidytree)

# CRAN libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tibble)
library(tidyr)
library(reshape2)
library(ape)
library(zoo)
library(treemap)
library(viridis)
library(patchwork)

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

## Filter all the data by the ID of the domain previously obtained
```R
# Filter orders (genus, class, family, order) by domain (virus, bacteria)
filtered_data <- filtered_data_1 %>%
  filter(taxid %in% bacteria_id$taxid)
```


## Taxonomic distribution across samples with relative abundance

```R
# Obtaine absolute abundance matrix
abundance_matrix <- filtered_data %>%
  group_by(sample, name) %>%
  summarise(total_percentage = sum(reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = total_percentage)

# Calculate relative abundance excluding NA
abundance_matrix_relative <- abundance_matrix %>%
  mutate(across(where(is.numeric), ~ (./ sum(.,na.rm = TRUE)) * 100, 
                .names = ".{col}")) %>% replace(is.na(.), 0) # Replace NA with 0 if any empty cells

abundance_matrix_relative_2 <- abundance_matrix_relative %>%
  select(name, starts_with("."))

```

## Preparing data for bar charts

```R
# Convert matrix to long format
abundance_long <- abundance_matrix_relative_2 %>%
  pivot_longer(cols = starts_with("."), names_to = "sample", values_to = "abundance")

# clean names
abundance_long$name <- trimws(abundance_long$name)
abundance_long$sample <- trimws(abundance_long$sample)

# convert to data.frame for better flexibility
abundance_long <- data.frame(
  name = abundance_long$name,
  abundance = abundance_long$abundance,
  sample = abundance_long$sample
)

clado_abundance <- abundance_long # for cladogram


# Verify the result
head(abundance_long)

# write taxa doc
write.csv(abundance_long, "taxa.csv")

# Add factor for comparison 
abundance_long_filt <- abundance_long %>%
  mutate(factor_abundance = cut(abundance, 
                                breaks = c(-Inf, 0.1, 1, Inf), 
                                labels = c("Factor 3", "Factor 2", "Factor 1")))
head(abundance_long_filt)

# Add dates
Sys.setlocale("LC_TIME", "C")

m <- paste0(".", sprintf("%02d", seq(1, 48)))
f <- c(
  "2023_2", "2023_2", "2023_2", "2023_2", "2023_3", "2023_3", "2023_2", "2023_1",
  "2023_1", "2023_1", "2023_10", "2023_2", "2023_10", rep("2024_4", 35))
dates <- data.frame(date = f, sample = m)
abundance_long <- merge(abundance_long, dates, by = "sample", all.x = TRUE)
abundance_long$date <- as.yearmon(abundance_long$date, format = "%Y_%m")
```

## Plot
``` R
# Filter abundance >= 2
abundance_long_filt <- abundance_long %>%
  filter(abundance >= 2)
head(abundance_long_filt)

# Colors
my_colors <- c(brewer.pal(12, "Paired"), brewer.pal(1, "Set1"))

# Plot
prop <- ggplot(abundance_long_filt, 
               aes(x = sample, y = abundance, fill = trimws(name))) +
  scale_fill_manual(values = my_colors)+
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  facet_grid(cols = vars(date), scales = "free", space = "free") +
  theme_minimal()+
  theme(strip.background = element_rect(fill = "grey"), 
        strip.placement = "outside",
        strip.text = element_text(face = "bold", size = 14),
        
        
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1,
                                   size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.position = "right") +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Taxon")

# Show the plot
print(prop)

# Filter abundance between 1 and 1.99
abundance_long_filt2 <- abundance_long %>%
  filter(abundance > 1 & abundance <= 2 )
head(abundance_long_filt2)

# Colors
my_colors2 <- c(brewer.pal(12, "Paired"), brewer.pal(9, "Set1"))

# Plot
prop2 <- ggplot(abundance_long_filt2, 
               aes(x = sample, y = abundance, fill = trimws(name))) +
  scale_fill_manual(values = my_colors2)+
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  facet_grid(cols = vars(date), scales = "free", space = "free") +
  theme_minimal()+
  theme(strip.background = element_rect(fill = "grey"), 
        strip.placement = "outside",
        strip.text = element_text(face = "bold", size = 14),
        
        
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1,
                                   size = 20),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.position = "right") +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Taxon")

# Show the plot
print(prop2)

```
## Cladogram

```R
# Filter abundance <= 1
abundance_long_filt3 <- clado_abundance %>%
  filter(abundance <= 1)
head(abundance_long_filt)

# Format wide
abundance_long_filt3 <- abundance_long_filt3 %>%
  pivot_wider(names_from = sample, values_from = abundance) %>% 
  replace(is.na(.), 0)

# Convert first column to row names
abundance_long_filt3$name <- trimws(abundance_long_filt3$name)
abundance_long_filt3 <- column_to_rownames(abundance_long_filt3, var = "name")

abundance_long_filt_db <- abundance_long_filt3 # for heatmap

# Distance matrix
dd <- dist(abundance_long_filt3, method = "euclidean")
hc <- hclust(dd, method = "average")

# Convert hc to phylo object
phylo_tree <- as.phylo(hc)

# Scale branch lengths using log10
phylo_tree <- compute.brlen(phylo_tree, fun = function(x) log10(x))

# Plot tree
p <- ggtree(phylo_tree, branch.length = "branch.length") +
  theme_tree2()+
  geom_tree(size = 1 , aes(color = branch.length)) + # branch color
  geom_tiplab(linesize = 1, size = 8) +    
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(size = 20),
        axis.ticks.x = element_line(size = 2)) +
  scale_color_gradient(name = "Log10 euclidean \n distance")+
  scale_x_ggtree()

# Plot heatmap
b3 <- gheatmap(p, abundance_long_filt_db, offset = 0.3, width = 0.5, 
         font.size = 4.5, hjust = 1, colnames_angle = 90) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
  ) +
  scale_y_continuous(expand=c(0, 10)) +
  scale_fill_viridis_c(option="D", direction = -1,
                       name="Relative \nabundance")
```


## Layout: combine plots

```R
layout <- (prop / prop2) | b3

# Adjust height and width proportions
final_plot <- layout + plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))


# Show result
print(final_plot)

# 55 X 45
```

