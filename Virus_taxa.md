# Obtain virus orders plots

## Set the working space
```R
# libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(tibble)
library(tidyr)
library(ggtree)
library(tidytree)
library(ape)
library(treeio)
library(taxize)
library(rentrez)
library(zoo)
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

# number of DNA viruses in the samples
sum(superkingdoms_taxid$name == "Viruses")

# Filter by viruses superkingdom the IDs
viruses_id <- superkingdoms_taxid %>%
  filter(name == "Viruses")
```


## Filter all the data by the ID of the domain previously obtained
```R
# Filter orders (genus, class, family, order) by domain (virus, bacteria)
filtered_data <- filtered_data_1 %>%
  filter(taxid %in% viruses_id$taxid)
```

## Taxonomic distribution across samples with relative abundance

```R
# Obtaine absolute abundance matrix
abundance_matrix <- filtered_data %>%
  group_by(sample, name) %>%
  summarise(total_percentage = sum(reads, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = total_percentage)

# Calculate relative abundance (percentage) excluding NA
abundance_matrix_relative <- abundance_matrix %>%
  mutate(across(where(is.numeric), ~ (./ sum(.,na.rm = TRUE)) * 100, .names = ".{col}")) %>%
  replace(is.na(.), 0) # Replace NA with 0 if any empty cells

abundance_matrix_relative_2 <- abundance_matrix_relative[c(1,50:97)] # clean the table

# Save taxa file
write.csv((abundance_matrix_relative_2), "taxa.csv")
```

## Generate heatmap

```R
# Melt the matrix specifying 'name' as the ID variable
melted_data <- melt(abundance_matrix_relative_2, id.vars = "name")

# Rename columns for clarity
colnames(melted_data) <- c("name", "sample", "value")
melted_data$name <- as.factor(trimws(melted_data$name))

# Add dates
m <- paste0(".", sprintf("%02d", seq(1, 48)))
f <- c(
  "2023_2", "2023_2", "2023_2", "2023_2", "2023_3", "2023_3", "2023_2", "2023_1",
  "2023_1", "2023_1", "2023_10", "2023_2", "2023_10", rep("2024_4", 35))
dates <- data.frame(date = f, sample = m)
melted_data <- merge(melted_data, dates, by = "sample", all.x = TRUE)
melted_data$date <- as.yearmon(melted_data$date, format = "%Y_%m")
head(melted_data)

# Custom order for clarity in the graph
unique(melted_data$name)  # Show unique values in name
setdiff(unique(melted_data$name), custom_order)  # Show values not in custom_order
custom_order <- c("Herpesvirales", "Tymovirales", "Ligamenvirales",
                  "Nidovirales", "Picornavirales", "Mononegavirales")
melted_data$name <- factor(melted_data$name, levels = custom_order) # change to custom order

# Plotting
taxon <- ggplot(melted_data, aes(x = name, y = sample, fill = value)) +
  geom_tile() +

  # Color scale
  scale_fill_gradientn(colors = c("white", "blue")) +

  # Theme customization
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0,
                               size = 14, face = "bold"),

    # Move x-axis title down using margin
    axis.title.x = element_text(face = "bold", size = 14, 
                                margin = margin(t = 1, r = 0, b = 10, l = 0)),  # Adjust t (top) to move down

    axis.title.y = element_text(face = "bold", size = 14),

    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.1, "cm"),
    legend.position = "bottom"
  ) +

  # Labels
  labs(x = "Taxon", y = "Sample", fill = "Relative abundance (%)")

# Show plot
print(taxon)
# Save as PDF at 12 x 10 inches
```


## Barplot graph

```R
# Define colors
custom_colors <- c("#1F77B4", "#FF7F0E", "#D62728", "#2CA02C", "#8C564B", "#7a0177")

# Plot
prop <- ggplot(filtered_data, 
               aes(x = sample, y = percentage, fill = trimws(name))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1,
                                   size = 12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),

        panel.grid.major = element_blank(),       
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.1, "cm"),
        legend.position = "bottom") +
  labs(x = "Sample", y = "Proportion (%)", fill = "Taxon")

# Show plot
print(prop)
# Save as PDF at 12 x 10 inches
```

## Cladogram

```R
# Convert first column to row names
abundance_matrix <- column_to_rownames(abundance_matrix, var = "name")

# Distance matrix
dd <- dist(scale(abundance_matrix), method = "euclidean")

# Cluster by ward.D2
hc <- hclust(dd, method = "ward.D2")

# Convert hc to phylo object
phylo_tree <- as.phylo(hc)

# Plot the cladogram
p <- ggtree(phylo_tree, layout = "rectangular") + 
  theme_tree() +
  coord_flip() + 
  geom_text2(aes(label = round(branch.length, 2)), 
             hjust = 0.5, size = 3, 
             vjust = 1) +
  scale_x_reverse()

# Combine plots in column 1 (p above taxon)
column_1 <- p / taxon + plot_layout(ncol = 1, heights = c(1,3.5)) # p above taxon

# Combine columns using patchwork
layout <- column_1 | prop  # Column 1 next to column 2

# Adjust height and width proportions
final_plot <- layout + plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

# Show result
print(final_plot)

# Size: 24 x 12
```
