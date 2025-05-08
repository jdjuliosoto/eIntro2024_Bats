# eIntro2024 Bats

## Microbiome Characterization of Bat Guano from Camino del Hierro, Arribes del Duero Natural Park
This repository contains the analysis pipeline used to characterize the microbiome composition of bat guano samples collected at Camino del Hierro , within the Arribes del Duero Natural Park.

To reproduce the results presented in our publication, follow the workflow in this order:

1) INSTALL.md – Setup environment and dependencies
2) Preprocessing.md – Quality control and preprocessing of raw sequencing reads
3) Tax_functional_profile.md – Taxonomic profiling using Centrifuge, Functional annotation with HUMAnN and Detection of virulence factors via BLAST against VFDB

Abundace graphs and diversity analyses were done with the kreports from Centrifuge

1) Virus_taxa.md
2) Bacteria_taxa.md
3) Diversity_analyses.md

Functional analyses were done with pathabundance_cpm.tsv from Humann3

1) Functional_analyses.R

