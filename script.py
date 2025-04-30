#!/usr/bin/env python3

# Filter taxa by score: score ≤ 300 and match lengths ≤ 60 bp

import pandas as pd
import os

# Define the folder where the input files are stored
input_folder = "/home/user/proyect/results/centrifuge/report"  # Change this to your actual input path
output_folder = "/home/user/proyect/results/centrifuge/filtered_results"  # Folder to save filtered files

# Iterate over numbers from 1 to 48
for i in range(1, 49):  
    # Build the input file name
    input_file = f"results_eICh24_{i}_bacteria.txt"
    input_path = os.path.join(input_folder, input_file)
    
    # Check if the file exists
    if not os.path.exists(input_path):
        print(f"File {input_file} does not exist. Skipping...")
        continue
    
    # Load the file
    try:
        df = pd.read_csv(input_path, sep="\t")
    except Exception as e:
        print(f"Could not load file {input_file}. Error: {e}")
        continue
    
    # Apply combined filters
    df_filtered = df[(df["score"] <= 300) & (df["hitLength"] <= 60)]
    
    # Build the output file name
    output_file = f"filtered_results_eICh24_{i}_bacteria.tsv"
    output_path = os.path.join(output_folder, output_file)
    
    # Export filtered data
    df_filtered.to_csv(output_path, sep="\t", index=False)
    print(f"Filtered file saved: {output_file}")
