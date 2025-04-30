#!/usr/bin/env python3

# Filter report.tsv files using taxIDs from the filtered results files

import pandas as pd
import os

# Define folders where files are stored
filtered_results_folder = "/home/user/proyect/results/centrifuge/filtered_results"  # Folder with filtered_results files
report_folder = "/home/user/proyect/results/centrifuge/report"  # Folder with report.tsv files
output_folder = "/home/user/proyect/results/centrifuge/filtered_reports"  # Folder to save filtered reports

# Iterate over numbers from 1 to 48
for i in range(1, 49):  
    # Build the filename for the filtered results file
    filtered_results_file = f"filtered_results_eICh24_{i}_bacteria.tsv"
    filtered_results_path = os.path.join(filtered_results_folder, filtered_results_file)
    
    # Check if the filtered file exists
    if not os.path.exists(filtered_results_path):
        print(f"Filtered file {filtered_results_file} does not exist. Skipping...")
        continue
    
    # Load the filtered results file to extract taxIDs
    try:
        filtered_results = pd.read_csv(filtered_results_path, sep="\t")
    except Exception as e:
        print(f"Could not load filtered file {filtered_results_file}. Error: {e}")
        continue
    
    # Extract unique taxIDs from the filtered results
    tax_ids = filtered_results['taxID'].unique()
    
    # Build the filename for the report file
    report_file = f"report_eICh24_{i}_bacteria.txt"
    report_path = os.path.join(report_folder, report_file)
    
    # Check if the report file exists
    if not os.path.exists(report_path):
        print(f"Report file {report_file} does not exist. Skipping...")
        continue
    
    # Load the report file
    try:
        report = pd.read_csv(report_path, sep="\t")
    except Exception as e:
        print(f"Could not load report file {report_file}. Error: {e}")
        continue
    
    # Filter the report using the extracted taxIDs
    report_filtered = report[report['taxID'].isin(tax_ids)]
    
    # Build output file name
    output_report_file = f"filtered_report_eICh24_{i}_bacteria.tsv"
    output_report_path = os.path.join(output_folder, output_report_file)
    
    # Export filtered report
    report_filtered.to_csv(output_report_path, sep="\t", index=False)
    print(f"Filtered report saved: {output_report_file}")
