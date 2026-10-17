#!/bin/bash
echo "=== Executing Per-Region PCR Coverage Analysis ==="

python3 - << 'EOF'
import os
import glob
import pandas as pd

input_dir = 'pcr_coverage'
output_file = 'pcr_failed_regions.csv'
clinical_threshold = 30.0

if not os.path.exists(input_dir):
    print(f"Error: Directory '{input_dir}' not found.")
    exit(1)

file_pattern = os.path.join(input_dir, '*.regions.bed')
bed_files = glob.glob(file_pattern)

if not bed_files:
    print(f"Error: No .regions.bed files found in '{input_dir}'.")
    exit(1)

failed_regions = []
total_regions_checked = 0

for file_path in bed_files:
    sample_name = os.path.basename(file_path).replace('.regions.bed', '')
    
    try:
        df = pd.read_csv(file_path, sep='\t', header=None, usecols=[0, 1, 2, 3, 4],
                         names=['chrom', 'start', 'end', 'region', 'depth'])
        
        total_regions_checked += len(df)
        low_cov = df[df['depth'] < clinical_threshold]
        
        for _, row in low_cov.iterrows():
            failed_regions.append({
                'PCR_Sample': sample_name,
                'Gene_Region': row['region'],
                'PCR_Depth': row['depth']
            })
    except Exception as e:
        print(f"Warning: Could not process {file_path}. Reason: {e}")

df_failed = pd.DataFrame(failed_regions)

if df_failed.empty:
    print(f"Analysis complete. No regions found below {clinical_threshold}X in PCR samples.")
else:
    df_failed = df_failed.sort_values(by=['Gene_Region', 'PCR_Depth'])
    df_failed.to_csv(output_file, index=False)
    
    print(f"--------------------------------------------------")
    print(f"Per-Region Analysis Results (< {clinical_threshold}X)")
    print(f"--------------------------------------------------")
    print(f"Total Regions Evaluated : {total_regions_checked}")
    print(f"Regions Below Threshold : {len(df_failed)}")
    print(f"\nBreakdown by Region:")
    print(df_failed['Gene_Region'].value_counts().to_string())
    print(f"--------------------------------------------------")
    print(f"Data exported to '{output_file}'.")
EOF
