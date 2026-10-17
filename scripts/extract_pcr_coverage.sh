#!/bin/bash
echo "=== Executing PCR Coverage Data Extraction (Multiple Thresholds) ==="

python3 - << 'EOF'
import os
import glob
import pandas as pd

input_dir = 'pcr_coverage'
output_file = 'pcr_coverage_summary.csv'
thresholds = ['10', '20', '30', '50', '100']

if not os.path.exists(input_dir):
    print(f"Error: Directory '{input_dir}' not found.")
    exit(1)

coverage_data = []
file_pattern = os.path.join(input_dir, '*.mosdepth.region.dist.txt')

for file_path in glob.glob(file_pattern):
    sample_name = os.path.basename(file_path).split('.mosdepth')[0]
    sample_data = {'PCR_Sample': sample_name}
    
    for t in thresholds:
        sample_data[f'{t}X'] = 0.0
        
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 3 and parts[0] == 'total':
                t_val = parts[1]
                if t_val in thresholds:
                    sample_data[f'{t_val}X'] = float(parts[2]) * 100

    coverage_data.append(sample_data)

if not coverage_data:
    print("Error: No valid coverage data could be extracted.")
    exit(1)

df_cov = pd.DataFrame(coverage_data)
df_cov.to_csv(output_file, index=False)
print(f"Extraction complete. Data exported to '{output_file}'.")
EOF
