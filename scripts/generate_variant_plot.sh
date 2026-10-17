#!/bin/bash
# ==============================================================================
# Script Name: generate_variant_plot.sh
# Description: Extracts QUAL and AF metrics from phased VCF files and generates
#              a comparative scatter plot for Probe-Capture vs LR-PCR ONT data.
# Requirements: bcftools, python3, pandas, matplotlib, seaborn
# ==============================================================================

echo "[INFO] Starting Variant Quality vs Allele Frequency analysis..."

# Step 1: Extract data from VCF files using bcftools
echo "[INFO] Extracting metrics from VCF files..."

# Empty existing files if they exist
> probe_qual_af.txt
> pcr_qual_af.txt

# Extract from Probe-Capture VCFs
for vcf in phased_results/ONT_UDI_*_phased.vcf.gz; do
    if [ -f "$vcf" ]; then
        bcftools query -f '%QUAL\t[%AF]\n' "$vcf" >> probe_qual_af.txt
    fi
done

# Extract from LR-PCR VCFs
for vcf in phased_results/barcode*_phased.vcf.gz; do
    if [ -f "$vcf" ]; then
        bcftools query -f '%QUAL\t[%AF]\n' "$vcf" >> pcr_qual_af.txt
    fi
done

echo "[INFO] Data extraction complete. Generating visualization..."

# Step 2: Create a temporary Python script to plot the data
cat << 'EOF' > temp_plot_qual_af.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def load_data(file_name, label):
    if not os.path.exists(file_name) or os.stat(file_name).st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(file_name, sep='\t', header=None, names=['QUAL', 'AF'])
        df['AF'] = pd.to_numeric(df['AF'], errors='coerce')
        df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce')
        df = df.dropna()
        df['Method'] = label
        return df
    except Exception as e:
        print(f"Error loading {file_name}: {e}")
        return pd.DataFrame()

# Load datasets
df_probe = load_data('probe_qual_af.txt', 'Probe-Capture (ONT)')
df_pcr = load_data('pcr_qual_af.txt', 'LR-PCR (ONT)')

# Combine data
df_all = pd.concat([df_probe, df_pcr])

if df_all.empty:
    print("Error: No data found to plot.")
    exit()

# Set up the plot
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df_all, x='AF', y='QUAL', hue='Method', alpha=0.6, 
                edgecolor=None, s=30, palette=['#1f77b4', '#d62728'])

# Add Clinical Threshold Lines
plt.axhline(y=30, color='gray', linestyle='--', linewidth=1.5, label='QUAL > 30 (Confidence Threshold)')
plt.axvline(x=0.2, color='green', linestyle=':', linewidth=1.5, label='AF = 0.2 (Low-Freq Boundary)')

# Formatting
plt.yscale('log')
plt.title('Comparative Variant Quality (QUAL) vs. Allele Frequency (AF)', fontsize=14, fontweight='bold')
plt.xlabel('Allele Frequency (AF)', fontsize=12)
plt.ylabel('Variant Quality Score (QUAL) - Log Scale', fontsize=12)
plt.xlim(0, 1.05)
plt.legend(loc='upper right')
plt.grid(True, which="both", ls="--", alpha=0.3)

# Save Plot
output_file = 'Figure_1_Variant_Quality_AF_Plot.png'
plt.tight_layout()
plt.savefig(output_file, dpi=300)
print(f"[SUCCESS] Plot saved as '{output_file}'")
EOF

# Step 3: Execute the Python script
python3 temp_plot_qual_af.py

# Step 4: Cleanup temporary files (Optional: comment out to keep raw data)
rm pcr_qual_af.txt probe_qual_af.txt temp_plot_qual_af.py

echo "[INFO] Analysis pipeline finished successfully."
