#!/bin/bash
echo "=== Executing Stratified Bland-Altman Analysis ==="

python3 - << 'END_XYZ'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

excel_file = 'faod_variant_metrics_report_pcr_vs_probe.xlsx'
if not os.path.exists(excel_file):
    print(f"Error: Target file '{excel_file}' not found. Execution aborted.")
    exit(1)

# Configure plot aesthetics for publication
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

# Extract clinically significant sample rows
df = pd.read_excel(excel_file, sheet_name='P_LP_Conf_Sample')

probe_df = df[df['Platform'] == 'Probe'][['Probe_Sample', 'Variant(g.)', 'Allele_Frequency']].rename(columns={'Allele_Frequency': 'AF_Probe'})
pcr_df = df[df['Platform'] == 'PCR'][['Probe_Sample', 'Variant(g.)', 'Allele_Frequency']].rename(columns={'Allele_Frequency': 'AF_PCR'})

merged_df = pd.merge(probe_df, pcr_df, on=['Probe_Sample', 'Variant(g.)']).dropna()

merged_df['Mean_AF'] = (merged_df['AF_Probe'] + merged_df['AF_PCR']) / 2
merged_df['Diff_AF'] = merged_df['AF_Probe'] - merged_df['AF_PCR']

# Stratify by intermediate vs. high allele frequency
merged_df['AF_Class'] = np.where(merged_df['Mean_AF'] < 0.80, 'Intermediate AF (<0.80)', 'High AF (\u22650.80)')

fig, ax = plt.subplots(figsize=(9, 6))
colors = {'Intermediate AF (<0.80)': '#E67E22', 'High AF (\u22650.80)': '#2980B9'}

for af_class, group in merged_df.groupby('AF_Class'):
    ax.scatter(group['Mean_AF'], group['Diff_AF'], label=af_class, s=90, 
               color=colors.get(af_class, '#7F8C8D'), alpha=0.85, edgecolor='black', linewidth=0.8)

mean_diff = merged_df['Diff_AF'].mean()
std_diff = merged_df['Diff_AF'].std()
loa_upper = mean_diff + 1.96 * std_diff
loa_lower = mean_diff - 1.96 * std_diff

ax.axhline(mean_diff, color='#2C3E50', linestyle='-', linewidth=2, label=f'Mean Bias ({mean_diff:+.3f})')
ax.axhline(loa_upper, color='#C0392B', linestyle='--', linewidth=1.5, label=f'+1.96 SD ({loa_upper:+.3f})')
ax.axhline(loa_lower, color='#C0392B', linestyle='--', linewidth=1.5, label=f'\u22121.96 SD ({loa_lower:+.3f})')

ax.set_xlabel('Mean Allele Frequency', fontsize=12, fontweight='bold')
ax.set_ylabel('Difference (Probe AF \u2212 PCR AF)', fontsize=12, fontweight='bold')
ax.set_title('Bland-Altman Agreement Stratified by AF Class', fontsize=13, fontweight='bold', pad=12)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=True)
ax.grid(True, linestyle=':', alpha=0.6)
ax.set_xlim(-0.05, 1.05)
plt.tight_layout()

output_img = 'Stratified_BlandAltman.png'
plt.savefig(output_img, dpi=300, bbox_inches='tight')

# Generate statistical summary report
output_text = f"""--------------------------------------------------
Bland-Altman Analysis Results
--------------------------------------------------
Overall Mean Bias (Difference) : {mean_diff:+.4f}
Standard Deviation of Diff     : {std_diff:.4f}
95% Limits of Agreement (LOA)  : [{loa_lower:+.4f}, {loa_upper:+.4f}]

Stratification by AF Class:
"""

for af_class, group in merged_df.groupby('AF_Class'):
    c_mean = group['Diff_AF'].mean()
    c_std = group['Diff_AF'].std()
    output_text += f" - {af_class} (n={len(group)}): Mean Bias = {c_mean:+.4f}, SD = {c_std:.4f}\n"
output_text += "--------------------------------------------------\n"

print(output_text)

with open("BlandAltman_Results.txt", "w") as f:
    f.write(output_text)

print(f"Analysis complete. Outputs successfully written to '{output_img}' and 'BlandAltman_Results.txt'.")
END_XYZ
