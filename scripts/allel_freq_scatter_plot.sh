#!/bin/bash
echo "=== Executing Allele Frequency Scatter Plot Analysis ==="

python3 - << 'END_XYZ'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os

excel_file = 'faod_variant_metrics_report_pcr_vs_probe.xlsx'
if not os.path.exists(excel_file):
    print(f"Error: Target file '{excel_file}' not found. Execution aborted.")
    exit(1)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

df = pd.read_excel(excel_file, sheet_name='P_LP_Conf_Sample')

probe_df = df[df['Platform'] == 'Probe'][['Probe_Sample', 'Variant(g.)', 'Allele_Frequency']].rename(columns={'Allele_Frequency': 'AF_Probe'})
pcr_df = df[df['Platform'] == 'PCR'][['Probe_Sample', 'Variant(g.)', 'Allele_Frequency']].rename(columns={'Allele_Frequency': 'AF_PCR'})

merged_df = pd.merge(probe_df, pcr_df, on=['Probe_Sample', 'Variant(g.)']).dropna()

x = merged_df['AF_Probe'].values
y = merged_df['AF_PCR'].values

r_val, p_val = pearsonr(x, y)

fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(x, y, s=90, color='#2C3E50', alpha=0.85, edgecolor='black', linewidth=0.8)
ax.plot([0, 1], [0, 1], color='#C0392B', linestyle='--', linewidth=1.5, label='Identity Line (y=x)')

ax.set_xlabel('Probe-based Allele Frequency', fontsize=12, fontweight='bold')
ax.set_ylabel('PCR-based Allele Frequency', fontsize=12, fontweight='bold')
ax.set_title('Allele Frequency Correlation', fontsize=13, fontweight='bold', pad=12)

# Format p-value for scientific reporting
p_val_str = "P < 0.001" if p_val < 0.001 else f"P = {p_val:.3f}"
textstr = f"Pearson r = {r_val:.4f}\n{p_val_str}"

props = dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='#BDC3C7', alpha=0.9)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11, verticalalignment='top', bbox=props)

ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower right', frameon=True)

plt.tight_layout()

output_img = 'AF_ScatterPlot.png'
plt.savefig(output_img, dpi=300, bbox_inches='tight')

print(f"Analysis complete. Scatter plot saved to '{output_img}'.")
END_XYZ
