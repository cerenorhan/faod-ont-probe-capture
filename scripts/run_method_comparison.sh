#!/bin/bash
echo "=== Executing Method Comparison: Coverage Retention Profiles ==="

python3 - << 'EOF'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

pcr_file = 'pcr_coverage_summary.csv'
probe_file = 'FAOD_ONT_report.xlsx'

if not os.path.exists(pcr_file) or not os.path.exists(probe_file):
    print("Error: Required input files not found. Ensure 'pcr_coverage_summary.csv' and the ONT Excel report are present.")
    exit(1)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

df_pcr = pd.read_csv(pcr_file)
pcr_cols = ['10X', '20X', '30X', '50X', '100X']

missing_pcr = [c for c in pcr_cols if c not in df_pcr.columns]
if missing_pcr:
    print(f"Error: Missing columns in PCR data: {', '.join(missing_pcr)}")
    exit(1)

pcr_data = df_pcr[pcr_cols].apply(pd.to_numeric, errors='coerce').dropna()
pcr_means = pcr_data.mean().values
pcr_stds = pcr_data.std().values

df_probe = pd.read_excel(probe_file, sheet_name='Summary', header=1)

if 'Demux Method' in df_probe.columns:
    df_probe = df_probe[df_probe['Demux Method'] == 'barbell']

probe_cols = ['10X (%)', '20X (%)', '30X (%)', '50X (%)', '100X (%)']
missing_probe = [c for c in probe_cols if c not in df_probe.columns]
if missing_probe:
    print(f"Error: Missing columns in Probe data: {', '.join(missing_probe)}")
    exit(1)

probe_data = df_probe[probe_cols].apply(pd.to_numeric, errors='coerce').dropna()
probe_means = probe_data.mean().values
probe_stds = probe_data.std().values

thresholds = ['10X', '20X', '30X', '50X', '100X']
x_pos = np.arange(len(thresholds))

fig, ax = plt.subplots(figsize=(8.5, 6))

# Lejantlardan (n=...) kısımları kaldırıldı
ax.plot(x_pos, pcr_means, marker='o', color='#E74C3C', linewidth=2.5, markersize=8, label='PCR Amplicon')
ax.fill_between(x_pos, pcr_means - pcr_stds, pcr_means + pcr_stds, color='#E74C3C', alpha=0.15)

ax.plot(x_pos, probe_means, marker='s', color='#2980B9', linewidth=2.5, markersize=8, label='Probe Capture ONT')
ax.fill_between(x_pos, probe_means - probe_stds, probe_means + probe_stds, color='#2980B9', alpha=0.15)

ax.set_xticks(x_pos)
ax.set_xticklabels(thresholds, fontweight='bold')
ax.set_ylabel('Target Fraction Covered (%)', fontweight='bold', fontsize=12)
ax.set_xlabel('Minimum Depth Threshold', fontweight='bold', fontsize=12)
ax.set_title('Coverage Retention Comparison Across Clinical Thresholds', fontweight='bold', fontsize=13, pad=15)
ax.set_ylim(0, 105)
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower left', frameon=True)

for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

plt.tight_layout()
output_img = 'Method_Comparison_Coverage.png'
plt.savefig(output_img, dpi=300, bbox_inches='tight')

output_text = f"""--------------------------------------------------
Method Comparison: Target Coverage Retention
--------------------------------------------------
PCR Amplicon (n={len(pcr_data)}) vs Probe Capture ONT (n={len(probe_data)}, Filter: barbell)

Mean Target Coverage (% ± SD):
Threshold | PCR Amplicon          | Probe Capture ONT
--------------------------------------------------
"""
for i, t in enumerate(thresholds):
    output_text += f"  >={t:<3}  | {pcr_means[i]:>6.2f}% ± {pcr_stds[i]:<5.2f} | {probe_means[i]:>6.2f}% ± {probe_stds[i]:<5.2f}\n"

output_text += "--------------------------------------------------\n"

print(output_text)
with open("Method_Comparison_Results.txt", "w") as f:
    f.write(output_text)

print(f"Analysis complete. Outputs saved to '{output_img}' and 'Method_Comparison_Results.txt'.")
EOF
