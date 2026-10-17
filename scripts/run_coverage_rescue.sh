#!/bin/bash
echo "=== Executing Coverage Rescue Analysis ==="

python3 - << 'END_XYZ'
import pandas as pd
import matplotlib.pyplot as plt
import os

metrics_file = 'faod_variant_metrics_report_pcr_vs_probe.xlsx'
report_file = 'FAOD_ONT_report.xlsx'

if not os.path.exists(metrics_file) or not os.path.exists(report_file):
    print(f"Error: Required input files not found. Ensure '{metrics_file}' and '{report_file}' are present.")
    exit(1)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

df_map = pd.read_excel(metrics_file, sheet_name='Sample_Mapping')
df_sum = pd.read_excel(report_file, sheet_name='Summary', header=1)

sample_col = 'Sample Name' if 'Sample Name' in df_sum.columns else 'Sample'
cov_col = '30X (%)'

df_sum_filtered = df_sum[[sample_col, cov_col]].copy()
df_sum_filtered[cov_col] = pd.to_numeric(df_sum_filtered[cov_col], errors='coerce')

merged = df_map.merge(df_sum_filtered, left_on='Probe_Sample', right_on=sample_col, how='left').rename(columns={cov_col: 'Probe_30X'})
merged = merged.merge(df_sum_filtered, left_on='PCR_Sample', right_on=sample_col, how='left').rename(columns={cov_col: 'PCR_30X'})
merged = merged.dropna(subset=['Probe_30X', 'PCR_30X'])

fig, ax = plt.subplots(figsize=(6.5, 7))

rescued_count = 0
total = len(merged)

for idx, row in merged.iterrows():
    pcr_val = row['PCR_30X']
    probe_val = row['Probe_30X']
    
    is_rescued = probe_val > pcr_val
    color = '#16A085' if is_rescued else '#7F8C8D'
    alpha = 0.85 if is_rescued else 0.4
    linewidth = 2.0 if is_rescued else 1.2
    
    if is_rescued:
        rescued_count += 1
        
    ax.plot([1, 2], [pcr_val, probe_val], color=color, alpha=alpha, linewidth=linewidth, marker='o', markersize=6)

ax.set_xticks([1, 2])
ax.set_xticklabels(['PCR Amplicon\nWorkflow', 'Probe Capture\nWorkflow (Rescued)'], fontweight='bold', fontsize=12)
ax.set_xlim(0.8, 2.2)
ax.set_ylabel('Target Genomic Regions \u226530X Coverage (%)', fontweight='bold', fontsize=12)
ax.set_title('Sequencing Depth Rescue at Clinical \u226530X Threshold', fontsize=13, fontweight='bold', pad=15)
ax.grid(axis='y', linestyle=':', alpha=0.6)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

stats_text = f"Total Paired Samples: {total}\nSamples with Rescued Depth: {rescued_count} ({rescued_count/total*100:.1f}%)"
props = dict(boxstyle='round,pad=0.6', facecolor='#F8F9FA', edgecolor='#BDC3C7', alpha=0.9)
ax.text(1.5, ax.get_ylim()[0] + (ax.get_ylim()[1]-ax.get_ylim()[0])*0.02, stats_text, ha='center', fontsize=11, bbox=props)

plt.tight_layout()

output_img = 'Coverage_Rescue_Plot.png'
plt.savefig(output_img, dpi=300, bbox_inches='tight')

output_text = f"""--------------------------------------------------
Coverage Rescue Analysis Results (>= 30X Threshold)
--------------------------------------------------
Total Paired Clinical Samples         : {total}
Samples with Improved Target Coverage : {rescued_count}
Percentage of Rescued Samples         : {rescued_count/total*100:.1f}%
--------------------------------------------------
"""

print(output_text)
with open("Coverage_Rescue_Results.txt", "w") as f:
    f.write(output_text)

print(f"Analysis complete. Outputs saved to '{output_img}' and 'Coverage_Rescue_Results.txt'.")
END_XYZ
