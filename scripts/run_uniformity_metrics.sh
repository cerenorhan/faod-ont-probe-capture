#!/bin/bash
echo "=== Executing Formal Depth Uniformity Analysis ==="

python3 - << 'EOF'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

report_file = 'FAOD_ONT_report.xlsx'

if not os.path.exists(report_file):
    print(f"Error: Required file '{report_file}' not found.")
    exit(1)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

df_sum = pd.read_excel(report_file, sheet_name='Summary', header=1)

if 'Demux Method' in df_sum.columns:
    df_sum = df_sum[df_sum['Demux Method'] == 'barbell']

required_cols = ['Fold 80', '1X (%)', '10X (%)', '20X (%)', '30X (%)', '50X (%)', '100X (%)']
missing_cols = [c for c in required_cols if c not in df_sum.columns]

if missing_cols:
    print(f"Error: Missing columns in the report: {', '.join(missing_cols)}")
    exit(1)

fold_80 = pd.to_numeric(df_sum['Fold 80'], errors='coerce').dropna()

thresholds = ['1X', '10X', '20X', '30X', '50X', '100X']
t_cols = [f'{t} (%)' for t in thresholds]
cov_data = df_sum[t_cols].apply(pd.to_numeric, errors='coerce').dropna()

cov_means = cov_data.mean()
cov_stds = cov_data.std()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

sns.boxplot(y=fold_80, ax=ax1, color='#BDC3C7', width=0.4, showfliers=False, boxprops=dict(alpha=0.7))
sns.stripplot(y=fold_80, ax=ax1, color='#2C3E50', alpha=0.7, size=6, jitter=True)

ax1.axhline(1.5, color='#E74C3C', linestyle='--', linewidth=1.5, label='Clinical Ideal (<1.5)')
ax1.set_title('A. Fold-80 Base Penalty Distribution', fontweight='bold', fontsize=13, pad=15)
ax1.set_ylabel('Fold-80 Penalty Score', fontweight='bold', fontsize=12)
ax1.legend(loc='upper right')
ax1.grid(axis='y', linestyle=':', alpha=0.6)
for spine in ['top', 'right', 'bottom']:
    ax1.spines[spine].set_visible(False)
ax1.set_xticks([])

x_pos = np.arange(len(thresholds))
ax2.plot(x_pos, cov_means, marker='o', color='#2980B9', linewidth=2.5, markersize=8)
ax2.fill_between(x_pos, cov_means - cov_stds, cov_means + cov_stds, color='#3498DB', alpha=0.2)

ax2.set_xticks(x_pos)
ax2.set_xticklabels(thresholds, fontweight='bold')
ax2.set_title('B. Target Coverage Uniformity Profile', fontweight='bold', fontsize=13, pad=15)
ax2.set_ylabel('Target Fraction Covered (%)', fontweight='bold', fontsize=12)
ax2.set_xlabel('Minimum Depth Threshold', fontweight='bold', fontsize=12)
ax2.set_ylim(0, 105)
ax2.grid(True, linestyle=':', alpha=0.6)
for spine in ['top', 'right']:
    ax2.spines[spine].set_visible(False)

plt.tight_layout()
output_img = 'Depth_Uniformity_Metrics.png'
plt.savefig(output_img, dpi=300, bbox_inches='tight')

output_text = f"""--------------------------------------------------
Formal Depth Uniformity Metrics (Filter: barbell)
--------------------------------------------------
Total Samples Analyzed : {len(fold_80)}

Fold-80 Base Penalty:
  Mean ± SD            : {fold_80.mean():.3f} ± {fold_80.std():.3f}
  Median               : {fold_80.median():.3f}
  Range (Min - Max)    : {fold_80.min():.3f} - {fold_80.max():.3f}

Coverage Retention across Thresholds (Mean %):
"""
for t in thresholds:
    t_col = f'{t} (%)'
    output_text += f"  \u2265 {t:<4}               : {cov_means[t_col]:.2f}%\n"
output_text += "--------------------------------------------------\n"

print(output_text)
with open("Uniformity_Metrics_Results.txt", "w") as f:
    f.write(output_text)

print(f"Analysis complete. Outputs saved to '{output_img}' and 'Uniformity_Metrics_Results.txt'.")
EOF

