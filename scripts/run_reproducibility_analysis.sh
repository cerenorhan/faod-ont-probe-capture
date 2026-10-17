#!/bin/bash
echo "=== Executing Reproducibility and Concordance Analysis ==="

python3 - << 'EOF'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import os
import warnings

# Suppress pandas HTML parsing warnings
warnings.filterwarnings('ignore', category=FutureWarning)

html_file = 'vcf_comparison_report.html'

if not os.path.exists(html_file):
    print(f"Error: Target file '{html_file}' not found.")
    exit(1)

# Parse all tables sequentially from the HTML DOM
try:
    tables = pd.read_html(html_file)
except Exception as e:
    print(f"Error parsing HTML tables: {e}")
    exit(1)

if len(tables) < 9:
    print("Error: Expected 9 tables in the HTML structure, found fewer. Check file integrity.")
    exit(1)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

# ---------------------------------------------------------
# 1. Pairwise Concordance Summary Table
# ---------------------------------------------------------
df_neg_ctrl_1_pair = tables[1].copy()
df_neg_ctrl_2_pair = tables[4].copy()
df_b68_pair = tables[7].copy()

df_b68_pair.insert(0, 'Cohort', 'Clinical Replicates (B68-B72)')
df_neg_ctrl_1_pair.insert(0, 'Cohort', 'Negative Control 1')
df_neg_ctrl_2_pair.insert(0, 'Cohort', 'Negative Control 2')

df_combined_pair = pd.concat([df_b68_pair, df_neg_ctrl_1_pair, df_neg_ctrl_2_pair], ignore_index=True)
csv_output = 'Pairwise_Concordance_Summary.csv'
df_combined_pair.to_csv(csv_output, index=False)

# ---------------------------------------------------------
# 2. AF Correlation Scatter Plot (B68 vs B72)
# ---------------------------------------------------------
df_b68_mat = tables[8].copy()

def extract_af(val):
    if pd.isna(val) or val == '-': return np.nan
    val_str = str(val)
    if 'AF=' in val_str:
        return float(val_str.split('AF=')[1])
    return np.nan

x_af = df_b68_mat['ONT_UDI_68'].apply(extract_af)
y_af = df_b68_mat['ONT_UDI_72'].apply(extract_af)

mask = ~x_af.isna() & ~y_af.isna()
x_valid = x_af[mask]
y_valid = y_af[mask]

fig1, ax1 = plt.subplots(figsize=(6, 6))

ax1.scatter(x_valid, y_valid, color='#2C3E50', s=70, alpha=0.85, edgecolor='black', linewidth=0.8)
ax1.plot([0, 1.05], [0, 1.05], color='#E74C3C', linestyle='--', linewidth=1.5, label='Identity Line (y=x)')

r_val, p_val = pearsonr(x_valid, y_valid)
p_str = "P < 0.001" if p_val < 0.001 else f"P = {p_val:.3f}"
stats_text = f"Pearson r = {r_val:.4f}\n{p_str}\nn = {len(x_valid)} variants"

props = dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='#BDC3C7', alpha=0.9)
ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes, fontsize=11, verticalalignment='top', bbox=props)

ax1.set_xlim(-0.02, 1.02)
ax1.set_ylim(-0.02, 1.02)
ax1.set_xlabel('ONT_UDI_68 Allele Frequency', fontweight='bold', fontsize=12)
ax1.set_ylabel('ONT_UDI_72 Allele Frequency', fontweight='bold', fontsize=12)
ax1.set_title('Intra-run Technical Replicate Concordance', fontweight='bold', fontsize=13, pad=15)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='lower right', frameon=True)

for spine in ['top', 'right']:
    ax1.spines[spine].set_visible(False)

plt.tight_layout()
scatter_output = 'AF_Correlation_Replicates.png'
fig1.savefig(scatter_output, dpi=300, bbox_inches='tight')

# ---------------------------------------------------------
# 3. Jaccard Similarity Heatmap Panel
# ---------------------------------------------------------
def clean_jaccard_matrix(df):
    df_clean = df.set_index(df.columns[0])
    for col in df_clean.columns:
        df_clean[col] = df_clean[col].astype(str).str.replace('%', '').astype(float)
    return df_clean

jac_neg_ctrl_1 = clean_jaccard_matrix(tables[0])
jac_neg_ctrl_2 = clean_jaccard_matrix(tables[3])
jac_b68 = clean_jaccard_matrix(tables[6])

fig2, axes = plt.subplots(1, 3, figsize=(16, 5), gridspec_kw={'width_ratios': [1, 2, 2]})

sns.heatmap(jac_b68, annot=True, fmt=".1f", cmap="Blues", ax=axes[0], cbar=False, vmin=90, vmax=100, linewidths=1, linecolor='white')
axes[0].set_title('A. Clinical Replicates (B68-B72)', fontweight='bold', pad=10)
axes[0].set_ylabel('')

sns.heatmap(jac_neg_ctrl_1, annot=True, fmt=".1f", cmap="Blues", ax=axes[1], cbar=False, vmin=90, vmax=100, linewidths=1, linecolor='white')
axes[1].set_title('B. Negative Control 1', fontweight='bold', pad=10)
axes[1].set_ylabel('')

sns.heatmap(jac_neg_ctrl_2, annot=True, fmt=".1f", cmap="Blues", ax=axes[2], cbar=True, vmin=90, vmax=100, linewidths=1, linecolor='white', cbar_kws={'label': 'Jaccard Similarity (%)'})
axes[2].set_title('C. Negative Control 2', fontweight='bold', pad=10)
axes[2].set_ylabel('')

for ax in axes:
    ax.tick_params(axis='x', rotation=45)
    ax.tick_params(axis='y', rotation=0)

plt.tight_layout()
heatmap_output = 'Jaccard_Similarity_Heatmaps.png'
fig2.savefig(heatmap_output, dpi=300, bbox_inches='tight')

print("--------------------------------------------------")
print("Reproducibility Extraction Completed Successfully")
print("--------------------------------------------------")
print(f"1. Tabular Summary exported to : {csv_output}")
print(f"2. AF Correlation Plot saved   : {scatter_output}")
print(f"3. Similarity Matrices saved   : {heatmap_output}")
print("--------------------------------------------------")
EOF
