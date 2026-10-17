#!/bin/bash
echo "=== Executing Per-Region Depth Rescue Analysis ==="

python3 - << 'EOF'
import os
import pandas as pd
import matplotlib.pyplot as plt

pcr_failed_file = 'pcr_failed_regions.csv'
map_file = 'faod_variant_metrics_report_pcr_vs_probe.xlsx'
probe_dir = 'probe_coverage'
clinical_threshold = 30.0

if not os.path.exists(pcr_failed_file) or not os.path.exists(probe_dir):
    print("Error: Required inputs ('pcr_failed_regions.csv' or 'probe_coverage' directory) not found.")
    exit(1)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

df_failed = pd.read_csv(pcr_failed_file)
df_map = pd.read_excel(map_file, sheet_name='Sample_Mapping')
sample_map = dict(zip(df_map['PCR_Sample'], df_map['Probe_Sample']))

results = []

for _, row in df_failed.iterrows():
    pcr_sample = row['PCR_Sample']
    region = row['Gene_Region']
    pcr_depth = row['PCR_Depth']
    
    probe_sample = sample_map.get(pcr_sample)
    if not probe_sample:
        continue
        
    probe_bed = os.path.join(probe_dir, f"{probe_sample}.regions.bed")
    
    if os.path.exists(probe_bed):
        try:
            df_probe = pd.read_csv(probe_bed, sep='\t', header=None, usecols=[3, 4], names=['region', 'depth'])
            probe_match = df_probe[df_probe['region'] == region]
            
            if not probe_match.empty:
                probe_depth = float(probe_match.iloc[0]['depth'])
                results.append({
                    'Label': f"{region} ({pcr_sample})",
                    'PCR_Depth': pcr_depth,
                    'Probe_Depth': probe_depth
                })
        except Exception as e:
            print(f"Warning: Failed to parse {probe_bed}. Reason: {e}")

if not results:
    print("Error: No matching Probe regions found for the failed PCR regions.")
    exit(1)

df_plot = pd.DataFrame(results)
df_plot = df_plot.sort_values(by='PCR_Depth', ascending=True).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(9, max(6, len(df_plot) * 0.4)))

for i, row in df_plot.iterrows():
    ax.plot([row['PCR_Depth'], row['Probe_Depth']], [i, i], color='#95A5A6', linewidth=2.5, zorder=1)
    
    ax.scatter(row['PCR_Depth'], i, color='#E74C3C', s=100, zorder=2, edgecolor='black', 
               label='PCR Amplicon (<30X)' if i == 0 else "")
    ax.scatter(row['Probe_Depth'], i, color='#27AE60', s=100, zorder=2, edgecolor='black', 
               label='Probe Capture (Rescued)' if i == 0 else "")

ax.axvline(clinical_threshold, color='#2C3E50', linestyle='--', linewidth=1.5, zorder=0, label=f'Clinical Threshold ({clinical_threshold}X)')

ax.set_yticks(range(len(df_plot)))
ax.set_yticklabels(df_plot['Label'], fontweight='bold')
ax.set_xlabel('Sequencing Depth (X)', fontweight='bold', fontsize=12)
ax.set_title('Per-Region Depth Rescue Analysis', fontweight='bold', fontsize=13, pad=15)
ax.grid(axis='x', linestyle=':', alpha=0.6)

for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

ax.legend(loc='lower right', frameon=True)

plt.tight_layout()
output_img = 'Per_Region_Depth_Rescue.png'
plt.savefig(output_img, dpi=300, bbox_inches='tight')

output_csv = 'Per_Region_Rescue_Data.csv'
df_plot.to_csv(output_csv, index=False)

print(f"--------------------------------------------------")
print(f"Per-Region Rescue Analysis Completed")
print(f"--------------------------------------------------")
print(f"Total Regions Rescued : {len(df_plot)}")
print(f"Outputs generated     : '{output_img}', '{output_csv}'")
print(f"--------------------------------------------------")
EOF
