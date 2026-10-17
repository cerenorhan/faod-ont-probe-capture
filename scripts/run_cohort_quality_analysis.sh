#!/bin/bash
echo "=== Executing Cohort Comparative Variant Quality vs. Allele Frequency Analysis ==="

python3 - << 'EOF'
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import glob

# Yeni oluşturduğunuz klasör yolları
probe_dir = 'probe_vcf'
pcr_dir = 'pcr_vcf'

def extract_metrics_from_folder(folder_path):
    qualities = []
    afs = []
    
    if not os.path.exists(folder_path):
        print(f"Warning: Folder '{folder_path}' not found. Skipping...")
        return afs, qualities

    # Klasördeki tüm .vcf dosyalarını bul
    vcf_files = glob.glob(os.path.join(folder_path, '*.vcf'))
    
    if not vcf_files:
        print(f"Warning: No .vcf files found in '{folder_path}'.")
        return afs, qualities

    print(f"Processing {len(vcf_files)} files in '{folder_path}'...")
    
    for vcf_file in vcf_files:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                    
                qual_str = parts[5]
                if qual_str == '.' or qual_str == '0':
                    continue
                    
                try:
                    qual = float(qual_str)
                except ValueError:
                    continue

                info = parts[7]
                fmt = parts[8].split(':')
                sample = parts[9].split(':')
                
                af_val = None
                
                if 'AF' in fmt:
                    af_idx = fmt.index('AF')
                    if len(sample) > af_idx and sample[af_idx] != '.':
                        try:
                            af_val = float(sample[af_idx].split(',')[0])
                        except ValueError:
                            pass
                            
                if af_val is None:
                    af_match = re.search(r'(?:^|;)AF=([\d\.]+)', info)
                    if af_match:
                        try:
                            af_val = float(af_match.group(1))
                        except ValueError:
                            pass
                
                if af_val is not None:
                    qualities.append(qual)
                    afs.append(af_val)
                    
    return afs, qualities

# Her iki klasörden de verileri topla
probe_afs, probe_quals = extract_metrics_from_folder(probe_dir)
pcr_afs, pcr_quals = extract_metrics_from_folder(pcr_dir)

if not probe_afs and not pcr_afs:
    print("Error: No valid metrics could be extracted from either folder.")
    exit(1)

# Görselleştirme ayarları
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

fig, ax = plt.subplots(figsize=(9, 6))

# Önce PCR'ı çiz (Alta kalsın)
if pcr_afs:
    ax.scatter(pcr_afs, pcr_quals, color='#E74C3C', alpha=0.5, s=40, edgecolor='none', label=f'PCR Amplicon (n={len(pcr_afs)} variants)')

# Sonra Probe'u çiz (Üste çıksın)
if probe_afs:
    ax.scatter(probe_afs, probe_quals, color='#2980B9', alpha=0.5, s=40, edgecolor='none', label=f'Probe Capture ONT (n={len(probe_afs)} variants)')

# Eşik çizgileri
ax.axvline(0.2, color='#7F8C8D', linestyle='--', linewidth=1.5, label='Low AF Threshold (0.2)')
ax.axhline(30, color='#27AE60', linestyle='--', linewidth=1.5, label='Quality Threshold (QUAL 30)')

ax.set_xlabel('Allele Frequency (AF)', fontweight='bold')

# Yüksek kaliteli varyantlar için log ölçeği kontrolü
max_qual = max(max(probe_quals) if probe_quals else 0, max(pcr_quals) if pcr_quals else 0)
if max_qual > 200:
    ax.set_yscale('log')
    ax.set_ylabel('Variant Quality Score (QUAL) - Log Scale', fontweight='bold')
else:
    ax.set_ylabel('Variant Quality Score (QUAL)', fontweight='bold')

ax.set_title('Cohort-Level Variant Reliability: PCR vs. Probe Capture', fontweight='bold', pad=15)
ax.grid(True, which="both", linestyle=':', alpha=0.5)

# Lejantı sağ alta yerleştir
ax.legend(frameon=True, loc='lower right')

for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

plt.tight_layout()
output_png = 'Cohort_Comparative_Variant_Reliability_Plot.png'
plt.savefig(output_png, dpi=300, bbox_inches='tight')

print("--------------------------------------------------")
print(f"Analysis complete. Generated: {output_png}")
print(f"Total variants extracted - Probe: {len(probe_afs)}, PCR: {len(pcr_afs)}")
print("--------------------------------------------------")
EOF
