#!/usr/bin/env python3
from pathlib import Path
import subprocess, itertools, csv

PCR_ROOT = Path("/home/nadirlab/Desktop/faod_pcr/FAOD_PCR_Analysis")
PROBE_ROOT = Path("/home/nadirlab/Desktop/CELEMICS_DATA/FAOD_PROBE_Analysis")
OUT = Path("/home/nadirlab/Desktop/faod_paper/faod_sonuclar_final/reanalysis_final")
OUT.mkdir(parents=True, exist_ok=True)

PCR = ["barcode93","barcode94"]
PROBE = ["barcode95","barcode96"]

def load_variants(vcf):
    cmd = ["bcftools","query","-f","%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n",str(vcf)]
    p = subprocess.run(cmd, check=True, text=True, capture_output=True)
    s=set()

    for line in p.stdout.splitlines():
        chrom,pos,ref,alt,gt = line.split("\t")[:5]

        if gt in {"0/0","./.",".|."}:
            continue

        for a in alt.split(","):
            s.add((chrom,int(pos),ref.upper(),a.upper()))

    return s

sets={}

for s in PCR:
    sets[s]=load_variants(PCR_ROOT/"norm_vcf"/f"{s}.norm.vcf.gz")

for s in PROBE:
    sets[s]=load_variants(PROBE_ROOT/"norm_vcf"/f"{s}.norm.vcf.gz")

rows=[]

for p,q in itertools.product(PCR,PROBE):
    a,b=sets[p],sets[q]
    inter=len(a & b)
    union=len(a | b)
    jac=inter/union if union else 0.0
    rows.append((p,q,len(a),len(b),inter,union,jac))

rows.sort(key=lambda x:x[-1], reverse=True)

out=OUT/"new_pair_fingerprint.tsv"

with out.open("w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["PCR","Probe","PCR_variants","Probe_variants","Shared","Union","Jaccard"])
    w.writerows(rows)

print()
print("Normalized-VCF fingerprint comparison")
print("PCR\tProbe\tPCR_n\tProbe_n\tShared\tUnion\tJaccard")

for r in rows:
    print("\t".join([*map(str,r[:-1]), f"{r[-1]:.6f}"]))

print()
print("Saved:", out)
