#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re
import math

PCR_ROOT = Path("/home/nadirlab/Desktop/faod_pcr/FAOD_PCR_Analysis")
PROBE_ROOT = Path("/home/nadirlab/Desktop/CELEMICS_DATA/FAOD_PROBE_Analysis")
PAPER_ROOT = Path("/home/nadirlab/Desktop/faod_paper/faod_sonuclar_final")
OUT = PAPER_ROOT / "reanalysis_final"
OUT.mkdir(parents=True, exist_ok=True)

MANIFEST = PAPER_ROOT / "master_sample_manifest.tsv"

STRICT_PLP = {
    "PATHOGENIC",
    "LIKELY_PATHOGENIC",
    "PATHOGENIC/LIKELY_PATHOGENIC"
}

def norm_class(x):
    if pd.isna(x):
        return ""

    s = str(x).strip().upper().replace(" ","_")

    return s

def is_strict_plp(x):
    s = norm_class(x)

    if "CONFLICT" in s:
        return False

    if "UNCERTAIN" in s:
        return False

    if s in {"VUS",""}:
        return False

    return s in STRICT_PLP

def parse_gvar(s):
    m = re.match(
        r'["\']?(chr[^:]+):(\d+)\s+([^>]+)>([^"\']+)["\']?',
        str(s).strip()
    )

    if not m:
        return ("",None,"","")

    return (
        m.group(1),
        int(m.group(2)),
        m.group(3).upper(),
        m.group(4).upper()
    )

def read_report(root,sample,method):
    f = root / "report_fix2" / f"{sample}.final.rsid.excel.tsv"

    if not f.exists():
        raise FileNotFoundError(f)

    df = pd.read_csv(
        f,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    required = [
        "Sample",
        "Gene",
        "Variant(g.)",
        "Zygosity",
        "Read_Depth",
        "Allele_Frequency",
        "ClinVar_Significance",
        "RefSeq_HGVSc",
        "RefSeq_HGVSp",
        "Technical_Status"
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise RuntimeError(
            f"{f}: missing columns {missing}"
        )

    df["Method"] = method
    df["Input_File"] = str(f)

    keys = df["Variant(g.)"].map(parse_gvar)

    df["CHROM"] = [k[0] for k in keys]
    df["POS"]   = [k[1] for k in keys]
    df["REF"]   = [k[2] for k in keys]
    df["ALT"]   = [k[3] for k in keys]

    df["Variant_Key"] = df.apply(
        lambda r: f"{r.CHROM}:{r.POS}:{r.REF}:{r.ALT}",
        axis=1
    )

    df["Strict_PLP"] = df["ClinVar_Significance"].map(is_strict_plp)

    return df

manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
    dtype=str
).fillna("")

primary = manifest[
    (manifest["role"] == "biological") &
    (manifest["pair_status"] == "confirmed")
].copy()

all_rows=[]

for _,m in primary.iterrows():

    datasets = [
        ("PCR",m["pcr_sample"],PCR_ROOT),
        ("Probe",m["probe_sample"],PROBE_ROOT),
    ]

    for method,sample,root in datasets:

        d = read_report(
            root,
            sample,
            method
        )

        d["patient_id"] = m["patient_id"]
        d["pair_status"] = m["pair_status"]

        all_rows.append(d)

allv = pd.concat(
    all_rows,
    ignore_index=True
)

plp = allv[
    allv["Strict_PLP"]
].copy()

plp = plp.sort_values(
    [
        "patient_id",
        "Method",
        "Variant_Key",
        "Read_Depth"
    ],
    ascending=[
        True,
        True,
        True,
        False
    ]
)

plp = plp.drop_duplicates(
    [
        "patient_id",
        "Method",
        "Variant_Key"
    ],
    keep="first"
)

pcr = plp[
    plp["Method"]=="PCR"
].copy()

probe = plp[
    plp["Method"]=="Probe"
].copy()

keep = [
    "patient_id",
    "Variant_Key",
    "Gene",
    "RefSeq_HGVSc",
    "RefSeq_HGVSp",
    "ClinVar_Significance",
    "Zygosity",
    "Read_Depth",
    "Allele_Frequency",
    "Technical_Status",
    "Sample"
]

pcr2 = pcr[keep].add_prefix("PCR_")
probe2 = probe[keep].add_prefix("Probe_")

pcr2 = pcr2.rename(
    columns={
        "PCR_patient_id":"patient_id",
        "PCR_Variant_Key":"Variant_Key"
    }
)

probe2 = probe2.rename(
    columns={
        "Probe_patient_id":"patient_id",
        "Probe_Variant_Key":"Variant_Key"
    }
)

cmp = pcr2.merge(
    probe2,
    on=[
        "patient_id",
        "Variant_Key"
    ],
    how="outer",
    indicator=True
)

cmp["Presence_Status"] = cmp["_merge"].map(
    {
        "both":"Matched",
        "left_only":"PCR_only",
        "right_only":"Probe_only"
    }
)

cmp["Genotype_Agreement"] = cmp.apply(
    lambda r:
    (
        r.get("PCR_Zygosity","")
        ==
        r.get("Probe_Zygosity","")
    )
    if r["Presence_Status"]=="Matched"
    else pd.NA,
    axis=1
)

def afnum(x):
    try:
        return float(str(x).split(",")[0])
    except:
        return math.nan

cmp["PCR_AF_num"] = cmp[
    "PCR_Allele_Frequency"
].map(afnum)

cmp["Probe_AF_num"] = cmp[
    "Probe_Allele_Frequency"
].map(afnum)

cmp[
    "AF_Difference_Probe_minus_PCR"
] = (
    cmp["Probe_AF_num"]
    -
    cmp["PCR_AF_num"]
)

summary = pd.DataFrame(
    {
        "Metric":[
            "Confirmed biological PCR-Probe pairs",
            "Strict P/LP PCR observations",
            "Strict P/LP Probe observations",
            "Matched strict P/LP sample-variant pairs",
            "PCR-only strict P/LP",
            "Probe-only strict P/LP",
            "Matched pairs with same zygosity",
            "Matched pairs with different zygosity",
        ],
        "Value":[
            len(primary),
            len(pcr),
            len(probe),
            (cmp["Presence_Status"]=="Matched").sum(),
            (cmp["Presence_Status"]=="PCR_only").sum(),
            (cmp["Presence_Status"]=="Probe_only").sum(),
            (
                (cmp["Presence_Status"]=="Matched")
                &
                (cmp["Genotype_Agreement"]==True)
            ).sum(),
            (
                (cmp["Presence_Status"]=="Matched")
                &
                (cmp["Genotype_Agreement"]==False)
            ).sum(),
        ]
    }
)

with pd.ExcelWriter(
    OUT/"strict_plp_pcr_probe_report.xlsx",
    engine="openpyxl"
) as xw:

    manifest.to_excel(
        xw,
        "Manifest",
        index=False
    )

    primary.to_excel(
        xw,
        "Primary_Pairs",
        index=False
    )

    allv.to_excel(
        xw,
        "All_Final_Variants",
        index=False
    )

    plp.to_excel(
        xw,
        "Strict_PLP_All",
        index=False
    )

    cmp.to_excel(
        xw,
        "PLP_Comparison",
        index=False
    )

    cmp[
        cmp["Presence_Status"]=="Matched"
    ].to_excel(
        xw,
        "Matched_PLP",
        index=False
    )

    cmp[
        cmp["Presence_Status"]=="PCR_only"
    ].to_excel(
        xw,
        "PCR_Only_PLP",
        index=False
    )

    cmp[
        cmp["Presence_Status"]=="Probe_only"
    ].to_excel(
        xw,
        "Probe_Only_PLP",
        index=False
    )

    summary.to_excel(
        xw,
        "Summary",
        index=False
    )

summary.to_csv(
    OUT/"strict_plp_summary.tsv",
    sep="\t",
    index=False
)

cmp.to_csv(
    OUT/"strict_plp_comparison.tsv",
    sep="\t",
    index=False
)

print()
print("STRICT P/LP PRIMARY COMPARISON")
print(summary.to_string(index=False))

print()
print(
    "Workbook:",
    OUT/"strict_plp_pcr_probe_report.xlsx"
)

print(
    "Comparison TSV:",
    OUT/"strict_plp_comparison.tsv"
)
