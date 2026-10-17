#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import subprocess
import re

PAPER = Path("/home/nadirlab/Desktop/faod_paper/faod_sonuclar_final")

PCR_ROOT = Path("/home/nadirlab/Desktop/faod_pcr/FAOD_PCR_Analysis")
PROBE_ROOT = Path("/home/nadirlab/Desktop/CELEMICS_DATA/FAOD_PROBE_Analysis")

INFILE = PAPER / "reanalysis_final" / "one_sided_plp_bam_audit.tsv"
OUTFILE = PAPER / "reanalysis_final" / "one_sided_plp_high_quality_audit.tsv"
OUTXLSX = PAPER / "reanalysis_final" / "one_sided_plp_high_quality_audit.xlsx"

# Two quality levels
QUALITY_SETTINGS = [
    ("Q20", 20, 20),
    ("Q30", 30, 20),
]

def parse_key(key):
    chrom, pos, ref, alt = str(key).split(":", 3)
    return chrom, int(pos), ref.upper(), alt.upper()

def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"

def bam_path(method, sample):
    root = PCR_ROOT if method == "PCR" else PROBE_ROOT
    return root / "bam" / f"{sample}.sorted.bam"

def parse_bases(bases, ref, alt):
    vtype = variant_type(ref, alt)

    snv_alt = 0
    ins = []
    dels = []

    i = 0

    while i < len(bases):
        c = bases[i]

        if c == "^":
            i += 2
            continue

        if c == "$":
            i += 1
            continue

        if c in ".,ACGTNacgtn":
            base = c.upper()

            if vtype == "SNV":
                if base == alt:
                    snv_alt += 1

            i += 1

            if i < len(bases) and bases[i] in "+-":
                sign = bases[i]
                i += 1

                m = re.match(r"(\d+)", bases[i:])

                if not m:
                    continue

                n = int(m.group(1))
                i += len(m.group(1))

                seq = bases[i:i+n].upper()
                i += n

                if sign == "+":
                    ins.append(seq)
                else:
                    dels.append(seq)

            continue

        if c in "*<>":
            i += 1
            continue

        i += 1

    if vtype == "SNV":
        return snv_alt

    if vtype == "INS":
        inserted = alt[len(ref):]
        return sum(x == inserted for x in ins)

    if vtype == "DEL":
        deleted = ref[len(alt):]
        return sum(x == deleted for x in dels)

    return None

def audit(method, sample, chrom, pos, ref, alt, mapq, baseq):
    bam = bam_path(method, sample)

    cmd = [
        "samtools", "mpileup",
        "-A",
        "-B",
        "-q", str(mapq),
        "-Q", str(baseq),
        "-d", "1000000",
        "-r", f"{chrom}:{pos}-{pos}",
        str(bam)
    ]

    p = subprocess.run(
        cmd,
        text=True,
        capture_output=True
    )

    if p.returncode != 0:
        return None, None, None

    line = p.stdout.strip()

    if not line:
        return 0, 0, 0.0

    fields = line.split("\t")

    depth = int(fields[3])
    bases = fields[4]

    alt_count = parse_bases(
        bases,
        ref,
        alt
    )

    if alt_count is None:
        return depth, None, None

    af = alt_count / depth if depth else 0.0

    return depth, alt_count, af

df = pd.read_csv(
    INFILE,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

rows = []

for _, r in df.iterrows():

    chrom, pos, ref, alt = parse_key(
        r["Variant_Key"]
    )

    result = {
        "patient_id": r["patient_id"],
        "Original_Status": r["Original_Status"],
        "Gene": r["Gene"],
        "HGVSc": r["HGVSc"],
        "HGVSp": r["HGVSp"],
        "Variant_Key": r["Variant_Key"],
        "Source_Method": r["Source_Method"],
        "Source_Sample": r["Source_Sample"],
        "Opposite_Method": r["Opposite_Method"],
        "Opposite_Sample": r["Opposite_Sample"],
        "Original_Opposite_DP_Q0":
            r["Opposite_BAM_DP"],
        "Original_Opposite_AF_Q0":
            r["Opposite_BAM_AF"]
    }

    for label, mapq, baseq in QUALITY_SETTINGS:

        sdp, salt, saf = audit(
            r["Source_Method"],
            r["Source_Sample"],
            chrom,
            pos,
            ref,
            alt,
            mapq,
            baseq
        )

        odp, oalt, oaf = audit(
            r["Opposite_Method"],
            r["Opposite_Sample"],
            chrom,
            pos,
            ref,
            alt,
            mapq,
            baseq
        )

        result[f"Source_DP_{label}"] = sdp
        result[f"Source_ALT_{label}"] = salt
        result[f"Source_AF_{label}"] = saf

        result[f"Opposite_DP_{label}"] = odp
        result[f"Opposite_ALT_{label}"] = oalt
        result[f"Opposite_AF_{label}"] = oaf

    rows.append(result)

out = pd.DataFrame(rows)

def interpretation(r):

    vtype = variant_type(
        *parse_key(r["Variant_Key"])[2:]
    )

    threshold = 0.20 if vtype == "SNV" else 0.15

    oaf = r["Opposite_AF_Q20"]
    odp = r["Opposite_DP_Q20"]

    if pd.isna(odp) or odp < 10:
        return "NON_EVALUABLE"

    if pd.isna(oaf):
        return "MANUAL_REVIEW"

    if oaf >= threshold:
        return "HIGH_QUALITY_SUPPORT_ABOVE_THRESHOLD"

    if oaf > 0:
        return "HIGH_QUALITY_LOW_LEVEL_SUPPORT"

    return "NO_HIGH_QUALITY_ALT_SUPPORT"

out["HQ_Interpretation"] = out.apply(
    interpretation,
    axis=1
)

out.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

with pd.ExcelWriter(
    OUTXLSX,
    engine="openpyxl"
) as writer:
    out.to_excel(
        writer,
        sheet_name="HQ_Audit",
        index=False
    )

print()
print("=" * 120)
print("HIGH-QUALITY ONE-SIDED P/LP AUDIT")
print("=" * 120)

cols = [
    "patient_id",
    "Original_Status",
    "Gene",
    "HGVSc",
    "Source_AF_Q20",
    "Opposite_DP_Q20",
    "Opposite_ALT_Q20",
    "Opposite_AF_Q20",
    "Opposite_AF_Q30",
    "HQ_Interpretation"
]

print(
    out[cols].to_string(
        index=False
    )
)

print()
print("Summary:")
print(
    out["HQ_Interpretation"]
    .value_counts(dropna=False)
    .to_string()
)

print()
print("Saved:")
print(OUTFILE)
print(OUTXLSX)
