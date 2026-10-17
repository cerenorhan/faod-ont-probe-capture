#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import subprocess
import re
import math

# ============================================================
# PATHS
# ============================================================

PAPER = Path("/home/nadirlab/Desktop/faod_paper/faod_sonuclar_final")

PCR_ROOT = Path("/home/nadirlab/Desktop/faod_pcr/FAOD_PCR_Analysis")
PROBE_ROOT = Path("/home/nadirlab/Desktop/CELEMICS_DATA/FAOD_PROBE_Analysis")

COMPARISON = PAPER / "reanalysis_final" / "strict_plp_comparison.tsv"
MANIFEST = PAPER / "master_sample_manifest.tsv"

OUTDIR = PAPER / "reanalysis_final"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_TSV = OUTDIR / "one_sided_plp_bam_audit.tsv"
OUT_XLSX = OUTDIR / "one_sided_plp_bam_audit.xlsx"


# ============================================================
# PARAMETERS
# Same primary technical thresholds used in the study
# ============================================================

MIN_DP = 10
MIN_AF_SNV = 0.20
MIN_AF_INDEL = 0.15


# ============================================================
# HELPERS
# ============================================================

def parse_variant_key(key):
    """
    Expected:
      chr1:75732852:G:A
      chr1:75732879:C:CT
      chr1:75728472:A:ATTAGGAT
    """
    parts = str(key).split(":", 3)

    if len(parts) != 4:
        raise ValueError(f"Cannot parse Variant_Key: {key}")

    chrom = parts[0]
    pos = int(parts[1])
    ref = parts[2].upper()
    alt = parts[3].upper()

    return chrom, pos, ref, alt


def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"

    if len(alt) > len(ref) and alt.startswith(ref):
        return "INS"

    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"

    if len(ref) == len(alt):
        return "MNV"

    return "COMPLEX"


def bam_path(method, sample):
    root = PCR_ROOT if method == "PCR" else PROBE_ROOT
    f = root / "bam" / f"{sample}.sorted.bam"

    if not f.exists():
        raise FileNotFoundError(f)

    return f


def vcf_path(method, sample):
    root = PCR_ROOT if method == "PCR" else PROBE_ROOT
    f = root / "norm_vcf" / f"{sample}.norm.vcf.gz"

    if not f.exists():
        raise FileNotFoundError(f)

    return f


def report_path(method, sample):
    root = PCR_ROOT if method == "PCR" else PROBE_ROOT
    return root / "report_fix2" / f"{sample}.final.rsid.excel.tsv"


# ============================================================
# PILEUP PARSER
# ============================================================

def parse_pileup_bases(bases, ref, alt):
    """
    Returns exact ALT-support count from samtools mpileup bases string.

    Handles:
      SNV
      left-anchored insertion
      left-anchored deletion

    MNV/complex variants are flagged unsupported here.
    """

    ref = ref.upper()
    alt = alt.upper()

    vtype = variant_type(ref, alt)

    snv_alt_count = 0
    ref_base_count = 0

    insertion_events = []
    deletion_events = []

    i = 0

    while i < len(bases):

        c = bases[i]

        # Start-of-read marker: ^ + mapping quality character
        if c == "^":
            i += 2
            continue

        # End-of-read marker
        if c == "$":
            i += 1
            continue

        # Reference match
        if c in ".,":
            ref_base_count += 1
            i += 1

            # insertion/deletion may immediately follow
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
                    insertion_events.append(seq)
                else:
                    deletion_events.append(seq)

            continue

        # Explicit mismatch base
        if c.upper() in {"A", "C", "G", "T", "N"}:

            base = c.upper()

            if vtype == "SNV" and base == alt:
                snv_alt_count += 1

            i += 1

            # mismatch can also be followed by indel event
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
                    insertion_events.append(seq)
                else:
                    deletion_events.append(seq)

            continue

        # deletion placeholder in subsequent deleted positions
        if c == "*":
            i += 1
            continue

        # reference skip
        if c in "<>":
            i += 1
            continue

        # unexpected character
        i += 1

    if vtype == "SNV":
        return snv_alt_count, ref_base_count, "SNV"

    elif vtype == "INS":
        inserted = alt[len(ref):]
        count = sum(
            1 for x in insertion_events
            if x == inserted
        )
        return count, ref_base_count, f"INS:{inserted}"

    elif vtype == "DEL":
        deleted = ref[len(alt):]
        count = sum(
            1 for x in deletion_events
            if x == deleted
        )
        return count, ref_base_count, f"DEL:{deleted}"

    else:
        return None, ref_base_count, vtype


def run_mpileup(method, sample, chrom, pos, ref, alt):

    bam = bam_path(method, sample)

    region = f"{chrom}:{pos}-{pos}"

    cmd = [
        "samtools",
        "mpileup",
        "-A",
        "-B",
        "-q", "0",
        "-Q", "0",
        "-d", "1000000",
        "-r", region,
        str(bam)
    ]

    p = subprocess.run(
        cmd,
        text=True,
        capture_output=True
    )

    if p.returncode != 0:
        return {
            "BAM_Depth": None,
            "BAM_ALT_Count": None,
            "BAM_ALT_AF": None,
            "BAM_REF_Count": None,
            "Pileup_Event": "ERROR",
            "Pileup_Error": p.stderr.strip()
        }

    line = p.stdout.strip()

    if not line:
        return {
            "BAM_Depth": 0,
            "BAM_ALT_Count": 0,
            "BAM_ALT_AF": 0.0,
            "BAM_REF_Count": 0,
            "Pileup_Event": "NO_COVERAGE",
            "Pileup_Error": ""
        }

    fields = line.split("\t")

    if len(fields) < 5:
        return {
            "BAM_Depth": None,
            "BAM_ALT_Count": None,
            "BAM_ALT_AF": None,
            "BAM_REF_Count": None,
            "Pileup_Event": "MALFORMED_PILEUP",
            "Pileup_Error": line
        }

    depth = int(fields[3])
    bases = fields[4]

    alt_count, ref_count, event = parse_pileup_bases(
        bases,
        ref,
        alt
    )

    if alt_count is None:
        af = None
    else:
        af = alt_count / depth if depth > 0 else 0.0

    return {
        "BAM_Depth": depth,
        "BAM_ALT_Count": alt_count,
        "BAM_ALT_AF": af,
        "BAM_REF_Count": ref_count,
        "Pileup_Event": event,
        "Pileup_Error": ""
    }


# ============================================================
# EXACT NORMALIZED VCF CHECK
# ============================================================

def check_exact_vcf(method, sample, chrom, pos, ref, alt):

    vcf = vcf_path(method, sample)

    cmd = [
        "bcftools",
        "query",
        "-r", f"{chrom}:{pos}-{pos}",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%DP]\t[%AD]\n",
        str(vcf)
    ]

    p = subprocess.run(
        cmd,
        text=True,
        capture_output=True
    )

    if p.returncode != 0:
        return {
            "VCF_Exact_Present": False,
            "VCF_GT": "",
            "VCF_DP": "",
            "VCF_AD": "",
            "VCF_Error": p.stderr.strip()
        }

    exact = False
    gt = ""
    dp = ""
    ad = ""

    for line in p.stdout.splitlines():

        fields = line.split("\t")

        if len(fields) < 4:
            continue

        c = fields[0]
        ppos = int(fields[1])
        r = fields[2].upper()
        alts = fields[3].upper().split(",")

        if (
            c == chrom
            and ppos == pos
            and r == ref
            and alt in alts
        ):
            exact = True

            if len(fields) > 4:
                gt = fields[4]

            if len(fields) > 5:
                dp = fields[5]

            if len(fields) > 6:
                ad = fields[6]

            break

    return {
        "VCF_Exact_Present": exact,
        "VCF_GT": gt,
        "VCF_DP": dp,
        "VCF_AD": ad,
        "VCF_Error": ""
    }


# ============================================================
# REPORT CHECK
# ============================================================

def check_report_locus(method, sample, chrom, pos, ref, alt):

    f = report_path(method, sample)

    if not f.exists():
        return {
            "Report_Exact_Present": False,
            "Report_Any_Same_Locus": False
        }

    df = pd.read_csv(
        f,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    exact = False
    same_locus = False

    pattern = re.compile(
        rf'^"?{re.escape(chrom)}:{pos}\s+'
    )

    for x in df["Variant(g.)"].astype(str):

        y = x.replace('"', '').strip()

        if pattern.search(y):
            same_locus = True

            expected = f"{chrom}:{pos} {ref}>{alt}"

            if y.upper() == expected.upper():
                exact = True

    return {
        "Report_Exact_Present": exact,
        "Report_Any_Same_Locus": same_locus
    }


# ============================================================
# AUDIT INTERPRETATION
# ============================================================

def interpret(depth, alt_count, af, ref, alt, vcf_present):

    if depth is None:
        return "ERROR"

    if depth < MIN_DP:
        return "NON_EVALUABLE_LOW_DEPTH"

    if alt_count is None:
        return "MANUAL_REVIEW_COMPLEX_VARIANT"

    vtype = variant_type(ref, alt)

    threshold = (
        MIN_AF_SNV
        if vtype == "SNV"
        else MIN_AF_INDEL
    )

    if vcf_present:
        return "EXACT_VARIANT_PRESENT_IN_OPPOSITE_VCF"

    if af >= threshold:
        return "BAM_SUPPORTED_ABOVE_AF_THRESHOLD"

    if alt_count > 0:
        return "LOW_LEVEL_BAM_SUPPORT_BELOW_AF_THRESHOLD"

    return "NO_ALT_SUPPORT_AT_ADEQUATE_DEPTH"


# ============================================================
# LOAD INPUTS
# ============================================================

cmp = pd.read_csv(
    COMPARISON,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False
)

# Only confirmed biological samples
manifest = manifest[
    (manifest["role"] == "biological")
    &
    (manifest["pair_status"] == "confirmed")
].copy()

mapping = manifest.set_index("patient_id")


one_sided = cmp[
    cmp["Presence_Status"].isin(
        ["PCR_only", "Probe_only"]
    )
].copy()

print()
print("One-sided strict P/LP observations:", len(one_sided))
print()


# ============================================================
# RUN AUDIT
# ============================================================

results = []

for _, row in one_sided.iterrows():

    patient = row["patient_id"]

    if patient not in mapping.index:
        raise RuntimeError(
            f"{patient} not found as confirmed biological patient in manifest"
        )

    chrom, pos, ref, alt = parse_variant_key(
        row["Variant_Key"]
    )

    if row["Presence_Status"] == "PCR_only":

        source_method = "PCR"
        opposite_method = "Probe"

        source_sample = row["PCR_Sample"]
        opposite_sample = mapping.loc[
            patient,
            "probe_sample"
        ]

        gene = row["PCR_Gene"]
        hgvsc = row["PCR_RefSeq_HGVSc"]
        hgvsp = row["PCR_RefSeq_HGVSp"]
        clinvar = row["PCR_ClinVar_Significance"]
        source_status = row["PCR_Technical_Status"]
        source_af = row["PCR_Allele_Frequency"]
        source_dp = row["PCR_Read_Depth"]

    else:

        source_method = "Probe"
        opposite_method = "PCR"

        source_sample = row["Probe_Sample"]
        opposite_sample = mapping.loc[
            patient,
            "pcr_sample"
        ]

        gene = row["Probe_Gene"]
        hgvsc = row["Probe_RefSeq_HGVSc"]
        hgvsp = row["Probe_RefSeq_HGVSp"]
        clinvar = row["Probe_ClinVar_Significance"]
        source_status = row["Probe_Technical_Status"]
        source_af = row["Probe_Allele_Frequency"]
        source_dp = row["Probe_Read_Depth"]

    print(
        f"{patient}: {row['Presence_Status']} "
        f"{gene} {hgvsc} "
        f"{source_sample} -> {opposite_sample}"
    )

    # Audit source BAM too
    source_bam = run_mpileup(
        source_method,
        source_sample,
        chrom,
        pos,
        ref,
        alt
    )

    # Audit opposite BAM
    opposite_bam = run_mpileup(
        opposite_method,
        opposite_sample,
        chrom,
        pos,
        ref,
        alt
    )

    # VCF exact variant presence
    opposite_vcf = check_exact_vcf(
        opposite_method,
        opposite_sample,
        chrom,
        pos,
        ref,
        alt
    )

    source_vcf = check_exact_vcf(
        source_method,
        source_sample,
        chrom,
        pos,
        ref,
        alt
    )

    # Report locus checks
    opposite_report = check_report_locus(
        opposite_method,
        opposite_sample,
        chrom,
        pos,
        ref,
        alt
    )

    audit_flag = interpret(
        opposite_bam["BAM_Depth"],
        opposite_bam["BAM_ALT_Count"],
        opposite_bam["BAM_ALT_AF"],
        ref,
        alt,
        opposite_vcf["VCF_Exact_Present"]
    )

    results.append({

        "patient_id": patient,

        "Original_Status": row["Presence_Status"],

        "Gene": gene,
        "HGVSc": hgvsc,
        "HGVSp": hgvsp,
        "ClinVar": clinvar,

        "Variant_Key": row["Variant_Key"],
        "Variant_Type": variant_type(ref, alt),

        "Source_Method": source_method,
        "Source_Sample": source_sample,
        "Source_Report_DP": source_dp,
        "Source_Report_AF": source_af,
        "Source_Technical_Status": source_status,

        "Source_BAM_DP": source_bam["BAM_Depth"],
        "Source_BAM_ALT_Count": source_bam["BAM_ALT_Count"],
        "Source_BAM_AF": source_bam["BAM_ALT_AF"],
        "Source_VCF_Exact_Present": source_vcf["VCF_Exact_Present"],

        "Opposite_Method": opposite_method,
        "Opposite_Sample": opposite_sample,

        "Opposite_BAM_DP": opposite_bam["BAM_Depth"],
        "Opposite_BAM_ALT_Count": opposite_bam["BAM_ALT_Count"],
        "Opposite_BAM_AF": opposite_bam["BAM_ALT_AF"],

        "Opposite_VCF_Exact_Present":
            opposite_vcf["VCF_Exact_Present"],

        "Opposite_VCF_GT":
            opposite_vcf["VCF_GT"],

        "Opposite_VCF_DP":
            opposite_vcf["VCF_DP"],

        "Opposite_VCF_AD":
            opposite_vcf["VCF_AD"],

        "Opposite_Report_Exact_Present":
            opposite_report["Report_Exact_Present"],

        "Opposite_Report_Any_Same_Locus":
            opposite_report["Report_Any_Same_Locus"],

        "Audit_Flag":
            audit_flag,

        "Pileup_Event":
            opposite_bam["Pileup_Event"]
    })


out = pd.DataFrame(results)

out.to_csv(
    OUT_TSV,
    sep="\t",
    index=False
)

with pd.ExcelWriter(
    OUT_XLSX,
    engine="openpyxl"
) as writer:

    out.to_excel(
        writer,
        sheet_name="One_Sided_Audit",
        index=False
    )

    for flag, group in out.groupby(
        "Audit_Flag",
        dropna=False
    ):
        sheet = str(flag)[:31]
        group.to_excel(
            writer,
            sheet_name=sheet,
            index=False
        )


# ============================================================
# TERMINAL SUMMARY
# ============================================================

print()
print("=" * 100)
print("ONE-SIDED STRICT P/LP BAM AUDIT")
print("=" * 100)

display_cols = [
    "patient_id",
    "Original_Status",
    "Gene",
    "HGVSc",
    "Source_Sample",
    "Opposite_Sample",
    "Opposite_BAM_DP",
    "Opposite_BAM_ALT_Count",
    "Opposite_BAM_AF",
    "Opposite_VCF_Exact_Present",
    "Audit_Flag"
]

print(
    out[display_cols].to_string(
        index=False
    )
)

print()
print("Audit categories:")
print(
    out["Audit_Flag"]
    .value_counts(dropna=False)
    .to_string()
)

print()
print("Saved:")
print(OUT_TSV)
print(OUT_XLSX)
