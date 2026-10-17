#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# ONT TARGETED-PANEL PIPELINE — CLINICAL STATUS v3.2 — FASTQ + UNFILTERED REPORT
#
# PURPOSE
#   Portable end-to-end ONT targeted-panel workflow:
#     POD5 -> Dorado SUP -> Barbell demultiplexing/trimming -> barcodeXX FASTQ
#     -> minimap2 -> mosdepth -> Clair3 -> ClinVar -> VEP -> final TSV reports.
#
# IMPORTANT PORTABILITY RULE
#   Put this script in the analysis directory and create a resources/ directory
#   beside it. The analysis itself is written into the current working directory.
#
# RECOMMENDED STRUCTURE
#   analysis/
#   ├── ont_panel_pipeline_pod5_dorado_barbell.sh
#   └── resources/
#       ├── hg38.fa
#       ├── hg38.fa.fai
#       ├── hg38.mmi
#       ├── clinvar.chr.vcf.gz
#       ├── clinvar.chr.vcf.gz.tbi
#       ├── targets.bed                 # project-specific target BED
#       └── validation_loci.tsv         # optional
#
# BASIC USAGE — POD5 INPUT
#   cd /path/to/analysis
#   POD5_DIR=/path/to/pod5 bash ./ont_panel_pipeline_pod5_dorado_barbell.sh
#
# PROJECT-SPECIFIC BED WITHOUT RENAMING
#   POD5_DIR=/path/to/pod5 \
#   TARGETS_BED="$PWD/resources/pku_group1_targets.bed" \
#   bash ./ont_panel_pipeline_pod5_dorado_barbell.sh
#
# IMPORTANT DEFAULTS FOR THE CURRENT R10.4.1 / E8.2 / 400 bps / 5 kHz RUN
#   DORADO_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.2.0
#   BARCODE_KIT=SQK-NBD114-96
#   CLAIR3_MODEL=r1041_e82_400bps_sup_v520
#
# BARBELL
#   Conservative demultiplexing is the default. To opt into Barbell --maximize:
#     BARBELL_MAXIMIZE=1 ... bash ./ont_panel_pipeline_pod5_dorado_barbell.sh
#
# FASTQ/BARCODE-DIRECTORY MODE (skip Dorado + Barbell)
#   INPUT_MODE=fastq bash ./ont_panel_pipeline_pod5_dorado_barbell.sh
#
# COMMON OVERRIDES
#   THREADS=16 SORT_MEM=2G VEP_CACHE=$HOME/vep_cache ...
#   DORADO_BIN=/path/to/dorado ...
#   BARCODE_KIT=SQK-NBD114-96 ...
#   DORADO_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.2.0 ...
#   CLAIR3_MODEL=r1041_e82_400bps_sup_v520 ...
###############################################################################

############################
# DIRECTORY VERIFICATION
############################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)" || {
  echo "[ERROR] Could not resolve SCRIPT_DIR" >&2
  exit 1
}

RUN_DIR="$(pwd -P)" || {
  echo "[ERROR] Current working directory is invalid" >&2
  exit 1
}

python3 - <<'PY' >/dev/null 2>&1 || {
import os
os.getcwd()
PY
  echo "[ERROR] Current working directory is invalid for Python. Please 'cd /' and re-enter the run folder, or reopen terminal." >&2
  exit 1
}

############################
# CONFIG
############################
THREADS="${THREADS:-10}"
SORT_MEM="${SORT_MEM:-2G}"

MIN_QUAL="${MIN_QUAL:-20}"
MIN_DP="${MIN_DP:-10}"
MIN_AF="${MIN_AF:-0.20}"
HET_MAX_AF="${HET_MAX_AF:-0.70}"

export MIN_QUAL MIN_DP MIN_AF HET_MAX_AF

RESOURCE_DIR="${RESOURCE_DIR:-$RUN_DIR/resources}"
TMPDIR_LOCAL="${TMPDIR_LOCAL:-$RUN_DIR/tmp}"
VEP_CACHE="${VEP_CACHE:-$HOME/vep_cache}"
CLAIR3_IMAGE="${CLAIR3_IMAGE:-hkubal/clair3}"
VEP_IMAGE="${VEP_IMAGE:-ensemblorg/ensembl-vep:release_115.0}"

REF_FA="${REF_FA:-$RESOURCE_DIR/hg38.fa}"
REF_MMI="${REF_MMI:-$RESOURCE_DIR/hg38.mmi}"
CLINVAR_CHR_VCF="${CLINVAR_CHR_VCF:-$RESOURCE_DIR/clinvar.chr.vcf.gz}"
CLINVAR_CHR_TBI="${CLINVAR_CHR_TBI:-$RESOURCE_DIR/clinvar.chr.vcf.gz.tbi}"
VALIDATION_LOCI="${VALIDATION_LOCI:-$RESOURCE_DIR/validation_loci.tsv}"

# Input mode:
#   pod5  = Dorado basecalling + Barbell demultiplexing, then analysis
#   fastq = barcodeXX directories already exist; skip preprocessing
INPUT_MODE="${INPUT_MODE:-fastq}"
POD5_DIR="${POD5_DIR:-}"

# Reproducible defaults for this chemistry/run. Override for another chemistry.
DORADO_MODEL="${DORADO_MODEL:-dna_r10.4.1_e8.2_400bps_sup@v5.2.0}"
DORADO_DEVICE="${DORADO_DEVICE:-cuda:0}"
BARCODE_KIT="${BARCODE_KIT:-SQK-NBD114-96}"
BARBELL_MAXIMIZE="${BARBELL_MAXIMIZE:-0}"
BARCODE_MIN_READS="${BARCODE_MIN_READS:-100}"
CLAIR3_MODEL="${CLAIR3_MODEL:-r1041_e82_400bps_hac_v410}"

# TARGETS_BED resolution:
#   1) explicit TARGETS_BED=...
#   2) resources/targets.bed
#   3) exactly one resources/*_targets.bed
# This lets each project carry its own target BED without editing this script.
TARGETS_BED="${TARGETS_BED:-}"
TARGETS_BED_PIPELINE="$RESOURCE_DIR/targets.pipeline.sorted.bed"

DORADO_OUTDIR="$RUN_DIR/basecalling"
DORADO_FASTQ="$DORADO_OUTDIR/dorado_sup.no_trim.fastq"
BARBELL_OUTDIR="$RUN_DIR/barbell_demux"
BARBELL_DONE="$BARBELL_OUTDIR/.barbell_complete"

mkdir -p \
  "$RESOURCE_DIR" \
  "$TMPDIR_LOCAL" \
  "$DORADO_OUTDIR" \
  "$BARBELL_OUTDIR" \
  "$RUN_DIR/report_fix2" \
  "$RUN_DIR/report_fix2/final_all" \
  "$RUN_DIR/logs" \
  "$RUN_DIR/merged_fastq" \
  "$RUN_DIR/nanoplot" \
  "$RUN_DIR/bam" \
  "$RUN_DIR/mosdepth_out" \
  "$RUN_DIR/clair3_out" \
  "$RUN_DIR/norm_vcf" \
  "$RUN_DIR/clinvar_vcf" \
  "$RUN_DIR/vep_ensembl" \
  "$RUN_DIR/vep_refseq"

############################
# LOGGING / ERRORS
############################
log(){
  echo "[$(date '+%F %T')] $*"
}

die(){
  echo "[ERROR] $*" >&2
  exit 1
}

trap 'echo "[ERROR] line=${LINENO} command=${BASH_COMMAND}" >&2' ERR

############################
# CHECKS
############################
need_cmd(){ command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"; }

need_cmd awk
need_cmd sed
need_cmd sort
need_cmd grep
need_cmd find
need_cmd bgzip
need_cmd tabix
need_cmd bcftools
need_cmd samtools
need_cmd minimap2
# Container engine: prefer Docker when available, otherwise Podman.
if command -v docker >/dev/null 2>&1; then
  CONTAINER_ENGINE="docker"
elif command -v podman >/dev/null 2>&1; then
  CONTAINER_ENGINE="podman"
else
  die "Neither docker nor podman was found; Clair3 and VEP require a container engine"
fi
log "Container engine: $CONTAINER_ENGINE"

# NanoPlot and mosdepth are intentionally optional.
HAVE_NANOPLOT=0
HAVE_MOSDEPTH=0
command -v NanoPlot >/dev/null 2>&1 && HAVE_NANOPLOT=1
command -v mosdepth >/dev/null 2>&1 && HAVE_MOSDEPTH=1

# Optional: skip NanoPlot QC if explicitly requested.
if [[ "${SKIP_NANOPLOT:-0}" == "1" ]]; then
  HAVE_NANOPLOT=0
fi

if [[ "$HAVE_NANOPLOT" -eq 0 ]]; then
  log "NanoPlot not found; NanoPlot QC will be skipped"
fi
if [[ "$HAVE_MOSDEPTH" -eq 0 ]]; then
  log "mosdepth not found; coverage QC will be skipped"
fi

[[ -s "$REF_FA" ]] || die "Missing reference FASTA: $REF_FA"
[[ -s "$REF_FA.fai" ]] || die "Missing FASTA index: $REF_FA.fai"
[[ -s "$REF_MMI" ]] || die "Missing minimap2 index: $REF_MMI"
[[ -s "$CLINVAR_CHR_VCF" ]] || die "Missing ClinVar chr VCF: $CLINVAR_CHR_VCF"
[[ -s "$CLINVAR_CHR_TBI" ]] || die "Missing ClinVar chr VCF index: $CLINVAR_CHR_TBI"
[[ -d "$VEP_CACHE/homo_sapiens/115_GRCh38" ]] || die "Missing Ensembl VEP cache: $VEP_CACHE/homo_sapiens/115_GRCh38"
[[ -d "$VEP_CACHE/homo_sapiens_refseq/115_GRCh38" ]] || die "Missing RefSeq VEP cache: $VEP_CACHE/homo_sapiens_refseq/115_GRCh38"

############################
# TARGET BED PREPARATION
############################
prepare_targets_bed(){
  local src="${TARGETS_BED:-}"
  local -a candidates=()

  if [[ -z "$src" ]]; then
    if [[ -s "$RESOURCE_DIR/targets.bed" ]]; then
      src="$RESOURCE_DIR/targets.bed"
    else
      shopt -s nullglob
      candidates=("$RESOURCE_DIR"/*_targets.bed)
      shopt -u nullglob
      if [[ "${#candidates[@]}" -eq 1 && -s "${candidates[0]}" ]]; then
        src="${candidates[0]}"
      elif [[ "${#candidates[@]}" -gt 1 ]]; then
        die "Multiple *_targets.bed files found in $RESOURCE_DIR. Set TARGETS_BED explicitly."
      fi
    fi
  fi

  [[ -n "$src" && -s "$src" ]] || die "Missing TARGETS_BED. Add resources/targets.bed, keep exactly one resources/*_targets.bed, or set TARGETS_BED=/path/to/file.bed"
  TARGETS_BED_RAW="$src"

  # Keep only valid BED data rows, then coordinate-sort. This protects downstream
  # tools from project BEDs that were written in exon order rather than genomic order.
  awk 'BEGIN{OFS="\t"}
       /^[[:space:]]*#/ || NF==0 {next}
       NF < 3 {print "[ERROR] Invalid BED row: "$0 > "/dev/stderr"; exit 2}
       $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ {print "[ERROR] Non-numeric BED coordinate: "$0 > "/dev/stderr"; exit 2}
       $2 < 0 || $3 <= $2 {print "[ERROR] Invalid BED interval: "$0 > "/dev/stderr"; exit 2}
       {print}' "$src" \
    | sort -k1,1V -k2,2n -k3,3n \
    > "$TARGETS_BED_PIPELINE"

  [[ -s "$TARGETS_BED_PIPELINE" ]] || die "Prepared TARGETS_BED is empty: $TARGETS_BED_PIPELINE"
  TARGETS_BED="$TARGETS_BED_PIPELINE"
  log "Target BED prepared: $TARGETS_BED"
}

prepare_targets_bed

############################
# DORADO + BARBELL PREPROCESSING
############################
find_working_dorado(){
  local candidates=()
  local c

  if [[ -n "${DORADO_BIN:-}" ]]; then
    candidates+=("$DORADO_BIN")
  fi

  c="$(command -v dorado 2>/dev/null || true)"
  [[ -n "$c" ]] && candidates+=("$c")

  # Local fallback known on this workstation; harmless on other machines.
  candidates+=("$HOME/Desktop/gokce/karnitin_SUP_reanalysis/dorado-2.1.2-linux-x64/bin/dorado")

  for c in "${candidates[@]}"; do
    [[ -n "$c" && -x "$c" ]] || continue
    if "$c" --version >/dev/null 2>&1; then
      DORADO_BIN="$c"
      return 0
    fi
  done

  return 1
}

run_dorado_basecalling(){
  [[ "$INPUT_MODE" == "pod5" ]] || return 0

  [[ -n "$POD5_DIR" ]] || {
    if [[ -d "$RUN_DIR/pod5" ]]; then
      POD5_DIR="$RUN_DIR/pod5"
    else
      die "INPUT_MODE=pod5 but POD5_DIR is not set. Example: POD5_DIR=/path/to/pod5 bash $0"
    fi
  }

  [[ -d "$POD5_DIR" ]] || die "POD5 directory not found: $POD5_DIR"
  find "$POD5_DIR" -maxdepth 1 -type f -name '*.pod5' -print -quit | grep -q . \
    || die "No .pod5 files found in: $POD5_DIR"

  find_working_dorado || die "No working Dorado executable found. Set DORADO_BIN=/path/to/dorado"
  log "Dorado executable: $DORADO_BIN"
  log "Dorado model: $DORADO_MODEL"
  log "Dorado device: $DORADO_DEVICE"

  if [[ -s "$DORADO_FASTQ" ]]; then
    log "Dorado FASTQ already exists; skipping basecalling: $DORADO_FASTQ"
    return 0
  fi

  log "Dorado SUP basecalling from POD5 (no trimming)"
  "$DORADO_BIN" basecaller "$DORADO_MODEL" "$POD5_DIR" \
    --device "$DORADO_DEVICE" \
    --no-trim \
    --emit-fastq \
    > "$DORADO_FASTQ" \
    2> "$RUN_DIR/logs/dorado_basecaller.log"

  [[ -s "$DORADO_FASTQ" ]] || die "Dorado finished without FASTQ output. Check logs/dorado_basecaller.log"
  log "Dorado FASTQ: $DORADO_FASTQ"
}

stage_barbell_fastqs(){
  # Barbell may emit NBxx (Native Barcoding), BCxx, or barcodeXX filenames.
  # Normalize these to barcodeXX directories so downstream Sample stays barcodeXX.
  # Very tiny demux files are usually stray assignments; by default require 100 reads.
  # Override with BARCODE_MIN_READS=0 to stage every non-empty barcode FASTQ.
  local f base digits n sample dest idx reads
  local staged=0
  local skipped=0
  local min_reads="${BARCODE_MIN_READS:-100}"

  while IFS= read -r -d '' f; do
    base="$(basename "$f")"
    digits=""

    if [[ "$base" =~ ^[Nn][Bb]0*([0-9]{1,3})([^0-9]|$) ]]; then
      digits="${BASH_REMATCH[1]}"
    elif [[ "$base" =~ ^[Bb][Cc]0*([0-9]{1,3})([^0-9]|$) ]]; then
      digits="${BASH_REMATCH[1]}"
    elif [[ "$base" =~ ^[Bb][Aa][Rr][Cc][Oo][Dd][Ee]0*([0-9]{1,3})([^0-9]|$) ]]; then
      digits="${BASH_REMATCH[1]}"
    else
      continue
    fi

    n=$((10#$digits))
    printf -v sample 'barcode%02d' "$n"

    # FASTQ has 4 lines/read. Count headers safely via line count / 4 so '@' in
    # quality strings cannot inflate the count. gzip files are streamed.
    case "$f" in
      *.fastq.gz|*.fq.gz)
        reads=$(( $(gzip -cd "$f" | wc -l) / 4 ))
        ;;
      *.fastq|*.fq)
        reads=$(( $(wc -l < "$f") / 4 ))
        ;;
      *)
        continue
        ;;
    esac

    if (( reads < min_reads )); then
      log "Skipping $base: ${reads} read(s) < BARCODE_MIN_READS=${min_reads}"
      skipped=$((skipped+1))
      continue
    fi

    mkdir -p "$RUN_DIR/$sample"
    idx=1
    case "$f" in
      *.fastq.gz|*.fq.gz) dest="$RUN_DIR/$sample/barbell_${idx}.fastq.gz" ;;
      *.fastq|*.fq)       dest="$RUN_DIR/$sample/barbell_${idx}.fastq" ;;
    esac
    while [[ -e "$dest" ]]; do
      idx=$((idx+1))
      case "$f" in
        *.fastq.gz|*.fq.gz) dest="$RUN_DIR/$sample/barbell_${idx}.fastq.gz" ;;
        *)                   dest="$RUN_DIR/$sample/barbell_${idx}.fastq" ;;
      esac
    done

    # Hard-link when possible (no extra disk use); copy only as a fallback.
    ln "$f" "$dest" 2>/dev/null || cp -p "$f" "$dest"
    log "Staged $base -> $sample (${reads} reads)"
    staged=$((staged+1))
  done < <(find "$BARBELL_OUTDIR" -type f \
      \( -iname '*.fastq' -o -iname '*.fq' -o -iname '*.fastq.gz' -o -iname '*.fq.gz' \) \
      -print0)

  [[ "$staged" -gt 0 ]] || die "Barbell completed but no NBxx/BCxx/barcodeXX FASTQ passed staging under $BARBELL_OUTDIR"
  log "Staged $staged Barbell FASTQ file(s) into barcodeXX directories; skipped $skipped below read threshold"
}

run_barbell_demux(){
  [[ "$INPUT_MODE" == "pod5" ]] || return 0
  need_cmd barbell
  [[ -s "$DORADO_FASTQ" ]] || die "Missing Dorado FASTQ for Barbell: $DORADO_FASTQ"

  if [[ -f "$BARBELL_DONE" ]]; then
    log "Barbell completion marker exists; reusing demultiplexed output"
  else
    log "Barbell demultiplexing/trimming with kit: $BARCODE_KIT"
    # Barbell notes that high sample counts may require a larger file-descriptor limit.
    ulimit -n 65535 2>/dev/null || true

    args=(kit -k "$BARCODE_KIT" -i "$DORADO_FASTQ" -o "$BARBELL_OUTDIR")
    if [[ "$BARBELL_MAXIMIZE" == "1" ]]; then
      args+=(--maximize)
      log "Barbell --maximize enabled"
    else
      log "Barbell conservative mode (no --maximize)"
    fi

    barbell "${args[@]}" > "$RUN_DIR/logs/barbell.log" 2>&1
    touch "$BARBELL_DONE"
  fi

  # Recreate barcodeXX staging only if no barcode directories currently exist.
  if ! find "$RUN_DIR" -maxdepth 1 -type d -name 'barcode[0-9]*' -print -quit | grep -q .; then
    stage_barbell_fastqs
  else
    log "barcodeXX directories already exist; staging step skipped"
  fi
}

case "$INPUT_MODE" in
  pod5)
    run_dorado_basecalling
    run_barbell_demux
    ;;
  fastq)
    log "INPUT_MODE=fastq: skipping Dorado and Barbell"
    ;;
  *)
    die "Unknown INPUT_MODE=$INPUT_MODE (allowed: pod5, fastq)"
    ;;
esac

############################
# SAMPLE DETECTION
############################
detect_samples(){
  find -L "$RUN_DIR" -maxdepth 1 -mindepth 1 -type d \
    -name 'barcode*' -printf '%f\n' \
    | grep -E '^barcode[0-9]+$' \
    | sort -V
}

SAMPLES="${SAMPLES_OVERRIDE:-$(detect_samples)}"
[[ -n "$SAMPLES" ]] || die "No FASTQ/FASTQ chunk files found in current directory"
log "Detected samples:"
printf '%s\n' "$SAMPLES"

############################
# HELPERS
############################
merge_fastq(){
  local s="$1"
  local src="$RUN_DIR/$s"
  local out="$RUN_DIR/merged_fastq/${s}.fastq"

  if [[ -s "$out" ]]; then
    log "Merged FASTQ exists: merged_fastq/${s}.fastq"
    return 0
  fi

  [[ -d "$src" ]] || die "Barcode directory not found: $src"

  log "Merging FASTQ for $s"

  mapfile -t files < <(
    find -L "$src" -maxdepth 1 -type f \
      \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) \
      -print | sort -V
  )

  [[ "${#files[@]}" -gt 0 ]] || die "No FASTQ files found in $src"

  : > "$out"

  for f in "${files[@]}"; do
    case "$f" in
      *.fastq.gz|*.fq.gz)
        gzip -cd "$f" >> "$out"
        ;;
      *.fastq|*.fq)
        cat "$f" >> "$out"
        ;;
      *)
        die "Unsupported FASTQ file: $f"
        ;;
    esac
  done

  [[ -s "$out" ]] || die "Merged FASTQ is empty for $s"

  log "Merged ${#files[@]} FASTQ file(s) for $s"
}

run_nanoplot(){
  local s="$1"

  if [[ "$HAVE_NANOPLOT" -eq 0 ]]; then
    log "Skipping NanoPlot for $s"
    return 0
  fi

  local fq="$RUN_DIR/merged_fastq/${s}.fastq"
  local outdir="$RUN_DIR/nanoplot/${s}"
  local logfile="$RUN_DIR/logs/${s}.nanoplot.log"

  mkdir -p "$outdir"

  if [[ -s "$outdir/NanoStats.txt" ]]; then
    log "NanoPlot already done for $s"
    return 0
  fi

  [[ -s "$fq" ]] || die "Missing merged FASTQ for NanoPlot: $fq"

  log "NanoPlot: $s"
  NanoPlot \
    --fastq "$fq" \
    --outdir "$outdir" \
    --threads "$THREADS" \
    --N50 \
    --loglength \
    > "$logfile" 2>&1
}

align_and_sort(){
  local s="$1"
  local fq="$RUN_DIR/merged_fastq/${s}.fastq"
  local bam="$RUN_DIR/bam/${s}.sorted.bam"

  if [[ -s "$bam" && -s "$bam.bai" ]]; then
    log "BAM already exists for $s"
    return 0
  fi

  log "Align + sort: $s"
  minimap2 -ax map-ont --secondary=no -t "$THREADS" "$REF_MMI" "$fq" \
    | samtools sort -@ "$THREADS" -m "$SORT_MEM" -T "$TMPDIR_LOCAL/${s}" -o "$bam" -

  samtools index -@ "$THREADS" "$bam"
}

run_mosdepth(){
  local s="$1"

  if [[ "$HAVE_MOSDEPTH" -eq 0 ]]; then
    log "Skipping mosdepth for $s"
    return 0
  fi

  local bam="$RUN_DIR/bam/${s}.sorted.bam"
  local prefix="$RUN_DIR/mosdepth_out/${s}"

  if [[ -s "${prefix}.mosdepth.summary.txt" || -s "${prefix}.summary.txt" ]]; then
    log "mosdepth already done for $s"
    return 0
  fi

  log "mosdepth: $s"

  if [[ -n "$TARGETS_BED" && -s "$TARGETS_BED" ]]; then
    log "mosdepth using BED: $TARGETS_BED"
    mosdepth --by "$TARGETS_BED" -t "$THREADS" "$prefix" "$bam" \
      > "$RUN_DIR/logs/${s}.mosdepth.log" 2>&1 || true
  else
    log "mosdepth without BED (genome-wide summary)"
    mosdepth -t "$THREADS" "$prefix" "$bam" \
      > "$RUN_DIR/logs/${s}.mosdepth.log" 2>&1 || true
  fi
}

run_clair3(){
  local s="$1"
  local bam_rel="bam/${s}.sorted.bam"
  local outdir_rel="clair3_out/${s}"
  local outdir_abs="$RUN_DIR/${outdir_rel}"
  local logfile="$RUN_DIR/logs/${s}.clair3.log"

  mkdir -p "$outdir_abs"

  if [[ -s "$outdir_abs/merge_output.vcf.gz" ]]; then
    log "Clair3 already done for $s"
    return 0
  fi

  log "Clair3: $s"

  local bed_arg=""
  if [[ -n "$TARGETS_BED" && -s "$TARGETS_BED" ]]; then
    bed_arg="--bed_fn=/data/resources/$(basename "$TARGETS_BED")"
    log "Clair3 restricted to BED: $TARGETS_BED"
  else
    log "Clair3 running without BED restriction"
  fi

  "$CONTAINER_ENGINE" run --rm \
    -v "$RUN_DIR":/data \
    "$CLAIR3_IMAGE" \
    bash --noprofile --norc -lc \
    "/opt/bin/run_clair3.sh \
      --bam_fn=/data/${bam_rel} \
      --ref_fn=/data/resources/hg38.fa \
      --threads=${THREADS} \
      --platform=ont \
      ${bed_arg} \
      --model_path=/opt/models/${CLAIR3_MODEL} \
      --output=/data/${outdir_rel}" \
    > "$logfile" 2>&1

  [[ -s "$outdir_abs/merge_output.vcf.gz" ]] || die "Clair3 finished without merge_output.vcf.gz for $s. Check logs/${s}.clair3.log"
}

normalize_vcf(){
  local s="$1"
  local in="$RUN_DIR/clair3_out/${s}/merge_output.vcf.gz"
  local out="$RUN_DIR/norm_vcf/${s}.norm.vcf.gz"

  [[ -s "$in" ]] || die "Missing Clair3 output VCF for $s: $in"

  if [[ -s "$out" && -s "$out.tbi" ]]; then
    log "Normalized VCF exists for $s"
    return 0
  fi

  log "Normalize VCF: $s"
  bcftools norm -f "$REF_FA" "$in" -Oz -o "$out"
  tabix -f -p vcf "$out"
}

annotate_clinvar(){
  local s="$1"
  local in="$RUN_DIR/norm_vcf/${s}.norm.vcf.gz"
  local out="$RUN_DIR/clinvar_vcf/${s}.annotated.clinvar.vcf.gz"
  local rebgzip="$RUN_DIR/clinvar_vcf/${s}.annotated.clinvar.rebgzip.vcf.gz"

  if [[ -s "$rebgzip" && -s "$rebgzip.tbi" ]]; then
    log "ClinVar annotated rebgzip VCF exists for $s"
    return 0
  fi

  log "ClinVar annotate: $s"
  bcftools annotate \
    -a "$CLINVAR_CHR_VCF" \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO \
    "$in" -Oz -o "$out"

  tabix -f -p vcf "$out"

  zcat "$out" | bgzip -c > "$rebgzip"
  tabix -f -p vcf "$rebgzip"
}

run_vep_ensembl(){
  local s="$1"
  local in_rel="clinvar_vcf/${s}.annotated.clinvar.rebgzip.vcf.gz"
  local out_rel="vep_ensembl/${s}.vep.hgvs.vcf.gz"
  local out_abs="$RUN_DIR/${out_rel}"
  local logfile="$RUN_DIR/logs/${s}.vep_ensembl.log"

  if [[ -s "$out_abs" && -s "$out_abs.tbi" ]]; then
    log "VEP Ensembl already done for $s"
    return 0
  fi

  log "VEP Ensembl: $s"
  "$CONTAINER_ENGINE" run --rm \
    -v "$RUN_DIR":/work \
    -v "$VEP_CACHE":/opt/vep/.vep \
    "$VEP_IMAGE" \
    vep \
      --format vcf \
      -i "/work/${in_rel}" \
      -o "/work/${out_rel}" \
      --vcf --compress_output bgzip \
      --cache --offline \
      --dir_cache /opt/vep/.vep \
      --assembly GRCh38 \
      --cache_version 115 \
      --fasta /work/resources/hg38.fa \
      --hgvs \
      --everything \
      --force_overwrite \
      > "$logfile" 2>&1

  tabix -f -p vcf "$out_abs"
}

run_vep_refseq(){
  local s="$1"
  local in_rel="clinvar_vcf/${s}.annotated.clinvar.rebgzip.vcf.gz"
  local out_rel="vep_refseq/${s}.vep.refseq.hgvs.vcf.gz"
  local out_abs="$RUN_DIR/${out_rel}"
  local logfile="$RUN_DIR/logs/${s}.vep_refseq.log"

  if [[ -s "$out_abs" && -s "$out_abs.tbi" ]]; then
    log "VEP RefSeq already done for $s"
    return 0
  fi

  log "VEP RefSeq: $s"
  "$CONTAINER_ENGINE" run --rm \
    -v "$RUN_DIR":/work \
    -v "$VEP_CACHE":/opt/vep/.vep \
    "$VEP_IMAGE" \
    vep \
      --format vcf \
      -i "/work/${in_rel}" \
      -o "/work/${out_rel}" \
      --vcf --compress_output bgzip \
      --cache --offline --refseq \
      --dir_cache /opt/vep/.vep \
      --assembly GRCh38 \
      --cache_version 115 \
      --fasta /work/resources/hg38.fa \
      --hgvs \
      --force_overwrite \
      > "$logfile" 2>&1

  tabix -f -p vcf "$out_abs"
}

############################
# MAIN ANALYSIS PIPELINE
############################
for s in $SAMPLES; do
  log "=== MAIN PIPELINE: $s ==="
  merge_fastq "$s"
  run_nanoplot "$s"
  align_and_sort "$s"
  run_mosdepth "$s"
  run_clair3 "$s"
  normalize_vcf "$s"
  annotate_clinvar "$s"
  run_vep_ensembl "$s"
  run_vep_refseq "$s"
done

############################
# REPORT SCRIPT GENERATION
############################
cat > "$RUN_DIR/report_fix2/make_ensembl_best.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
vcf="$1"
sample="$(basename "$vcf" .vep.hgvs.vcf.gz)"
bcftools +split-vep "$vcf" -d -A tab \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%HGVSc\t%HGVSp\t%Existing_variation\n' \
| awk -F'\t' 'BEGIN{OFS="\t"}
{
  key=$1":"$2":"$3":"$4
  score=0
  if($8!="." && $8!="") score+=100
  if($9!="." && $9!="") score+=50
  if($7=="YES") score+=10
  if($10!="." && $10!="") score+=2
  if($11!="." && $11!="") score+=1
  if(!(key in best) || score>best[key]){ best[key]=score; line[key]=$0 }
}
END{ for(k in line) print line[k] }' \
| sort -k1,1V -k2,2n > "report_fix2/${sample}.ensembl.best.tsv"
BASH
chmod +x "$RUN_DIR/report_fix2/make_ensembl_best.sh"

cat > "$RUN_DIR/report_fix2/make_refseq_best.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
vcf="$1"
sample="$(basename "$vcf" .vep.refseq.hgvs.vcf.gz)"
bcftools +split-vep "$vcf" -d -A tab \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%HGVSc\t%HGVSp\n' \
| awk -F'\t' 'BEGIN{OFS="\t"}
{
  key=$1":"$2":"$3":"$4
  np="."; prot="."
  hgvsp=$10; gsub(/%3D/,"=",hgvsp)

  if(hgvsp ~ /^NP_[0-9]+\.[0-9]+:p\./){
      split(hgvsp,a,":")
      np=a[1]
      prot=a[2]
  }

  feature=$6
  mane=$7
  mane_plus=$8
  hgvsc=$9

  score=0

  # Highest priority: RefSeq transcript explicitly matching MANE Select
  if(mane!="." && mane!="" && feature==mane) score+=1000

  # Second priority: MANE Plus Clinical match
  if(mane_plus!="." && mane_plus!="" && feature==mane_plus) score+=500

  # Then any MANE-tagged transcript
  if(mane!="." && mane!="") score+=200
  if(mane_plus!="." && mane_plus!="") score+=100

  if(feature ~ /^NM_/) score+=50
  if(hgvsc!="." && hgvsc!="") score+=10
  if(np!=".") score+=5
  if(hgvsp!="." && hgvsp!="") score+=2

  if(!(key in best) || score>best[key]){
      best[key]=score
      line[key]=$1 OFS $2 OFS $3 OFS $4 OFS $5 OFS feature OFS np OFS prot OFS hgvsc OFS hgvsp
  }
}
END{ for(k in line) print line[k] }' \
| sort -k1,1V -k2,2n > "report_fix2/${sample}.refseq.best.tsv"
BASH
chmod +x "$RUN_DIR/report_fix2/make_refseq_best.sh"

cat > "$RUN_DIR/report_fix2/make_core_tsv.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
vcf="$1"
sample="$(basename "$vcf" .annotated.clinvar.rebgzip.vcf.gz)"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%GQ\t%DP\t%AF\t%AD]\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t%INFO/GENEINFO\n' "$vcf" > "report_fix2/${sample}.core.tsv"
BASH
chmod +x "$RUN_DIR/report_fix2/make_core_tsv.sh"

cat > "$RUN_DIR/report_fix2/join_core_ensembl.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
sample="$1"
awk -F'\t' 'BEGIN{OFS="\t"}
FNR==NR{
  key=$1":"$2":"$3":"$4
  egene[key]=$5; enst[key]=$6; canonical[key]=$7; mane[key]=$8; mane_plus[key]=$9; mane_hgvsc[key]=$10; mane_hgvsp[key]=$11; rsvar[key]=$12
  next
}
{
  chr=$1; pos=$2; ref=$3; alt=$4; qual=$5; gt=$6; gq=$7; dp=$8; af=$9; ad=$10; clnsig=$11; rev=$12; dis=$13; geneinfo=$14
  key=chr":"pos":"ref":"alt
  gene=egene[key]
  if(gene=="" || gene=="."){
    if(geneinfo!="." && geneinfo!=""){ split(geneinfo,a,":"); gene=a[1] } else gene="."
  }
  print chr,pos,ref,alt,qual,gt,gq,dp,af,ad,clnsig,rev,dis,gene,enst[key],canonical[key],mane[key],mane_plus[key],mane_hgvsc[key],mane_hgvsp[key],rsvar[key],geneinfo
}' "report_fix2/${sample}.ensembl.best.tsv" "report_fix2/${sample}.core.tsv" > "report_fix2/${sample}.core_ensembl.tsv"
BASH
chmod +x "$RUN_DIR/report_fix2/join_core_ensembl.sh"

cat > "$RUN_DIR/report_fix2/join_with_refseq_and_finalize.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
sample="$1"
allout="report_fix2/${sample}.all_technical.tsv"
awk -v sample="$sample" -v allout="$allout" -F'\t' 'BEGIN{OFS="\t"}
FNR==NR{
  key=$1":"$2":"$3":"$4
  rgene[key]=$5; rnm[key]=$6; rnp[key]=$7; rprot[key]=$8; rhgvsc[key]=$9; rhgvsp[key]=$10
  next
}
BEGINFILE{
  if(FILENAME ~ /core_ensembl.tsv$/){
    header="Sample\tGene\tRefSeq_NM\tRefSeq_NP\tProtein_Change\tVariant(g.)\tGenotype\tZygosity\tRead_Depth\tAllele_Frequency\tAllele_Depth\tGenotype_Quality\tVariant_Quality\tMANE_SELECT\tMANE_PLUS_CLINICAL\tMANE_HGVSc\tMANE_HGVSp\tCANONICAL\tClinVar_Significance\tClinVar_ReviewStatus\tClinVar_Disease\tGeneInfo\trsID\tRefSeq_HGVSc\tRefSeq_HGVSp\tTechnical_Status"
    print header
    print header > allout
  }
}
function is_plp(sig, s){
  s=tolower(sig)
  if(s ~ /conflict/) return 0
  return (s ~ /pathogenic/ || s ~ /likely_pathogenic/)
}
{
  chr=$1; pos=$2; ref=$3; alt=$4; qual=$5; gt=$6; gq=$7; dp=$8; af=$9; ad=$10; clnsig=$11; rev=$12; dis=$13; gene=$14; enst=$15; canonical=$16; mane=$17; mane_plus=$18; mane_hgvsc=$19; mane_hgvsp=$20; rsid=$21; geneinfo=$22
  key=chr":"pos":"ref":"alt
  gene2=rgene[key]; if(gene=="" || gene==".") gene=gene2
  sub(/:.*/,"",gene); sub(/\|.*/,"",gene)
  nm=rnm[key]; np=rnp[key]; prot=rprot[key]; ref_hgvsc=rhgvsc[key]; ref_hgvsp=rhgvsp[key]
  if(nm=="") nm="."
  if(np=="") np="."
  if(prot=="") prot="."
  if(ref_hgvsc=="") ref_hgvsc="."
  if(ref_hgvsp=="") ref_hgvsp="."

  status="FAIL_MISSING_METRICS"; z="."
  if(qual!="." && qual!="" && dp!="." && dp!="" && af!="." && af!=""){
    qual_num=qual+0; dp_num=dp+0; af_num=af+0
    if(af_num >= ENVIRON["MIN_AF"] && af_num <= ENVIRON["HET_MAX_AF"]) z="heterozygote"
    else if(af_num > ENVIRON["HET_MAX_AF"]) z="homozygote"

    if(dp_num < ENVIRON["MIN_DP"]) status="FAIL_LOW_DP"
    else if(af_num < ENVIRON["MIN_AF"]) status="FAIL_LOW_AF"
    else if(qual_num >= ENVIRON["MIN_QUAL"]) status="PASS_STANDARD"
    else if(is_plp(clnsig)) status="PASS_LOW_QUAL_DP_AF_SUPPORTED_CLINVAR_PLP"
    else status="PASS_LOW_QUAL_DP_AF_SUPPORTED"
  }

  varg=chr ":" pos " " ref ">" alt
  row=sample OFS gene OFS nm OFS np OFS prot OFS varg OFS "GT=" gt OFS z OFS dp OFS af OFS ad OFS gq OFS qual OFS mane OFS mane_plus OFS mane_hgvsc OFS mane_hgvsp OFS canonical OFS clnsig OFS rev OFS dis OFS geneinfo OFS rsid OFS ref_hgvsc OFS ref_hgvsp OFS status
  print row > allout
  if(status !~ /^FAIL_/) print row
}' "report_fix2/${sample}.refseq.best.tsv" "report_fix2/${sample}.core_ensembl.tsv" > "report_fix2/${sample}.final.tsv"
BASH
chmod +x "$RUN_DIR/report_fix2/join_with_refseq_and_finalize.sh"

cat > "$RUN_DIR/report_fix2/add_rsid.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
sample="$1"
bcftools +split-vep "vep_ensembl/${sample}.vep.hgvs.vcf.gz" -d -A tab \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%Existing_variation\n' \
| awk -F'\t' 'BEGIN{OFS="\t"}
{
  key=$1":"$2":"$3":"$4
  rs="."
  if($5 ~ /rs[0-9]+/){ n=split($5,a,"&"); for(i=1;i<=n;i++) if(a[i] ~ /^rs[0-9]+$/){ rs=a[i]; break } }
  print key,rs
}' | sort -u > "report_fix2/${sample}.rsid_map.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
FNR==NR{ rs[$1]=$2; next }
NR==1{ print; next }
{
  split($6,a,/[: >]/); key=a[1]":"a[2]":"a[3]":"a[4]
  if((key in rs) && ($23=="." || $23=="")) $23=rs[key]
  print
}' "report_fix2/${sample}.rsid_map.tsv" "report_fix2/${sample}.final.tsv" > "report_fix2/${sample}.final.rsid.tsv"
BASH
chmod +x "$RUN_DIR/report_fix2/add_rsid.sh"

############################
# PER-SAMPLE REPORTS
############################
for s in $SAMPLES; do
  log "=== REPORT PIPELINE: $s ==="
  bash "$RUN_DIR/report_fix2/make_core_tsv.sh" "$RUN_DIR/clinvar_vcf/${s}.annotated.clinvar.rebgzip.vcf.gz"
  bash "$RUN_DIR/report_fix2/make_ensembl_best.sh" "$RUN_DIR/vep_ensembl/${s}.vep.hgvs.vcf.gz"
  bash "$RUN_DIR/report_fix2/make_refseq_best.sh" "$RUN_DIR/vep_refseq/${s}.vep.refseq.hgvs.vcf.gz"
  bash "$RUN_DIR/report_fix2/join_core_ensembl.sh" "$s"
  bash "$RUN_DIR/report_fix2/join_with_refseq_and_finalize.sh" "$s"
  bash "$RUN_DIR/report_fix2/add_rsid.sh" "$s"

  # Unfiltered caller-emitted report: no QUAL / DP / AF thresholds are applied here.
  # This retains every variant emitted by Clair3 within the target BED, while keeping
  # Technical_Status only as an audit annotation. It does NOT rescue variants that
  # Clair3 never emitted; those remain in the separate BAM-pileup safety-net output.
  awk -F'\t' 'BEGIN{OFS="\t"}
    FNR==NR{rs[$1]=$2; next}
    NR==1{print; next}
    {
      split($6,a,/[: >]/); key=a[1]":"a[2]":"a[3]":"a[4]
      if((key in rs) && ($23=="." || $23=="")) $23=rs[key]
      print
    }' \
    "$RUN_DIR/report_fix2/${s}.rsid_map.tsv" \
    "$RUN_DIR/report_fix2/${s}.all_technical.tsv" \
    > "$RUN_DIR/report_fix2/${s}.unfiltered.rsid.tsv"

  awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {$6="\""$6"\""; $7="\""$7"\""; print}' \
    "$RUN_DIR/report_fix2/${s}.final.rsid.tsv" > "$RUN_DIR/report_fix2/${s}.final.rsid.excel.tsv"

  awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {$6="\""$6"\""; $7="\""$7"\""; print}' \
    "$RUN_DIR/report_fix2/${s}.unfiltered.rsid.tsv" > "$RUN_DIR/report_fix2/${s}.unfiltered.rsid.excel.tsv"
done

############################
# MERGED OUTPUTS
############################
FIRST_SAMPLE="$(for s in $SAMPLES; do printf '%s\n' "$s"; break; done)"
[[ -n "$FIRST_SAMPLE" ]] || die "No samples detected for merge"

head -n 1 "$RUN_DIR/report_fix2/${FIRST_SAMPLE}.final.rsid.tsv" > "$RUN_DIR/report_fix2/total_final.rsid.tsv"
for s in $SAMPLES; do
  tail -n +2 "$RUN_DIR/report_fix2/${s}.final.rsid.tsv" >> "$RUN_DIR/report_fix2/total_final.rsid.tsv"
done

awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {$6="\""$6"\""; $7="\""$7"\""; print}' \
  "$RUN_DIR/report_fix2/total_final.rsid.tsv" > "$RUN_DIR/report_fix2/total_final.rsid.excel.tsv"

# Technical audit: includes FAIL calls; routine final TSV excludes them.
head -n 1 "$RUN_DIR/report_fix2/${FIRST_SAMPLE}.all_technical.tsv" > "$RUN_DIR/report_fix2/total_all_technical.tsv"
for s in $SAMPLES; do
  tail -n +2 "$RUN_DIR/report_fix2/${s}.all_technical.tsv" >> "$RUN_DIR/report_fix2/total_all_technical.tsv"
done

# Completely unfiltered post-caller report: every Clair3-emitted variant in the target BED,
# irrespective of QUAL, DP, AF, or Technical_Status. This is the closest analogue to
# manually scanning a raw caller VCF while retaining our ClinVar/VEP annotations.
head -n 1 "$RUN_DIR/report_fix2/${FIRST_SAMPLE}.unfiltered.rsid.tsv" > "$RUN_DIR/report_fix2/total_unfiltered.rsid.tsv"
for s in $SAMPLES; do
  tail -n +2 "$RUN_DIR/report_fix2/${s}.unfiltered.rsid.tsv" >> "$RUN_DIR/report_fix2/total_unfiltered.rsid.tsv"
done

awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {$6="\""$6"\""; $7="\""$7"\""; print}' \
  "$RUN_DIR/report_fix2/total_unfiltered.rsid.tsv" > "$RUN_DIR/report_fix2/total_unfiltered.rsid.excel.tsv"

# Research-only convenience output for technically supported ClinVar VUS calls.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} tolower($19) ~ /uncertain/ {print}' \
  "$RUN_DIR/report_fix2/total_final.rsid.tsv" > "$RUN_DIR/report_fix2/research_vus.tsv"

# Safety net for known ClinVar P/LP SNVs that Clair3 did not emit.
#
# v3.2 SAFETY-NET V2
#
# A ClinVar P/LP SNV that was not emitted by Clair3 is checked directly
# against the BAM. If DP/AF support passes the assay thresholds, the variant:
#
#   1. remains explicitly marked as caller-missed,
#   2. is annotated independently by VEP Ensembl + VEP RefSeq,
#   3. receives MANE-first transcript selection,
#   4. is added to the supported per-sample final report,
#   5. is included in total_final and total_complete.
#
# The original Clair3 VCF is NEVER modified.
# Provenance remains:
#   Record_Source = BAM_PILEUP_SAFETY_NET
#   Caller_Emitted = NO
#   Technical_Status = SUPPORTED_CALLER_MISSED_CLINVAR_PLP_SNV

python3 - "$RUN_DIR" "$CLINVAR_CHR_VCF" "$TARGETS_BED" "$MIN_DP" "$MIN_AF" $SAMPLES <<'PYPLP'
import csv
import subprocess
import sys
from pathlib import Path

run_dir = Path(sys.argv[1])
clinvar = Path(sys.argv[2])
target_bed = Path(sys.argv[3])

min_dp = int(float(sys.argv[4]))
min_af = float(sys.argv[5])

samples = sys.argv[6:]

out = run_dir / "report_fix2" / "clinvar_plp_caller_missed_snv.tsv"
ref_fa = run_dir / "resources" / "hg38.fa"


def is_plp(sig):

    s = (sig or "").lower()

    if "conflict" in s:
        return False

    return (
        "pathogenic" in s
        or "likely_pathogenic" in s
    )


def raw_has_variant(sample, chrom, pos, ref, alt):

    vcf = (
        run_dir
        / "clair3_out"
        / sample
        / "merge_output.vcf.gz"
    )

    if not vcf.exists():
        return False

    p = subprocess.run(
        [
            "bcftools", "query",
            "-r", f"{chrom}:{pos}-{pos}",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
            str(vcf)
        ],
        capture_output=True,
        text=True,
        check=False
    )

    for line in p.stdout.splitlines():

        f = line.split("\t")

        if len(f) < 4:
            continue

        if (
            f[0] == chrom
            and f[1] == str(pos)
            and f[2] == ref
            and alt in f[3].split(",")
        ):
            return True

    return False


def snv_pileup(sample, chrom, pos, alt):

    bam = run_dir / "bam" / f"{sample}.sorted.bam"

    if not bam.exists():
        return 0, 0, 0, 0.0

    p = subprocess.run(
        [
            "samtools", "mpileup",
            "-aa",
            "-B",
            "-Q", "0",
            "-q", "0",
            "-f", str(ref_fa),
            "-r", f"{chrom}:{pos}-{pos}",
            str(bam)
        ],
        capture_output=True,
        text=True,
        check=False
    )

    if not p.stdout.strip():
        return 0, 0, 0, 0.0

    f = p.stdout.rstrip("\n").split("\t")

    if len(f) < 5:
        return 0, 0, 0, 0.0

    mp_depth = int(f[3])
    bases = f[4]

    ref_count = 0
    alt_count = 0

    i = 0

    while i < len(bases):

        c = bases[i]

        if c == "^":
            i += 2
            continue

        if c == "$":
            i += 1
            continue

        if c in "+-":

            j = i + 1

            while j < len(bases) and bases[j].isdigit():
                j += 1

            n = int(bases[i+1:j] or 0)

            i = j + n
            continue

        if c in ".,":
            ref_count += 1

        elif c.upper() == alt.upper():
            alt_count += 1

        i += 1

    informative = ref_count + alt_count

    af = (
        alt_count / informative
        if informative
        else 0.0
    )

    return (
        mp_depth,
        informative,
        alt_count,
        af
    )


q = subprocess.run(
    [
        "bcftools", "query",
        "-R", str(target_bed),
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\t%ID\t"
        "%INFO/CLNSIG\t%INFO/CLNREVSTAT\t"
        "%INFO/CLNDN\t%INFO/GENEINFO\n",
        str(clinvar)
    ],
    capture_output=True,
    text=True,
    check=True
)


candidates = []

for line in q.stdout.splitlines():

    f = line.split("\t")

    if len(f) < 9:
        continue

    (
        chrom,
        pos,
        ref,
        alts,
        vid,
        sig,
        rev,
        disease,
        geneinfo
    ) = f[:9]

    if not is_plp(sig):
        continue

    for alt in alts.split(","):

        # Current safety-net is deliberately SNV-only.
        if (
            len(ref) == 1
            and len(alt) == 1
            and alt.upper() in {"A","C","G","T"}
        ):

            candidates.append(
                (
                    chrom,
                    int(pos),
                    ref,
                    alt,
                    vid,
                    sig,
                    rev,
                    disease,
                    geneinfo
                )
            )


with out.open("w", newline="") as fh:

    w = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow([
        "Sample",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "ClinVar_ID",
        "ClinVar_Significance",
        "ClinVar_ReviewStatus",
        "ClinVar_Disease",
        "GeneInfo",
        "Mpileup_Depth",
        "Informative_Depth",
        "Alt_Reads",
        "Alt_AF",
        "Technical_Status"
    ])

    for sample in samples:

        for (
            chrom,
            pos,
            ref,
            alt,
            vid,
            sig,
            rev,
            disease,
            geneinfo
        ) in candidates:

            if raw_has_variant(
                sample,
                chrom,
                pos,
                ref,
                alt
            ):
                continue

            (
                mpdp,
                idp,
                alt_reads,
                af
            ) = snv_pileup(
                sample,
                chrom,
                pos,
                alt
            )

            if (
                idp >= min_dp
                and af >= min_af
            ):

                w.writerow([
                    sample,
                    chrom,
                    pos,
                    ref,
                    alt,
                    vid,
                    sig,
                    rev,
                    disease,
                    geneinfo,
                    mpdp,
                    idp,
                    alt_reads,
                    f"{af:.4f}",
                    "SUPPORTED_CALLER_MISSED_CLINVAR_PLP_SNV"
                ])
PYPLP


###############################################################################
# SAFETY-NET ANNOTATION
###############################################################################

SAFETY_DIR="$RUN_DIR/safetynet_annotation"

rm -rf "$SAFETY_DIR"
mkdir -p "$SAFETY_DIR"

SAFETY_MISSED="$RUN_DIR/report_fix2/clinvar_plp_caller_missed_snv.tsv"

N_SAFETY="$(
  awk 'NR>1{n++} END{print n+0}' "$SAFETY_MISSED"
)"

log "Caller-missed supported ClinVar P/LP SNVs: $N_SAFETY"


if [[ "$N_SAFETY" -gt 0 ]]; then

  log "Annotating caller-missed ClinVar P/LP SNVs"

  SAFETY_UNIQUE="$SAFETY_DIR/safetynet_unique.tsv"
  SAFETY_VCF="$SAFETY_DIR/safetynet.vcf"
  SAFETY_VCFGZ="$SAFETY_DIR/safetynet.vcf.gz"

  SAFETY_VEP_ENS="$SAFETY_DIR/safetynet.vep.ensembl.vcf.gz"
  SAFETY_VEP_REF="$SAFETY_DIR/safetynet.vep.refseq.vcf.gz"


  # ----------------------------------------------------------
  # Unique genomic variants
  # ----------------------------------------------------------

  awk -F'\t' '
    BEGIN{OFS="\t"}
    NR==1{next}
    !seen[$2 FS $3 FS $4 FS $5]++{
      print $2,$3,$4,$5
    }
  ' "$SAFETY_MISSED" \
  | sort -k1,1V -k2,2n \
  > "$SAFETY_UNIQUE"


  # ----------------------------------------------------------
  # Minimal valid GRCh38 VCF with contig definitions
  # ----------------------------------------------------------

  {
    echo '##fileformat=VCFv4.2'
    echo '##reference=GRCh38'

    awk '
      BEGIN{OFS=""}
      {
        print "##contig=<ID=" $1 ",length=" $2 ">"
      }
    ' "$REF_FA.fai"

    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'

    awk -F'\t' '
      BEGIN{OFS="\t"}
      {
        print $1,$2,".",$3,$4,".","PASS","."
      }
    ' "$SAFETY_UNIQUE"

  } > "$SAFETY_VCF"


  bcftools norm \
    -f "$REF_FA" \
    -m -any \
    "$SAFETY_VCF" \
    -Oz \
    -o "$SAFETY_VCFGZ"

  tabix -f -p vcf "$SAFETY_VCFGZ"


  # ----------------------------------------------------------
  # VEP Ensembl
  # ----------------------------------------------------------

  "$CONTAINER_ENGINE" run --rm \
    -v "$RUN_DIR":/work \
    -v "$VEP_CACHE":/opt/vep/.vep \
    "$VEP_IMAGE" \
    vep \
      --format vcf \
      -i /work/safetynet_annotation/safetynet.vcf.gz \
      -o /work/safetynet_annotation/safetynet.vep.ensembl.vcf.gz \
      --vcf \
      --compress_output bgzip \
      --cache \
      --offline \
      --dir_cache /opt/vep/.vep \
      --assembly GRCh38 \
      --cache_version 115 \
      --fasta /work/resources/hg38.fa \
      --hgvs \
      --everything \
      --force_overwrite

  tabix -f -p vcf "$SAFETY_VEP_ENS"


  # ----------------------------------------------------------
  # VEP RefSeq
  # ----------------------------------------------------------

  "$CONTAINER_ENGINE" run --rm \
    -v "$RUN_DIR":/work \
    -v "$VEP_CACHE":/opt/vep/.vep \
    "$VEP_IMAGE" \
    vep \
      --format vcf \
      -i /work/safetynet_annotation/safetynet.vcf.gz \
      -o /work/safetynet_annotation/safetynet.vep.refseq.vcf.gz \
      --vcf \
      --compress_output bgzip \
      --cache \
      --offline \
      --refseq \
      --dir_cache /opt/vep/.vep \
      --assembly GRCh38 \
      --cache_version 115 \
      --fasta /work/resources/hg38.fa \
      --hgvs \
      --everything \
      --force_overwrite

  tabix -f -p vcf "$SAFETY_VEP_REF"


  # ----------------------------------------------------------
  # MANE-first Ensembl record
  # ----------------------------------------------------------

  bcftools +split-vep \
    "$SAFETY_VEP_ENS" \
    -d \
    -A tab \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%CANONICAL\t%MANE_SELECT\t%MANE_PLUS_CLINICAL\t%HGVSc\t%HGVSp\t%Existing_variation\n' \
  | awk -F'\t' '
    BEGIN{OFS="\t"}
    {
      key=$1":"$2":"$3":"$4

      score=0

      if($8!="." && $8!="") score+=1000
      if($9!="." && $9!="") score+=500
      if($7=="YES") score+=100
      if($10!="." && $10!="") score+=10
      if($11!="." && $11!="") score+=5

      if(
        !(key in best)
        || score > best[key]
      ){
        best[key]=score
        line[key]=$0
      }
    }

    END{
      for(k in line)
        print line[k]
    }
  ' \
  | sort -k1,1V -k2,2n \
  > "$SAFETY_DIR/ensembl.best.tsv"


  # ----------------------------------------------------------
  # All RefSeq consequences
  # ----------------------------------------------------------

  bcftools +split-vep \
    "$SAFETY_VEP_REF" \
    -d \
    -A tab \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%HGVSc\t%HGVSp\n' \
    > "$SAFETY_DIR/refseq.all.tsv"


  # ----------------------------------------------------------
  # RefSeq MANE-first selection
  # ----------------------------------------------------------

  python3 - \
    "$SAFETY_DIR/ensembl.best.tsv" \
    "$SAFETY_DIR/refseq.all.tsv" \
    "$SAFETY_DIR/refseq.best.tsv" <<'PYREFBEST'

import csv
import sys

ens_file = sys.argv[1]
ref_file = sys.argv[2]
out_file = sys.argv[3]


mane = {}


with open(
    ens_file,
    newline=""
) as fh:

    for row in csv.reader(
        fh,
        delimiter="\t"
    ):

        if len(row) < 12:
            continue

        key = tuple(row[:4])

        mane_select = (
            row[7]
            if len(row) > 7
            else "."
        )

        mane_plus = (
            row[8]
            if len(row) > 8
            else "."
        )

        mane[key] = (
            mane_select,
            mane_plus
        )


best = {}


with open(
    ref_file,
    newline=""
) as fh:

    for row in csv.reader(
        fh,
        delimiter="\t"
    ):

        if len(row) < 8:
            continue

        (
            chrom,
            pos,
            ref,
            alt,
            gene,
            feature,
            hgvsc,
            hgvsp
        ) = row[:8]

        key = (
            chrom,
            pos,
            ref,
            alt
        )

        ms, mp = mane.get(
            key,
            (".",".")
        )

        score = 0

        if (
            ms not in ("",".")
            and feature == ms
        ):
            score += 10000

        if (
            mp not in ("",".")
            and feature == mp
        ):
            score += 5000

        if feature.startswith("NM_"):
            score += 100

        if hgvsc not in ("","."):
            score += 20

        if hgvsp not in ("","."):
            score += 10

        np = "."
        protein = "."

        if ":" in hgvsp:

            a,b = hgvsp.split(
                ":",
                1
            )

            if a.startswith("NP_"):
                np = a
                protein = b

        candidate = [
            chrom,
            pos,
            ref,
            alt,
            gene,
            feature,
            np,
            protein,
            hgvsc,
            hgvsp
        ]

        if (
            key not in best
            or score > best[key][0]
        ):
            best[key] = (
                score,
                candidate
            )


with open(
    out_file,
    "w",
    newline=""
) as oh:

    w = csv.writer(
        oh,
        delimiter="\t",
        lineterminator="\n"
    )

    for key in sorted(
        best,
        key=lambda x: (
            x[0],
            int(x[1])
        )
    ):
        w.writerow(
            best[key][1]
        )

PYREFBEST


  ###########################################################################
  # Merge annotated safety-net variants into per-sample final reports
  ###########################################################################

  python3 - "$RUN_DIR" <<'PYMERGERESCUE'

import csv
import re
import sys
from pathlib import Path


run_dir = Path(sys.argv[1])

report_dir = (
    run_dir
    / "report_fix2"
)

safety_dir = (
    run_dir
    / "safetynet_annotation"
)

missed = (
    report_dir
    / "clinvar_plp_caller_missed_snv.tsv"
)

ens_file = (
    safety_dir
    / "ensembl.best.tsv"
)

ref_file = (
    safety_dir
    / "refseq.best.tsv"
)


def clean(x):

    x = (
        str(x).strip()
        if x is not None
        else "."
    )

    return (
        x
        if x not in {
            "",
            ".",
            "NA",
            "N/A",
            "None",
            "nan"
        }
        else "."
    )


def first_rsid(x):

    m = re.search(
        r"\brs\d+\b",
        x or ""
    )

    return (
        m.group(0)
        if m
        else "."
    )


ens = {}


with ens_file.open(
    newline=""
) as fh:

    for r in csv.reader(
        fh,
        delimiter="\t"
    ):

        if len(r) < 12:
            continue

        key = tuple(r[:4])

        ens[key] = {
            "Gene":clean(r[4]),
            "CANONICAL":clean(r[6]),
            "MANE_SELECT":clean(r[7]),
            "MANE_PLUS_CLINICAL":clean(r[8]),
            "MANE_HGVSc":clean(r[9]),
            "MANE_HGVSp":clean(r[10]),
            "Existing_variation":clean(r[11]),
        }


refseq = {}


with ref_file.open(
    newline=""
) as fh:

    for r in csv.reader(
        fh,
        delimiter="\t"
    ):

        if len(r) < 10:
            continue

        key = tuple(r[:4])

        refseq[key] = {
            "Gene":clean(r[4]),
            "RefSeq_NM":clean(r[5]),
            "RefSeq_NP":clean(r[6]),
            "Protein_Change":clean(r[7]),
            "RefSeq_HGVSc":clean(r[8]),
            "RefSeq_HGVSp":clean(r[9]),
        }


rescues = []


with missed.open(
    newline=""
) as fh:

    for m in csv.DictReader(
        fh,
        delimiter="\t"
    ):

        key = (
            m["CHROM"],
            m["POS"],
            m["REF"],
            m["ALT"]
        )

        e = ens.get(key)
        r = refseq.get(key)

        if not e or not r:

            print(
                "[WARN] Safety-net annotation missing:",
                m.get("Sample"),
                key,
                file=sys.stderr
            )

            continue

        af = m.get(
            "Alt_AF",
            "."
        )

        idp = m.get(
            "Informative_Depth",
            "."
        )

        row = {
            "Sample":m.get("Sample","."),
            "Gene":(
                r["Gene"]
                if r["Gene"] != "."
                else e["Gene"]
            ),
            "RefSeq_NM":r["RefSeq_NM"],
            "RefSeq_NP":r["RefSeq_NP"],
            "Protein_Change":r["Protein_Change"],
            "Variant(g.)":(
                f'{m["CHROM"]}:{m["POS"]} '
                f'{m["REF"]}>{m["ALT"]}'
            ),
            "Genotype":".",
            "Zygosity":(
                "heterozygote"
                if (
                    af not in ("",".")
                    and float(af) <= 0.70
                )
                else "homozygote"
            ),
            "Read_Depth":idp,
            "Allele_Frequency":af,
            "Allele_Depth":".",
            "Genotype_Quality":".",
            "Variant_Quality":".",
            "MANE_SELECT":e["MANE_SELECT"],
            "MANE_PLUS_CLINICAL":e["MANE_PLUS_CLINICAL"],
            "MANE_HGVSc":e["MANE_HGVSc"],
            "MANE_HGVSp":e["MANE_HGVSp"],
            "CANONICAL":e["CANONICAL"],
            "ClinVar_Significance":m.get(
                "ClinVar_Significance",
                "."
            ),
            "ClinVar_ReviewStatus":m.get(
                "ClinVar_ReviewStatus",
                "."
            ),
            "ClinVar_Disease":m.get(
                "ClinVar_Disease",
                "."
            ),
            "GeneInfo":m.get(
                "GeneInfo",
                "."
            ),
            "rsID":first_rsid(
                e["Existing_variation"]
            ),
            "RefSeq_HGVSc":r["RefSeq_HGVSc"],
            "RefSeq_HGVSp":r["RefSeq_HGVSp"],
            "Technical_Status":
                "SUPPORTED_CALLER_MISSED_CLINVAR_PLP_SNV",
        }

        rescues.append(row)


# Dedicated audit table
if rescues:

    out = (
        report_dir
        / "safetynet_rescued_annotated.tsv"
    )

    fields = list(rescues[0].keys())

    with out.open(
        "w",
        newline=""
    ) as oh:

        w = csv.DictWriter(
            oh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()
        w.writerows(rescues)


# Add/update per-sample final.rsid.tsv
for rescue in rescues:

    sample = rescue["Sample"]

    path = (
        report_dir
        / f"{sample}.final.rsid.tsv"
    )

    if not path.exists():
        continue

    with path.open(
        newline=""
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        fields = reader.fieldnames
        rows = list(reader)

    found = False

    for row in rows:

        if (
            row.get("Variant(g.)")
            == rescue["Variant(g.)"]
        ):

            for col in fields:

                if col in rescue:
                    row[col] = rescue[col]

            found = True
            break

    if not found:

        rows.append({
            col:rescue.get(
                col,
                "."
            )
            for col in fields
        })

    with path.open(
        "w",
        newline=""
    ) as oh:

        w = csv.DictWriter(
            oh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()
        w.writerows(rows)


print(
    f"[INFO] Annotated safety-net rows merged: "
    f"{len(rescues)}"
)

PYMERGERESCUE


  ###########################################################################
  # Rebuild per-sample Excel copies after rescue merge
  ###########################################################################

  for s in $SAMPLES; do

    awk -F'\t' '
      BEGIN{OFS="\t"}
      NR==1{
        print
        next
      }
      {
        $6="\""$6"\""
        $7="\""$7"\""
        print
      }
    ' \
      "$RUN_DIR/report_fix2/${s}.final.rsid.tsv" \
      > "$RUN_DIR/report_fix2/${s}.final.rsid.excel.tsv"

  done


  ###########################################################################
  # Rebuild total_final after safety-net rescue
  ###########################################################################

  head -n 1 \
    "$RUN_DIR/report_fix2/${FIRST_SAMPLE}.final.rsid.tsv" \
    > "$RUN_DIR/report_fix2/total_final.rsid.tsv"

  for s in $SAMPLES; do

    tail -n +2 \
      "$RUN_DIR/report_fix2/${s}.final.rsid.tsv" \
      >> "$RUN_DIR/report_fix2/total_final.rsid.tsv"

  done


  awk -F'\t' '
    BEGIN{OFS="\t"}
    NR==1{
      print
      next
    }
    {
      $6="\""$6"\""
      $7="\""$7"\""
      print
    }
  ' \
    "$RUN_DIR/report_fix2/total_final.rsid.tsv" \
    > "$RUN_DIR/report_fix2/total_final.rsid.excel.tsv"

fi


###############################################################################
# UNIFIED COMPLETE REPORT
###############################################################################

python3 - "$RUN_DIR" <<'PYUNIFIED'

import csv
import re
import sys
from pathlib import Path


run_dir = Path(sys.argv[1])

report_dir = (
    run_dir
    / "report_fix2"
)

primary = (
    report_dir
    / "total_all_technical.tsv"
)

missed = (
    report_dir
    / "clinvar_plp_caller_missed_snv.tsv"
)

annotated_rescue = (
    report_dir
    / "safetynet_rescued_annotated.tsv"
)

out = (
    report_dir
    / "total_complete_report.tsv"
)


base_cols = [
    "Sample",
    "Gene",
    "RefSeq_NM",
    "RefSeq_NP",
    "Protein_Change",
    "Variant(g.)",
    "Genotype",
    "Zygosity",
    "Read_Depth",
    "Allele_Frequency",
    "Allele_Depth",
    "Genotype_Quality",
    "Variant_Quality",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "MANE_HGVSc",
    "MANE_HGVSp",
    "CANONICAL",
    "ClinVar_Significance",
    "ClinVar_ReviewStatus",
    "ClinVar_Disease",
    "GeneInfo",
    "rsID",
    "RefSeq_HGVSc",
    "RefSeq_HGVSp",
    "Technical_Status"
]


extra_cols = [
    "Record_Source",
    "Report_Category",
    "Routine_Reportable",
    "Research_VUS_Flag",
    "Caller_Emitted",
    "ClinVar_ID",
    "Mpileup_Depth",
    "Informative_Depth",
    "Alt_Reads",
    "Pileup_AF"
]


cols = base_cols + extra_cols


def vus(sig):

    return (
        "uncertain"
        in (sig or "").lower()
    )


# pileup metadata indexed by sample+variant
pileup = {}


if missed.exists():

    with missed.open(
        newline=""
    ) as fh:

        for m in csv.DictReader(
            fh,
            delimiter="\t"
        ):

            var = (
                f'{m.get("CHROM",".")}:'
                f'{m.get("POS",".")} '
                f'{m.get("REF",".")}>'
                f'{m.get("ALT",".")}'
            )

            pileup[
                (
                    m.get("Sample","."),
                    var
                )
            ] = m


with out.open(
    "w",
    newline=""
) as oh:

    w = csv.DictWriter(
        oh,
        fieldnames=cols,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore"
    )

    w.writeheader()


    # --------------------------------------------------------
    # Clair3 caller-emitted variants
    # --------------------------------------------------------

    if primary.exists():

        with primary.open(
            newline=""
        ) as fh:

            for r in csv.DictReader(
                fh,
                delimiter="\t"
            ):

                status = (
                    r.get(
                        "Technical_Status"
                    )
                    or ""
                ).strip()

                is_fail = status.startswith(
                    "FAIL_"
                )

                is_vus = (
                    vus(
                        r.get(
                            "ClinVar_Significance",
                            ""
                        )
                    )
                    and not is_fail
                )

                if is_fail:

                    category = "FAIL_QC"
                    routine = "NO"

                elif is_vus:

                    category = "RESEARCH_VUS"
                    routine = "NO"

                else:

                    category = "ROUTINE_SUPPORTED"
                    routine = "YES"

                r.update({
                    "Record_Source":"CLAIR3",
                    "Report_Category":category,
                    "Routine_Reportable":routine,
                    "Research_VUS_Flag":(
                        "YES"
                        if is_vus
                        else "NO"
                    ),
                    "Caller_Emitted":"YES",
                    "ClinVar_ID":".",
                    "Mpileup_Depth":".",
                    "Informative_Depth":".",
                    "Alt_Reads":".",
                    "Pileup_AF":".",
                })

                w.writerow({
                    c:r.get(c,".")
                    for c in cols
                })


    # --------------------------------------------------------
    # Fully annotated caller-missed rescue variants
    # --------------------------------------------------------

    if annotated_rescue.exists():

        with annotated_rescue.open(
            newline=""
        ) as fh:

            for r in csv.DictReader(
                fh,
                delimiter="\t"
            ):

                key = (
                    r.get("Sample","."),
                    r.get("Variant(g.)",".")
                )

                m = pileup.get(
                    key,
                    {}
                )

                r.update({
                    "Record_Source":
                        "BAM_PILEUP_SAFETY_NET",

                    "Report_Category":
                        "ROUTINE_SUPPORTED",

                    "Routine_Reportable":
                        "YES",

                    "Research_VUS_Flag":
                        "NO",

                    "Caller_Emitted":
                        "NO",

                    "ClinVar_ID":
                        m.get(
                            "ClinVar_ID",
                            "."
                        ),

                    "Mpileup_Depth":
                        m.get(
                            "Mpileup_Depth",
                            "."
                        ),

                    "Informative_Depth":
                        m.get(
                            "Informative_Depth",
                            "."
                        ),

                    "Alt_Reads":
                        m.get(
                            "Alt_Reads",
                            "."
                        ),

                    "Pileup_AF":
                        m.get(
                            "Alt_AF",
                            "."
                        ),
                })

                w.writerow({
                    c:r.get(c,".")
                    for c in cols
                })

PYUNIFIED


# Excel-safe copy of the unified report.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {$6="\""$6"\""; $7="\""$7"\""; print}' \
  "$RUN_DIR/report_fix2/total_complete_report.tsv" > "$RUN_DIR/report_fix2/total_complete_report.excel.tsv"

log "DONE"
log "Input mode: $INPUT_MODE"
log "Target BED used: $TARGETS_BED"
log "Clair3 model: $CLAIR3_MODEL"
log "Per-sample final TSV: report_fix2/<sample>.final.rsid.tsv"
log "Per-sample Excel-safe TSV: report_fix2/<sample>.final.rsid.excel.tsv"
log "Merged final TSV: report_fix2/total_final.rsid.tsv"
log "Merged Excel-safe TSV: report_fix2/total_final.rsid.excel.tsv"
log "Technical audit TSV: report_fix2/total_all_technical.tsv"
log "Unfiltered caller-emitted TSV: report_fix2/total_unfiltered.rsid.tsv"
log "Unfiltered Excel-safe TSV: report_fix2/total_unfiltered.rsid.excel.tsv"
log "Research VUS TSV: report_fix2/research_vus.tsv"
log "Caller-missed ClinVar P/LP SNV safety net: report_fix2/clinvar_plp_caller_missed_snv.tsv"
log "Unified complete report: report_fix2/total_complete_report.tsv"
log "Unified Excel-safe report: report_fix2/total_complete_report.excel.tsv"

###############################################################################
# VALIDATION LOCI — BAM PILEUP QC
###############################################################################

VALIDATION_OUT="$RUN_DIR/report_fix2/validation_pileup.tsv"

if [[ -s "$VALIDATION_LOCI" ]]; then

  log "Running known-validation-locus BAM pileup QC"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "Sample" \
    "Gene" \
    "Chrom" \
    "Pos" \
    "Ref" \
    "Alt" \
    "Expected_Zygosity" \
    "Expected_HGVS" \
    "Depth" \
    "Ref_Reads" \
    "Alt_Reads" \
    "Alt_AF" \
    "Pileup_Status" \
    > "$VALIDATION_OUT"

  tail -n +2 "$VALIDATION_LOCI" |
  while IFS=$'\t' read -r sample gene chrom pos ref alt expected_z expected_hgvs; do

    found=0
    for s in $SAMPLES; do
      if [[ "$s" == "$sample" ]]; then
        found=1
        break
      fi
    done

    [[ "$found" -eq 1 ]] || continue

    bam="$RUN_DIR/bam/${sample}.sorted.bam"

    if [[ ! -s "$bam" ]]; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\tNO_BAM\n' \
        "$sample" "$gene" "$chrom" "$pos" "$ref" "$alt" \
        "$expected_z" "$expected_hgvs" >> "$VALIDATION_OUT"
      continue
    fi

    line="$(
      samtools mpileup \
        -aa \
        -B \
        -Q 0 \
        -q 0 \
        -f "$REF_FA" \
        -r "${chrom}:${pos}-${pos}" \
        "$bam" 2>/dev/null |
      head -n1
    )"

    if [[ -z "$line" ]]; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\tNO_COVERAGE\n' \
        "$sample" "$gene" "$chrom" "$pos" "$ref" "$alt" \
        "$expected_z" "$expected_hgvs" >> "$VALIDATION_OUT"
      continue
    fi

    depth="$(printf '%s\n' "$line" | cut -f4)"
    bases="$(printf '%s\n' "$line" | cut -f5)"

    result="$(
      python3 - "$bases" "$alt" "$depth" <<'PY'
import re
import sys

bases = sys.argv[1]
alt = sys.argv[2].upper()
reported_depth = int(sys.argv[3])

bases = re.sub(r'\^.', '', bases)
bases = bases.replace('$', '')

ref_count = 0
alt_count = 0

i = 0
while i < len(bases):
    c = bases[i]

    if c in '.,':
        ref_count += 1
        i += 1
        continue

    if c.upper() == alt:
        alt_count += 1
        i += 1
        continue

    if c in '+-':
        j = i + 1
        num = ''

        while j < len(bases) and bases[j].isdigit():
            num += bases[j]
            j += 1

        if num:
            n = int(num)
            i = j + n
        else:
            i += 1
        continue

    i += 1

informative = ref_count + alt_count
af = alt_count / informative if informative else 0.0

print(f"{reported_depth}\t{ref_count}\t{alt_count}\t{af:.4f}")
PY
    )"

    IFS=$'\t' read -r depth2 ref_reads alt_reads alt_af <<< "$result"

    status="$(
      python3 - "$depth2" "$alt_af" "$MIN_DP" "$MIN_AF" "$HET_MAX_AF" <<'PY'
import sys

dp = int(sys.argv[1])
af = float(sys.argv[2])
min_dp = int(float(sys.argv[3]))
min_af = float(sys.argv[4])
het_max = float(sys.argv[5])

if dp < min_dp:
    print("LOW_DEPTH")
elif af < min_af:
    print("NO_ALT_SUPPORT")
elif af <= het_max:
    print("ALT_SUPPORTED_HET_RANGE")
else:
    print("ALT_SUPPORTED_HOM_RANGE")
PY
    )"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$sample" \
      "$gene" \
      "$chrom" \
      "$pos" \
      "$ref" \
      "$alt" \
      "$expected_z" \
      "$expected_hgvs" \
      "$depth2" \
      "$ref_reads" \
      "$alt_reads" \
      "$alt_af" \
      "$status" \
      >> "$VALIDATION_OUT"

  done

  log "Validation pileup report: report_fix2/validation_pileup.tsv"
fi

###############################################################################
# VALIDATION CONCORDANCE CLASSIFICATION
###############################################################################

VALIDATION_PILEUP="$RUN_DIR/report_fix2/validation_pileup.tsv"
VALIDATION_CLASSIFIED="$RUN_DIR/report_fix2/validation_concordance.tsv"

if [[ -s "$VALIDATION_PILEUP" ]]; then

  log "Classifying validation concordance against final caller reports"

  python3 - \
    "$VALIDATION_PILEUP" \
    "$RUN_DIR/report_fix2" \
    "$VALIDATION_CLASSIFIED" \
    "$MIN_DP" \
    "$MIN_AF" \
    "$HET_MAX_AF" <<'PY'

import csv
import sys
from pathlib import Path

pileup_file = Path(sys.argv[1])
report_dir = Path(sys.argv[2])
out_file = Path(sys.argv[3])

min_dp = int(float(sys.argv[4]))
min_af = float(sys.argv[5])
het_max = float(sys.argv[6])


def normalize_zygosity(x):
    x = (x or "").strip().lower()

    if x in {"heterozygote", "heterozygous", "het"}:
        return "heterozygote"

    if x in {"homozygote", "homozygous", "hom"}:
        return "homozygote"

    return x


# ------------------------------------------------------------
# Read final caller reports
# key:
# (sample, chrom, pos, ref, alt)
# ------------------------------------------------------------

caller = {}

for report in sorted(report_dir.glob("barcode*.final.rsid.tsv")):

    with report.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            sample = row.get("Sample", "").strip()
            variant = row.get("Variant(g.)", "").strip()

            # Expected format:
            # chr5:132378438 G>C

            try:
                loc, change = variant.split(None, 1)
                chrom, pos = loc.split(":", 1)
                ref, alt = change.split(">", 1)
            except ValueError:
                continue

            key = (
                sample,
                chrom,
                pos,
                ref,
                alt,
            )

            caller[key] = {
                "zygosity": normalize_zygosity(
                    row.get("Zygosity", "")
                ),
                "dp": row.get("Read_Depth", ""),
                "af": row.get("Allele_Frequency", ""),
            }


# ------------------------------------------------------------
# Compare known loci with caller output
# ------------------------------------------------------------

with pileup_file.open(newline="") as inp, \
     out_file.open("w", newline="") as out:

    reader = csv.DictReader(inp, delimiter="\t")

    fields = list(reader.fieldnames) + [
        "Pileup_Zygosity",
        "Caller_Called",
        "Caller_Zygosity",
        "Caller_DP",
        "Caller_AF",
        "Final_Status",
        "Interpretation",
    ]

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for row in reader:

        sample = row["Sample"]
        chrom = row["Chrom"]
        pos = row["Pos"]
        ref = row["Ref"]
        alt = row["Alt"]

        expected = normalize_zygosity(
            row["Expected_Zygosity"]
        )

        dp = int(float(row["Depth"]))
        af = float(row["Alt_AF"])

        # Infer genotype independently from pileup AF
        if dp < min_dp:
            pileup_z = "not_assessable"
        elif af < min_af:
            pileup_z = "no_alt_support"
        elif af <= het_max:
            pileup_z = "heterozygote"
        else:
            pileup_z = "homozygote"

        key = (
            sample,
            chrom,
            pos,
            ref,
            alt,
        )

        call = caller.get(key)

        if call is None:
            caller_called = "NO"
            caller_z = "."
            caller_dp = "."
            caller_af = "."
        else:
            caller_called = "YES"
            caller_z = call["zygosity"]
            caller_dp = call["dp"]
            caller_af = call["af"]

        # ----------------------------------------------------
        # Final classification
        # ----------------------------------------------------

        if dp < min_dp:

            status = "LOW_COVERAGE_NOT_ASSESSABLE"
            interpretation = (
                "Coverage is below the minimum threshold; "
                "absence or genotype cannot be interpreted reliably."
            )

        elif af < min_af:

            status = "NO_ALT_SUPPORT"
            interpretation = (
                "Coverage is sufficient but alternative-allele support "
                "is below the reporting threshold."
            )

        elif call is None:

            status = "CALLER_MISSED_VARIANT"
            interpretation = (
                "Alternative reads support the expected variant at "
                "sufficient depth and allele fraction, but the variant "
                "was not emitted in the final caller report."
            )

        elif pileup_z != expected:

            status = "DISCORDANT_GENOTYPE"
            interpretation = (
                "The variant is detected, but pileup-based allele "
                "fraction does not match the expected zygosity."
            )

        elif caller_z != expected:

            status = "DISCORDANT_GENOTYPE"
            interpretation = (
                "The variant is called, but caller/report zygosity "
                "does not match the expected zygosity."
            )

        else:

            status = "CONCORDANT"
            interpretation = (
                "Expected variant and zygosity are supported by both "
                "pileup evidence and the final caller report."
            )

        row["Pileup_Zygosity"] = pileup_z
        row["Caller_Called"] = caller_called
        row["Caller_Zygosity"] = caller_z
        row["Caller_DP"] = caller_dp
        row["Caller_AF"] = caller_af
        row["Final_Status"] = status
        row["Interpretation"] = interpretation

        writer.writerow(row)

PY

  log "Validation concordance report: report_fix2/validation_concordance.tsv"

fi

###############################################################################
# VALIDATION CONCORDANCE V2
# Distinguishes:
#   1) raw Clair3 call
#   2) final-report pass
#   3) BAM pileup support
###############################################################################

VALIDATION_V2_OUT="$RUN_DIR/report_fix2/validation_concordance_v2.tsv"

if [[ -s "$RUN_DIR/report_fix2/validation_pileup.tsv" ]]; then

  log "Building validation concordance v2"

  python3 - \
    "$RUN_DIR/report_fix2/validation_pileup.tsv" \
    "$RUN_DIR" \
    "$VALIDATION_V2_OUT" \
    "$MIN_DP" \
    "$MIN_AF" \
    "$HET_MAX_AF" <<'PY'

import csv
import subprocess
import sys
from pathlib import Path

pileup_file = Path(sys.argv[1])
run_dir = Path(sys.argv[2])
out_file = Path(sys.argv[3])

min_dp = int(float(sys.argv[4]))
min_af = float(sys.argv[5])
het_max = float(sys.argv[6])


def norm_z(x):
    x = (x or "").strip().lower()
    if x in {"het", "heterozygous", "heterozygote"}:
        return "heterozygote"
    if x in {"hom", "homozygous", "homozygote"}:
        return "homozygote"
    return x


def gt_to_z(gt):
    gt = gt.replace("|", "/")
    if gt in {"0/1", "1/0"}:
        return "heterozygote"
    if gt == "1/1":
        return "homozygote"
    return "."


def query_raw_vcf(sample, chrom, pos, ref, alt):
    vcf = run_dir / "clair3_out" / sample / "merge_output.vcf.gz"

    if not vcf.exists():
        return None

    cmd = [
        "bcftools", "query",
        "-r", f"{chrom}:{pos}-{pos}",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%DP\t%AF]\n",
        str(vcf),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )

    for line in result.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue

        c, p, r, a, qual, gt, dp, af = f[:8]

        # ALT may theoretically contain multiple alleles
        alts = a.split(",")

        if (
            c == chrom
            and p == str(pos)
            and r == ref
            and alt in alts
        ):
            return {
                "qual": qual,
                "gt": gt,
                "dp": dp,
                "af": af,
                "zygosity": gt_to_z(gt),
            }

    return None


def final_report_call(sample, chrom, pos, ref, alt):
    report = run_dir / "report_fix2" / f"{sample}.final.rsid.tsv"

    if not report.exists():
        return None

    with report.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            variant = row.get("Variant(g.)", "").strip()

            try:
                loc, change = variant.split(None, 1)
                c, p = loc.split(":", 1)
                r, a = change.split(">", 1)
            except ValueError:
                continue

            if (
                c == chrom
                and p == str(pos)
                and r == ref
                and a == alt
            ):
                return row

    return None


with pileup_file.open(newline="") as inp, \
     out_file.open("w", newline="") as out:

    reader = csv.DictReader(inp, delimiter="\t")

    extra = [
        "Raw_Clair3_Called",
        "Raw_Clair3_GT",
        "Raw_Clair3_Zygosity",
        "Raw_Clair3_DP",
        "Raw_Clair3_AF",
        "Raw_Clair3_QUAL",
        "Passed_Final_Filter",
        "Final_Status_V2",
        "Interpretation_V2",
    ]

    writer = csv.DictWriter(
        out,
        fieldnames=list(reader.fieldnames) + extra,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for row in reader:

        sample = row["Sample"]
        chrom = row["Chrom"]
        pos = row["Pos"]
        ref = row["Ref"]
        alt = row["Alt"]

        expected = norm_z(row["Expected_Zygosity"])

        pileup_dp = int(float(row["Depth"]))
        pileup_af = float(row["Alt_AF"])

        raw = query_raw_vcf(
            sample, chrom, pos, ref, alt
        )

        passed = final_report_call(
            sample, chrom, pos, ref, alt
        )

        if raw:
            raw_called = "YES"
            raw_gt = raw["gt"]
            raw_z = raw["zygosity"]
            raw_dp = raw["dp"]
            raw_af = raw["af"]
            raw_qual = raw["qual"]
        else:
            raw_called = "NO"
            raw_gt = "."
            raw_z = "."
            raw_dp = "."
            raw_af = "."
            raw_qual = "."

        passed_final = "YES" if passed else "NO"

        # ------------------------------------------------
        # Classification
        # ------------------------------------------------

        if pileup_dp < min_dp:

            status = "LOW_COVERAGE_NOT_ASSESSABLE"
            interpretation = (
                "Read depth is below the validation threshold; "
                "absence or genotype cannot be assessed reliably."
            )

        elif raw is None and pileup_af >= min_af:

            status = "CALLER_MISSED_VARIANT"
            interpretation = (
                "The BAM contains sufficient alternative-read support, "
                "but Clair3 did not emit the variant in its raw VCF."
            )

        elif raw is not None and passed is None:

            status = "CALLER_CALLED_FILTERED_OUT"
            interpretation = (
                "Clair3 emitted the variant, but the call did not pass "
                "the final reporting thresholds."
            )

        elif raw is not None and passed is not None:

            report_z = norm_z(
                passed.get("Zygosity", "")
            )

            if report_z == expected:
                status = "CONCORDANT"
                interpretation = (
                    "The expected variant was called, passed the final "
                    "filters, and the reported zygosity matches the "
                    "expected result."
                )
            else:
                status = "DISCORDANT_GENOTYPE"
                interpretation = (
                    "The variant was called and passed final filters, "
                    "but reported zygosity differs from the expected result."
                )

        else:

            status = "NO_ALT_SUPPORT"
            interpretation = (
                "No reportable alternative-allele support was observed."
            )

        row["Raw_Clair3_Called"] = raw_called
        row["Raw_Clair3_GT"] = raw_gt
        row["Raw_Clair3_Zygosity"] = raw_z
        row["Raw_Clair3_DP"] = raw_dp
        row["Raw_Clair3_AF"] = raw_af
        row["Raw_Clair3_QUAL"] = raw_qual
        row["Passed_Final_Filter"] = passed_final
        row["Final_Status_V2"] = status
        row["Interpretation_V2"] = interpretation

        writer.writerow(row)

PY

  log "Validation concordance v2: report_fix2/validation_concordance_v2.tsv"

fi
