#!/bin/bash

# ==========================================
# DİZİN VE YOL ATAMALARI
# ==========================================
BASE_DIR="/mnt/c/Users/corha/Desktop/faod_sonuclar_final"
cd "$BASE_DIR" || exit

FASTQ_DIR="/mnt/c/Users/corha/Desktop/celemic_new_demux/demux_trimmed_merged"
BAM_DIR="$FASTQ_DIR/bam"
MOSDEPTH_DIR="$BASE_DIR/probe_coverage"
TARGET_BED="/mnt/c/Users/corha/Desktop/celemic_new_demux/demux_trimmed_merged/hedef_paneliniz.bed" 
OUT_CSV="$BASE_DIR/ONT_Panel_QC_Sonuclari.csv"

echo "Sample,Output_Mb,Total_Reads,Mean_Length,Q20_%,Q30_%,Alignment_Rate_%,On_Target_Rate_%,Mean_Depth,Median_Depth,Cov_>=30x_%,Cov_>=100x_%" > "$OUT_CSV"

echo "Esnek arama ile analiz başlatılıyor..."

for fq in "$FASTQ_DIR"/ONT_UDI_*.fastq; do
    
    [ -e "$fq" ] || continue
    SAMPLE_ID=$(basename "$fq" .fastq)
    echo "İşleniyor: $SAMPLE_ID..."

    # 1. FASTQ METRİKLERİ (SeqKit kolon isimlerinden dinamik okuma ile Q20/Q30 hatası giderildi)
    SEQ_STATS=$(seqkit stats -T -a "$fq")
    TOTAL_READS=$(echo "$SEQ_STATS" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="num_seqs") c=i} NR>1{print $c}')
    YIELD_MB=$(echo "$SEQ_STATS" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sum_len") c=i} NR>1{printf "%.2f", $c/1000000}')
    MEAN_LENGTH=$(echo "$SEQ_STATS" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="avg_len") c=i} NR>1{print $c}')
    Q20_PCT=$(echo "$SEQ_STATS" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Q20(%)") c=i} NR>1{print $c}')
    Q30_PCT=$(echo "$SEQ_STATS" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Q30(%)") c=i} NR>1{print $c}')

    # 2. ALIGNMENT METRİKLERİ (Esnek dosya arama)
    # İçinde sample ID geçen ilk .bam uzantılı dosyayı yakalar (Örn: ONT_UDI_57.sorted.bam)
    BAM_FILE=$(find "$BAM_DIR" -maxdepth 1 -type f -name "${SAMPLE_ID}*.bam" | head -n 1)
    
    if [ -n "$BAM_FILE" ] && [ -f "$BAM_FILE" ]; then
        TOTAL_MAPPED_PRIMARY=$(samtools view -c -F 2308 "$BAM_FILE")
        TOTAL_INPUT_PRIMARY=$(samtools view -c -F 2304 "$BAM_FILE")
        
        if [ "$TOTAL_INPUT_PRIMARY" -gt 0 ]; then
            ALIGNMENT_RATE=$(echo "scale=2; ($TOTAL_MAPPED_PRIMARY / $TOTAL_INPUT_PRIMARY) * 100" | bc)
        else
            ALIGNMENT_RATE="0"
        fi
        
        ON_TARGET_PRIMARY=$(samtools view -c -F 2308 -L "$TARGET_BED" "$BAM_FILE")
        if [ "$TOTAL_MAPPED_PRIMARY" -gt 0 ]; then
            ON_TARGET_RATE=$(echo "scale=2; ($ON_TARGET_PRIMARY / $TOTAL_MAPPED_PRIMARY) * 100" | bc)
        else
            ON_TARGET_RATE="0"
        fi
    else
        ALIGNMENT_RATE="BAM_YOK"
        ON_TARGET_RATE="BAM_YOK"
    fi

    # 3. DERİNLİK METRİKLERİ (Esnek dosya arama)
    SUMMARY_FILE=$(find "$MOSDEPTH_DIR" -maxdepth 1 -type f -name "${SAMPLE_ID}*.mosdepth.summary.txt" | head -n 1)

    if [ -n "$SUMMARY_FILE" ] && [ -f "$SUMMARY_FILE" ]; then
        # Summary dosyasının ön ekini alıp diğer dosyaları o isme göre buluruz
        MOSDEPTH_PREFIX="${SUMMARY_FILE%.mosdepth.summary.txt}"
        DIST_FILE="${MOSDEPTH_PREFIX}.mosdepth.global.dist.txt"
        THRESHOLDS_FILE="${MOSDEPTH_PREFIX}.thresholds.bed.gz"

        MEAN_DEPTH=$(awk '$1 == "total" {print $4}' "$SUMMARY_FILE")
        TARGET_SIZE=$(awk '{sum+=($3-$2)} END {print sum}' "$TARGET_BED")
        
        BASES_30X=$(zcat "$THRESHOLDS_FILE" | awk '{sum+=$5} END {print sum}')
        BASES_100X=$(zcat "$THRESHOLDS_FILE" | awk '{sum+=$6} END {print sum}')
        
        COV_30X_PCT=$(echo "scale=2; ($BASES_30X / $TARGET_SIZE) * 100" | bc)
        COV_100X_PCT=$(echo "scale=2; ($BASES_100X / $TARGET_SIZE) * 100" | bc)
        
        MEDIAN_DEPTH=$(python3 -c "
import pandas as pd
try:
    df = pd.read_csv('$DIST_FILE', sep='\t', names=['chrom', 'depth', 'fraction'])
    df_target = df[df['chrom'] == 'total']
    median = df_target[df_target['fraction'] <= 0.5]['depth'].min()
    print(int(median))
except:
    print('Hata')
")
    else:
        MEAN_DEPTH="MOSDEPTH_YOK"
        MEDIAN_DEPTH="MOSDEPTH_YOK"
        COV_30X_PCT="MOSDEPTH_YOK"
        COV_100X_PCT="MOSDEPTH_YOK"
    fi

    # 4. SONUÇ YAZDIRMA
    echo "$SAMPLE_ID,$YIELD_MB,$TOTAL_READS,$MEAN_LENGTH,$Q20_PCT,$Q30_PCT,$ALIGNMENT_RATE,$ON_TARGET_RATE,$MEAN_DEPTH,$MEDIAN_DEPTH,$COV_30X_PCT,$COV_100X_PCT" >> "$OUT_CSV"

done

echo "Tüm metrikler başarıyla çekildi. Dosya hazır: $OUT_CSV"
