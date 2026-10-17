#!/bin/bash

# ==========================================
# 1. DOSYA YOLLARI
# ==========================================
BASE_DIR="/mnt/c/Users/corha/Desktop/faod_sonuclar_final"
TARGET_BED="$BASE_DIR/targets.bed"
OUT_CSV="$BASE_DIR/Concordance_Metrics.csv"

# Senin belirttiğin doğru klasör yolları
PCR_BASE="/mnt/c/Users/corha/Desktop/PCR_Probeta_Calisilan/clair3_out"
PROBE_BASE="/mnt/c/Users/corha/Desktop/celemic_new_demux/demux_trimmed_merged/clair3_out"

# ==========================================
# 2. SENİN VERDİĞİN DOĞRU EŞLEŞTİRME LİSTESİ
# ==========================================
declare -A PAIRS=(
    ["ONT_UDI_57"]="barcode01"
    ["ONT_UDI_58"]="barcode09"
    ["ONT_UDI_59"]="barcode15"
    ["ONT_UDI_60"]="barcode18"
    ["ONT_UDI_61"]="barcode43"
    ["ONT_UDI_62"]="barcode31"
    ["ONT_UDI_63"]="barcode33"
    ["ONT_UDI_64"]="barcode36"
    ["ONT_UDI_65"]="barcode37"
    ["ONT_UDI_66"]="barcode23"
    ["ONT_UDI_67"]="barcode29"
    ["ONT_UDI_68"]="barcode30"
)

# ==========================================
# 3. ANALİZ
# ==========================================
echo "Sample_Probe,Sample_PCR,Total_PCR(Truth),Total_Probe(Query),TP(Concordant),FN(PCR_Only),FP(Probe_Only),Sensitivity_%,Precision_%,F1_Score" > "$OUT_CSV"

for probe_id in "${!PAIRS[@]}"; do
    pcr_id="${PAIRS[$probe_id]}"
    
    # Senin belirttiğin dosya yolları
    pcr_vcf="$PCR_BASE/$pcr_id/merge_output.vcf.gz"
    probe_vcf="$PROBE_BASE/$probe_id/merge_output.vcf.gz"

    if [[ -f "$pcr_vcf" && -f "$probe_vcf" ]]; then
        echo "İşleniyor: Probe:$probe_id <-> PCR:$pcr_id"
        
        TMP_ISEC="$BASE_DIR/tmp_$probe_id"
        mkdir -p "$TMP_ISEC"
        
        # Kesistirme (PCR 'Truth' olarak baz alındı)
        bcftools isec -p "$TMP_ISEC" -T "$TARGET_BED" "$pcr_vcf" "$probe_vcf" 2>/dev/null
        
        FN=$(bcftools view -H "$TMP_ISEC/0000.vcf" 2>/dev/null | wc -l)
        FP=$(bcftools view -H "$TMP_ISEC/0001.vcf" 2>/dev/null | wc -l)
        TP=$(bcftools view -H "$TMP_ISEC/0002.vcf" 2>/dev/null | wc -l)
        
        # Metrikler
        SENS=$(echo "scale=2; if(($TP+$FN)>0) ($TP/($TP+$FN))*100 else 0" | bc)
        PREC=$(echo "scale=2; if(($TP+$FP)>0) ($TP/($TP+$FP))*100 else 0" | bc)
        F1=$(echo "scale=4; if(($TP+$FP+$FN)>0) (2*$TP)/(2*$TP+$FP+$FN) else 0" | bc)
        
        echo "$probe_id,$pcr_id,$((TP+FN)),$((TP+FP)),$TP,$FN,$FP,$SENS,$PREC,$F1" >> "$OUT_CSV"
        rm -rf "$TMP_ISEC"
    else
        echo "HATA: Dosya bulunamadı -> Probe:$probe_id veya PCR:$pcr_id"
    fi
done

echo "Analiz tamamlandı. Sonuç tablosu: $OUT_CSV"
