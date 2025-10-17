#!/usr/bin/env bash
set -euo pipefail

# =======================================================================
# Processing ENCODE ChIP-Seq for TF binding from unfiltered mapping files
# =======================================================================
# Usage:
#   pipeline.sh [-m] [-c] <encode.narrowPeak> <rep1.bam> <rep2.bam> <ctrl1.bam> [ctrl2.bam]
#
# Options:
#   -m    merge ctrl1.bam and ctrl2.bam into one file, meaningful if they represent portions of
#         a single control experiment.
#
#   -c    CRISPR MODE: use 0.01 as FDR threshold for peak calling with MACS2, instead of the
#         default (0.05).

ARE_CTRLS_MERGED=false
if [ "${1:-}" = "-m" ]; then
    ARE_CTRLS_MERGED=true
    shift
fi

FDR_THRESHOLD=0.05
if [ "${1:-}" = "-c" ]; then
    FDR_THRESHOLD=0.01
    shift
fi

if [ $# -lt 5 ]; then
    echo "Usage: $0 [-m] [-c] <encode.narrowPeak> <rep1.bam> <rep2.bam> <ctrl1.bam> [ctrl2.bam]"
    exit 1
fi

ENCODE="$1"
REP1="$2"
REP2="$3"
CTRL1="$4"
CTRL2="${5:-}"

# Merge controls if requested
if $ARE_CTRLS_MERGED && [ -n "$CTRL2" ]; then
    samtools merge -f ctrl1.bam "$CTRL1" "$CTRL2"
    CTRL1="ctrl1.bam"
    echo "Controls have been merged into ctrl1.bam"
fi


# Filtering
BAMS=("$REP1" "$REP2" "$CTRL1")
[ -n "$CTRL2" ] && BAMS+=("$CTRL2")

FILTERED_NAMES=("filtered1" "filtered2" "filtered_ctrl1")
[ -n "$CTRL2" ] && FILTERED_NAMES+=("filtered_ctrl2")

echo "Filtering mapping files..."
for i in "${!BAMS[@]}"; do
    samtools view -bq 1 "$BAMS[$i]" > "${FILTERED_NAMES[$i]}.bam"
done


# Mapping QC
for i in "${!BAMS[@]}"; do
    echo "$BAMS[$i] mapping statistics:"
    samtools flagstat "$BAMS[$i]"

    READS_AMT=$(samtools view -c "$BAMS[$i]")
    echo "    Total reads: $READS_AMT"

    MAPPING_AMT=$(samtools view -c -F 4 "$BAMS[$i]")
    echo "    Mapping reads: $MAPPING_AMT"

    UMAPPING_AMT=$(samtools view -c "${FILTERED_NAMES[$i]}.bam")
    echo "    Uniquely mapping reads: $UMAPPING_AMT"

    MMAPPING_AMT=$((MAPPING_AMT - UMAPPING_AMT))
    echo "    Multi-mapping reads (mapping but not uniquely): $MMAPPING_AMT"
done


# Peak calling
echo "Calling peaks..."
if [ -n "$CTRL2" ]; then
    macs2 callpeak -t f1.bam -c fc1.bam -q "$FDR_THRESHOLD" -g hs -n REP1
    macs2 callpeak -t f2.bam -c fc2.bam -q "$FDR_THRESHOLD" -g hs -n REP2
    macs2 callpeak -t f1.bam f2.bam -c fc1.bam fc2.bam -q "$FDR_THRESHOLD" -g hs -n MERGE
else
    macs2 callpeak -t f1.bam -c fc1.bam -q "$FDR_THRESHOLD" -g hs -n REP1
    macs2 callpeak -t f2.bam -c fc1.bam -q "$FDR_THRESHOLD" -g hs -n REP2
    macs2 callpeak -t f1.bam f2.bam -c fc1.bam -q "$FDR_THRESHOLD" -g hs -n MERGE
fi


# R plots
R < REP1_model.r --vanilla
R < REP2_model.r --vanilla
R < MERGE_model.r --vanilla
echo "R plots are available as pdfs"


# Computing overlaps
rep1_count=$(wc -l < REP1_peaks.narrowPeak)
rep2_count=$(wc -l < REP2_peaks.narrowPeak)
merge_count=$(wc -l < MERGE_peaks.narrowPeak)

if [ "$rep1_count" -le "$rep2_count" ]; then
    REP_COMP_NAME="REP1_in_REP2"
    bedtools intersect -a REP1_peaks.narrowPeak -b REP2_peaks.narrowPeak -u > "${REP_COMP_NAME}.narrowPeak"
    
    summitProx=$(bedtools closest -a REP1_summits.bed -b REP2_summits.bed -d | awk '$NF <= 100' | wc -l)
    overlaps=$(wc -l < "${REP_COMP_NAME}.narrowPeak")
    fraction=$(awk "BEGIN {printf \"%.3f\", $overlaps/$rep1_count}")

    echo "REP1 (${rep1_count} peaks) vs REP2 (${rep2_count} peaks):"
else
    REP_COMP_NAME="REP2_in_REP1"
    bedtools intersect -a REP2_peaks.narrowPeak -b REP1_peaks.narrowPeak -u > "${REP_COMP_NAME}.narrowPeak"
    
    summitProx=$(bedtools closest -a REP2_summits.bed -b REP1_summits.bed -d | awk '$NF <= 100' | wc -l)
    overlaps=$(wc -l < "${REP_COMP_NAME}.narrowPeak")
    fraction=$(awk "BEGIN {printf \"%.3f\", $overlaps/$rep2_count}")

    echo "REP2 (${rep2_count} peaks) vs REP1 (${rep1_count} peaks):"
fi
echo "    Summit proximity: $summitProx"
echo "    Overlapping peaks: $overlaps"
echo "    Fraction of overlapping peaks: $fraction"

bedtools intersect -a REP1_peaks.narrowPeak -b MERGE_peaks.narrowPeak -u > REP1_in_MERGE.narrowPeak    
summitProx=$(bedtools closest -a REP1_summits.bed -b MERGE_summits.bed -d | awk '$NF <= 100' | wc -l)
overlaps=$(wc -l < REP1_in_MERGE.narrowPeak)
fraction=$(awk "BEGIN {printf \"%.3f\", $overlaps/$rep1_count}")
echo "REP1 (${rep1_count} peaks) vs MERGE (${merge_count} peaks):"
echo "    Summit proximity: $summitProx"
echo "    Overlapping peaks: $overlaps"
echo "    Fraction of overlapping peaks: $fraction"

bedtools intersect -a REP2_peaks.narrowPeak -b MERGE_peaks.narrowPeak -u > REP2_in_MERGE.narrowPeak    
summitProx=$(bedtools closest -a REP2_summits.bed -b MERGE_summits.bed -d | awk '$NF <= 100' | wc -l)
overlaps=$(wc -l < REP2_in_MERGE.narrowPeak)
fraction=$(awk "BEGIN {printf \"%.3f\", $overlaps/$rep2_count}")
echo "REP2 (${rep2_count} peaks) vs MERGE (${merge_count} peaks):"
echo "    Summit proximity: $summitProx"
echo "    Overlapping peaks: $overlaps"
echo "    Fraction of overlapping peaks: $fraction"


# Comparison with the ENCODE results
bedtools sort -i "$ENCODE" > ENCODE_peaks.narrowPeak

NAMES=("REP1" "REP2" "MERGE" $REP_COMP_NAME "REP1_in_MERGE" "REP2_in_MERGE")
for NAME in "${NAMES[@]}"; do
    echo "Jaccard index for $NAME vs $ENCODE:"
    bedtools jaccard -a "${NAME}_peaks.narrowPeak" -b ENCODE_peaks.narrowPeak
done


echo "All done!"