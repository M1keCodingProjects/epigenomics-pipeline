#!/usr/bin/env bash
set -euo pipefail

# =======================================================================
# Processing ENCODE ChIP-Seq for TF binding from unfiltered mapping files
# =======================================================================
# Usage:
#   pipeline.sh [-h] [-m] [-c] [-o output.txt] <encode.narrowPeak> <blacklist.bed> <rep1.bam> <rep2.bam> <ctrl1.bam> [ctrl2.bam]
#
# Options:
#   -h    print this help message and exit.
#
#   -m    merge ctrl1.bam and ctrl2.bam into one file, meaningful if they represent portions of
#         a single control experiment.
#
#   -c    CRISPR MODE: use 0.01 as FDR threshold for peak calling with MACS2, instead of the
#         default (0.05).
#
#   -o    use the provided name for the output file containing the quality control parameters
#         of the analysis, instead of the default (qc_params.txt)

HELP_MSG="Usage: $0 [-h] [-m] [-c] [-o output.txt] <encode.narrowPeak> <blacklist.bed> <rep1.bam> <rep2.bam> <ctrl1.bam> [ctrl2.bam]"

ARE_CTRLS_MERGED=false
FDR_THRESHOLD=0.05
OUTPUT="./qc_params.txt"
while getopts "hmco:" opt; do
    case ${opt} in
        h )
            echo $HELP_MSG
            echo
            echo " Options:"
            echo "   -h    print this help message and exit."
            echo
            echo "   -m    merge ctrl1.bam and ctrl2.bam into one file, meaningful if they represent portions of"
            echo "         a single control experiment."
            echo
            echo "   -c    CRISPR MODE: use 0.01 as FDR threshold for peak calling with MACS2, instead of the"
            echo "         default (0.05)."
            echo
            echo "   -o    use the provided name for the output file containing the quality control parameters"
            echo "         of the analysis, instead of the default (qc_params.txt)"
            exit 1
            ;;
        m )
            ARE_CTRLS_MERGED=true
            ;;
        c )
            FDR_THRESHOLD=0.01
            ;;
        o )
            OUTPUT="$OPTARG"
            ;;
        \? )
            echo $HELP_MSG
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

if [ $# -lt 5 ]; then
    echo $HELP_MSG
    exit 1
fi

ENCODE="$1"
BLACKLIST="$2"
REP1="$3"
REP2="$4"
CTRL1="$5"
CTRL2="${6:-}"

> "$OUTPUT"
echo "Full analysis output available in $OUTPUT."

# Merge controls if requested
if "$ARE_CTRLS_MERGED" && [ -n "$CTRL2" ]; then
    echo
    echo "Merging control files..."
    
    samtools merge -f -@ 8 ctrl.bam "$CTRL1" "$CTRL2"
    
    CTRL1="ctrl.bam"
    unset CTRL2
    
    echo "Controls have been merged into ctrl.bam." | tee -a "$OUTPUT"
else
    echo
    echo "Controls have not been merged." | tee -a "$OUTPUT"
fi

# Filtering
BAMS=("$REP1" "$REP2" "$CTRL1")
[ -n "${CTRL2:-}" ] && BAMS+=("$CTRL2")

FILTERED_NAMES=("filtered1" "filtered2" "filtered_ctrl1")
[ -n "${CTRL2:-}" ] && FILTERED_NAMES+=("filtered_ctrl2")

echo
echo "Filtering mapping files..."
for i in "${!BAMS[@]}"; do
    samtools view -bq 1 -@ 8 "${BAMS[$i]}" > "${FILTERED_NAMES[$i]}.bam"
done


# Mapping QC
for i in "${!BAMS[@]}"; do
    echo | tee -a "$OUTPUT"
    echo "${BAMS[$i]} mapping statistics:" | tee -a "$OUTPUT"

    READS_AMT=$(samtools view -c -@ 8 "${BAMS[$i]}")
    echo "  Total reads: $READS_AMT" | tee -a "$OUTPUT"

    MAPPING_AMT=$(samtools view -c -F 4 -@ 8 "${BAMS[$i]}")
    echo "  Mapping reads: $MAPPING_AMT" | tee -a "$OUTPUT"

    UMAPPING_AMT=$(samtools view -c -@ 8 "${FILTERED_NAMES[$i]}.bam")
    echo "  Uniquely mapping reads: $UMAPPING_AMT" | tee -a "$OUTPUT"

    MMAPPING_AMT=$((MAPPING_AMT - UMAPPING_AMT))
    echo "  Multi-mapping reads (mapping but not uniquely): $MMAPPING_AMT" | tee -a "$OUTPUT"

    MAPPING_F=$(awk "BEGIN {printf \"%.3f\", $MAPPING_AMT/$READS_AMT*100}")
    echo "  Mapping reads as a percentage of total reads: ${MAPPING_F}%" | tee -a "$OUTPUT"

    MMAPPING_F=$(awk "BEGIN {printf \"%.3f\", $MMAPPING_AMT/$MAPPING_AMT*100}")
    echo "  Multi-mapping reads as a percentage of mapping reads: ${MMAPPING_F}%" | tee -a "$OUTPUT"

    echo
    samtools flagstat -@ 8 "${BAMS[$i]}"
    echo
done


# Peak calling
echo "Calling peaks..."
if [ -n "${CTRL2:-}" ]; then
    echo "Calling peaks on ${REP1}..."
    macs2 callpeak -t filtered1.bam -c filtered_ctrl1.bam -q "$FDR_THRESHOLD" -g hs -B --SPMR -n REP1 > macs2_REP1_log.txt 2>&1 &

    echo "Calling peaks on ${REP2}..."
    macs2 callpeak -t filtered2.bam -c filtered_ctrl2.bam -q "$FDR_THRESHOLD" -g hs -B --SPMR -n REP2 > macs2_REP2_log.txt 2>&1 &
    
    echo "Calling peaks on the 2 replicates merged..."
    macs2 callpeak -t filtered1.bam filtered2.bam -c filtered_ctrl1.bam filtered_ctrl2.bam -q "$FDR_THRESHOLD" -g hs -B --SPMR -n MERGE > macs2_MERGE_log.txt 2>&1 &
else
    echo "Calling peaks on ${REP1}..."
    macs2 callpeak -t filtered1.bam -c filtered_ctrl1.bam -q "$FDR_THRESHOLD" -g hs -B --SPMR -n REP1 > macs2_REP1_log.txt 2>&1 &
    
    echo "Calling peaks on ${REP2}..."
    macs2 callpeak -t filtered2.bam -c filtered_ctrl1.bam -q "$FDR_THRESHOLD" -g hs -B --SPMR -n REP2 > macs2_REP2_log.txt 2>&1 &
    
    echo "Calling peaks on the 2 replicates merged..."
    macs2 callpeak -t filtered1.bam filtered2.bam -c filtered_ctrl1.bam -q "$FDR_THRESHOLD" -g hs -B --SPMR -n MERGE > macs2_MERGE_log.txt 2>&1 &
fi

echo "All peak calling jobs started, wait a while..."
wait

echo
echo "Check the outcome and logs of the peak calling in the macs2_*_log.txt files."

# R plots
echo
echo "Plotting..."
R < REP1_model.r --vanilla > /dev/null
R < REP2_model.r --vanilla > /dev/null
R < MERGE_model.r --vanilla > /dev/null

echo
echo "R plots are available as pdfs."

# Computing overlaps
echo
echo "Computing overlaps..."
REP1_COUNT=$(wc -l < REP1_peaks.narrowPeak)
REP2_COUNT=$(wc -l < REP2_peaks.narrowPeak)
if [ "$REP1_COUNT" -le "$REP2_COUNT" ]; then
    SMALL="REP1"
    BIG="REP2"
    SMALL_COUNT="$REP1_COUNT"
    BIG_COUNT="$REP2_COUNT"
else
    SMALL="REP2"
    BIG="REP1"
    SMALL_COUNT="$REP2_COUNT"
    BIG_COUNT="$REP1_COUNT"
fi

REP_COMP_NAME="${SMALL}_in_${BIG}"
bedtools intersect -a ${SMALL}_peaks.narrowPeak -b ${BIG}_peaks.narrowPeak -u > "${REP_COMP_NAME}_peaks.narrowPeak"

SUMMIT_PROX=$(bedtools closest -a ${SMALL}_summits.bed -b ${BIG}_summits.bed -d | awk '$NF <= 100' | wc -l)
OVERLAPS=$(wc -l < "${REP_COMP_NAME}_peaks.narrowPeak")
FRACTION=$(awk "BEGIN {printf \"%.3f\", $OVERLAPS/$SMALL_COUNT}")
PROXIMITY_FRAC=$(awk "BEGIN {printf \"%.3f\", $SUMMIT_PROX/$SMALL_COUNT}")

echo | tee -a "$OUTPUT"
echo "${SMALL} (${SMALL_COUNT} peaks) in ${BIG} (${BIG_COUNT} peaks):" | tee -a "$OUTPUT"
echo "  Summit proximity: $SUMMIT_PROX" | tee -a "$OUTPUT"
echo "  Overlapping peaks: $OVERLAPS" | tee -a "$OUTPUT"
echo "  Fraction of overlapping peaks: $FRACTION" | tee -a "$OUTPUT"
echo "  Summit proximity ratio: $PROXIMITY_FRAC" | tee -a "$OUTPUT"


MERGE_COUNT=$(wc -l < MERGE_peaks.narrowPeak)
OVERLAPS=$(bedtools intersect -a REP1_peaks.narrowPeak -b MERGE_peaks.narrowPeak -u | wc -l)
FRACTION=$(awk "BEGIN {printf \"%.3f\", $OVERLAPS/$REP1_COUNT}")

echo | tee -a "$OUTPUT"
echo "REP1 (${REP1_COUNT} peaks) in MERGE (${MERGE_COUNT} peaks):" | tee -a "$OUTPUT"
echo "  Overlapping peaks: $OVERLAPS" | tee -a "$OUTPUT"
echo "  Fraction of overlapping peaks: $FRACTION" | tee -a "$OUTPUT"


OVERLAPS=$(bedtools intersect -a REP2_peaks.narrowPeak -b MERGE_peaks.narrowPeak -u | wc -l)
FRACTION=$(awk "BEGIN {printf \"%.3f\", $OVERLAPS/$REP2_COUNT}")

echo | tee -a "$OUTPUT"
echo "REP2 (${REP2_COUNT} peaks) in MERGE (${MERGE_COUNT} peaks):" | tee -a "$OUTPUT"
echo "  Overlapping peaks: $OVERLAPS" | tee -a "$OUTPUT"
echo "  Fraction of overlapping peaks: $FRACTION" | tee -a "$OUTPUT"


# Comparison with the ENCODE results
echo
echo "Comparing with ENCODE results..."
bedtools sort -i "$ENCODE" > ENCODE_peaks.narrowPeak

ENCODE_COUNT=$(wc -l < ENCODE_peaks.narrowPeak)
echo | tee -a "$OUTPUT"
echo "ENCODE reports ${ENCODE_COUNT} peaks."

NAMES=("REP1" "REP2" "MERGE" "$REP_COMP_NAME")
for NAME in "${NAMES[@]}"; do
    echo | tee -a "$OUTPUT"
    echo "Jaccard index for ${NAME} vs ENCODE:" | tee -a "$OUTPUT"
    JACCARD_INDEX=$(bedtools jaccard -a "${NAME}_peaks.narrowPeak" -b ENCODE_peaks.narrowPeak)
    echo "$JACCARD_INDEX" | column -t | tee -a "$OUTPUT"
    
    # Extract actual JI:
    JACCARD_INDEX=$(echo "$JACCARD_INDEX" | awk "NR==2 {print \$3}")
    if [ "$NAME" = "MERGE" ]; then
        MERGE_JI="$JACCARD_INDEX"
    elif [ "$NAME" = "$REP_COMP_NAME" ]; then
        INTERSECT_JI="$JACCARD_INDEX"
    fi
done


# Determine final peak set:
echo | tee -a "$OUTPUT"
echo "Filtering..."
if awk "BEGIN {exit !($INTERSECT_JI > $MERGE_JI)}"; then
    FINAL_SET="$REP_COMP_NAME"
    bedtools intersect -a ${SMALL}_summits.bed -b ${FINAL_SET}_peaks.narrowPeak -u > ${FINAL_SET}_summits.bed
else
    FINAL_SET="MERGE"
fi


# Filter out black-listed regions:
bedtools intersect -a ${FINAL_SET}_peaks.narrowPeak -b "${BLACKLIST}" -v > ${FINAL_SET}_filtered_peaks.narrowPeak
bedtools intersect -a ${FINAL_SET}_summits.bed -b ${FINAL_SET}_filtered_peaks.narrowPeak -u > ${FINAL_SET}_filtered_summits.bed

FINAL_COUNT=$(wc -l < ${FINAL_SET}_peaks.narrowPeak)
FILTERED_COUNT=$(wc -l < ${FINAL_SET}_filtered_peaks.narrowPeak)
echo "$((FINAL_COUNT - FILTERED_COUNT)) regions were black-listed and have been removed." | tee -a "$OUTPUT"

echo "The final peak set is ${FINAL_SET}_filtered_peaks.narrowPeak, you can find the corresponding summits in ${FINAL_SET}_filtered_summits.bed." | tee -a "$OUTPUT"


# Generating boxplots:
echo
echo "Generating boxplots..."

# Check if the script already exists:
if [ ! -f "boxplots.R" ]; then
    cat << "EOF" > boxplots.R
args <- commandArgs(trailingOnly=TRUE)
overlapData <- read.table(args[1])[ ,9]
uniqueData  <- read.table(args[2])[ ,9]

pdf("boxplots.pdf")
boxplot(overlapData, uniqueData, names = c("ENCODE Overlap", "Not in ENCODE"), ylab = "-log10(q-value)")
dev.off()
EOF
fi

ENCODE_INTERSECT="ENCODE_and_${FINAL_SET}_filtered_peaks.narrowPeak"
FINAL_SET_NO_ENCODE="${FINAL_SET}_not_in_ENCODE_filtered_peaks.narrowPeak"

bedtools intersect -a "${FINAL_SET}_filtered_peaks.narrowPeak" -b ENCODE_peaks.narrowPeak -u > $ENCODE_INTERSECT
bedtools intersect -a "${FINAL_SET}_filtered_peaks.narrowPeak" -b ENCODE_peaks.narrowPeak -v > $FINAL_SET_NO_ENCODE

Rscript boxplots.R $ENCODE_INTERSECT $FINAL_SET_NO_ENCODE
echo "Boxplots are available in boxplots.pdf." | tee -a "$OUTPUT"

echo
echo "All done!"