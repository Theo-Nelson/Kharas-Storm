#!/usr/bin/env bash
# Unified modkit pileup for BCL2 + MYC + SP1 from the Kharas_Subset_MYC_BCL2_SP1 BAMs.
# Produces one consistent run across all three genes (modkit 0.5.0, bioconda osx-arm64).
# Thresholds match the original run_pileup_BCL2_MYC.sh.
set -euo pipefail
MODKIT="${MODKIT:-modkit}"                  # macOS-arm64 build: bioconda ont-modkit 0.5.0
REF="${REF:-ref_chr8_12_18.fa}"             # chr8+chr12+chr18 extracted from hg38.fa.gz (+ .fai)
BAMDIR="${BAMDIR:-.}"                       # dir containing *_SP1.bam (contains MYC, BCL2, SP1)
OUT="${OUT:-pileup_out_v2}"
# regions derived from hg38.ncbiRefSeq_V2.gtf (min start / max end of transcript+exon features)
REGIONS="BCL2:chr18:63123346-63320090 MYC:chr8:127735434-127742951 SP1:chr12:53380176-53416446"
mkdir -p "$OUT"
for bam in "$BAMDIR"/*_SP1.bam; do
  base=$(basename "$bam" .bam); od="$OUT/$base"; mkdir -p "$od"
  for spec in $REGIONS; do
    g="${spec%%:*}"; reg="${spec#*:}"
    "$MODKIT" pileup "$bam" "$od/${base}_${g}_filtered_mod.bed" \
      --ref "$REF" --region "$reg" \
      --filter-threshold A:0.8 --mod-thresholds a:0.99 \
      --max-depth 1000 -t 6 --with-header --interval-size 100000 \
      --log-filepath "$od/${base}_${g}_log.txt"
  done
done
