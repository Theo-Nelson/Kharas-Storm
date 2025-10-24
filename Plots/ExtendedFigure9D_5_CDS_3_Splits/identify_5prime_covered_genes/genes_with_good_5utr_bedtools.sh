#!/usr/bin/env bash
set -euo pipefail

############################################
# CONFIG
############################################
# 5'UTR per-gene BED *segments* (multiple rows per gene allowed; BED6)
# Produced from your GTF with the helper you ran.
UTR5_SEG_BED="hg38.utr5_by_gene.sorted.bed"

# Base directory containing sample subdirs with *_filtered_mod.bed
BASE="/athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2"

# Coverage threshold and fraction of 5'UTR required
COV=20
FRAC=0.6667

# Max concurrent jobs (tune to your storage)
N_JOBS=48

OUTDIR="genes_with_good_5utr_bedtools_v2"
mkdir -p "$OUTDIR" "$OUTDIR/tmp"

# Samples (12)
SAMPLES=(
"16-Jun-25_XuejingYang_RNADirect_Sample10_no_sample_id_20250616_1418_2B_PBE86098_1adb3d7b_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample11_no_sample_id_20250616_1419_2C_PBE87545_8e394b62_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample12_no_sample_id_20250616_1420_2D_PBE87459_e98ccd9b_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample1_no_sample_id_20250616_1413_1A_PBE87470_f8113319_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample2_no_sample_id_20250616_1413_1B_PBE86077_49f05d8e_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample4_no_sample_id_20250616_1414_1D_PBE87294_0ee51c63_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample5_no_sample_id_20250616_1414_1E_PBE87462_0fc36d5d_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample7_no_sample_id_20250616_1415_1G_PBE82517_0e771aa9_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample8_no_sample_id_20250616_1416_1H_PBE87488_fd92caf3_pod5"
"16-Jun-25_XuejingYang_RNADirect_Sample9_no_sample_id_20250616_1417_2A_PBE88213_ea08403b_pod5"
)

# 8 modification IDs present in your dirs
MODS=(17596 17802 19227 19228 19229 69426 a m)

############################################
# PRECOMPUTE: total 5'UTR length per gene (from segment BED)
############################################
LEN_TSV="${OUTDIR}/utr5_total_length_per_gene.tsv"
if [[ ! -s "$LEN_TSV" ]]; then
  # Sum length across all 5'UTR segments per gene (no strand needed for length)
  awk 'BEGIN{OFS="\t"} $0!~/^#/ {len[$4]+=($3-$2)} END{for (g in len) print g, len[g]}' "$UTR5_SEG_BED" \
    | sort -k1,1 > "$LEN_TSV"
fi

############################################
# HELPER: run up to N_JOBS in parallel
############################################
joblim() {
  local njobs
  while true; do
    njobs=$(jobs -pr | wc -l)
    (( njobs < N_JOBS )) && break
    sleep 0.5
  done
}

############################################
# STEP 1: Build merged cov>=COV intervals per sample+mod (NO resort)
############################################
build_cov_bed() {
  local sample="$1" mod="$2"
  local inbed="${BASE}/${sample}/${sample}_${mod}_filtered_mod.bed"
  local outbed="${OUTDIR}/${sample}.${mod}.cov${COV}.merged.sorted.bed"

  if [[ -s "$outbed" ]]; then
    echo "[SKIP] $outbed"
    return 0
  fi
  if [[ ! -s "$inbed" ]]; then
    echo "[WARN] missing $inbed"; : > "$outbed"; return 0
  fi

  # Stream filter col10 (valid_coverage) >= COV, then merge.
  # Assumes *_filtered_mod.bed is already sorted by (chrom,start).
  # bedtools merge requires sorted input; we pass -i - (stdin).
  awk -v FS='\t' -v OFS='\t' -v C="$COV" 'NR>1 && NF>=10 && $10+0 >= C {print $1,$2,$3}' "$inbed" \
    | bedtools merge -i - \
    > "$outbed"
}
export -f build_cov_bed
export BASE OUTDIR COV

echo "[INFO] Building merged coverage intervals (>=${COV})"
for mod in "${MODS[@]}"; do
  for s in "${SAMPLES[@]}"; do
    joblim
    bash -c "build_cov_bed \"$s\" \"$mod\"" &
  done
done
wait

############################################
# STEP 2: Per-sample per-gene fraction covered over the UTR5 union
# Method: intersect -wo (UTR segments vs merged cov intervals), sum overlap per gene,
# then divide by total UTR length from LEN_TSV.
############################################
coverage_per_gene() {
  local sample="$1" mod="$2"
  local covbed="${OUTDIR}/${sample}.${mod}.cov${COV}.merged.sorted.bed"
  local out="${OUTDIR}/${sample}.${mod}.utr5_fraction.tsv"

  if [[ -s "$out" ]]; then
    echo "[SKIP] $out"; return 0
  fi
  if [[ ! -s "$covbed" ]]; then
    echo "[WARN] missing covbed $covbed"; : > "$out"; return 0
  fi

  # When no coverage intervals, everything is 0; just emit 0 for all genes
  if [[ ! -s "$covbed" || ! -s "$UTR5_SEG_BED" ]]; then
    awk 'BEGIN{OFS="\t"} {print $1, 0}' "$LEN_TSV" > "$out"
    return 0
  fi

  # Intersect segments with cov intervals, sum overlap per gene (col 4 of A), then join with total length.
  # -sorted speeds up if both inputs are grouped by chrom then start.
  bedtools intersect -a "$UTR5_SEG_BED" -b "$covbed" -wo -sorted \
    | awk 'BEGIN{OFS="\t"} {gene=$4; olap=$NF; sum[gene]+=olap} END{for(g in sum) print g, sum[g]}' \
    | sort -k1,1 \
    | join -t $'\t' -a 2 -e 0 -o 2.1,1.2 -1 1 -2 1 - "$LEN_TSV" \
    | awk -v OFS="\t" '{gene=$1; olap=$2+0; tot=$3+0; frac=(tot>0?olap/tot:0); print gene, frac}' \
    > "$out"

  # Ensure all genes appear (even with 0), if you want a dense table:
  # join -t $'\t' -a 2 -a 1 -e 0 -o 2.1,1.2 -1 1 -2 1 <(sort -k1,1 "$out") <(awk '{print $1,0}' "$LEN_TSV" | sort -k1,1) | sort -u -k1,1 > "$out.tmp" && mv "$out.tmp" "$out"
}
export -f coverage_per_gene
export UTR5_SEG_BED LEN_TSV

echo "[INFO] Computing per-gene UTR5 coverage fraction per sample/mod"
for mod in "${MODS[@]}"; do
  for s in "${SAMPLES[@]}"; do
    joblim
    bash -c "coverage_per_gene \"$s\" \"$mod\"" &
  done
done
wait

############################################
# STEP 3: For each mod, keep genes with fraction >= FRAC in ALL 12 samples
############################################
echo "[INFO] Selecting genes passing in all 12 samples (FRAC >= ${FRAC})"
for mod in "${MODS[@]}"; do
  # Start with genes that pass in the first sample:
  first="${OUTDIR}/${SAMPLES[0]}.${mod}.utr5_fraction.tsv"
  if [[ ! -s "$first" ]]; then
    echo "[WARN] missing $first; making empty pass list for mod $mod"
    : > "${OUTDIR}/genes_pass_all12_${mod}.list"
    continue
  fi

  awk -v F="$FRAC" '$2+0 >= F {print $1}' "$first" | sort -u > "${OUTDIR}/_pass.${mod}.list"

  # Intersect with the pass-lists from the remaining samples
  for ((i=1;i<${#SAMPLES[@]};i++)); do
    s="${SAMPLES[$i]}"
    tab="${OUTDIR}/${s}.${mod}.utr5_fraction.tsv"
    if [[ ! -s "$tab" ]]; then
      echo "[WARN] missing $tab; emptying pass list for mod $mod"
      : > "${OUTDIR}/_pass.${mod}.list"
      break
    fi
    awk -v F="$FRAC" '$2+0 >= F {print $1}' "$tab" | sort -u > "${OUTDIR}/_pass.${mod}.list.$i"
    comm -12 <(sort -u "${OUTDIR}/_pass.${mod}.list") <(sort -u "${OUTDIR}/_pass.${mod}.list.$i") \
      > "${OUTDIR}/_pass.${mod}.next"
    mv "${OUTDIR}/_pass.${mod}.next" "${OUTDIR}/_pass.${mod}.list"
    rm -f "${OUTDIR}/_pass.${mod}.list.$i"
  done

  mv "${OUTDIR}/_pass.${mod}.list" "${OUTDIR}/genes_pass_all12_${mod}.list"
  echo "[OK] ${OUTDIR}/genes_pass_all12_${mod}.list"
done

echo "[DONE] Outputs in $OUTDIR
- ${OUTDIR}/*.cov${COV}.merged.sorted.bed (per sample+mod merged cov>=${COV})
- ${OUTDIR}/*.utr5_fraction.tsv (per sample+mod: gene<TAB>fraction)
- ${OUTDIR}/genes_pass_all12_<mod>.list (8 files; genes passing in all 12 samples)"

