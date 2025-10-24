#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
UTR5_BED="hg38.utr5_by_gene.sorted.bed"     # per-gene 5'UTR BED6
OUTDIR="genes_with_good_5utr_bedtools_v2"   # where cov beds live & outputs go
COV=20                                      # part of cov bed filenames
MIN_FRAC="0"                                # keep genes with fraction > MIN_FRAC (0 = any coverage)

# 12 samples
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

# 8 modification IDs (file stems)
MODS=(17596 17802 19227 19228 19229 69426 a m)

# --- helpers ---
is_primary() {
  # keep only primary chroms: chr1-22, chrX, chrY, chrM
  awk 'BEGIN{FS=OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/{print}'
}

sort_bed() {
  # coordinate sort: chr then start
  sort -k1,1 -k2,2n
}

canon_and_sort() {
  # stdin -> filter to primary -> sort -> stdout
  is_primary | sort_bed
}

# --- checks ---
if [[ ! -s "$UTR5_BED" ]]; then
  echo "[ERROR] Missing $UTR5_BED" >&2
  exit 1
fi
if [[ ! -d "$OUTDIR" ]]; then
  echo "[ERROR] Missing dir $OUTDIR" >&2
  exit 1
fi

# Build a canonical, sorted UTR5 file once (primary chroms only)
UTR5_CANON="${OUTDIR}/_utr5.primary.sorted.bed"
if [[ ! -s "$UTR5_CANON" ]]; then
  echo "[INFO] making primary-only, sorted UTR5: $UTR5_CANON"
  canon_and_sort < "$UTR5_BED" > "$UTR5_CANON"
fi

echo "[INFO] Recomputing 5'UTR fractions with canonical chroms + consistent sorting..."
for mod in "${MODS[@]}"; do
  echo "[INFO]  -> mod ${mod}"
  # per-sample in parallel
  for s in "${SAMPLES[@]}"; do
    (
      B_RAW="${OUTDIR}/${s}.${mod}.cov${COV}.merged.sorted.bed"
      B_CANON="${OUTDIR}/${s}.${mod}.cov${COV}.merged.primary.sorted.bed"
      OUT_TSV="${OUTDIR}/${s}.${mod}.utr5_fraction.tsv"
      OUT_LIST="${OUTDIR}/${s}.${mod}.genes_frac_gt${MIN_FRAC}.list"

      if [[ ! -s "$B_RAW" ]]; then
        echo "[WARN] missing/empty: $B_RAW"
        : > "$OUT_TSV"; : > "$OUT_LIST"
        exit 0
      fi

      # filter alts & resort (even if already sorted)
      if [[ ! -s "$B_CANON" || "$B_CANON" -ot "$B_RAW" ]]; then
        canon_and_sort < "$B_RAW" > "$B_CANON"
      fi

      # run coverage (sweeping algo) now that both inputs share the same primary-only ordering
      if ! bedtools coverage -a "$UTR5_CANON" -b "$B_CANON" -sorted \
            | awk -v OFS="\t" '{print $4, $NF}' > "$OUT_TSV" ; then
        echo "[WARN] -sorted failed for $s $mod; retrying without -sorted (slower)"
        bedtools coverage -a "$UTR5_CANON" -b "$B_CANON" \
            | awk -v OFS="\t" '{print $4, $NF}' > "$OUT_TSV"
      fi

      awk -v thr="$MIN_FRAC" '$2>thr{print $1}' "$OUT_TSV" | sort -u > "$OUT_LIST"
      echo "[OK] $OUT_TSV"
    ) &
  done
  wait

  # Intersect genes across all 12 samples for this mod
  echo "[INFO]  -> intersecting gene lists across all 12 for mod ${mod}"
  first="${OUTDIR}/${SAMPLES[0]}.${mod}.genes_frac_gt${MIN_FRAC}.list"
  if [[ ! -s "$first" ]]; then
    echo "[WARN]   empty/missing first list for ${mod}; intersection will be empty."
    : > "${OUTDIR}/genes_anycov_all12_${mod}.list"
    continue
  fi
  sort -u "$first" > "${OUTDIR}/.tmp.${mod}.list"
  for ((i=1;i<${#SAMPLES[@]};i++)); do
    s="${SAMPLES[$i]}"
    lst="${OUTDIR}/${s}.${mod}.genes_frac_gt${MIN_FRAC}.list"
    if [[ ! -s "$lst" ]]; then
      echo "[WARN]   missing list for $s ${mod}; intersection becomes empty."
      : > "${OUTDIR}/.tmp.${mod}.list"
      break
    fi
    sort -u "$lst" > "${OUTDIR}/.tmp.${mod}.list.$i"
    comm -12 "${OUTDIR}/.tmp.${mod}.list" "${OUTDIR}/.tmp.${mod}.list.$i" > "${OUTDIR}/.tmp.${mod}.list.next"
    mv "${OUTDIR}/.tmp.${mod}.list.next" "${OUTDIR}/.tmp.${mod}.list"
    rm -f "${OUTDIR}/.tmp.${mod}.list.$i"
  done
  mv "${OUTDIR}/.tmp.${mod}.list" "${OUTDIR}/genes_anycov_all12_${mod}.list"
  echo "[OK] ${OUTDIR}/genes_anycov_all12_${mod}.list"
done

echo "[DONE] Updated per-sample fractions and all-12 intersections in $OUTDIR"

