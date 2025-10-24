#!/usr/bin/env bash
set -euo pipefail

# Directory that already contains the per-sample per-mod fraction files
OUTDIR="genes_with_good_5utr_bedtools_v2"

# 12 samples (same order you used before)
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

# 8 modifications
MODS=(17596 17802 19227 19228 19229 69426 a m)

# Threshold to consider a gene as "covered at all" in a sample.
# You asked to "filter for those > 0", so we keep genes with fraction > 0 in every sample.
MIN_FRAC_PER_SAMPLE=0

mkdir -p "$OUTDIR/tmp"

for mod in "${MODS[@]}"; do
  echo "[INFO] processing mod ${mod}"

  # Collect file list for this mod and ensure they exist
  files=()
  for s in "${SAMPLES[@]}"; do
    f="${OUTDIR}/${s}.${mod}.utr5_fraction.tsv"
    if [[ ! -s "$f" ]]; then
      echo "[WARN] missing ${f}; will produce an empty list for ${mod}"
      files=()
      break
    fi
    # Ensure sorted by gene for any downstream debugging/joins (not strictly required for AWK method)
    sort -k1,1 "$f" > "${OUTDIR}/tmp/${s}.${mod}.sorted.tsv"
    files+=("${OUTDIR}/tmp/${s}.${mod}.sorted.tsv")
  done

  out="${OUTDIR}/genes_pass_all12_${mod}.tsv"
  : > "$out"

  if (( ${#files[@]} == 0 )); then
    echo -e "gene\tmin_fraction_across_12" > "$out"
    echo "[OK] ${out} (empty due to missing inputs)"
    continue
  fi

  # AWK across all 12: keep genes present in *every* file with fraction > MIN_FRAC_PER_SAMPLE
  # and write the minimum fraction across the 12 (2nd column).
  awk -v n="${#files[@]}" -v thr="$MIN_FRAC_PER_SAMPLE" '
    {
      g=$1; f=$2+0;
      # track only if this sample shows > thr
      if (f>thr) {
        cnt[g]++; sum[g]+=f;
        if (!(g in min) || f<min[g]) min[g]=f;
      }
      # also store a flag we saw this gene (for completeness, not strictly needed)
      seen[g]=1;
    }
    END {
      print "gene\tmin_fraction_across_12";
      for (g in seen) {
        if (cnt[g]==n) {
          printf "%s\t%.6f\n", g, min[g];
        }
      }
    }
  ' "${files[@]}" | sort -k1,1 > "$out"

  echo "[OK] wrote ${out} ($(wc -l < "$out") lines incl. header)"
done

# Cleanup temp if you like:
# rm -rf "$OUTDIR/tmp"

echo "[DONE] Lists with per-gene min fraction (2nd column) are in ${OUTDIR}/genes_pass_all12_*.tsv"

