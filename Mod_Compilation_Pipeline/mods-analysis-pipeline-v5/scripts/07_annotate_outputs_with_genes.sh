#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash scripts/07_annotate_outputs_with_genes.sh <ROOT_DIR> <SAMPLE_LIST> <REF_GTF>
#
# Example:
#   bash scripts/07_annotate_outputs_with_genes.sh \
#     /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/yale_brain/NEW_MODEL/comp1 \
#     /athena/masonlab/scratch/users/thn4005/scripts/mods-analysis-pipeline-v5/config/yale_brain.txt \
#     /athena/masonlab/scratch/users/thn4005/scripts/references_v2/references/human/hg38.ncbiRefSeq.gtf
#
# Notes:
# - SAMPLE_LIST is accepted to mirror the pipeline signature, but not used here.
# - Requires: bedtools. Optionally uses gtf2bed if available; otherwise falls back to awk parsing.

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <ROOT_DIR> <SAMPLE_LIST> <REF_GTF>" >&2
  exit 1
fi

ROOT_DIR="$(readlink -f "$1")"
SAMPLE_LIST="$2"   # not used; accepted for signature parity
REF_GTF="$(readlink -f "$3")"

OUT_DIR="${ROOT_DIR}/outputs"
MERGED_DIR="${OUT_DIR}/merged_matrix"
Z_DIR="${OUT_DIR}/zscore_filtered"

MERGED_ANN_DIR="${OUT_DIR}/merged_matrix_annotated"
Z_ANN_DIR="${OUT_DIR}/zscore_filtered_annotated"
TMP_DIR="${OUT_DIR}/__annot_tmp"

mkdir -p "$MERGED_ANN_DIR" "$Z_ANN_DIR" "$TMP_DIR"

BEDTOOLS="$(command -v bedtools || true)"
if [[ -z "$BEDTOOLS" ]]; then
  echo "[ERROR] bedtools not found in PATH." >&2
  exit 1
fi

GTF2BED_BIN="${GTF2BED:-$(command -v gtf2bed || true)}"

EXONS_BED="${TMP_DIR}/exons_gene_name.bed"

echo "[INFO] Building exon BED with gene names from: ${REF_GTF}"
if [[ -n "$GTF2BED_BIN" ]]; then
  # Use gtf2bed if available (preferred)
  cat "$REF_GTF" | "$GTF2BED_BIN" - | grep -w "exon" | \
  awk '
  BEGIN{OFS="\t"}
  {
    gene=""; 
    if (match($0, /gene_name "[^"]+"/, a)) { gene=a[0]; gsub(/gene_name "|"/, "", gene) }
    else if (match($0, /gene_id "[^"]+"/, b)) { gene=b[0]; gsub(/gene_id "|"/, "", gene) }
    else { gene="NA" }
    print $1, $2, $3, gene, ".", $6
  }' > "$EXONS_BED"
else
  # Fallback: awk-only GTF -> BED (0-based start)
  awk '
  BEGIN{OFS="\t"}
  $0 !~ /^#/ && $3=="exon" {
    chrom=$1; start=$4-1; end=$5; strand=$7;
    gene="";
    if (match($0, /gene_name "[^"]+"/, a)) { gene=a[0]; gsub(/gene_name "|"/, "", gene) }
    else if (match($0, /gene_id "[^"]+"/, b)) { gene=b[0]; gsub(/gene_id "|"/, "", gene) }
    else { gene="NA" }
    print chrom, start, end, gene, ".", strand
  }' "$REF_GTF" > "$EXONS_BED"
fi

if [[ ! -s "$EXONS_BED" ]]; then
  echo "[ERROR] Failed to create exon BED with gene names." >&2
  exit 1
fi
echo "[INFO] Exon BED written: $EXONS_BED"

annotate_tsv () {
  local in_tsv="$1"
  local out_tsv="$2"
  local base tmp_sites tmp_ix tmp_map

  base="$(basename "$in_tsv")"
  tmp_sites="${TMP_DIR}/${base}.sites.bed"
  tmp_ix="${TMP_DIR}/${base}.intersect.bed"
  tmp_map="${TMP_DIR}/${base}.site2gene.tsv"

  echo "[INFO] Annotating $(basename "$in_tsv")"

  if [[ ! -s "$in_tsv" ]]; then
    echo "[WARN] Empty input: $in_tsv (skipping)"
    return 0
  fi

  # Build sites BED (chrom, start, end) from the TSV (skip header)
  tail -n +2 "$in_tsv" | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' > "$tmp_sites"

  # Intersect with exons and extract (chrom,start,end,gene); dedup & aggregate per site
  if [[ -s "$tmp_sites" ]]; then
    $BEDTOOLS intersect -wa -wb -a "$tmp_sites" -b "$EXONS_BED" > "$tmp_ix" || true

    if [[ -s "$tmp_ix" ]]; then
      # Pull gene name (4th col of B-file → field NF-2), dedup, then aggregate per site
      awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$(NF-2)}' "$tmp_ix" | sort -u | \
      awk 'BEGIN{OFS="\t"}
           {k=$1"\t"$2"\t"$3; g[k]=(k in g)?g[k]";"$4:$4}
           END{for (k in g) print k, g[k]}' > "$tmp_map"
    else
      : > "$tmp_map"
    fi
  else
    : > "$tmp_map"
  fi

  # Append gene_name column (or "intergenic" if no mapping)
  awk -v OFS="\t" -v FS="\t" '
    FNR==NR { key=$1"\t"$2"\t"$3; map[key]=$4; next }
    FNR==1  { print $0, "gene_name"; next }
            { key=$1"\t"$2"\t"$3; g=(key in map?map[key]:"intergenic"); print $0, g }
  ' "$tmp_map" "$in_tsv" > "$out_tsv"

  echo "[✅] Wrote: $out_tsv"
}

# -------- Annotate merged matrices --------
if [[ -d "$MERGED_DIR" ]]; then
  shopt -s nullglob
  mm_files=( "$MERGED_DIR"/*_merged_union_matrix_wide.tsv "$MERGED_DIR"/*_merged_union_matrix_long.tsv "$MERGED_DIR"/merged_union_matrix_long.tsv )
  if (( ${#mm_files[@]} )); then
    for f in "${mm_files[@]}"; do
      out="${MERGED_ANN_DIR}/$(basename "${f%.tsv}")_annotated.tsv"
      annotate_tsv "$f" "$out"
    done
  else
    echo "[INFO] No merged matrix TSVs found in ${MERGED_DIR} (skipping)."
  fi
  shopt -u nullglob
else
  echo "[INFO] ${MERGED_DIR} not found (skipping merged matrix annotation)."
fi

# -------- Annotate z-score filtered outputs --------
if [[ -d "$Z_DIR" ]]; then
  shopt -s nullglob
  z_files=( "$Z_DIR"/*.tsv )
  if (( ${#z_files[@]} )); then
    for f in "${z_files[@]}"; do
      out="${Z_ANN_DIR}/$(basename "${f%.tsv}")_annotated.tsv"
      annotate_tsv "$f" "$out"
    done
  else
    echo "[INFO] No TSVs found in ${Z_DIR} (skipping)."
  fi
  shopt -u nullglob
else
  echo "[INFO] ${Z_DIR} not found (skipping z-score annotation)."
fi

echo "[DONE] Annotations complete."

