#!/usr/bin/env bash
set -euo pipefail

# -----------------------
# Config (edit as needed)
# -----------------------
MODKIT_BIN="${MODKIT_BIN:-modkit}"                 # path to modkit in your env
REF_FASTA="${REFERENCE_GENOME:-hg38.fa}"  # REQUIRED: set to your hg38 FASTA
GTF="${GTF:-hg38.ncbiRefSeq_V2.gtf}"               # your GTF in this folder

# Thresholds/settings
REF_BASE="${REF:-A}"            # fundamental base for filter-threshold (A for RNA m6A runs)
REF_THR="${REF_THR:-0.8}"       # e.g. 0.8  -> --filter-threshold A:0.8
MOD_CODE="${MOD:-a}"            # 'a' == m6A (A+a)
MOD_THR="${MOD_THR:-0.99}"      # e.g. 0.99 -> --mod-thresholds a:0.99
MAX_DEPTH="${MAX_DEPTH:-1000}"
THREADS="${THREADS:-8}"         # set to 64 if you’ve got the cores
INTERVAL_SIZE="${INTERVAL_SIZE:-100000}"  # --interval-size
OUTDIR="${OUTDIR:-pileup_out}"

# Toggle: compress outputs if bgzip exists
BGZIP="${BGZIP:-$(command -v bgzip || true)}"

# -----------------------
# Checks
# -----------------------
if ! command -v "$MODKIT_BIN" >/dev/null 2>&1; then
  echo "ERROR: modkit not found at MODKIT_BIN='$MODKIT_BIN'." >&2
  exit 1
fi
if [ ! -f "$REF_FASTA" ]; then
  echo "ERROR: REF_FASTA not found: $REF_FASTA" >&2
  exit 1
fi
if [ ! -f "$GTF" ]; then
  echo "ERROR: GTF not found: $GTF" >&2
  exit 1
fi

# Optional: re-index BAMs if needed
if command -v samtools >/dev/null 2>&1; then
  for b in *.bam; do
    [ -f "$b" ] || continue
    bai="${b}.bai"
    if [ ! -f "$bai" ] || [ "$bai" -ot "$b" ]; then
      echo "[info] (re)indexing $b ..."
      samtools index -@ 4 -b "$b"
    fi
  done
fi

mkdir -p "$OUTDIR"

# -----------------------
# Get gene regions (1-based inclusive)
# -----------------------
get_region () {
  local gene="$1"
  awk -v gene="$gene" 'BEGIN{chr=""; s=0; e=0}
    $0 ~ "gene_name \""gene"\"" && $3 ~ /^(transcript|exon)$/ {
      if (chr==""){chr=$1; s=$4; e=$5}
      if ($4 < s) s=$4;
      if ($5 > e) e=$5;
    }
    END{
      if (chr=="") { exit 2 }
      printf "%s:%d-%d\n", chr, s, e
    }' "$GTF"
}

BCL2_REGION="$(get_region BCL2 || true)"
MYC_REGION="$(get_region MYC || true)"

if [ -z "$BCL2_REGION" ] || [ -z "$MYC_REGION" ]; then
  echo "ERROR: could not extract BCL2 or MYC region from $GTF" >&2
  echo "  BCL2='$BCL2_REGION'  MYC='$MYC_REGION'" >&2
  exit 1
fi

echo "[info] Regions:"
echo "  BCL2: $BCL2_REGION"
echo "  MYC : $MYC_REGION"

# -----------------------
# Run modkit pileup
# -----------------------
run_one () {
  local bam="$1"
  local region="$2"
  local gene="$3"

  local base="$(basename "$bam" .bam)"
  local sample_outdir="$OUTDIR/$base"
  mkdir -p "$sample_outdir"

  local log="${sample_outdir}/${base}_${gene}_filtered_mod_bed_log.txt"
  local out="${sample_outdir}/${base}_${gene}_filtered_mod.bed"

  echo "[run] $base :: $gene -> $out"

  if [ -n "$BGZIP" ]; then
    # write compressed if bgzip is available
    out="${out}.gz"
    "$MODKIT_BIN" pileup "$bam" /dev/stdout \
      --ref "$REF_FASTA" \
      --region "$region" \
      --filter-threshold "${REF_BASE}:${REF_THR}" \
      --mod-thresholds "${MOD_CODE}:${MOD_THR}" \
      --log-filepath "$log" \
      --max-depth "$MAX_DEPTH" \
      -t "$THREADS" \
      --with-header \
      --interval-size "$INTERVAL_SIZE" \
    | "$BGZIP" -c > "$out"
  else
    "$MODKIT_BIN" pileup "$bam" "$out" \
      --ref "$REF_FASTA" \
      --region "$region" \
      --filter-threshold "${REF_BASE}:${REF_THR}" \
      --mod-thresholds "${MOD_CODE}:${MOD_THR}" \
      --log-filepath "$log" \
      --max-depth "$MAX_DEPTH" \
      -t "$THREADS" \
      --with-header \
      --interval-size "$INTERVAL_SIZE"
  fi
}

export MODKIT_BIN REF_FASTA REF_BASE REF_THR MOD_CODE MOD_THR MAX_DEPTH THREADS INTERVAL_SIZE OUTDIR BGZIP
export -f run_one

# Loop all BAMs in current dir
shopt -s nullglob
BAMS=(*.bam)
if [ ${#BAMS[@]} -eq 0 ]; then
  echo "ERROR: no *.bam files found in $(pwd)" >&2
  exit 1
fi

for bam in "${BAMS[@]}"; do
  run_one "$bam" "$BCL2_REGION" "BCL2"
  run_one "$bam" "$MYC_REGION"  "MYC"
done

echo "[done] Outputs under: $OUTDIR/"

