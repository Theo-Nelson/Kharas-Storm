#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=union_mod_sites
#SBATCH --output=logs/union_sites_%j.out
#SBATCH --error=logs/union_sites_%j.err

set -euo pipefail

# ============================
# CONFIGURATION SECTION
# ============================
MOD_CODES=(a m 17596 17802 19227 19228 19229 69426)
ROOT_DIR=$1
FILTERED_BEDS_DIR="${ROOT_DIR}/outputs/filtered_beds"
UNION_OUTPUT_DIR="${ROOT_DIR}/outputs/union_beds/unions"

mkdir -p "$UNION_OUTPUT_DIR"

# ============================
# GENERATE UNION FILES
# ============================
echo "📌 Generating union sites for each mod code..."

for mod in "${MOD_CODES[@]}"; do
    echo "🔄 Processing mod code: $mod"
    out_union="$UNION_OUTPUT_DIR/all_samples_${mod}_union_sites.bed"

    find "$FILTERED_BEDS_DIR" -name "*_${mod}_filtered_passonly.bed" \
        -exec grep -v '^#' {} \; \
        | awk -F '\t' '{print $1"\t"$2"\t"$3}' \
        | sort -k1,1 -k2,2n -k3,3n | uniq > "$out_union"

    echo "  ✅ Union file written: $out_union"
done

echo "✅ All union BEDs generated."

