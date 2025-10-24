#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=intersect_union_mods
#SBATCH --output=logs/intersect_union_%j.out
#SBATCH --error=logs/intersect_union_%j.err

# ============================
# CONFIGURATION SECTION
# ============================
MOD_CODES=(a m 17596 17802 19227 19228 19229 69426)
ROOT_DIR=$1  # Assumes script is run from base directory
UNION_DIR="$ROOT_DIR/outputs/union_beds/unions"
FILTERED_BEDS_DIR="$ROOT_DIR/outputs/filtered_beds"
SAMPLE_DIRS=( $(< $2) )
INTERSECT_OUT_DIR="$ROOT_DIR/outputs/union_beds/by_sample"

mkdir -p "$INTERSECT_OUT_DIR"

source ~/.bashrc
micromamba activate samtools_env

# ============================
# INTERSECTION LOGIC
# ============================
echo "🔗 Intersecting union sites with original modkit BEDs..."

for mod in "${MOD_CODES[@]}"; do
    union_bed="$UNION_DIR/all_samples_${mod}_union_sites.bed"

    for sample_path in "${SAMPLE_DIRS[@]}"; do
        sample_name=$(basename "$sample_path")
        in_bed="$sample_path/${sample_name}_${mod}_filtered_mod.bed"

        if [[ ! -f "$in_bed" ]]; then
            echo "  ⚠️ Missing BED: $in_bed"
            continue
        fi

        out_bed="$INTERSECT_OUT_DIR/${sample_name}_${mod}_union_sites_with_attrs.bed"

        echo "  📎 Intersecting $in_bed with $union_bed"

        awk -v mod="$mod" 'BEGIN{OFS="\t"} /^#/ {print; next} $4 == mod {print}' "$in_bed" > tmp_filtered_${mod}.bed

        bedtools intersect -a tmp_filtered_${mod}.bed -b "$union_bed" -wa -header > "$out_bed"
        rm tmp_filtered_${mod}.bed

        echo "  ✅ Wrote: $out_bed"
    done
done

echo "🏁 All union intersection jobs complete."

