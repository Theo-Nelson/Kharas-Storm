#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=modkit_full_pipeline
#SBATCH --output=logs/run_full_pipeline_%j.out
#SBATCH --error=logs/run_full_pipeline_%j.err

# ============================
# CONFIG
# ============================
MOD_CODES=(a m 17596 17802 19227 19228 19229 69426)
SAMPLE_LIST=$2
if [ -z "$1" ]; then
  echo "❌ Error: Must provide ROOT_DIR as the first argument."
  exit 1
fi
ROOT_DIR=$1
FILTER_SCRIPT=scripts/filter_bed_sites.py

FILTERED_DIR=$ROOT_DIR/outputs/filtered_beds
UNION_DIR=$ROOT_DIR/outputs/union_beds/unions
INTERSECT_DIR=$ROOT_DIR/outputs/union_beds/by_sample
MERGED_MATRIX=$ROOT_DIR/outputs/merged_matrix/merged_union_matrix_long.tsv
ZFILTER_DIR=$ROOT_DIR/outputs/zscore_filtered
ZTHRESH=3.0

source ~/.bashrc

# ============================
# STEP 1: Filter modkit BEDs
# ============================
echo "🚦 Step 1: Filtering modkit BEDs"
python3 scripts/01_filter_beds.py \
  --sample_dirs $(cat $SAMPLE_LIST) \
  --mod_codes ${MOD_CODES[@]} \
  --output_dir $FILTERED_DIR \
  --filter_script $FILTER_SCRIPT

# ============================
# STEP 2: Generate union sites
# ============================
echo "🚦 Step 2: Generating union BEDs"
UNION_JOBID=$(sbatch --parsable scripts/02_generate_union_beds.sbatch $ROOT_DIR)
echo "  📋 Submitted union job as $UNION_JOBID"

# ============================
# STEP 3: Intersect with union
# ============================
echo "🚦 Step 3: Intersecting with union sites"
INTERSECT_JOBID=$(sbatch --parsable --dependency=afterok:$UNION_JOBID scripts/03_intersect_with_union.sbatch $ROOT_DIR $SAMPLE_LIST)
echo "  📋 Submitted intersection job as $INTERSECT_JOBID (after $UNION_JOBID)"

# ============================
# STEP 4 & 5: Merge and z-filter
# ============================
echo "🚦 Step 4/5: Merging and z-score filtering (after intersection)"
MERGE_JOBID=$(sbatch --parsable --dependency=afterok:$INTERSECT_JOBID <<EOF
#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=merge_zfilter
#SBATCH --output=logs/merge_zfilter_%j.out
#SBATCH --error=logs/merge_zfilter_%j.err

source ~/.bashrc
micromamba activate samtools_env

echo "🧬 Running merge to long-format matrix..."
python3 scripts/04_merge_to_matrix.py --input_dir $INTERSECT_DIR --output_dir $ROOT_DIR/outputs/merged_matrix

echo "📉 Running z-score filtering..."
python3 scripts/05_filter_by_zscore.py \
  --input_dir $ROOT_DIR/outputs/merged_matrix \
  --output_dir $ZFILTER_DIR \
  --z_thresh $ZTHRESH

python3 scripts/06_summarize_pipeline_stats.py \
  --root_dir $ROOT_DIR \
  --sample_list $SAMPLE_LIST \
  --mod_codes a m 17596 17802 19227 19228 19229 69426 \
  --z_thresh $ZTHRESH

micromamba deactivate
EOF
)

echo "✅ Pipeline complete."

