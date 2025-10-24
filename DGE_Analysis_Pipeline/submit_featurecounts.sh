#!/bin/bash

# === CONFIGURATION ===
BASE_DIR="/athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2"
REF_GTF="/athena/masonlab/scratch/users/thn4005/scripts/references_v2/references/human/hg38.ncbiRefSeq.gtf"
MODEL="directRNA"

# === OUTPUT LOGS ===
LOG_DIR="${BASE_DIR}/featurecounts_logs"
mkdir -p "$LOG_DIR"

# === SCAN FOR BAM FILES AND SUBMIT JOBS ===
find "$BASE_DIR" -type f -name "*_aligned.sorted.bam" | while read -r bam_path; do
    sample_dir=$(dirname "$bam_path")
    bam_file=$(basename "$bam_path")
    sample_name="${bam_file%%_aligned.sorted.bam}"
    counts_dir="${sample_dir}/counts"
    mkdir -p "$counts_dir"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=fc_${sample_name}
#SBATCH --output=${LOG_DIR}/${sample_name}_featurecounts.out
#SBATCH --error=${LOG_DIR}/${sample_name}_featurecounts.err
#SBATCH --partition=scu-cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Activate micromamba environment
source ~/.bashrc
micromamba activate subread_env

# Run featureCounts
featureCounts -L \
  -T \${SLURM_CPUS_PER_TASK} \
  -a "$REF_GTF" \
  -o "${counts_dir}/${sample_name}_${MODEL}_counts.txt" \
  "$bam_path"

micromamba deactivate
EOF

    echo "✅ Submitted job for $sample_name"
done

