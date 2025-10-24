#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --job-name=modalign

# Usage:
# sbatch mods-alignment-from-folder.sbatch <reference_genome> <input_folder> 

REFERENCE_GENOME=$1
INPUT_FOLDER=$2

# Derive sample name from folder
SAMPLE_NAME=$(basename "$INPUT_FOLDER")
SAMPLE_OUTPUT_FOLDER="$INPUT_FOLDER"  # Output will go to same folder

# Locate input BAM
BASECALLED_BAM=$(find "$INPUT_FOLDER" -name "calls_*.bam" | head -n 1)
if [ ! -f "$BASECALLED_BAM" ]; then
    echo "[ERROR] No calls_*.bam file found in $INPUT_FOLDER"
    exit 1
fi

# Define output paths
FASTQ="${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_multi_mod_basecalled.fastq"
SAM="${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_aligned.sam"
BAM="${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_aligned.bam"
SORTED_BAM="${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_aligned.sorted.bam"

# Load environment
source ~/.bashrc

# Step 1: Convert BAM to FASTQ
micromamba activate samtools_env
echo "[INFO] Converting BAM to FASTQ..."
samtools fastq -T "*" "$BASECALLED_BAM" > "$FASTQ"
micromamba deactivate

# Step 2: Align FASTQ to reference genome using minimap2
micromamba activate minimap2_env
echo "[INFO] Aligning FASTQ to genome with minimap2..."
minimap2 -ax splice -t 32 -uf -k14 -y "$REFERENCE_GENOME" "$FASTQ" > "$SAM"
micromamba deactivate

# Step 3: Convert SAM to sorted BAM and index
micromamba activate samtools_env
echo "[INFO] Converting SAM to BAM, sorting, and indexing..."
samtools view -bS "$SAM" > "$BAM"
samtools sort "$BAM" -o "$SORTED_BAM"
samtools index "$SORTED_BAM"
rm "$SAM" "$BAM"
micromamba deactivate

echo "[INFO] Alignment complete. Output saved to: $SORTED_BAM"

# Step 4: Generate and submit modkit jobs
echo "[INFO] Generating and submitting modkit jobs..."

MODKIT_BIN="/athena/masonlab/scratch/users/thn4005/programs/dist_modkit_v0.5.0_5120ef7/modkit"
MODS=("17596" "a" "m" "17802" "69426" "19228" "19229" "19227")
REF_BASES=("A" "A" "C" "T" "A" "C" "G" "T")  # Matching reference bases

for i in "${!MODS[@]}"; do
    MOD="${MODS[$i]}"
    REF="${REF_BASES[$i]}"
    JOB_SCRIPT="${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_${MOD}_modkit.sbatch"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --job-name=modkit_${SAMPLE_NAME}_${MOD}

source ~/.bashrc

echo "[INFO] Running modkit pileup and summary for ${MOD}..."

${MODKIT_BIN} pileup ${SORTED_BAM} ${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_${MOD}_filtered_mod.bed \\
    --ref ${REFERENCE_GENOME} \\
    --filter-threshold ${REF}:.8 \\
    --mod-thresholds ${MOD}:0.99 \\
    --log-filepath ${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_${MOD}_filtered_mod_bed_log \\
    --max-depth 1000 \\
    -t 64 \\
    --with-header \\
    --interval-size 100000

${MODKIT_BIN} summary ${SORTED_BAM} \\
    --filter-threshold ${REF}:.8 \\
    --mod-thresholds ${MOD}:0.99 \\
    > ${SAMPLE_OUTPUT_FOLDER}/${SAMPLE_NAME}_${MOD}_filtered_mod_summary.txt
EOF

    sbatch "$JOB_SCRIPT"
done

