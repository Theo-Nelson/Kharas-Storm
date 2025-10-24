#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --job-name=modkit_16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5_m

source ~/.bashrc

echo "[INFO] Running modkit pileup and summary for m..."

/athena/masonlab/scratch/users/thn4005/programs/dist_modkit_v0.5.0_5120ef7/modkit pileup /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5_aligned.sorted.bam /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5_m_filtered_mod.bed \
    --ref /athena/masonlab/scratch/users/thn4005/scripts/references_v2/references/human/hg38.fa \
    --filter-threshold C:.8 \
    --mod-thresholds m:0.99 \
    --log-filepath /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5_m_filtered_mod_bed_log \
    --max-depth 1000 \
    -t 64 \
    --with-header \
    --interval-size 100000

/athena/masonlab/scratch/users/thn4005/programs/dist_modkit_v0.5.0_5120ef7/modkit summary /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5_aligned.sorted.bam \
    --filter-threshold C:.8 \
    --mod-thresholds m:0.99 \
    > /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5/16-Jun-25_XuejingYang_RNADirect_Sample6_no_sample_id_20250616_1415_1F_PBE87531_88e60ff5_pod5_m_filtered_mod_summary.txt
