#!/bin/bash
#SBATCH --partition=scu-cpu
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --job-name=modkit_16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5_17596

source ~/.bashrc

echo "[INFO] Running modkit pileup and summary for 17596..."

/athena/masonlab/scratch/users/thn4005/programs/dist_modkit_v0.5.0_5120ef7/modkit pileup /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5_aligned.sorted.bam /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5_17596_filtered_mod.bed \
    --ref /athena/masonlab/scratch/users/thn4005/scripts/references_v2/references/human/hg38.fa \
    --filter-threshold A:.8 \
    --mod-thresholds 17596:0.99 \
    --log-filepath /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5_17596_filtered_mod_bed_log \
    --max-depth 1000 \
    -t 64 \
    --with-header \
    --interval-size 100000

/athena/masonlab/scratch/users/thn4005/programs/dist_modkit_v0.5.0_5120ef7/modkit summary /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5_aligned.sorted.bam \
    --filter-threshold A:.8 \
    --mod-thresholds 17596:0.99 \
    > /athena/masonlab/scratch/users/thn4005/scripts/projects_v2/basecalled_samples/v2/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5/16-Jun-25_XuejingYang_RNADirect_Sample3_no_sample_id_20250616_1414_1C_PBE87494_3ba06475_pod5_17596_filtered_mod_summary.txt
