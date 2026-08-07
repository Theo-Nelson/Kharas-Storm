# Direct RNA Analysis to Examine METTL3 Inhibitor Resistance

This repository contains pipelines for alignment, base modification calling, and differential expression analysis of **Nanopore direct RNA sequencing data** from THP-1 leukemia cells treated with METTL3 inhibitors and/or YTHDC1 overexpression.

**Figure 1. Workflow for Nanopore Direct-RNA Sequencing Analysis in Kharas-Storm Repository.**  
![KharasLab_FlowChart_GitHub](KharasLab_FlowChart_GitHub.png)

All pipelines were developed in the **Mason Laboratory (WCM)**.

---

## Study Overview

Nanopore direct RNA sequencing was used to characterize **m⁶A and other base modifications** in THP-1 cell lines under four experimental conditions:

| Condition | Replicates | Description |
|------------|-------------|-------------|
| THP1 GFP DMSO 48 h | 3 | Control (untreated) |
| THP1 GFP STM2457 15 µM 48 h | 3 | METTL3 inhibitor treatment |
| THP1 GFP-YTHDC1 DMSO 48 h | 3 | YTHDC1 overexpression |
| THP1 GFP-YTHDC1 STM2457 15 µM 48 h | 3 | Combined YTHDC1 + METTL3i |

Sequencing was performed on **Oxford Nanopore PromethION** (Study ERP182419, ENA project PRJEB100992, submission 2025-10-29). All analyses use the same naming convention as the raw sample metadata (as documented in `THP1_GFP-YTHDC1 OE_STM_48H_05202025.xlsx`).

---

## Repository Structure

```
Kharas-Storm/
├── Alignments_Pipeline/                # Alignment and Modkit Summarization of multi-mod Dorado-basecalled reads
│   ├── mods-alignment-from-folder.sbatch
│   └── Modkit_Scripts_Logs/
│       └── [sample folders + sbatch + logs]
│
├── Mod_Compilation_Pipeline/           # Modkit post-processing and z-score filtering
│   ├── run_full_pipeline.sbatch
│   ├── scripts/
│   │   ├── 01_filter_beds.py
│   │   ├── 02_generate_union_beds.sbatch
│   │   ├── 03_intersect_with_union.sbatch
│   │   ├── 04_merge_to_matrix.py
│   │   ├── 05_filter_by_zscore.py
│   │   ├── 06_summarize_pipeline_stats.py
│   │   ├── 07_annotate_outputs_with_genes.sh
│   │   └── filter_bed_sites.py
│   └── outputs/
│       ├── filtered_beds/
│       ├── union_beds/
│       ├── merged_matrix/
│       ├── zscore_filtered/
│       └── annotated TSVs per modification
│
├── DGE_Analysis_Pipeline/              # SARTools + DESeq2 RNA-seq differential analysis
│   ├── Kharas_SarTools_DESeq2_Analysis/
│   │   ├── target_kharas.txt
│   │   ├── raw/ (featureCounts output)
│   │   └── DESeq2 for Kharas THP1 DC1_Report.html
│   └── SARToolsScript.R
│
├── Plots/                              # Figure generation scripts and rendered panels
│   ├── Figure5B_and_ExtendedFigure9A_m6a_distribution/
│   │   ├── new_Fig1.py / new_Fig1_V2.py
│   │   └── Plots/[per-mod-code bars, ecdf, upset]
│   ├── ExtendedFigure_11B/             # Revision panel: nanopore m⁶A + YTHDC1 iCLIP enrichment
│   │   ├── README.md                   # Panel-by-panel description + validation notes
│   │   ├── SuppFig9C_combined.{svg,pdf,png}
│   │   ├── panels/                     # BCL2 exon 1, BCL2 exon 2, MYC, SP1 (svg/pdf/png)
│   │   └── scripts/
│   │       ├── run_pileup_all_genes.sh
│   │       ├── plot_m6A_with_iclip_enrichment.py
│   │       └── ref_chr8_12_18.fa.fai
│   └── OLD_ExtendedFigure9C_MYC_BCL2_Dot_Plot/   # Original (pre-revision) MYC/BCL2 dot plots
│       ├── run_pileup_BCL2_MYC.sh
│       ├── plot_m6A_fraction_points_with_coverage.py
│       ├── m6A_plots_points_cov_exon/
│       └── data/ + dist_modkit_v0.5.1_8fa79e3/   # superseded copies, now in Data/ and Software/
│
├── Data/                               # Gene-subset inputs for the figure scripts
│   ├── Kharas_Subset_MYC_BCL2_SP1/     # 12 MYC/BCL2/SP1-subset BAMs + .bai
│   │   └── pileup_out_v2/              # Unified modkit pileup (BCL2, MYC, SP1 per sample)
│   └── high_confidence_MYC_BCL2_SP1_enrichment_plus_percent_modified_heatmap_values.tsv
│
└── Software/                           # Vendored analysis binaries
    └── dist_modkit_v0.5.1_8fa79e3/     # ONT modkit 0.5.1 (Linux x86-64) + docs + models
```

---

## Pipeline Summaries

### 1. **Alignments Pipeline**

Automates alignment of **basecalled multi-mod BAMs** to the reference genome and submits per-modification `modkit` jobs.

**Key script:**  
`mods-alignment-from-folder.sbatch`

**Steps**
1. Convert `calls_*.bam` → FASTQ via `samtools fastq`
2. Align to reference using `minimap2 -ax splice -uf -k14`
3. Sort + index BAMs with `samtools`
4. Run `modkit pileup` + `modkit summary` per modification code  
   (A/m⁶A → `a`, m⁵C → `m`, inosine/m⁶A → `17596`, etc.)

**Usage**
```bash
sbatch mods-alignment-from-folder.sbatch <reference_genome.fa> <input_folder>
```

---

### 2. **Mod Compilation Pipeline**

Aggregates, filters, and annotates **per-base modification calls** across samples.

**Primary driver:**  
`run_full_pipeline.sbatch`

**Major stages**
1. **Filter BEDs:** apply coverage + stoichiometry criteria (`filter_bed_sites.py`)  
2. **Generate union sites:** merge all samples per modification  
3. **Intersect union × sample BEDs:** retain shared sites with attributes  
4. **Merge to long/wide matrices:** `04_merge_to_matrix.py`  
5. **Compute z-scores:** per-site, per-sample modification ratios (`05_filter_by_zscore.py`)  
6. **Summarize + annotate:** site statistics and gene mapping via `bedtools intersect`

**Run example**
```bash
sbatch run_full_pipeline.sbatch /path/to/project_root sample_list.txt
```

**Output Highlights**
- `*_merged_union_matrix_wide__with_z.tsv`  
  — per-site wide matrix with z-scores across all samples  
- `*_z_gt3.0_annotated.tsv`  
  — filtered sites where any z > 3  
- `stats/*.txt` — summary counts per mod type  

---

### 3. **Differential Gene Expression (DGE) Pipeline**

Performs DESeq2 analysis of gene-level counts (e.g. featureCounts) with **SARTools** reporting.

**Script:**  
`SARToolsScript.R`

**Key Parameters**
```r
projectName <- "DESeq2 for Kharas THP1 DC1"
condRef <- "THP1GFPDMSO48h"
alpha <- 0.05
pAdjustMethod <- "BH"
typeTrans <- "VST"
```

**Run in R**
```r
source("SARToolsScript.R")
```

Generates a full HTML report with QC plots, PCA, MA plots, and differential gene tables.

---

### 4. **Figure Generation (`Plots/`)**

#### `Figure5B_and_ExtendedFigure9A_m6a_distribution/`
Distribution plots (bar charts, ECDFs, UpSet plots) of modified sites per modification code
(`a`, `m`, `17596`, `17802`, `19227`, `19228`, `19229`, `69426`), generated by `new_Fig1.py` /
`new_Fig1_V2.py` from the Mod Compilation Pipeline outputs.

#### `ExtendedFigure_11B/` — nanopore m⁶A + YTHDC1 iCLIP enrichment
Revision panels for **Extended Figure 11B** (files retain the earlier `SuppFig9C` naming). Each
gene panel keeps the **original nanopore rendering unchanged** and adds aligned **YTHDC1 iCLIP
enrichment** tracks below it, addressing YTHDC1 binding at retained m⁶A sites under METTL3
inhibition (STM2457).

Per gene, top → bottom: fraction m⁶A (per-sample dots, 4 groups) · gene model · number of called
bases (coverage) · iCLIP enrichment (−4..+5 window) · iCLIP enrichment (+1..+5 window). Filled
points = GFP-YTHDC1 Vehicle, open points = GFP-YTHDC1 STM2457; **purple = retained m⁶A sites,
grey = non-retained**.

**Scripts**
- `scripts/run_pileup_all_genes.sh` — single unified `modkit pileup` over BCL2, MYC and SP1 so all
  three genes are mutually consistent. Thresholds match the original run
  (`--filter-threshold A:0.8 --mod-thresholds a:0.99 --max-depth 1000 --interval-size 100000`).
  Regions from `hg38.ncbiRefSeq_V2.gtf`: BCL2 `chr18:63123346-63320090`,
  MYC `chr8:127735434-127742951`, SP1 `chr12:53380176-53416446`.
- `scripts/plot_m6A_with_iclip_enrichment.py` — draws the panels, reusing the original
  `plot_m6A_fraction_points_with_coverage.py` module for the nanopore section.

```bash
BAMDIR=Data/Kharas_Subset_MYC_BCL2_SP1 REF=ref_chr8_12_18.fa OUT=pileup_out_v2 \
  bash Plots/ExtendedFigure_11B/scripts/run_pileup_all_genes.sh
```

The reference (`ref_chr8_12_18.fa`, chr8 + chr12 + chr18 extracted from `hg38.fa.gz`) is not
committed — only its `.fai` is, so the exact contig set and lengths can be verified.

Validation details (agreement with the enrichment table, the effect of unifying the pileup run, and
one known outlier site) are documented in `Plots/ExtendedFigure_11B/README.md`.

#### `OLD_ExtendedFigure9C_MYC_BCL2_Dot_Plot/`
The original, pre-revision MYC/BCL2 dot plots (`run_pileup_BCL2_MYC.sh`,
`plot_m6A_fraction_points_with_coverage.py`, rendered panels in `m6A_plots_points_cov_exon/`),
retained for provenance. The plotting module here is imported directly by the Extended Figure 11B
script.

---

## Data and Software Directories

### `Data/`
Gene-subset inputs used by the figure scripts (the full alignments are in ENA — see
[Data Availability](#data-availability)):

- `Kharas_Subset_MYC_BCL2_SP1/` — 12 per-sample BAMs (+ `.bai`) subset to the MYC, BCL2 and SP1 loci.
- `Kharas_Subset_MYC_BCL2_SP1/pileup_out_v2/` — the unified `modkit pileup` output, one folder per
  sample with `*_{BCL2,MYC,SP1}_filtered_mod.bed` and matching `_log.txt`.
- `high_confidence_MYC_BCL2_SP1_enrichment_plus_percent_modified_heatmap_values.tsv` — 66
  high-confidence m⁶A sites with per-group mean % modified and YTHDC1 iCLIP fold-enrichment for
  three windows (`minus4`, `plus1_plus5`, `binding_minus4_plus5`), plus retained/non-retained status.
  This table supplies the enrichment values plotted in Extended Figure 11B.

### `Software/`
`dist_modkit_v0.5.1_8fa79e3/` — the vendored ONT **modkit 0.5.1** distribution (binary, docs, and
the `r1041_e82_400bps_hac_v5.2.0@v0.1.0` model) used by the alignment and compilation pipelines.
Note that the bundled binary is
**Linux x86-64**; the Extended Figure 11B pileups were re-run on macOS with the bioconda
`ont-modkit` 0.5.0 build using identical thresholds.

---

## Dependencies

Software versions used for both alignment and modification pipelines:

- `samtools=1.21` (bioconda), `htslib=1.21`
- `bedtools=2.31.1` (bioconda)
- `python=3.12`
- `pandas`, `numpy` (>=2.3)
- `minimap2` (>=2.28; separate env recommended)
- `modkit=0.5.x` (ONT prebuilt binary, vendored in `Software/dist_modkit_v0.5.1_8fa79e3/`)
- `matplotlib` (figure scripts in `Plots/`)
- plus: `bedtools`, `curl`, `graphviz`, `imagemagick` (for QA/report figures)

All pipelines are SLURM-based and were run on the **Cornell Cayuga Cluster** nodes.

---

## Data Availability

Raw direct RNA sequencing data are deposited under:

| Accession | Instrument | Study | Submission Date | 
|------------|-------------|--------|----------------|
| **PRJEB100992 / ERP182419** | ONT PromethION | “Nuclear Condensates Drive Adaptive Resistance to METTL3 Inhibition” | 2025‑10‑29 |

Sample metadata (RNA concentrations and QC ratios) for 12 THP-1 libraries analyzed here are provided in `THP1_GFP-YTHDC1 OE_STM_48H_05202025.xlsx`.

Locus-subset BAMs (MYC, BCL2, SP1) and their `modkit` pileups are included in `Data/` so the
Extended Figure 11B panels can be regenerated without downloading the full alignments.

Associated reference genomes used within this analysis are available here: https://zenodo.org/records/21842348

---

## Citation

**Citation TBD**  

---

For questions, please contact **Theo Nelson** (thn4005@med.cornell.edu).


