# Supplementary Figure 9C (revision) — nanopore m⁶A + YTHDC1 iCLIP enrichment

Revision panels that keep the **original nanopore section unchanged** and add an aligned
**YTHDC1 iCLIP enrichment** track below each gene, addressing the reviewers' questions about
YTHDC1 binding at retained m⁶A sites and under METTL3 inhibition (STM2457).

## Files
- `SuppFig9C_combined.{svg,pdf,png}` — one stacked figure: BCL2 exon 1, BCL2 exon 2, MYC, SP1.
- `panels/SuppFig9C_{BCL2_exon1,BCL2_exon2,MYC,SP1}_iclip.{svg,pdf,png}` — individual gene panels.

## Panel layout (per gene, top → bottom)
1. **Fraction m⁶A** — per-sample nanopore dots, 4 groups (GFP Vehicle/STM2457, GFP-YTHDC1 Vehicle/STM2457).
2. Gene model (exons navy; 5′UTR grey, CDS green, 3′UTR brown).
3. **Number of Called Bases** — nanopore coverage track.
4. **iCLIP enrichment (−4..+5 window)** — YTHDC1 binding window.
5. **iCLIP enrichment (+1..+5 window)**.

Panels 1–3 are the **exact original published rendering** (same pipeline/script/data). Panels 4–5 are new.

## iCLIP enrichment encoding
Per m⁶A site: **DC1 (GFP-YTHDC1) Vehicle = filled**, **DC1 STM2457 = open**, joined by a thin
line showing the Vehicle→STM shift. **Retained m⁶A sites = purple, non-retained = grey.**
Dashed line at enrichment = 1 (no enrichment). Values are fold-enrichment ratios from
`high_confidence_MYC_BCL2_SP1_enrichment_plus_percent_modified_heatmap_values.tsv`
(columns `binding_minus4_plus5_dc1_{dmso,stm}_mean_enrichment` and
`plus1_plus5_dc1_{dmso,stm}_mean_enrichment`).

## How it was generated
All three genes come from a **single unified modkit pileup run**, so BCL2, MYC and SP1 are mutually
consistent (`scripts/run_pileup_all_genes.sh`):

- **modkit 0.5.0** (bioconda osx-arm64 build; the bundled binary is Linux-x86_64 and the original
  pipeline used 0.5.1), run on the `Kharas_Subset_MYC_BCL2_SP1/` BAMs.
- Thresholds match `run_pileup_BCL2_MYC.sh`: `--filter-threshold A:0.8 --mod-thresholds a:0.99`,
  `--max-depth 1000`, `--interval-size 100000`.
- Regions from `hg38.ncbiRefSeq_V2.gtf`: BCL2 `chr18:63123346-63320090`,
  MYC `chr8:127735434-127742951`, SP1 `chr12:53380176-53416446`.
- Reference: chr8+chr12+chr18 extracted from `hg38.fa.gz` (`ref_chr8_12_18.fa`, see `.fai`).
- Panels drawn by `plot_m6A_with_iclip_enrichment.py`, which reuses the original
  `plot_m6A_fraction_points_with_coverage.py` for the nanopore panels.

### Validation vs. the data table (`*_mean_percent_modified` columns)
All 66 sites matched on genomic position (what the enrichment-track alignment depends on).
Group-mean %-modified agreement using the unified run:

| gene | comparisons | max diff | notes |
|------|------------:|---------:|-------|
| BCL2 |         180 | 0.59 pp  | 178/180 exact; one outlier site (below) |
| MYC  |          40 | 4.3e-14  | exact (machine precision) |
| SP1  |          44 | 4.3e-14  | exact (machine precision) |

MYC and SP1 reproduce the table to machine precision, confirming the reference, thresholds,
sample→group mapping and coordinate convention are correct.

**Effect of unifying the run.** Versus the previously supplied `pileup_out/`: BCL2 panels are
**pixel-identical** (0 pixels changed — the published BCL2 rendering is untouched); MYC shifted
slightly (0.37% of pixels) and now matches the data table exactly, where the supplied pileup was off
by up to 1.04 pp; SP1 unchanged.

**Known outlier — BCL2 `m6A_198596` (chr18:63,319,841).** The table's two *DMSO* values at this one
site disagree with the pileup (YTHDC1 DMSO 53.85 vs 53.26; GFP DMSO 50.90 vs 50.393); the STM groups
match exactly. Both the supplied `pileup_out/` and this fresh run give the same value, so the **table
is the outlier** at this site, not the pileup. Does not affect the figure — the table is used only for
enrichment values and retained/non-retained status; its %-modified columns are never plotted.

**Scope of these checks.** They validate *stoichiometry only*. The **coverage tracks** have no
independent reference (the table has no coverage column), and the **per-replicate spread** is
unvalidated (the table stores only group means).
- **Scripts:** `plot_m6A_with_iclip_enrichment.py` (reuses the original module for the nanopore
  panels and adds the enrichment panels + combined layout). `run_sp1_pileup.sh` documents the SP1 pileup.

## Note on the enrichment normalization
The precise definition of the enrichment ratio (e.g. IP/input per window) is not in this repo — it
was inferred from the column names and value distribution (centered near 1). Confirm the exact
definition from the iCLIP analysis for the figure legend / methods text.
