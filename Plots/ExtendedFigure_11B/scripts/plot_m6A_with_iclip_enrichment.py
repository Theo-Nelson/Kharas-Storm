#!/usr/bin/env python3
"""
Supplementary Figure 9C (revision): original nanopore m6A panels + YTHDC1 iCLIP enrichment.

Reproduces the published nanopore panels EXACTLY by reusing the original pipeline module
(plot_m6A_fraction_points_with_coverage.py) for data prep, colors, gene model, and ticks,
then adds an aligned iCLIP enrichment sub-panel (DC1 DMSO vs DC1 STM) beneath the coverage
track, sharing the genomic x-axis.

Nanopore section (Fraction m6A dots + gene model + coverage) is byte-for-byte the same
rendering as the original; only a new bottom panel is added.
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import plot_m6A_fraction_points_with_coverage as P

# ---- iCLIP enrichment window columns in the data table ----
WINDOWS = {
    "binding": ("binding_minus4_plus5_dc1_dmso_mean_enrichment",
                "binding_minus4_plus5_dc1_stm_mean_enrichment",
                "iCLIP enrichment\n(−4..+5 window)"),
    "minus4":  ("minus4_dc1_dmso_mean_enrichment",
                "minus4_dc1_stm_mean_enrichment",
                "iCLIP enrichment\n(−4 window)"),
    "plus1_plus5": ("plus1_plus5_dc1_dmso_mean_enrichment",
                    "plus1_plus5_dc1_stm_mean_enrichment",
                    "iCLIP enrichment\n(+1..+5 window)"),
}

GROUP_ORDER = ["GFP DMSO", "GFP STM", "YTHDC1 DMSO", "YTHDC1 STM"]

RETAINED_COLOR = "#6a3d9a"   # purple  - retained m6A sites
NONRET_COLOR   = "#9e9e9e"   # grey    - non-retained m6A sites


def load_enrichment(tsv, gene):
    df = pd.read_csv(tsv, sep="\t")
    g = df[df["gene_name"] == gene].copy()
    # 1-based genomic position of the m6A A (matches nanopore pos = chromStart+1 = end)
    g["pos1"] = g["start"].astype(int) + 1
    return g


def draw_enrichment(ax, enr, start_bp, end_bp, win):
    """Plot DC1 DMSO (filled) vs DC1 STM (open) enrichment for one window per m6A site.
    Retained sites in purple, non-retained in grey. Baseline (=1) dashed."""
    sub = enr[(enr["pos1"] >= start_bp) & (enr["pos1"] <= end_bp)]
    dcol, scol, label = WINDOWS[win]
    ymax = 1.0
    for _, r in sub.iterrows():
        x = r["pos1"]
        col = RETAINED_COLOR if bool(r["retained"]) else NONRET_COLOR
        ydm, yst = float(r[dcol]), float(r[scol])
        ymax = max(ymax, ydm, yst)
        # connector shows the DMSO -> STM shift at this site
        ax.plot([x, x], [ydm, yst], color=col, lw=0.6, alpha=0.55, zorder=1)
        ax.scatter([x], [ydm], s=26, facecolors=col, edgecolors=col,
                   linewidths=0.8, alpha=0.95, zorder=3)          # DMSO = filled
        ax.scatter([x], [yst], s=26, facecolors="white", edgecolors=col,
                   linewidths=1.0, alpha=0.95, zorder=3)          # STM  = open
    ax.axhline(1.0, ls=":", lw=1.0, color="0.45", zorder=0)       # no-enrichment baseline
    ax.set_ylim(0, ymax * 1.12)
    ax.set_ylabel(label, fontsize=8)
    ax.yaxis.set_label_coords(-0.058, 0.5)


def enrichment_legend(ax):
    handles = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor=RETAINED_COLOR,
               markeredgecolor=RETAINED_COLOR, markersize=7, label='Retained m6A'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=NONRET_COLOR,
               markeredgecolor=NONRET_COLOR, markersize=7, label='Non-retained m6A'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='0.25',
               markeredgecolor='0.25', markersize=7, label='DC1 Vehicle (filled)'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='white',
               markeredgecolor='0.25', markersize=7, label='DC1 STM2457 (open)'),
    ]
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)


# ---- per-block sub-panel heights (inches) ----
H_MAIN, H_EXON, H_COV, H_ENR = 2.00, 0.55, 1.35, 1.15
G_ME, G_EC, G_CE, G_EE = 0.03, 0.30, 0.28, 0.16   # inner gaps
BOTTOM_LABELS = 0.95                              # genomic tick labels under each block
TOP_TITLE = 0.35
INTER_BLOCK = 0.55
X_IN, W_IN = 1.25, 8.6                            # plot x offset / width (inches)


def block_height(n_enr):
    return (H_MAIN + G_ME + H_EXON + G_EC + H_COV + G_CE
            + n_enr * H_ENR + (n_enr - 1) * G_EE)


def build_context(pass_df, enr, gene, chrom, start_bp, end_bp,
                  exon_1b, utr5_1b, cds_1b, utr3_1b, title):
    region = pass_df[(pass_df["chrom"] == chrom) &
                     (pass_df["pos"] >= start_bp) & (pass_df["pos"] <= end_bp)].copy()
    region["total"] = region["modified"] + region["unmodified"]
    region = region[region["total"] > 0]
    region["frac"] = np.clip(region["modified"] / region["total"], 0.0, 1.0)
    cov = (P._COV[(P._COV["chrom"] == chrom) &
                  (P._COV["pos"] >= start_bp) & (P._COV["pos"] <= end_bp)].sort_values("pos"))
    return dict(region=region, enr=enr, gene=gene, chrom=chrom,
                start=start_bp, end=end_bp,
                x_cov=cov["pos"].to_numpy(), y_cov=cov["total_called_bases"].to_numpy(),
                exon=exon_1b, utr5=utr5_1b, cds=cds_1b, utr3=utr3_1b, title=title)


def draw_block(fig, W, H, y_top_in, ctx, windows, point_size, max_xticks,
               show_group_legend=True, show_enr_legend=True):
    """Draw one gene's 5-panel block into the figure. y_top_in = inches from TOP."""
    def rect(y_top, h):                    # y_top in inches-from-top -> axes rect
        y0 = H - (y_top + h)
        return [X_IN / W, y0 / H, W_IN / W, h / H]

    y = y_top_in
    ax_main = fig.add_axes(rect(y, H_MAIN));                       y += H_MAIN + G_ME
    ax_exon = fig.add_axes(rect(y, H_EXON), sharex=ax_main);       y += H_EXON + G_EC
    ax_cov  = fig.add_axes(rect(y, H_COV),  sharex=ax_main);       y += H_COV + G_CE
    enr_axes = []
    for k in range(len(windows)):
        enr_axes.append(fig.add_axes(rect(y, H_ENR), sharex=ax_main))
        y += H_ENR + (G_EE if k < len(windows) - 1 else 0)

    for (sample, grp), df_sg in ctx["region"].groupby(["sample", "group"], sort=False):
        color = P.GROUP_COLORS.get(grp, P.GROUP_COLORS["Unlabeled"])
        ax_main.scatter(df_sg["pos"].to_numpy(), df_sg["frac"].to_numpy(),
                        s=point_size, alpha=0.90, edgecolors="none", color=color)
    ax_main.set_ylim(0.0, 1.10)
    ax_main.set_ylabel("Fraction m6A")
    ax_main.set_title(f"{ctx['gene']} — {ctx['title']}", fontsize=11)
    ax_main.tick_params(axis='x', labelbottom=False)
    if show_group_legend:
        present = set(ctx["region"]["group"].unique())
        handles = [Line2D([0], [0], marker='o', linestyle='none',
                          markerfacecolor=P.GROUP_COLORS[g], markeredgecolor='none',
                          markersize=6, label=g)
                   for g in GROUP_ORDER if g in present]
        ax_main.legend(handles=handles, ncol=1, fontsize=8.5, frameon=False,
                       loc="upper left", bbox_to_anchor=(1.015, 1.0), borderaxespad=0.0)

    P.draw_exon_and_regions(ax_exon, ctx["exon"], ctx["utr5"], ctx["cds"], ctx["utr3"])
    ax_exon.tick_params(axis='x', labelbottom=False)

    ax_cov.plot(ctx["x_cov"], ctx["y_cov"], linewidth=1.0, color="0.25")
    ax_cov.set_ylabel("Number of\nCalled Bases", fontsize=9)
    ax_cov.set_ylim(bottom=0)
    ax_cov.tick_params(axis='x', labelbottom=False)

    for i, (ax_e, win) in enumerate(zip(enr_axes, windows)):
        draw_enrichment(ax_e, ctx["enr"], ctx["start"], ctx["end"], win)
        if i == 0 and show_enr_legend:
            enrichment_legend(ax_e)
        if i == len(enr_axes) - 1:
            ax_e.set_xlabel("Genomic position")
            P.set_genomic_ticks(ax_e, ctx["chrom"], ctx["start"], ctx["end"], max_xticks)
        else:
            ax_e.tick_params(axis='x', labelbottom=False)
    # ensure shared x-limits are exactly the region
    ax_main.set_xlim(ctx["start"], ctx["end"])


def save_single(ctx, windows, outdir, point_size=8.0, max_xticks=8):
    n = len(windows)
    H = TOP_TITLE + block_height(n) + BOTTOM_LABELS
    W = 12.0
    fig = plt.figure(figsize=(W, H))
    draw_block(fig, W, H, TOP_TITLE, ctx, windows, point_size, max_xticks)
    tag = f"{ctx['gene']}_{ctx['title'].replace(' ','_').replace(':','-').replace('(','').replace(')','')}_iclip"
    out = os.path.join(outdir, tag)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out + "." + ext, dpi=300)
    plt.close(fig)
    print(f"[OK] {tag}")


def save_combined(contexts, windows, outdir, name="SuppFig9C_iclip_combined",
                  point_size=8.0, max_xticks=8):
    n = len(windows)
    bh = block_height(n)
    per = TOP_TITLE + bh + BOTTOM_LABELS
    W = 12.0
    H = len(contexts) * per + (len(contexts) - 1) * INTER_BLOCK
    fig = plt.figure(figsize=(W, H))
    y = 0.0
    for ctx in contexts:
        draw_block(fig, W, H, y + TOP_TITLE, ctx, windows, point_size, max_xticks)
        y += per + INTER_BLOCK
    out = os.path.join(outdir, name)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out + "." + ext, dpi=200)
    plt.close(fig)
    print(f"[OK] combined -> {out}.[png/pdf/svg]  ({len(contexts)} blocks, {H:.1f}in tall)")


def gene_regions(pass_df, gtf, gene, per_exon):
    """Yield (start, end, exon_1b, utr5_1b, cds_1b, utr3_1b, title, chrom) tuples."""
    feats = P.parse_gtf_features_for_gene(gtf, gene)
    tx_id, tx = P.pick_best_transcript(feats, preferred_tx=None)
    if not tx:
        print(f"[WARN] no transcript for {gene}"); return
    chrom = tx["chrom"]
    bed_chroms = set(pass_df["chrom"].unique())
    target = chrom if chrom in bed_chroms else (chrom.replace("chr","") if chrom.startswith("chr") else "chr"+chrom)
    exons_0b = P.merge_intervals([(s, e) for (s, e, _) in tx["exons"]])
    gstart = min(s for s, _ in exons_0b) + 1
    gend = max(e for _, e in exons_0b)
    exons_1b = P.intervals_to_1based(exons_0b)
    utr5 = P.intervals_to_1based(P.merge_intervals(tx["utr5"]))
    cds  = P.intervals_to_1based(P.merge_intervals(tx["cds"]))
    utr3 = P.intervals_to_1based(P.merge_intervals(tx["utr3"]))
    yield (target, gstart, gend, exons_1b, utr5, cds, utr3, f"{tx_id} (all exons)")
    if per_exon:
        for i, (s, e, exnum) in enumerate(sorted(tx["exons"], key=lambda x: (x[0], x[1]))):
            s1, e1 = s + 1, e
            def clip(ivls):
                out = []
                for a, b in P.merge_intervals(ivls):
                    aa, bb = max(a, s), min(b, e)
                    if aa < bb: out.append((aa, bb))
                return P.intervals_to_1based(out)
            yield (target, s1, e1, [(s1, e1)], clip(tx["utr5"]), clip(tx["cds"]),
                   clip(tx["utr3"]), f"{tx_id} exon {exnum if exnum is not None else i+1}")


def prep_pipeline(gtf, input_dir, glob_pat, z_thresh=3.0, min_cov=20, cov_agg="max"):
    raw_df, pass_df = P.read_beds_raw_and_passonly(input_dir, glob_pat)
    cov_df = P.build_coverage_selectedmods(input_dir, glob_pat, cov_agg)
    ztable = P.compute_z_from_raw(raw_df)
    pass_df = pass_df.rename(columns={"score": "score_pass"}).merge(
        ztable, on=["chrom", "pos", "sample"], how="left")
    z_nonan = ztable.dropna(subset=["z"])
    zpos = z_nonan.groupby(["chrom", "pos"])["z"].max().reset_index(name="zmax")
    zpos = zpos[zpos["zmax"] > z_thresh][["chrom", "pos"]]
    pass_df["total"] = pass_df["modified"] + pass_df["unmodified"]
    pass_df = pass_df.merge(zpos, on=["chrom", "pos"], how="inner")
    pass_df = pass_df[pass_df["total"] > min_cov].copy()
    P._COV = cov_df
    return pass_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--input-dir", default="pileup_out")
    ap.add_argument("--glob", default="*/*_filtered_mod.bed*")
    ap.add_argument("--enrichment-tsv", required=True)
    ap.add_argument("--genes", nargs="+", default=["BCL2", "MYC"])
    ap.add_argument("--windows", nargs="+", default=["binding", "plus1_plus5"],
                    choices=list(WINDOWS.keys()))
    ap.add_argument("--per-exon", action="store_true")
    ap.add_argument("--regions", nargs="+", default=None,
                    help="Ordered block spec for combined fig, e.g. BCL2:exon1 BCL2:exon2 MYC:all SP1:all")
    ap.add_argument("--combined", action="store_true", help="Also emit one stacked combined figure")
    ap.add_argument("--combined-name", default="SuppFig9C_iclip_combined")
    ap.add_argument("--outdir", default="iclip_out")
    ap.add_argument("--point-size", type=float, default=8.0)
    ap.add_argument("--max-xticks", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pass_df = prep_pipeline(args.gtf, args.input_dir, args.glob)

    # collect all region contexts per gene
    contexts_by_key = {}
    for gene in args.genes:
        enr = load_enrichment(args.enrichment_tsv, gene)
        for (chrom, s, e, ex1b, u5, cds, u3, title) in gene_regions(pass_df, args.gtf, gene, True):
            ctx = build_context(pass_df, enr, gene, chrom, s, e, ex1b, u5, cds, u3, title)
            key_all = f"{gene}:all"
            if "(all exons)" in title:
                contexts_by_key[key_all] = ctx
            else:
                exnum = title.split("exon")[-1].strip()
                contexts_by_key[f"{gene}:exon{exnum}"] = ctx
            # single-panel outputs (all exons always; per exon if requested)
            if ("(all exons)" in title) or args.per_exon:
                save_single(ctx, args.windows, args.outdir, args.point_size, args.max_xticks)

    if args.combined:
        if args.regions:
            order = args.regions
        else:
            order = [k for k in contexts_by_key if k.endswith(":all")]
        ctxs = [contexts_by_key[k] for k in order if k in contexts_by_key]
        missing = [k for k in order if k not in contexts_by_key]
        if missing:
            print(f"[WARN] combined: missing region keys {missing}; available={list(contexts_by_key)}")
        if ctxs:
            save_combined(ctxs, args.windows, args.outdir, args.combined_name,
                          args.point_size, args.max_xticks)


if __name__ == "__main__":
    main()
