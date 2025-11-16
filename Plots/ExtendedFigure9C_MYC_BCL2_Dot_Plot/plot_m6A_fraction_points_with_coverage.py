#!/usr/bin/env python3
"""
m6A plots (pass-only + z-from-RAW + coverage), with UTR/CDS bands and exports.

Key behavior (this revision):
- Saves PNG, PDF, and SVG.
- Position indexing uses chromStart (BED 0-based); we store/plot as 1-based via pos = chromStart + 1.
- Group-mean % labels removed (hollow mean markers optional via --label-group-means).
- Coverage track ("Number of Called Bases") is built from RAW 'score' for ONLY mod codes
  [19227, 19228, 19229, 69426].
- Coverage is first computed PER-SAMPLE, then aggregated across samples using --cov-agg
  (max | mean | sum). Default = max (prevents inflation from multi-sample summation).
- De-duplication guard: within each sample, for (chrom,pos,mod_code) take MAX(score)
  before per-sample summation (avoids double counting if shards/dupes exist).
- 5′UTR band is GREY.

Outputs (after ALL filters):
  passonly_sites.tsv  (per-sample, with z from RAW, mod/fail ratios from RAW, region)
  group_means.tsv
  summary_stats.txt
"""

import argparse, os, glob, re, sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -------- bedMethyl columns (modkit bedMethyl) --------
BED_COLS = [
    "chrom","start","end","mod_code","score","strand","start_compat","end_compat",
    "rgb","valid_cov","frac_mod","n_mod","n_can","n_other","n_del","n_fail","n_diff","n_nocall"
]

# -------- Sample grouping (by number in path) --------
NUMBER_TO_GROUP = {
    1:"GFP DMSO",2:"GFP DMSO",3:"GFP DMSO",
    4:"GFP STM",5:"GFP STM",6:"GFP STM",
    7:"YTHDC1 DMSO",8:"YTHDC1 DMSO",9:"YTHDC1 DMSO",
    10:"YTHDC1 STM",11:"YTHDC1 STM",12:"YTHDC1 STM",
}
GROUP_COLORS = {
    "GFP DMSO":"#1f77b4","GFP STM":"#ff7f0e","YTHDC1 DMSO":"#2ca02c","YTHDC1 STM":"#d62728","Unlabeled":"#7f7f7f",
}

# ============================ helpers ============================
def sample_number_from_path(path: str) -> int | None:
    for pat in (r"Sample\s*([0-9]+)", r"Sample([0-9]+)"):
        m = re.search(pat, path, flags=re.IGNORECASE)
        if m: return int(m.group(1))
    return None

def sample_label_from_path(path: str) -> str:
    parent = os.path.basename(os.path.dirname(path))
    if parent and parent != os.path.basename(path): return parent
    base = os.path.basename(path)
    return re.sub(r"\.bed(\.gz)?$", "", base, flags=re.IGNORECASE)

def apply_passonly_filters(df: pd.DataFrame) -> pd.DataFrame:
    for col in ["valid_cov","n_mod","n_fail","n_diff","score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    mask = (df["n_diff"] <= 3 * df["valid_cov"]) & (df["n_mod"] > df["n_fail"] + 1)
    return df.loc[mask].copy()

def _read_parts(input_dir: str, pattern: str, passonly: bool, mod_code_filter=None) -> pd.DataFrame:
    """
    Read bedMethyl; return per (chrom,pos,sample,group) sums.

    mod_code_filter:
      - None => keep ALL mod codes
      - str  => keep only that code
      - Iterable[str] => keep only codes in the set

    De-duplication guard: for a given (chrom,pos,mod_code) WITHIN a sample,
    keep MAX(score) before aggregating within-sample.
    """
    paths = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not paths:
        raise SystemExit(f"No bedMethyl files matched: {input_dir}/{pattern}")

    # Normalize filter to a set of strings or None
    filt_set = None
    if mod_code_filter is None:
        filt_set = None
    elif isinstance(mod_code_filter, str):
        filt_set = {mod_code_filter}
    else:
        filt_set = set(str(x) for x in mod_code_filter)

    parts = []
    for p in paths:
        df = pd.read_csv(
            p, sep="\t", header=None, names=BED_COLS,
            comment="#", usecols=list(range(18)),
            dtype={"chrom": str, "mod_code": str}, low_memory=False
        )
        if filt_set is not None:
            df = df[df["mod_code"].isin(filt_set)]
        if df.empty:
            continue
        if passonly:
            df = apply_passonly_filters(df)
            if df.empty:
                continue

        sample_num = sample_number_from_path(p)
        group = NUMBER_TO_GROUP.get(sample_num, "Unlabeled")
        sample_label = sample_label_from_path(p)

        # === Position indexing from chromStart ===
        # BED is 0-based start, 1-based end. We store/plot 1-based coordinate = chromStart + 1.
        df["pos"] = df["start"].astype(int) + 1

        # De-duplication within a sample: keep max(score) per (chrom,pos,mod_code)
        df = (df.groupby(["chrom","pos","mod_code"], as_index=False)
                .agg(score=("score","max"),
                     n_mod=("n_mod","sum"), n_can=("n_can","sum"), n_fail=("n_fail","sum")))

        tmp = df[["chrom","pos","n_mod","n_can","n_fail","score"]].copy()
        tmp.rename(columns={"n_mod":"modified","n_can":"unmodified","n_fail":"fail"}, inplace=True)
        tmp["sample"] = sample_label
        tmp["group"]  = group
        parts.append(tmp)

    if not parts:
        tag = "after pass-only filters" if passonly else "from raw"
        tag2 = " (filtered mod codes)" if filt_set is not None else " (all mods)"
        raise SystemExit(f"No rows found {tag}{tag2}.")
    # Collapse across mod codes WITHIN a sample at a position
    out = (pd.concat(parts, ignore_index=True)
             .groupby(["chrom","pos","sample","group"], as_index=False)[["modified","unmodified","fail","score"]].sum())
    return out

def read_beds_raw_and_passonly(input_dir: str, pattern: str):
    # Only m6A ("a") for z/plotting
    raw_a   = _read_parts(input_dir, pattern, passonly=False, mod_code_filter="a")
    ponly_a = _read_parts(input_dir, pattern, passonly=True,  mod_code_filter="a")
    return raw_a, ponly_a

def build_coverage_selectedmods(input_dir: str, pattern: str, cov_agg: str) -> pd.DataFrame:
    """
    Coverage from RAW for ONLY codes [19227, 19228, 19229, 69426].
    Steps:
      1) Read & de-duplicate within-sample (max score per mod code at pos = chromStart+1).
      2) Sum across the four mod codes to get per-sample coverage at (chrom,pos).
      3) Aggregate across samples at each (chrom,pos) using cov_agg:
         - 'max'  : max per-sample coverage (default; prevents inflation)
         - 'mean' : average per-sample coverage
         - 'sum'  : sum across samples (cohort-wide view)
    Returns: chrom, pos, total_called_bases
    """
    codes = ["19227", "19228", "19229", "69426"]
    raw_sel = _read_parts(input_dir, pattern, passonly=False, mod_code_filter=codes)

    # Per-sample coverage at (chrom,pos) = sum of score across selected codes
    per_sample = (raw_sel.groupby(["chrom","pos","sample"], as_index=False)["score"]
                        .sum()
                        .rename(columns={"score":"per_sample_called"}))

    # Aggregate across samples
    if cov_agg == "sum":
        cov = per_sample.groupby(["chrom","pos"], as_index=False)["per_sample_called"].sum()
    elif cov_agg == "mean":
        cov = per_sample.groupby(["chrom","pos"], as_index=False)["per_sample_called"].mean()
    else:  # "max" (default)
        cov = per_sample.groupby(["chrom","pos"], as_index=False)["per_sample_called"].max()

    cov.rename(columns={"per_sample_called":"total_called_bases"}, inplace=True)
    return cov

# -------- z from RAW (pipeline-consistent) --------
def compute_z_from_raw(raw_df: pd.DataFrame) -> pd.DataFrame:
    d = raw_df.copy()
    for col in ["modified","fail","score"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["mod_ratio"]  = d["modified"] / d["score"].replace(0, np.nan)
    d["fail_ratio"] = d["fail"]     / d["score"].replace(0, np.nan)

    out = []
    for (chrom, pos), g in d.groupby(["chrom","pos"], as_index=False):
        mu  = g["fail_ratio"].mean(skipna=True)
        sd  = g["fail_ratio"].std(ddof=0, skipna=True)
        sd  = np.nan if (pd.isna(sd) or sd == 0) else sd
        gg  = g.copy()
        gg["z"] = (gg["mod_ratio"] - mu) / sd if not pd.isna(sd) else np.nan
        out.append(gg[["chrom","pos","sample","mod_ratio","fail_ratio","z","score"]])
    zt = pd.concat(out, ignore_index=True)
    zt.rename(columns={"score":"score_raw"}, inplace=True)
    return zt

# -------- GTF parsing & annotation --------
def parse_gtf_features_for_gene(gtf_path: str, gene: str):
    by_tx = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            a = line.rstrip("\n").split("\t")
            if len(a) < 9: continue
            chrom, source, feat, s, e, score, strand, frame, attrs = a
            if f'gene_name "{gene}"' not in attrs: continue

            tx_m = re.search(r'transcript_id "([^"]+)"', attrs)
            exnum_m = re.search(r'exon_number "([^"]+)"', attrs)
            tx_id = tx_m.group(1) if tx_m else None
            if not tx_id: continue

            s0 = int(s) - 1  # 0-based
            e0 = int(e)      # half-open
            rec = by_tx.setdefault(tx_id, {"chrom":chrom,"strand":strand,"exons":[], "cds":[], "utr5":[], "utr3":[]})
            if feat == "exon":
                exon_num = int(exnum_m.group(1)) if exnum_m else None
                rec["exons"].append((s0,e0,exon_num))
            elif feat == "CDS":
                rec["cds"].append((s0,e0))
            elif feat in ("5UTR","five_prime_UTR"):
                rec["utr5"].append((s0,e0))
            elif feat in ("3UTR","three_prime_UTR"):
                rec["utr3"].append((s0,e0))
    return {"by_tx": by_tx}

def pick_best_transcript(feats, preferred_tx=None):
    by_tx = feats["by_tx"]
    if not by_tx:
        return None, None
    if preferred_tx and preferred_tx in by_tx:
        return preferred_tx, by_tx[preferred_tx]
    def cds_len(t): return sum(e-s for (s,e) in t["cds"]) if t["cds"] else 0
    best_tx = max(by_tx.items(), key=lambda kv: (cds_len(kv[1]), sum(e-s for (s,e,_) in kv[1]["exons"])))[0]
    return best_tx, by_tx[best_tx]

def merge_intervals(ivls):
    if not ivls: return []
    ivls = sorted(ivls, key=lambda x:(x[0],x[1]))
    merged = [list(ivls[0])]
    for s,e in ivls[1:]:
        if s > merged[-1][1]:
            merged.append([s,e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [(s,e) for s,e in merged]

def intervals_to_1based(ivls): return [(s+1,e) for s,e in ivls]

def label_region_for_pos(pos1, exons1, utr5_1b, cds_1b, utr3_1b):
    def contains(ivls, p):
        for s,e in ivls:
            if s <= p <= e: return True
        return False
    if contains(utr5_1b, pos1): return "5UTR"
    if contains(cds_1b,  pos1): return "CDS"
    if contains(utr3_1b, pos1): return "3UTR"
    if contains(exons1,  pos1): return "exon_other"
    return "non_exonic"

# -------- Axis helpers --------
def set_genomic_ticks(ax, chrom: str, start: int, end: int, max_ticks: int):
    n = max(2, min(max_ticks, 20))
    ticks = np.linspace(start, end, num=n, dtype=int)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{chrom}:{v:,}" for v in ticks], rotation=90, ha="right", fontsize=9)

def draw_exon_and_regions(ax, exon_intervals, utr5, cds, utr3):
    ax.set_ylim(0, 1); ax.set_yticks([])
    for sp in ("top","right","left","bottom"):
        if sp != "bottom": ax.spines[sp].set_visible(False)
    color_map = {"5UTR":"#7f7f7f", "CDS":"#2ca02c", "3UTR":"#8c564b"}  # grey, green, brown
    def draw_band(ivls, label):
        for s,e in ivls:
            ax.add_patch(plt.Rectangle((s, 0.10), e - s, 0.30, alpha=0.3, color=color_map[label]))
    draw_band(exon_intervals, "CDS")
    draw_band(utr5, "5UTR"); draw_band(cds, "CDS"); draw_band(utr3, "3UTR")
    for (s,e) in exon_intervals:
        ax.add_patch(plt.Rectangle((s, 0.55), e - s, 0.35, color="navy"))
    ax.text(0.0, 1.05, "Exons (navy) • 5′UTR (grey) • CDS (green) • 3′UTR (brown)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=8)

# -------- Group means --------
def compute_group_means(sub: pd.DataFrame, weighted: bool) -> pd.DataFrame:
    df = sub.copy()
    df["total"] = df["modified"] + df["unmodified"]
    df = df[df["total"] > 0]
    if weighted:
        gm = (df.groupby(["pos","group"], as_index=False)
                .agg(mod_sum=("modified","sum"), tot_sum=("total","sum")))
        gm["mean_fraction"] = np.clip(gm["mod_sum"] / gm["tot_sum"], 0.0, 1.0)
        gm["mean_type"] = "weighted"
        return gm[["pos","group","mean_fraction","mean_type"]]
    else:
        df["frac"] = np.clip(df["modified"] / df["total"], 0.0, 1.0)
        gm = (df.groupby(["pos","group"], as_index=False)
                .agg(mean_fraction=("frac","mean")))
        gm["mean_fraction"] = np.clip(gm["mean_fraction"], 0.0, 1.0)
        gm["mean_type"] = "unweighted"
        return gm

# ============================ plotting & export ============================
def plot_region(
    sub: pd.DataFrame, gene: str, chrom: str, start_bp: int, end_bp: int,
    exon_ivl_1based, utr5_1b, cds_1b, utr3_1b,
    outdir: str, title_suffix: str,
    max_xticks: int, point_size: float,
    label_group_means: bool, group_mean_weighted: bool,
    export_rows: list, export_group_rows: list, export_summaries: list,
    cov_df: pd.DataFrame
):
    region = sub[(sub["pos"] >= start_bp) & (sub["pos"] <= end_bp)].copy()
    if region.empty:
        print(f"[WARN] {gene} {chrom}:{start_bp}-{end_bp} has no points after filters")
        return

    region["total"] = region["modified"] + region["unmodified"]
    region = region[region["total"] > 0]
    region["frac"]  = np.clip(region["modified"] / region["total"], 0.0, 1.0)

    # region label
    region["region"] = [label_region_for_pos(p, exon_ivl_1based, utr5_1b, cds_1b, utr3_1b)
                        for p in region["pos"].to_numpy()]

    # export rows
    export_rows.append(
        region.assign(gene=gene, chrom=chrom)[
            ["gene","chrom","pos","sample","group",
             "modified","unmodified","total",
             "score_raw","score_pass","mod_ratio","fail_ratio","z","frac","region"]
        ]
    )

    # coverage (already aggregated across samples according to --cov-agg)
    cov_region = (cov_df[(cov_df["chrom"] == chrom) &
                         (cov_df["pos"]   >= start_bp) &
                         (cov_df["pos"]   <= end_bp)]
                  .sort_values("pos"))
    x_cov = cov_region["pos"].to_numpy()
    y_cov = cov_region["total_called_bases"].to_numpy()

    # group means (no numeric labels)
    gm = compute_group_means(region[["pos","group","modified","unmodified","frac"]].copy(),
                             weighted=group_mean_weighted)
    export_group_rows.append(gm.assign(gene=gene, chrom=chrom)[["gene","chrom","pos","group","mean_fraction","mean_type"]])

    # summaries
    group_summary  = (region.groupby("group", as_index=False)
                      .agg(n_sites=("pos","nunique"), mean_frac=("frac","mean"),
                           median_frac=("frac","median"), mean_cov=("total","mean")))
    sample_summary = (region.groupby(["sample","group"], as_index=False)
                      .agg(n_sites=("pos","nunique"), mean_frac=("frac","mean"), mean_cov=("total","mean")))
    export_summaries.append({
        "gene": f"{gene} {title_suffix}",
        "total_passing_sites": int(region["pos"].nunique()),
        "groups": group_summary,
        "samples": sample_summary
    })

    # figure
    fig = plt.figure(figsize=(12, 6.0))
    ax_main = plt.axes([0.10, 0.64, 0.72, 0.28])
    ax_exon = plt.axes([0.10, 0.45, 0.72, 0.12], sharex=ax_main)
    ax_cov  = plt.axes([0.10, 0.16, 0.72, 0.20], sharex=ax_main)

    seen_group = set()
    for (sample, grp), df_sg in region.groupby(["sample","group"], sort=False):
        color = GROUP_COLORS.get(grp, GROUP_COLORS["Unlabeled"])
        ax_main.scatter(df_sg["pos"].to_numpy(), df_sg["frac"].to_numpy(),
                        s=point_size, alpha=0.90, edgecolors="none",
                        color=color, label=(grp if grp not in seen_group else None))
        seen_group.add(grp)

    if label_group_means:
        for _, df_g in gm.groupby("group", sort=False):
            ax_main.scatter(df_g["pos"].to_numpy(), df_g["mean_fraction"].to_numpy(),
                            s=point_size*2.2, facecolors="none", edgecolors="black", linewidths=1.0)

    ax_main.set_ylim(0.0, 1.10)
    ax_main.set_ylabel("Fraction m6A")
    ax_main.set_title(f"{gene} — {title_suffix} • pass-only + z + cov filters")
    ax_main.tick_params(axis='x', labelbottom=False)
    ax_main.legend(ncol=1, fontsize=9, frameon=False, loc="upper left",
                   bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

    draw_exon_and_regions(ax_exon, exon_ivl_1based, utr5_1b, cds_1b, utr3_1b)
    ax_exon.tick_params(axis='x', labelbottom=False)

    ax_cov.plot(x_cov, y_cov, linewidth=1.0, color="0.25")
    ax_cov.set_ylabel("Number of Called Bases", fontsize=9)
    ax_cov.set_xlabel("Genomic position")
    ax_cov.set_ylim(bottom=0)
    set_genomic_ticks(ax_cov, chrom, start_bp, end_bp, max_xticks)

    tag = f"{gene}_{title_suffix.replace(' ','_').replace(':','-')}"
    out = os.path.join(outdir, f"{tag}_m6A_points_cov_passonly_z_cov")
    plt.savefig(out + ".png", dpi=300, bbox_inches="tight")
    plt.savefig(out + ".pdf", dpi=300, bbox_inches="tight")
    plt.savefig(out + ".svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {tag} → {out}.png / .pdf / .svg")

def write_exports(outdir: str, pass_rows: list, group_rows: list, summaries: list):
    os.makedirs(outdir, exist_ok=True)
    if pass_rows:
        df_pass = pd.concat(pass_rows, ignore_index=True)
        df_pass = df_pass.sort_values(["gene","chrom","pos","group","sample"])
        df_pass.to_csv(os.path.join(outdir, "passonly_sites.tsv"), sep="\t", index=False)
        print(f"[OK] wrote {os.path.join(outdir,'passonly_sites.tsv')} (rows={len(df_pass)})")

    if group_rows:
        df_gm = pd.concat(group_rows, ignore_index=True)
        df_gm = df_gm.sort_values(["gene","chrom","pos","group"])
        df_gm.to_csv(os.path.join(outdir, "group_means.tsv"), sep="\t", index=False)
        print(f"[OK] wrote {os.path.join(outdir,'group_means.tsv')} (rows={len(df_gm)})")

    summ_path = os.path.join(outdir, "summary_stats.txt")
    with open(summ_path, "w") as fh:
        for s in summaries:
            fh.write(f"=== {s['gene']} ===\n")
            fh.write(f"total_passing_sites\t{s['total_passing_sites']}\n")
            fh.write("\n# per-group summary\n")
            s["groups"][["group","n_sites","mean_frac","median_frac","mean_cov"]]\
                .to_csv(fh, sep="\t", index=False)
            fh.write("\n# per-sample summary\n")
            s["samples"][["sample","group","n_sites","mean_frac","mean_cov"]]\
                .to_csv(fh, sep="\t", index=False)
            fh.write("\n\n")
    print(f"[OK] wrote {summ_path}")

# ============================ main ============================
def main():
    ap = argparse.ArgumentParser(description="m6A plots with pass-only + z (from raw) + coverage filters, UTR/CDS bands, per-exon, and region export.")
    ap.add_argument("--gtf", required=True, help="GTF with exon/CDS/UTR features")
    ap.add_argument("--input-dir", default="pileup_out", help="Directory with modkit bedMethyl outputs")
    ap.add_argument("--glob", default="*/*_filtered_mod.bed*", help="Glob under input-dir to collect files")
    ap.add_argument("--genes", nargs="+", default=["BCL2","MYC"], help="Genes to plot")
    ap.add_argument("--preferred-transcript", default=None, help="Transcript ID to use for exon/UTR/CDS (e.g., NM_002467.6)")
    ap.add_argument("--z-thresh", type=float, default=3.0, help="Keep sites where ANY sample z > this threshold (z computed from RAW)")
    ap.add_argument("--min-cov", type=int, default=20, help="Minimum per-sample coverage (n_mod+n_can)")
    ap.add_argument("--max-xticks", type=int, default=8, help="Max genomic tick labels")
    ap.add_argument("--point-size", type=float, default=8.0, help="Dot size")
    ap.add_argument("--outdir", default="m6A_plots_points_cov", help="Output directory")
    ap.add_argument("--label-group-means", action="store_true", help="Show hollow markers at group means (no numeric labels)")
    ap.add_argument("--group-mean-weighted", action="store_true", help="Use coverage-weighted means")
    ap.add_argument("--per-exon", action="store_true", help="Also emit one figure per exon")
    ap.add_argument("--cov-agg", choices=["max","mean","sum"], default="max",
                    help="Aggregate per-sample coverage across samples (default: max)")
    args = ap.parse_args()

    # 1) RAW (for z; m6A only) and PASS-ONLY (plot/export; m6A only)
    raw_df, pass_df = read_beds_raw_and_passonly(args.input_dir, args.glob)

    # 1b) RAW coverage for selected mod codes with configurable aggregation across samples
    cov_df = build_coverage_selectedmods(args.input_dir, args.glob, args.cov_agg)

    # 2) z from RAW, merge onto PASS-ONLY (keep both scores explicitly)
    ztable = compute_z_from_raw(raw_df)  # chrom,pos,sample,mod_ratio,fail_ratio,z,score_raw
    pass_df = pass_df.rename(columns={"score":"score_pass"})
    pass_df = pass_df.merge(ztable, on=["chrom","pos","sample"], how="left")

    # 3) Site gating by z
    z_nonan = ztable.dropna(subset=["z"])
    if z_nonan.empty:
        raise SystemExit("All z are NaN across all positions; cannot gate by z. (Check scores / fail ratios.)")
    zpos = (z_nonan.groupby(["chrom","pos"])["z"].max().reset_index(name="zmax"))
    zpos = zpos[zpos["zmax"] > args.z_thresh][["chrom","pos"]]
    if zpos.empty:
        raise SystemExit(f"No sites passed z > {args.z_thresh}.")
    pass_df["total"] = pass_df["modified"] + pass_df["unmodified"]
    pass_df = pass_df.merge(zpos, on=["chrom","pos"], how="inner")
    pass_df = pass_df[pass_df["total"] > args.min_cov].copy()
    if pass_df.empty:
        raise SystemExit(f"No points left after min coverage > {args.min_cov}.")

    pass_rows, group_rows, summaries = [], [], []

    # 5) Plot per gene
    for gene in args.genes:
        feats = parse_gtf_features_for_gene(args.gtf, gene)
        tx_id, tx = pick_best_transcript(feats, preferred_tx=args.preferred_transcript)
        if not tx:
            print(f"[WARN] No transcript features for {gene}; skipping.")
            continue

        chrom = tx["chrom"]
        bed_chroms = set(pass_df["chrom"].unique())
        target_chrom = chrom if chrom in bed_chroms else (chrom.replace("chr","") if chrom.startswith("chr") else "chr"+chrom)

        exons_0b = merge_intervals([(s,e) for (s,e,_) in tx["exons"]])
        utr5_0b  = merge_intervals(tx["utr5"])
        cds_0b   = merge_intervals(tx["cds"])
        utr3_0b  = merge_intervals(tx["utr3"])

        gstart = min(s for s,_ in exons_0b) + 1  # 1-based
        gend   = max(e for _,e in exons_0b)      # 1-based end already

        exons_1b = intervals_to_1based(exons_0b)
        utr5_1b  = intervals_to_1based(utr5_0b)
        cds_1b   = intervals_to_1based(cds_0b)
        utr3_1b  = intervals_to_1based(utr3_0b)

        sub_gene = pass_df[(pass_df["chrom"] == target_chrom) & (pass_df["pos"] >= gstart) & (pass_df["pos"] <= gend)].copy()
        if sub_gene.empty:
            print(f"[WARN] {gene} region has no points after filters")
            continue

        # Whole-transcript
        plot_region(
            sub=sub_gene, gene=gene, chrom=target_chrom, start_bp=gstart, end_bp=gend,
            exon_ivl_1based=exons_1b, utr5_1b=utr5_1b, cds_1b=cds_1b, utr3_1b=utr3_1b,
            outdir=args.outdir, title_suffix=f"{tx_id} (all exons)",
            max_xticks=args.max_xticks, point_size=args.point_size,
            label_group_means=args.label_group_means, group_mean_weighted=args.group_mean_weighted,
            export_rows=pass_rows, export_group_rows=group_rows, export_summaries=summaries,
            cov_df=cov_df
        )

        # Per-exon
        if args.per_exon:
            for i,(s,e,exnum) in enumerate(sorted(tx["exons"], key=lambda x:(x[0],x[1]))):
                s1, e1 = s+1, e  # 1-based exon range
                exon_label = f"{tx_id} exon {exnum if exnum is not None else i+1}"
                def clip(ivls):
                    out=[]
                    for a,b in merge_intervals(ivls):
                        aa,bb = max(a, s), min(b, e)
                        if aa<bb: out.append((aa,bb))
                    return intervals_to_1based(out)
                plot_region(
                    sub=sub_gene, gene=gene, chrom=target_chrom, start_bp=s1, end_bp=e1,
                    exon_ivl_1based=[(s1,e1)], utr5_1b=clip(tx["utr5"]), cds_1b=clip(tx["cds"]), utr3_1b=clip(tx["utr3"]),
                    outdir=args.outdir, title_suffix=exon_label,
                    max_xticks=args.max_xticks, point_size=args.point_size,
                    label_group_means=args.label_group_means, group_mean_weighted=args.group_mean_weighted,
                    export_rows=pass_rows, export_group_rows=group_rows, export_summaries=summaries,
                    cov_df=cov_df
                )

    # 6) write tables
    write_exports(args.outdir, pass_rows, group_rows, summaries)

if __name__ == "__main__":
    main()

