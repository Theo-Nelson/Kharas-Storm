#!/usr/bin/env python3
"""
Per-group ECDFs of gene-level max stoichiometry (percent_modified).

Outputs PNG, PDF, and SVG for each plot.

Pipeline per mod table:
1) Per site, compute per-group mean percent_modified, using only replicates with score >= 20.
2) For each gene and group, take the max of that group's site-means across all sites for that gene.
3) Plot ECDFs for three gene sets (all / up_rescued / down_rescued) per group.
   - Combined plot: 12 lines (3 sets x 4 groups)
   - Separate plots: 1 per group, 3 lines each
4) Write a TSV with gene-level maxima per group.

Usage:
  python plot_ecdf_max_stoich_per_group.py GT3_DIR UP_TXT DOWN_TXT [--outdir OUTDIR]
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Tuple, Set

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------- Config ---------
# Sample -> group mapping (1-based sample indices)
GROUPS = {
    "GFP_DMSO":       [1, 2, 3],
    "GFP_STM15uM":    [4, 5, 6],
    "GFPYTHDC1_DMSO": [7, 8, 9],
    "GFPYTHDC1_STM":  [10, 11, 12],
}
SCORE_MIN = 20.0  # replicates must have score >= 20 to be included

# New line colors (intentionally avoid tab:blue/tab:orange/tab:green/tab:pink from the bar plot)
GROUP_COLORS: Dict[str, str] = {
    "GFP_DMSO":       "tab:red",
    "GFP_STM15uM":    "tab:purple",
    "GFPYTHDC1_DMSO": "tab:brown",
    "GFPYTHDC1_STM":  "tab:gray",
}
# Line styles per gene set (consistent across groups)
SET_STYLES: Dict[str, str] = {
    "all": "solid",
    "up_rescued": "dashed",
    "down_rescued": "dotted",
}

# --------- Helpers ---------
def read_gene_list(path: str) -> Set[str]:
    s = set()
    with open(path, "r") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                s.add(re.sub(r"\s+", "", g).upper())
    return s

def infer_mod_from_filename(path: str) -> str:
    base = os.path.basename(path)
    m = re.match(r"^([^_]+)_merged_union_matrix_wide__z_gt3\.0_annotated\.tsv$", base)
    return m.group(1) if m else base

def find_sample_num(col: str) -> int:
    m = re.search(r"Sample\s*(\d+)", col, flags=re.IGNORECASE)
    return int(m.group(1)) if m else -1

def pick_gene_col(df: pd.DataFrame) -> str:
    cands = [c for c in df.columns if c.lower() in ("gene_name", "gene", "genename")]
    return cands[0] if cands else df.columns[-1]

def ecdf(y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    y = y[np.isfinite(y)]
    if y.size == 0:
        return np.array([]), np.array([])
    xs = np.sort(y)
    ys = np.arange(1, xs.size + 1) / xs.size
    return xs, ys

def savefig_all(basepath_no_ext: str):
    """Save current matplotlib figure as .png, .pdf, .svg using the same basename."""
    # PNG (raster): use dpi for nice resolution
    plt.savefig(f"{basepath_no_ext}.png", dpi=200)
    # Vector formats
    plt.savefig(f"{basepath_no_ext}.pdf")
    plt.savefig(f"{basepath_no_ext}.svg")

# --------- Core per-table work ---------
def process_one_table(tsv_path: str, up_set: Set[str], down_set: Set[str], outdir: str):
    mod = infer_mod_from_filename(tsv_path)
    print(f"[INFO] Processing {mod} -> {os.path.basename(tsv_path)}")

    # read
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str, engine="python")
    except Exception as e:
        print(f"[WARN] Failed to read {tsv_path}: {e}", file=sys.stderr)
        return

    gene_col = pick_gene_col(df)
    df = df[df[gene_col].notna()].copy()
    df["_gene_upper"] = df[gene_col].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    # drop intergenic
    df = df[df["_gene_upper"] != "INTERGENIC"].copy()
    if df.empty:
        print(f"[WARN] {mod}: no non-intergenic rows, skipping.")
        return

    # locate percent_modified and score columns
    pct_cols = [c for c in df.columns if c.endswith("__percent_modified")]
    score_cols = [c for c in df.columns if c.endswith("__score")]
    if not pct_cols or not score_cols:
        # relaxed fallback if headers are slightly malformed
        pct_cols = [c for c in df.columns if "percent_modified" in c]
        score_cols = [c for c in df.columns if re.search(r"(^|_)score($|__)", c, flags=re.IGNORECASE)]

    # map sample->column
    sample_to_pct, sample_to_score = {}, {}
    for c in pct_cols:
        sn = find_sample_num(c)
        if 1 <= sn <= 12:
            sample_to_pct[sn] = c
    for c in score_cols:
        sn = find_sample_num(c)
        if 1 <= sn <= 12:
            sample_to_score[sn] = c

    # cast numeric
    for sn, c in sample_to_pct.items():
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for sn, c in sample_to_score.items():
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # per-site per-group mean (only replicates with score>=20)
    group_avg_cols: List[str] = []
    for gname, samples in GROUPS.items():
        vals = []
        for sn in samples:
            pc = sample_to_pct.get(sn)
            sc = sample_to_score.get(sn)
            if pc is None or sc is None:
                vals.append(pd.Series(np.nan, index=df.index))
            else:
                ok = df[sc] >= SCORE_MIN
                vals.append(df[pc].where(ok, np.nan))
        grp_avg = pd.concat(vals, axis=1).mean(axis=1, skipna=True)
        out_c = f"group_avg__{gname}"
        df[out_c] = grp_avg
        group_avg_cols.append(out_c)

    # For each gene and each group, take max across sites (rows)
    records = []
    for gname in GROUPS.keys():
        col = f"group_avg__{gname}"
        if col not in df.columns:
            continue
        tmp = df[["_gene_upper", col]].rename(columns={"_gene_upper": "gene", col: "stoich"})
        # collapse per gene: max stoich across sites (rows)
        gene_max = tmp.groupby("gene", as_index=False)["stoich"].max()
        gene_max["group"] = gname
        records.append(gene_max)

    if not records:
        print(f"[WARN] {mod}: no group averages present, skipping.")
        return

    gene_max_all = pd.concat(records, ignore_index=True)
    # annotate set membership
    gene_max_all["set"] = "all"
    gene_max_all.loc[gene_max_all["gene"].isin(up_set), "set"] = "up_rescued"
    gene_max_all.loc[gene_max_all["gene"].isin(down_set), "set"] = "down_rescued"

    # write TSV with gene-level maxima per group
    os.makedirs(outdir, exist_ok=True)
    tsv_by_group = os.path.join(outdir, f"{mod}__gene_max_stoich__by_group.tsv")
    gene_max_all.to_csv(tsv_by_group, sep="\t", index=False)

    # ---- plots ----
    # Combined 12-line ECDF
    base_comb = os.path.join(outdir, f"{mod}__ecdf_max_stoich__12lines")
    plt.figure(figsize=(7.5, 6.0))

    ordered_groups = list(GROUPS.keys())

    for gname in ordered_groups:
        sub = gene_max_all[gene_max_all["group"] == gname]
        if sub.empty:
            continue
        color = GROUP_COLORS.get(gname, "black")

        # All
        x_all, y_all = ecdf(sub["stoich"].to_numpy(dtype=float))
        if x_all.size:
            plt.plot(x_all, y_all, label=f"{gname}: all (n={x_all.size})",
                     linewidth=2, color=color, linestyle=SET_STYLES["all"])
        # Up
        sub_up = sub[sub["set"] == "up_rescued"]["stoich"].to_numpy(dtype=float)
        x_up, y_up = ecdf(sub_up)
        if x_up.size:
            plt.plot(x_up, y_up, label=f"{gname}: up (n={x_up.size})",
                     linewidth=2, color=color, linestyle=SET_STYLES["up_rescued"])
        # Down
        sub_dn = sub[sub["set"] == "down_rescued"]["stoich"].to_numpy(dtype=float)
        x_dn, y_dn = ecdf(sub_dn)
        if x_dn.size:
            plt.plot(x_dn, y_dn, label=f"{gname}: down (n={x_dn.size})",
                     linewidth=2, color=color, linestyle=SET_STYLES["down_rescued"])

    plt.xlabel("Gene-level max stoichiometry per group (percent_modified)")
    plt.ylabel("ECDF")
    plt.title(f"{mod}: ECDFs by group & gene set\n(per-site group means; score ≥ {SCORE_MIN})")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    savefig_all(base_comb)
    plt.close()
    print(f"[OK] {mod}: wrote combined plots -> {base_comb}.png/.pdf/.svg")

    # Separate per-group plots (3 lines each)
    for gname in ordered_groups:
        sub = gene_max_all[gene_max_all["group"] == gname]
        if sub.empty:
            continue

        base_grp = os.path.join(outdir, f"{mod}__ecdf_max_stoich__{gname}")
        plt.figure(figsize=(6.5, 5.0))

        color = GROUP_COLORS.get(gname, "black")

        # all
        x_all, y_all = ecdf(sub["stoich"].to_numpy(dtype=float))
        if x_all.size:
            plt.plot(x_all, y_all, label=f"All (n={x_all.size})",
                     linewidth=2, color=color, linestyle=SET_STYLES["all"])
        # up
        x_up, y_up = ecdf(sub.loc[sub["set"]=="up_rescued","stoich"].to_numpy(dtype=float))
        if x_up.size:
            plt.plot(x_up, y_up, label=f"up_rescued (n={x_up.size})",
                     linewidth=2, color=color, linestyle=SET_STYLES["up_rescued"])
        # down
        x_dn, y_dn = ecdf(sub.loc[sub["set"]=="down_rescued","stoich"].to_numpy(dtype=float))
        if x_dn.size:
            plt.plot(x_dn, y_dn, label=f"down_rescued (n={x_dn.size})",
                     linewidth=2, color=color, linestyle=SET_STYLES["down_rescued"])

        plt.xlabel("Gene-level max stoichiometry (percent_modified)")
        plt.ylabel("ECDF")
        plt.title(f"{mod} — {gname}\n(per-site group means; score ≥ {SCORE_MIN})")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        savefig_all(base_grp)
        plt.close()
        print(f"[OK] {mod}: wrote group plots -> {base_grp}.png/.pdf/.svg")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gt3_dir", help="Directory containing *_z_gt3.0_annotated.tsv files")
    ap.add_argument("up_list", help="Path to up_rescued.txt")
    ap.add_argument("down_list", help="Path to down_rescued.txt")
    ap.add_argument("--outdir", default="ecdf_maxstoich_by_gene_pergroup", help="Output directory")
    args = ap.parse_args()

    up_set = read_gene_list(args.up_list)
    down_set = read_gene_list(args.down_list)

    gt3_files = sorted(
        os.path.join(args.gt3_dir, f)
        for f in os.listdir(args.gt3_dir)
        if f.endswith("_z_gt3.0_annotated.tsv")
    )
    if not gt3_files:
        print(f"[ERROR] No *_z_gt3.0_annotated.tsv files found in {args.gt3_dir}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    for tsv in gt3_files:
        process_one_table(tsv, up_set, down_set, args.outdir)

    print("\nDone.")

if __name__ == "__main__":
    main()

