#!/usr/bin/env python3
import os
import argparse
import re
from glob import glob

import pandas as pd
import numpy as np

"""
Reads per-mod WIDE matrices from script 04 and appends per-sample z-score columns.

Per-sample z-score at each site:
z_{sample} = ( (count_modified/score)_{sample} - mean_over_samples(count_fail/score) ) / std_over_samples(count_fail/score)

- Uses each sample's BED score as the denominator (as in your previous pipeline).
- If std_over_samples == 0 for a site, z is NaN for all samples at that site.
- Writes:
  * {mod}_merged_union_matrix_wide__with_z.tsv
  * {mod}_merged_union_matrix_wide__z_gt{z}.tsv   (rows where ANY per-sample z > threshold), if --z_thresh is given
"""

COUNT_MOD_SUFFIX = "__count_modified"
COUNT_FAIL_SUFFIX = "__count_fail"
SCORE_SUFFIX = "__score"

def infer_samples(columns):
    """
    Infer sample names from columns like '{sample}__count_modified' or '{sample}__score'.
    """
    samples = set()
    for c in columns:
        if c.endswith(COUNT_MOD_SUFFIX):
            samples.add(c[: -len(COUNT_MOD_SUFFIX)])
        elif c.endswith(SCORE_SUFFIX):
            samples.add(c[: -len(SCORE_SUFFIX)])
    return sorted(samples)

def add_z_scores(df):
    samples = infer_samples(df.columns)

    # Build per-sample ratios
    mod_ratio_cols = {}
    fail_ratio_cols = {}

    for s in samples:
        cm_col = f"{s}{COUNT_MOD_SUFFIX}"
        cf_col = f"{s}{COUNT_FAIL_SUFFIX}"
        sc_col = f"{s}{SCORE_SUFFIX}"

        # ensure numeric
        for col in [cm_col, cf_col, sc_col]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
            else:
                df[col] = np.nan

        # avoid div by 0
        denom = df[sc_col].replace(0, np.nan)

        mod_ratio = df[cm_col] / denom
        fail_ratio = df[cf_col] / denom

        mod_ratio_cols[s] = mod_ratio
        fail_ratio_cols[s] = fail_ratio

    # DataFrames of ratios (index aligned with df)
    mod_ratio_df = pd.DataFrame(mod_ratio_cols)
    fail_ratio_df = pd.DataFrame(fail_ratio_cols)

    # Cross-sample stats of fail_ratio
    mean_fail = fail_ratio_df.mean(axis=1, skipna=True)
    std_fail = fail_ratio_df.std(axis=1, ddof=0, skipna=True)  # population std; use ddof=0 to be stable with few samples
    std_fail = std_fail.replace(0, np.nan)

    # Per-sample z = (mod_ratio - mean_fail) / std_fail
    z_cols = {}
    for s in samples:
        z = (mod_ratio_df[s] - mean_fail) / std_fail
        z_cols[s] = z.replace([np.inf, -np.inf], np.nan)

    z_df = pd.DataFrame({f"z__{s}": z_cols[s] for s in samples})
    out = pd.concat([df, z_df], axis=1)
    return out, [f"z__{s}" for s in samples]

def process_one(in_path, out_dir, z_thresh=None):
    base = os.path.basename(in_path)
    mod = base.split("_merged_union_matrix_wide.tsv")[0]

    df = pd.read_csv(in_path, sep="\t", dtype=str)  # read as str first, then coerce numerics selectively
    # Coerce numeric columns that look numeric
    for c in df.columns:
        if c in {"start","end"}:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        elif "__" in c:
            # try numeric; if fails we'll leave as string (e.g., strand/color)
            df[c] = pd.to_numeric(df[c], errors="ignore")

    with_z, z_cols = add_z_scores(df)

    out_full = os.path.join(out_dir, f"{mod}_merged_union_matrix_wide__with_z.tsv")
    with_z.to_csv(out_full, sep="\t", index=False)
    print(f"✅ Wrote {out_full} (rows={len(with_z)})")

    if z_thresh is not None:
        # keep rows where any sample's z > threshold
        mask = False
        for zc in z_cols:
            # coerce to numeric in case mixed types
            with_z[zc] = pd.to_numeric(with_z[zc], errors="coerce")
            mask = mask | (with_z[zc] > z_thresh)
        filtered = with_z[mask].copy()
        out_filt = os.path.join(out_dir, f"{mod}_merged_union_matrix_wide__z_gt{z_thresh}.tsv")
        filtered.to_csv(out_filt, sep="\t", index=False)
        print(f"✅ Wrote {out_filt} (rows={len(filtered)})")

def main(input_dir, output_dir, z_thresh):
    os.makedirs(output_dir, exist_ok=True)
    files = sorted([
        p for p in glob(os.path.join(input_dir, "*_merged_union_matrix_wide.tsv"))
        if os.path.isfile(p)
    ])
    if not files:
        raise FileNotFoundError(f"No *_merged_union_matrix_wide.tsv files found under {input_dir}")

    for f in files:
        process_one(f, output_dir, z_thresh=z_thresh)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Append per-sample z-score columns to per-mod WIDE matrices; optionally write a thresholded file.")
    ap.add_argument("--input_dir", required=True, help="Directory containing per-mod *_merged_union_matrix_wide.tsv")
    ap.add_argument("--output_dir", required=True, help="Directory for outputs")
    ap.add_argument("--z_thresh", type=float, default=None, help="Optional z threshold; if set, also emits *_z_gt{z}.tsv")
    args = ap.parse_args()
    main(args.input_dir, args.output_dir, args.z_thresh)

