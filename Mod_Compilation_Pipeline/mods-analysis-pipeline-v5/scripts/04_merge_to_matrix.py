#!/usr/bin/env python3
import os
import re
import argparse
from glob import glob
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np

"""
Input files (by_sample):  <sample>_<mod>_union_sites_with_attrs.bed

Header example:
#chrom  chromStart  chromEnd  name  score  strand  thickStart  thickEnd  color
valid_coverage  percent_modified  count_modified  count_canonical
count_other_mod count_delete  count_fail  count_diff  count_nocall

Output (per-mod):
{mod}_merged_union_matrix_wide.tsv
Columns:
chrom start end mod
{SAMPLE}__score
{SAMPLE}__strand
{SAMPLE}__thickStart
{SAMPLE}__thickEnd
{SAMPLE}__color
{SAMPLE}__valid_coverage
{SAMPLE}__percent_modified
{SAMPLE}__count_modified
{SAMPLE}__count_canonical
{SAMPLE}__count_other_mod
{SAMPLE}__count_delete
{SAMPLE}__count_fail
{SAMPLE}__count_diff
{SAMPLE}__count_nocall
...
"""

FALLBACK_MOD_CODES = ["a", "m", "17596", "17802", "19227", "19228", "19229", "69426"]
FILENAME_RE = re.compile(r"(.+)_([A-Za-z0-9]+)_union_sites_with_attrs\.bed$")

# Fields that are per-sample and should become columns; keep a sensible order
SAMPLE_FIELDS_ORDER = [
    "score", "strand", "thickStart", "thickEnd", "color",
    "valid_coverage", "percent_modified",
    "count_modified", "count_canonical", "count_other_mod",
    "count_delete", "count_fail", "count_diff", "count_nocall",
]

NUMERIC_FIELDS = {
    "score", "valid_coverage", "percent_modified",
    "count_modified", "count_canonical", "count_other_mod",
    "count_delete", "count_fail", "count_diff", "count_nocall",
    "thickStart", "thickEnd", "chromStart", "chromEnd"
}
STRING_FIELDS = {"strand", "color"}

def to_number(val):
    if val is None:
        return np.nan
    try:
        # percent_modified may be float; others often int
        if isinstance(val, str) and "." in val:
            return float(val)
        return float(val)
    except Exception:
        return np.nan

def parse_header(line):
    cols = line.lstrip("#").strip().split("\t")
    return cols

def read_bed(path, expected_mod=None):
    rows = []
    header = None
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                header = parse_header(line)
                continue
            if line.startswith("chrom"):  # some files print non-# header
                header = line.strip().split("\t")
                continue
            if header is None:
                raise ValueError(f"No header detected in {path}")

            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < len(header):
                # pad if needed (defensive)
                parts = parts + [""] * (len(header) - len(parts))

            rec = dict(zip(header, parts))

            # normalize keys we use downstream
            chrom = rec.get("chrom", rec.get("#chrom"))
            start = rec.get("chromStart")
            end   = rec.get("chromEnd")
            mod   = rec.get("name")
            if expected_mod and mod != expected_mod:
                # keep anyway, but we trust file name more
                mod = expected_mod

            # coerce numerics
            start_i = int(float(start)) if start is not None and start != "" else None
            end_i   = int(float(end)) if end is not None and end != "" else None

            # stash typed record
            typed = {
                "chrom": chrom,
                "start": start_i,
                "end": end_i,
                "mod": mod,
            }
            for k, v in rec.items():
                if k in {"chrom", "#chrom", "chromStart", "chromEnd", "name"}:
                    continue
                if k in NUMERIC_FIELDS:
                    typed[k] = to_number(v)
                elif k in STRING_FIELDS:
                    typed[k] = v if v != "" else "."
                else:
                    # keep unknowns as strings
                    typed[k] = v
            rows.append(typed)
    return rows

def main(input_dir, output_dir, mod_codes):
    os.makedirs(output_dir, exist_ok=True)

    bed_paths = sorted(glob(os.path.join(input_dir, "*_union_sites_with_attrs.bed")))
    if not bed_paths:
        raise FileNotFoundError(f"No *_union_sites_with_attrs.bed found in {input_dir}")

    per_mod_sites = defaultdict(lambda: OrderedDict())   # key: (chrom,start,end,mod) -> row dict
    samples = set()

    for bed in bed_paths:
        base = os.path.basename(bed)
        m = FILENAME_RE.match(base)
        if not m:
            # fallback parse; try to split on last two underscores
            stem = base.replace("_union_sites_with_attrs.bed", "")
            parts = stem.split("_")
            sample_name = "_".join(parts[:-1]) if len(parts) > 1 else stem
            mod_from_name = parts[-1] if len(parts) > 1 else None
        else:
            sample_name, mod_from_name = m.group(1), m.group(2)

        samples.add(sample_name)
        rows = read_bed(bed, expected_mod=mod_from_name)

        for r in rows:
            key = (r["chrom"], r["start"], r["end"], r["mod"])
            if key not in per_mod_sites[r["mod"]]:
                per_mod_sites[r["mod"]][key] = {
                    "chrom": r["chrom"],
                    "start": r["start"],
                    "end": r["end"],
                    "mod": r["mod"],
                }
            # attach per-sample fields
            for fld in SAMPLE_FIELDS_ORDER:
                if fld in r:
                    per_mod_sites[r["mod"]][key][f"{sample_name}__{fld}"] = r[fld]

    # Build per-mod DataFrames
    for mod in sorted(set(list(per_mod_sites.keys()) + list(mod_codes))):
        if mod not in per_mod_sites or not per_mod_sites[mod]:
            # skip mods with no rows
            continue

        rows_list = list(per_mod_sites[mod].values())
        df = pd.DataFrame(rows_list)

        # make sure all sample columns exist
        all_cols = []
        for sample in sorted(samples):
            for fld in SAMPLE_FIELDS_ORDER:
                all_cols.append(f"{sample}__{fld}")
        for c in all_cols:
            if c not in df.columns:
                df[c] = np.nan

        # fill NaNs with sensible defaults
        for sample in samples:
            # numeric
            for fld in {"score","valid_coverage","percent_modified","count_modified",
                        "count_canonical","count_other_mod","count_delete","count_fail",
                        "count_diff","count_nocall","thickStart","thickEnd"}:
                col = f"{sample}__{fld}"
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors="coerce")
                    df[col] = df[col].fillna(0)
            # string
            for fld in {"strand","color"}:
                col = f"{sample}__{fld}"
                if col in df.columns:
                    df[col] = df[col].fillna(".")

        # order columns
        base = ["chrom", "start", "end", "mod"]
        per_sample_cols = []
        for sample in sorted(samples):
            for fld in SAMPLE_FIELDS_ORDER:
                per_sample_cols.append(f"{sample}__{fld}")
        df = df[base + per_sample_cols].sort_values(["chrom","start","end"]).reset_index(drop=True)

        out_path = os.path.join(output_dir, f"{mod}_merged_union_matrix_wide.tsv")
        df.to_csv(out_path, sep="\t", index=False)
        print(f"✅ Wrote {out_path}  (rows={len(df)}, samples={len(samples)})")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Merge per-sample union intersections into per-mod WIDE matrices with one column per attribute per sample.")
    ap.add_argument("--input_dir", required=True, help="Directory containing *_union_sites_with_attrs.bed")
    ap.add_argument("--output_dir", required=True, help="Directory to write per-mod *_merged_union_matrix_wide.tsv")
    ap.add_argument("--mod_codes", nargs="+", default=FALLBACK_MOD_CODES, help="Modification codes to track")
    args = ap.parse_args()
    main(args.input_dir, args.output_dir, args.mod_codes)

