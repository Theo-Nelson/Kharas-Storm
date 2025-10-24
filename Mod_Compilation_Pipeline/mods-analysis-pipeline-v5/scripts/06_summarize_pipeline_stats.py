#!/usr/bin/env python3
import os
import argparse
from typing import List, Tuple

HEADER_PREFIXES = ("#", "chrom")  # lines starting with these are not counted

def read_sample_list(sample_list_path: str) -> List[Tuple[str, str]]:
    """Return list of (sample_name, sample_dir) from a file of directories."""
    samples = []
    with open(sample_list_path, "r") as f:
        for line in f:
            p = line.strip()
            if not p:
                continue
            samples.append((os.path.basename(p.rstrip("/")), p.rstrip("/")))
    return samples

def count_sites(bed_path: str) -> int:
    """Count data rows (non-header) in a BED-like file. Returns 0 if missing."""
    if not os.path.isfile(bed_path):
        return 0
    n = 0
    with open(bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(HEADER_PREFIXES):
                continue
            n += 1
    return n

def find_z_filtered_file(zdir: str, mod: str, z: float) -> str:
    """
    Prefer wide-output from revised step 5:
      {mod}_merged_union_matrix_wide__z_gt{z}.tsv
    Fallback to earlier per-mod long format:
      {mod}_zscore_gt{z}.tsv
    Returns path or '' if not found.
    """
    cand1 = os.path.join(zdir, f"{mod}_merged_union_matrix_wide__z_gt{z}.tsv")
    if os.path.isfile(cand1):
        return cand1
    cand2 = os.path.join(zdir, f"{mod}_zscore_gt{z}.tsv")
    if os.path.isfile(cand2):
        return cand2
    return ""

def summarize_mod(
    mod: str,
    samples: List[Tuple[str, str]],
    root_dir: str,
    z_thresh: float,
    out_dir: str,
) -> str:
    """
    Build a text summary for one modification and write it to {out_dir}/{mod}_stats.txt.
    Returns the path to the written file.
    """
    filtered_beds_dir = os.path.join(root_dir, "outputs", "filtered_beds")
    unions_dir        = os.path.join(root_dir, "outputs", "union_beds", "unions")
    by_sample_dir     = os.path.join(root_dir, "outputs", "union_beds", "by_sample")
    zfilt_dir         = os.path.join(root_dir, "outputs", "zscore_filtered")

    lines = []
    lines.append(f"=== Modification: {mod} ===")
    lines.append("Per-sample counts")
    lines.append("sample\tfiltered_mod_sites\tpassonly_sites\tafter_union_sites")

    # Per-sample counts
    for sample_name, sample_dir in samples:
        pre_path      = os.path.join(sample_dir, f"{sample_name}_{mod}_filtered_mod.bed")
        passonly_path = os.path.join(filtered_beds_dir, f"{sample_name}_{mod}_filtered_passonly.bed")
        union_path    = os.path.join(by_sample_dir, f"{sample_name}_{mod}_union_sites_with_attrs.bed")

        n_pre      = count_sites(pre_path)
        n_pass     = count_sites(passonly_path)
        n_union    = count_sites(union_path)

        lines.append(f"{sample_name}\t{n_pre}\t{n_pass}\t{n_union}")

    # Union sites (global)
    union_file = os.path.join(unions_dir, f"all_samples_{mod}_union_sites.bed")
    n_union_global = count_sites(union_file)
    lines.append("")
    lines.append(f"union_sites_total\t{n_union_global}")

    # Z-score pass (global; ANY sample z > threshold)
    z_file = find_z_filtered_file(zfilt_dir, mod, z_thresh)
    if z_file and os.path.isfile(z_file):
        n_zpass = count_sites(z_file)  # these are TSV rows; header count logic still works
        lines.append(f"zscore_pass_any(z>{z_thresh})\t{n_zpass}")
        lines.append(f"zscore_pass_file\t{os.path.basename(z_file)}")
    else:
        lines.append(f"zscore_pass_any(z>{z_thresh})\t0")
        lines.append("zscore_pass_file\t<missing>")

    # Write file
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{mod}_stats.txt")
    with open(out_path, "w") as out:
        out.write("\n".join(lines) + "\n")
    return out_path

def main():
    ap = argparse.ArgumentParser(
        description="Summarize per-modification stats from the completed modkit pipeline."
    )
    ap.add_argument("--root_dir", required=True, help="Pipeline root directory (same as in run_full_pipeline.sbatch)")
    ap.add_argument("--sample_list", required=True, help="Path to file listing sample directories (same as pipeline arg)")
    ap.add_argument("--mod_codes", nargs="+", default=["a","m","17596","17802","19227","19228","19229","69426"],
                    help="Modification codes to summarize")
    ap.add_argument("--z_thresh", type=float, default=3.0, help="Z-score threshold used in step 5 (default: 3.0)")
    ap.add_argument("--output_dir", default=None, help="Directory for stats text files (default: {root}/outputs/stats)")
    args = ap.parse_args()

    root_dir = os.path.abspath(args.root_dir)
    out_dir = args.output_dir or os.path.join(root_dir, "outputs", "stats")

    samples = read_sample_list(args.sample_list)
    if not samples:
        raise SystemExit(f"No samples found in {args.sample_list}")

    written = []
    for mod in args.mod_codes:
        p = summarize_mod(mod, samples, root_dir, args.z_thresh, out_dir)
        print(f"✅ Wrote {p}")
        written.append(p)

    print(f"\nDone. Wrote {len(written)} stats files to: {out_dir}")

if __name__ == "__main__":
    main()

