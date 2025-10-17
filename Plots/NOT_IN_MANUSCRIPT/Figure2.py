#!/usr/bin/env python3
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# Paths & I/O pattern
# =========================
BASE_DIR = "/Users/theo/Downloads/kharas_zscore_filtered_annotated/zscore_filtered_annotated"
PATTERN = "*_merged_union_matrix_wide__z_gt3.0_annotated.tsv"   # only gt3 files
OUT_DIR = os.path.join(BASE_DIR, "plots_site_level_comparisons")
os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# Sample → Group mappings
# =========================
SAMPLE_TO_GROUP = {
    # THP1 GFP DMSO 48h
    "1": "DMSO", "2": "DMSO", "3": "DMSO",
    # THP1 GFP STM 15uM 48h
    "4": "STM", "5": "STM", "6": "STM",
    # THP1 GFP-YTHDC1 DMSO 48h
    "7": "YTHDC1 DMSO", "8": "YTHDC1 DMSO", "9": "YTHDC1 DMSO",
    # THP1 GFP-YTHDC1 STM 15uM 48h
    "10": "YTHDC1 STM", "11": "YTHDC1 STM", "12": "YTHDC1 STM",
}
GROUPS = ["DMSO", "STM", "YTHDC1 DMSO", "YTHDC1 STM"]

# Friendly names for titles
MOD_FRIENDLY = {
    "17596": "inosine",
    "a": "m6A",
    "m": "m5C",
    "17802": "pseU",
    "69426": "2’OmeA",
    "19228": "2’OmeC",
    "19229": "2’OmeG",
    "19227": "2’OmeT",
}

# Colorblind-friendly palette
GROUP_COLORS = {
    "DMSO": "#0072B2",         # blue
    "STM": "#E69F00",          # orange
    "YTHDC1 DMSO": "#009E73",  # teal/green
    "YTHDC1 STM": "#CC79A7",   # magenta
}

# Coverage threshold and constants
MIN_VALID = 20
TINY = 1e-9

# Clipping just for the "% change in percent-modified vs control" plot
REL_PCT_CLIP_MIN = -200.0   # percent
REL_PCT_CLIP_MAX =  200.0   # percent

# =========================
# Helpers
# =========================
def extract_sample_prefixes(columns):
    """
    Map each per-sample 'prefix' (up to ...pod5__) to:
      num, mod_col, valid_col
    Accepts __valid_coverage or falls back to __score.
    """
    stems = {}
    for col in columns:
        if not (col.endswith("__count_modified") or col.endswith("__valid_coverage") or col.endswith("__score")):
            continue
        m_prefix = re.match(r"^(.*?pod5)__", col)
        if not m_prefix:
            continue
        prefix = m_prefix.group(1)
        m_num = re.search(r"Sample(\d+)", prefix)
        if not m_num:
            continue
        num = m_num.group(1)

        if prefix not in stems:
            stems[prefix] = {"num": num, "mod_col": None, "valid_col": None}
        if col.endswith("__count_modified"):
            stems[prefix]["mod_col"] = col
        elif col.endswith("__valid_coverage"):
            stems[prefix]["valid_col"] = col
        elif col.endswith("__score") and stems[prefix]["valid_col"] is None:
            stems[prefix]["valid_col"] = col

    return {k: v for k, v in stems.items() if v["mod_col"] and v["valid_col"]}

def site_id_series(df):
    return (df["chrom"].astype(str) + ":" +
            df["start"].astype(str) + ":" +
            df["end"].astype(str) + ":" +
            df["mod"].astype(str))

def ensure_dirs(*paths):
    for p in paths:
        os.makedirs(p, exist_ok=True)

def save_fig_both(path_base):
    plt.savefig(path_base + ".png", dpi=300, bbox_inches="tight")
    plt.savefig(path_base + ".svg", bbox_inches="tight")

def control_for_pair(a, b):
    """
    UPDATED CONTROL PRIORITY:
      1) DMSO if present (ALWAYS)
      2) else YTHDC1 DMSO if present
      3) else (STM vs YTHDC1 STM): control = STM
      4) else fall back to the first argument as control
    Returns (ctrl, other)
    """
    if "DMSO" in (a, b):
        ctrl = "DMSO"
        other = b if a == "DMSO" else a
        return (ctrl, other)
    if "YTHDC1 DMSO" in (a, b):
        ctrl = "YTHDC1 DMSO"
        other = b if a == "YTHDC1 DMSO" else a
        return (ctrl, other)
    if {a, b} == {"STM", "YTHDC1 STM"}:
        return ("STM", "YTHDC1 STM")
    return (a, b)

def pooled_group_metrics(df, stems_map):
    """
    Pooled across replicates per group:
      <group>__Nmod = sum count_modified
      <group>__Nvalid = sum valid_coverage (or score)
      <group>__Pct = Nmod / Nvalid
    """
    df = df.copy()
    df["_site_id"] = site_id_series(df)
    group_sums = {}

    for prefix, meta in stems_map.items():
        grp = SAMPLE_TO_GROUP.get(meta["num"])
        if not grp:
            continue
        nmod = pd.to_numeric(df[meta["mod_col"]], errors="coerce").fillna(0)
        nval = pd.to_numeric(df[meta["valid_col"]], errors="coerce").fillna(0)
        chunk = pd.DataFrame({"_site_id": df["_site_id"], f"{grp}__Nmod": nmod, f"{grp}__Nvalid": nval})
        g = chunk.groupby("_site_id").sum()
        if grp not in group_sums:
            group_sums[grp] = g
        else:
            group_sums[grp] = group_sums[grp].add(g, fill_value=0)

    if not group_sums:
        return pd.DataFrame()

    all_sites = set().union(*(g.index for g in group_sums.values()))
    out = pd.DataFrame(index=sorted(all_sites))

    for grp, g in group_sums.items():
        out = out.join(g, how="left")
        out[f"{grp}__Pct"] = out[f"{grp}__Nmod"] / out[f"{grp}__Nvalid"].replace(0, np.nan)

    return out

def comparison_label(ctrl, other):
    # filenames are always "other_vs_ctrl"
    return f"{other}_vs_{ctrl}"

# =========================
# Main
# =========================
tsv_files = sorted(glob.glob(os.path.join(BASE_DIR, PATTERN)))
if not tsv_files:
    raise SystemExit(f"No files matched {os.path.join(BASE_DIR, PATTERN)}")

# all 6 comparisons with updated control priority
pairs = []
for i in range(len(GROUPS)):
    for j in range(i+1, len(GROUPS)):
        pairs.append(control_for_pair(GROUPS[i], GROUPS[j]))

filter_log_lines = [f"Coverage filter per comparison: require Nvalid ≥ {MIN_VALID} in BOTH groups\n",
                    "Control priority: DMSO > YTHDC1 DMSO > (STM in {STM,YTHDC1 STM}) > first-arg fallback\n"]

for path in tsv_files:
    mod_code = os.path.basename(path).split("_merged_union_matrix_wide__z_gt3.0_annotated.tsv")[0]
    friendly = MOD_FRIENDLY.get(mod_code, mod_code)

    mod_dir = os.path.join(OUT_DIR, mod_code)
    plots_dir = os.path.join(mod_dir, "plots")
    tables_dir = os.path.join(mod_dir, "tables")
    ensure_dirs(mod_dir, plots_dir, tables_dir)

    # read minimal columns
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    need = ["chrom", "start", "end", "mod"]
    need += [c for c in header if c.endswith(("__count_modified", "__valid_coverage", "__score"))]
    df = pd.read_csv(path, sep="\t", usecols=need, low_memory=False)

    stems = extract_sample_prefixes(df.columns)
    if not stems:
        print(f"[{mod_code}] No sample columns detected, skipping.")
        continue

    gdf = pooled_group_metrics(df, stems)
    if gdf.empty:
        print(f"[{mod_code}] No pooled metrics, skipping.")
        continue

    for ctrl, other in pairs:
        base = comparison_label(ctrl, other)  # other_vs_ctrl for filenames
        lbl = f"{friendly}: {other} vs {ctrl}"

        needed = [f"{ctrl}__Nmod", f"{ctrl}__Nvalid", f"{ctrl}__Pct",
                  f"{other}__Nmod", f"{other}__Nvalid", f"{other}__Pct"]
        if any(c not in gdf.columns for c in needed):
            continue

        sub = gdf[needed].copy()
        n_all = sub.shape[0]

        # coverage filter in BOTH groups
        mask_cov = (sub[f"{ctrl}__Nvalid"] >= MIN_VALID) & (sub[f"{other}__Nvalid"] >= MIN_VALID)
        sub = sub.loc[mask_cov].copy()
        n_pass = sub.shape[0]

        # values
        ctrl_nmod = sub[f"{ctrl}__Nmod"].astype(float)
        othr_nmod = sub[f"{other}__Nmod"].astype(float)

        ctrl_pct = sub[f"{ctrl}__Pct"].astype(float).fillna(0.0)   # fractions 0..1
        othr_pct = sub[f"{other}__Pct"].astype(float).fillna(0.0)

        # (1) % change in # modified molecules, relative to control
        pct_change_nmod = 100.0 * (othr_nmod - ctrl_nmod) / np.maximum(ctrl_nmod, 1.0)

        # (2) Δ percent-modified in percentage points
        delta_percent_modified_pp = 100.0 * (othr_pct - ctrl_pct)

        # (3) % change in percent-modified relative to control fraction (PER-SITE GUARD)
        # Guard per site: require Pct_ctrl ≥ 10 / Nvalid_ctrl_site
        ctrl_site_cov = sub[f"{ctrl}__Nvalid"].astype(float)
        ctrl_frac_thresholds = 10.0 / np.maximum(ctrl_site_cov, 1.0)  # per-site threshold (fraction units)
        rel_mask = ctrl_pct >= ctrl_frac_thresholds

        pct_change_percent_modified_vs_ctrl = np.full_like(ctrl_pct, np.nan, dtype=float)
        pct_change_percent_modified_vs_ctrl[rel_mask] = 100.0 * (
            (othr_pct[rel_mask] - ctrl_pct[rel_mask]) / np.maximum(ctrl_pct[rel_mask], TINY)
        )

        # Save tidy table (unclipped; NaN where guard fails)
        tidy = pd.DataFrame({
            "site_id": sub.index,  # chrom:start:end:mod
            f"{ctrl}__Nmod": ctrl_nmod,
            f"{ctrl}__Nvalid": sub[f"{ctrl}__Nvalid"].astype(float),
            f"{ctrl}__Pct": ctrl_pct,
            f"{other}__Nmod": othr_nmod,
            f"{other}__Nvalid": sub[f"{other}__Nvalid"].astype(float),
            f"{other}__Pct": othr_pct,
            "pct_change_Nmod_vs_ctrl": pct_change_nmod,
            "delta_percent_modified_pp": delta_percent_modified_pp,
            "pct_change_percent_modified_vs_ctrl": pct_change_percent_modified_vs_ctrl,
            "per_site_ctrl_frac_threshold": ctrl_frac_thresholds,
        })
        tidy_path = os.path.join(tables_dir, f"{mod_code}__{base}__site_metrics.tsv")
        tidy.to_csv(tidy_path, sep="\t", index=False)

        # filtering + guard stats
        n_rel_kept = int(np.isfinite(pct_change_percent_modified_vs_ctrl).sum())
        filter_log_lines.append(
            f"[{mod_code}] {base}: total_sites={n_all}, after_cov(both ≥{MIN_VALID})={n_pass}, "
            f"for rel %change plot kept {n_rel_kept} sites (per-site guard)"
        )

        color = GROUP_COLORS.get(other, "#666666")

        # 1) Histogram: % change in # modified molecules
        plt.figure(figsize=(7.6, 4.8))
        plt.hist(pct_change_nmod, bins=60, edgecolor="black", alpha=0.85, color=color)
        plt.axvline(0, linestyle="--", linewidth=1.5, color="#444444")
        plt.xlabel(f"% change in # modified molecules vs {ctrl}")
        plt.ylabel("Count of sites")
        plt.title(f"{lbl} — % change in # modified molecules")
        plt.tight_layout()
        save_fig_both(os.path.join(plots_dir, f"{mod_code}__{base}__hist_pct_change_Nmod"))
        plt.close()

        # 2) Scatter: Nmod_other vs Nmod_control
        plt.figure(figsize=(6.6, 6.0))
        plt.scatter(ctrl_nmod, othr_nmod, s=8, alpha=0.6, edgecolors="none", c=color)
        max_xy = float(np.nanmax([ctrl_nmod.max(), othr_nmod.max(), 1.0]))
        plt.plot([0, max_xy], [0, max_xy], linestyle="--", color="#444444", linewidth=1.2)
        plt.xlabel(f"{ctrl}: # modified molecules (pooled)")
        plt.ylabel(f"{other}: # modified molecules (pooled)")
        plt.title(f"{lbl} — modified molecules")
        plt.tight_layout()
        save_fig_both(os.path.join(plots_dir, f"{mod_code}__{base}__scatter_Nmod_vs_ctrl"))
        plt.close()

        # 3) Histogram: Δ percent-modified (percentage points)
        plt.figure(figsize=(7.6, 4.8))
        plt.hist(delta_percent_modified_pp, bins=60, edgecolor="black", alpha=0.85, color=color)
        plt.axvline(0, linestyle="--", linewidth=1.5, color="#444444")
        plt.xlabel(f"Δ percent-modified (pp), {other} − {ctrl}")
        plt.ylabel("Count of sites")
        plt.title(f"{lbl} — change in percent-modified (pp)")
        plt.tight_layout()
        save_fig_both(os.path.join(plots_dir, f"{mod_code}__{base}__hist_delta_percent_modified_pp"))
        plt.close()

        # 4) Histogram: % change in percent-modified vs control (PER-SITE GUARDED + CLIPPED)
        plot_vals = pct_change_percent_modified_vs_ctrl[np.isfinite(pct_change_percent_modified_vs_ctrl)]
        plot_vals = np.clip(plot_vals, REL_PCT_CLIP_MIN, REL_PCT_CLIP_MAX)

        plt.figure(figsize=(7.6, 4.8))
        plt.hist(plot_vals, bins=60, edgecolor="black", alpha=0.85, color=color)
        plt.axvline(0, linestyle="--", linewidth=1.5, color="#444444")
        plt.xlabel(f"% change in percent-modified vs {ctrl} (per-site guarded; clipped to {REL_PCT_CLIP_MIN:g}..{REL_PCT_CLIP_MAX:g}%)")
        plt.ylabel("Count of sites")
        plt.title(f"{lbl} — % change in percent-modified (guarded & clipped)")
        plt.tight_layout()
        save_fig_both(os.path.join(plots_dir, f"{mod_code}__{base}__hist_pct_change_percent_modified"))
        plt.close()

# Write filtering report
report_path = os.path.join(OUT_DIR, "filtering_stats.txt")
with open(report_path, "w") as fh:
    fh.write("\n".join(filter_log_lines) + "\n")

print("Wrote filtering report:", report_path)
print("All per-mod tables/plots under:", OUT_DIR)

