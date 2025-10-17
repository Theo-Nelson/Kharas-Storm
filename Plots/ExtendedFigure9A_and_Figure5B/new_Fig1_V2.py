#!/usr/bin/env python3
import os
import re
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Silence a harmless FutureWarning inside upsetplot
warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot.data")

# UpSetPlot (required)
try:
    from upsetplot import from_memberships, UpSet
except Exception as e:
    raise SystemExit(
        "This script requires 'upsetplot'. Install with:\n"
        "  pip install upsetplot\n\n"
        f"Import error: {e}"
    )

# =========================
# Paths & I/O pattern
# =========================
BASE_DIR = "/Users/theo/Downloads/kharas_zscore_filtered_annotated/zscore_filtered_annotated"
PATTERN = "*_merged_union_matrix_wide__z_gt3.0_annotated.tsv"   # gt3 files only
OUT_DIR = os.path.join(BASE_DIR, "plots_gt3_sites_pass_filter")
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
GROUP_TO_REPS = {
    "DMSO": ["1", "2", "3"],
    "STM": ["4", "5", "6"],
    "YTHDC1 DMSO": ["7", "8", "9"],
    "YTHDC1 STM": ["10", "11", "12"],
}
GROUP_ORDER = ["DMSO", "STM", "YTHDC1 DMSO", "YTHDC1 STM"]

# Human-friendly mod titles
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

# Colorblind-friendly (Okabe–Ito-ish) for 4 groups
GROUP_COLORS = {
    "DMSO": "#0072B2",         # blue
    "STM": "#E69F00",          # orange
    "YTHDC1 DMSO": "#009E73",  # green/teal
    "YTHDC1 STM": "#CC79A7",   # magenta
}

STD_DDOF = 1  # sample SD across replicates

# =========================
# Helpers
# =========================
def extract_sample_prefixes(columns):
    """
    Build mapping per sample 'prefix' (the per-sample column stem) with:
      - num: sample number as string (e.g., '10')
      - mod_col: ...__count_modified
      - fail_col: ...__count_fail
      - pct_col:  ...__percent_modified (optional)
    Works with long headers like "...Sample10...pod5__count_modified".
    """
    prefixes = {}
    for col in columns:
        if not (col.endswith("__count_modified") or col.endswith("__count_fail") or col.endswith("__percent_modified")):
            continue
        m_prefix = re.match(r"^(.*?pod5)__", col)  # up to "...pod5__"
        if not m_prefix:
            continue
        prefix = m_prefix.group(1)
        m_num = re.search(r"Sample(\d+)", prefix)
        if not m_num:
            continue
        num = m_num.group(1)  # keep as string (1..12)

        if prefix not in prefixes:
            prefixes[prefix] = {"num": num, "mod_col": None, "fail_col": None, "pct_col": None}

        if col.endswith("__count_modified"):
            prefixes[prefix]["mod_col"] = col
        elif col.endswith("__count_fail"):
            prefixes[prefix]["fail_col"] = col
        elif col.endswith("__percent_modified"):
            prefixes[prefix]["pct_col"] = col

    # require both mod and fail
    return {k: v for k, v in prefixes.items() if v["mod_col"] is not None and v["fail_col"] is not None}

def safe_to_numeric(s):
    return pd.to_numeric(s, errors="coerce")

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

def ecdf(vals):
    if len(vals) == 0:
        return np.array([]), np.array([])
    v = np.sort(np.asarray(vals, dtype=float))
    y = np.arange(1, len(v) + 1) / len(v)
    return v, y

# =========================
# Main
# =========================
tsv_files = sorted(glob.glob(os.path.join(BASE_DIR, PATTERN)))
if not tsv_files:
    raise SystemExit(f"No files matched {os.path.join(BASE_DIR, PATTERN)}")

summary_rows = []

for path in tsv_files:
    mod_code = os.path.basename(path).split("_merged_union_matrix_wide__z_gt3.0_annotated.tsv")[0]
    friendly = MOD_FRIENDLY.get(mod_code, mod_code)

    # Output subdirs for this mod
    mod_dir   = os.path.join(OUT_DIR, mod_code)
    bars_dir  = os.path.join(mod_dir, "bars")
    ecdf_dir  = os.path.join(mod_dir, "ecdf")
    upset_dir = os.path.join(mod_dir, "upset")
    ensure_dirs(mod_dir, bars_dir, ecdf_dir, upset_dir)

    # Read only needed columns
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    needed_cols = ["chrom", "start", "end", "mod"]
    needed_cols += [c for c in header if c.endswith(("__count_modified", "__count_fail", "__percent_modified"))]

    df = pd.read_csv(path, sep="\t", usecols=needed_cols, low_memory=False)
    df["_site_id"] = site_id_series(df)

    pref_map = extract_sample_prefixes(df.columns)
    if not pref_map:
        print(f"[{mod_code}] No complete (mod/fail) column pairs found; skipping.")
        continue

    # ----- Per-sample pass info -----
    per_sample_counts = {}   # prefix -> {sample_num, n_sites_pass}
    passed_sets = {}         # sample_num -> set of site_ids passing
    passed_pct  = {}         # sample_num -> array of percent_modified for passed rows

    for prefix, meta in pref_map.items():
        s_num = meta["num"]
        mvals = safe_to_numeric(df[meta["mod_col"]]).fillna(0)
        fvals = safe_to_numeric(df[meta["fail_col"]]).fillna(0)
        mask  = (mvals > (fvals + 1))  # pass rule

        per_sample_counts[prefix] = {"sample_num": s_num, "n_sites_pass": int(mask.sum())}
        passed_sets[s_num] = set(df.loc[mask, "_site_id"].values)

        if meta["pct_col"] is not None:
            passed_pct[s_num] = safe_to_numeric(df.loc[mask, meta["pct_col"]]).dropna().values
        else:
            passed_pct[s_num] = np.array([])

    # ----- Bar plots: mean ± SD per group -----
    group_to_vals = {g: [] for g in GROUP_ORDER}
    for rec in per_sample_counts.values():
        grp = SAMPLE_TO_GROUP.get(rec["sample_num"])
        if grp:
            group_to_vals[grp].append(rec["n_sites_pass"])

    means, stds, cols = [], [], []
    for grp in GROUP_ORDER:
        vals = group_to_vals.get(grp, [])
        if len(vals) == 0:
            mean_v, std_v = 0.0, 0.0
        elif len(vals) == 1:
            mean_v, std_v = float(vals[0]), 0.0
        else:
            arr = np.array(vals, dtype=float)
            mean_v = float(arr.mean())
            std_v  = float(arr.std(ddof=STD_DDOF))
        means.append(mean_v)
        stds.append(std_v)
        cols.append(GROUP_COLORS[grp])
        summary_rows.append({
            "mod": mod_code,
            "mod_friendly": friendly,
            "group": grp,
            "mean_sites_pass": mean_v,
            "std_sites_pass": std_v,
            "replicate_values": vals,
        })

    plt.figure(figsize=(7.6, 4.8))
    x = np.arange(len(GROUP_ORDER))
    plt.bar(x, means, yerr=stds, capsize=6, color=cols, edgecolor="black")
    plt.xticks(x, GROUP_ORDER, rotation=0)
    plt.ylabel("Number of Sites")
    plt.title(f"Number of {friendly} Sites")
    plt.tight_layout()
    save_fig_both(os.path.join(bars_dir, f"{mod_code}_sites_pass_gt3_bar_mean_sd"))
    plt.close()

    # ----- ECDFs: one line per replicate; same color & single legend per group -----
    plt.figure(figsize=(7.6, 4.8))
    for grp in GROUP_ORDER:
        reps = GROUP_TO_REPS[grp]
        color = GROUP_COLORS[grp]
        for i, r in enumerate(reps):
            vals = passed_pct.get(r, np.array([]))
            xv, yv = ecdf(vals)
            if len(xv) == 0:
                continue
            # Label only first replicate line to have a single legend entry per group
            label = grp if i == 0 else "_nolegend_"
            # Use slight transparency so multiple replicate lines are visible
            plt.plot(xv, yv, label=label, linewidth=1.8, alpha=0.85, color=color)
    plt.xlabel("Percent modified")
    plt.ylabel("ECDF")
    plt.title(f"Distribution of {friendly} Sites for the passed sites from each condition")
    handles, labels = plt.gca().get_legend_handles_labels()
    # Keep one handle per unique label (they're already unique: one per group)
    plt.legend(handles, labels, frameon=False)
    plt.tight_layout()
    save_fig_both(os.path.join(ecdf_dir, f"{mod_code}_sites_pass_gt3_ecdf_percent_modified_replicate_lines"))
    plt.close()

    # ----- UpSet (within-group): replicate overlap (3-way) -----
    for grp, reps in GROUP_TO_REPS.items():
        rep_sets = [passed_sets.get(r, set()) for r in reps]
        union = set().union(*rep_sets)
        memberships = []
        for site in union:
            cats = []
            for r, s in zip(reps, rep_sets):
                if site in s:
                    cats.append(f"Sample {r}")
            memberships.append(cats if cats else [])
        series_within = from_memberships(memberships)

        plt.figure(figsize=(9, 5))
        UpSet(series_within, subset_size="count", show_counts=True).plot()
        plt.suptitle(f"{friendly}: Replicate overlap (passed sites) — {grp}", y=1.02)
        plt.tight_layout()
        out_base = os.path.join(upset_dir, f"{mod_code}_upset_within_{grp.replace(' ', '_')}")
        save_fig_both(out_base)
        plt.close()

    # ----- UpSet (4-way across groups): using replicate INTERSECTIONS within each group -----
    group_intersections = {}
    for grp, reps in GROUP_TO_REPS.items():
        sets_ = [passed_sets.get(r, set()) for r in reps]
        if len(sets_) == 0:
            group_intersections[grp] = set()
        else:
            inter = sets_[0]
            for s in sets_[1:]:
                inter = inter & s
            group_intersections[grp] = inter

    labels4 = GROUP_ORDER[:]  # fixed order
    sets4 = [group_intersections.get(g, set()) for g in labels4]
    union4 = set().union(*sets4)

    memberships4 = []
    for site in union4:
        cats = [name for name, s in zip(labels4, sets4) if site in s]
        memberships4.append(cats if cats else [])
    series4 = from_memberships(memberships4)

    plt.figure(figsize=(10, 6))
    UpSet(series4, subset_size="count", show_counts=True).plot()
    plt.suptitle(f"{friendly}: 4-way UpSet of replicate-intersection sites", y=1.02)
    plt.tight_layout()
    save_fig_both(os.path.join(upset_dir, f"{mod_code}_upset_4way_groups"))
    plt.close()

# ----- Save summary table for bars -----
summary_df = pd.DataFrame(summary_rows)
summary_csv = os.path.join(OUT_DIR, "summary_sites_pass_gt3_by_group.csv")
summary_df.to_csv(summary_csv, index=False)

print("Wrote summary CSV:", summary_csv)
print("All figures saved under:", OUT_DIR)

