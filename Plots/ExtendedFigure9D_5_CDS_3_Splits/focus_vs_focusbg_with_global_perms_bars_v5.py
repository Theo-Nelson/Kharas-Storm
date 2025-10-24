#!/usr/bin/env python3
import os, re, glob, argparse, math, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

# ----------------------------- constants -----------------------------
REGIONS = ["5UTR","CDS","3UTR"]
COND_ORDER = [
    ("GFP_DMSO",     [1,2,3],    "DMSO",     "#d62728"),  # red
    ("GFP_STM15uM",  [4,5,6],    "STM",      "#9467bd"),  # purple
    ("DC1_DMSO",     [7,8,9],    "DC1 DMSO", "#8c564b"),  # brown
    ("DC1_STM15uM",  [10,11,12], "DC1 STM",  "#7f7f7f"),  # grey
]

# ----------------------------- GTF parsing & helpers -----------------------------
def parse_gtf_canonical(gtf_path):
    cols = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
    rows = []
    with open(gtf_path, "r") as fh:
        for ln in fh:
            if not ln or ln.startswith("#"): continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) != 9: continue
            rows.append(dict(zip(cols, parts)))
    gtf = pd.DataFrame(rows)
    if gtf.empty:
        raise ValueError("Parsed GTF is empty.")

    def attr_get(s, key):
        m = re.search(rf'{key} "([^"]+)"', s)
        return m.group(1) if m else None

    gtf["gene_id"] = gtf["attribute"].apply(lambda s: attr_get(s, "gene_id"))
    gtf["gene_name"] = gtf["attribute"].apply(lambda s: attr_get(s, "gene_name"))
    gtf["transcript_id"] = gtf["attribute"].apply(lambda s: attr_get(s, "transcript_id"))
    for c in ("start","end"):
        gtf[c] = pd.to_numeric(gtf[c], errors="coerce").astype("Int64")

    exon_like = gtf[gtf["feature"].isin(["exon","CDS","UTR","five_prime_UTR","three_prime_UTR"])]
    if exon_like.empty: exon_like = gtf
    tx_span = exon_like.groupby(["gene_id","transcript_id","seqname","strand"], as_index=False).agg(
        tx_start=("start","min"), tx_end=("end","max")
    )

    cds = gtf[gtf["feature"]=="CDS"]
    if cds.empty:
        warnings.warn("No CDS in GTF; cannot assign UTR/CDS.")
        return {}, {}

    cds_agg = cds.groupby(["gene_id","transcript_id"], as_index=False).agg(
        cds_start=("start","min"), cds_end=("end","max")
    )

    canon = pd.merge(tx_span, cds_agg, on=["gene_id","transcript_id"], how="inner")
    if canon.empty: return {}, {}

    canon["tx_len"] = canon["tx_end"] - canon["tx_start"] + 1
    canon.sort_values(["gene_id","tx_len"], ascending=[True, False], inplace=True)
    canon = canon.groupby("gene_id", as_index=False).head(1)

    g_names = gtf.dropna(subset=["gene_id","gene_name"]).groupby("gene_id")["gene_name"].agg(lambda x: x.iloc[0])
    canon = canon.merge(g_names.rename("gene_name"), on="gene_id", how="left")

    rec_by_id, rec_by_name = {}, {}
    for _, r in canon.iterrows():
        rec = {
            "chrom": r["seqname"], "strand": r["strand"],
            "cds_start": int(r["cds_start"]), "cds_end": int(r["cds_end"]),
            "gene_start": int(r["tx_start"]), "gene_end": int(r["tx_end"]),
        }
        rec_by_id[r["gene_id"]] = rec
        gnm = r.get("gene_name")
        if isinstance(gnm, str) and gnm:
            rec_by_name.setdefault(gnm, rec)
    return rec_by_id, rec_by_name

def assign_region(pos, strand, gs, cs, ce, ge):
    p = int(pos); gs=int(gs); cs=int(cs); ce=int(ce); ge=int(ge)
    if p < gs or p > ge: return "intergenic"
    if cs <= p <= ce: return "CDS"
    if strand == "+":
        if gs <= p < cs: return "5UTR"
        if ce < p <= ge: return "3UTR"
    else:
        if gs <= p < cs: return "3UTR"
        if ce < p <= ge: return "5UTR"
    return "intergenic"

def pick_gene_column(df):
    for c in ["gene","gene_name","gene_id","name","Gene","GeneName"]:
        if c in df.columns: return c
    return None

def pick_center(df):
    if "site_center" in df.columns: return pd.to_numeric(df["site_center"], errors="coerce")
    if "pos" in df.columns: return pd.to_numeric(df["pos"], errors="coerce")
    if "position" in df.columns: return pd.to_numeric(df["position"], errors="coerce")
    if {"start","end"}.issubset(df.columns):
        return (pd.to_numeric(df["start"], errors="coerce") + pd.to_numeric(df["end"], errors="coerce"))/2.0
    raise ValueError("No position column (need site_center / pos / position / start+end).")

# ----------------------------- sample column discovery -----------------------------
def discover_sample_cols(df):
    sc_map = defaultdict(list)  # score cols per Sample#
    st_map = defaultdict(list)  # stoich/fraction cols per Sample#
    for c in df.columns:
        mS = re.search(r"Sample(\d+)", c)
        if not mS: 
            continue
        s = int(mS.group(1))
        lc = c.lower()
        if lc.endswith("__score") or lc.endswith("_score"):
            sc_map[s].append(c)
        if "percent_modified" in lc or "stoich" in lc or "fraction_modified" in lc:
            st_map[s].append(c)
    for k in list(sc_map.keys()):
        sc_map[k] = sorted(sc_map[k])
    for k in list(st_map.keys()):
        st_map[k] = sorted(st_map[k])
    return sc_map, st_map

def group_cols_for_samples(sc_map, st_map, samples):
    score_cols, stoich_cols = [], []
    for s in samples:
        score_cols += sc_map.get(s, [])
        stoich_cols += st_map.get(s, [])
    return score_cols, stoich_cols

# ----------------------------- odds ratio math -----------------------------
def or_ci_from_counts(a,b,c,d, pseudo=0.5):
    if min(a,b,c,d) == 0:
        a,b,c,d = a+pseudo, b+pseudo, c+pseudo, d+pseudo
    orat = (a*d)/(b*c)
    lnOR = math.log(orat)
    se = math.sqrt(1.0/a + 1.0/b + 1.0/c + 1.0/d)
    lo, hi = math.exp(lnOR-1.96*se), math.exp(lnOR+1.96*se)
    return orat, lnOR, se, lo, hi

# ----------------------------- per-condition summarization -----------------------------
def summarize_condition(df, samples, focus, up_set, down_set, score_min, rec_by_id, rec_by_name):
    sc_map, st_map = discover_sample_cols(df)
    sc_cols, st_cols = group_cols_for_samples(sc_map, st_map, samples)
    if not sc_cols:
        return None, None

    work = df.copy()

    # max score across selected samples
    max_sc = work[sc_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
    work = work[max_sc >= score_min].copy()
    if work.empty:
        return None, None

    gcol = pick_gene_column(work)
    if gcol is None:
        raise ValueError("No gene column found in matrix.")

    # drop explicitly intergenic
    mask_inter = work[gcol].astype(str).str.lower().eq("intergenic")
    if mask_inter.any():
        work = work[~mask_inter].copy()
    if work.empty:
        return None, None

    # map canonical spans
    recs = []
    for g in work[gcol].astype(str):
        r = rec_by_name.get(g) or rec_by_id.get(g)
        recs.append(r)
    work["_rec"] = recs
    work = work[work["_rec"].notnull()].copy()
    if work.empty:
        return None, None

    work["_cs"] = [r["cds_start"] for r in work["_rec"]]
    work["_ce"] = [r["cds_end"] for r in work["_rec"]]
    work["_gs"] = [r["gene_start"] for r in work["_rec"]]
    work["_ge"] = [r["gene_end"] for r in work["_rec"]]
    work["_strand"] = [r["strand"] for r in work["_rec"]]
    centers = pick_center(work)

    regs = []
    for p,s,gs,cs,ce,ge in zip(centers, work["_strand"], work["_gs"], work["_cs"], work["_ce"], work["_ge"]):
        try:
            regs.append(assign_region(p, s, gs, cs, ce, ge))
        except Exception:
            regs.append("intergenic")
    work["_region"] = regs
    work = work[work["_region"].isin(REGIONS)].copy()
    if work.empty: 
        return None, None

    # within-focus only
    focus_mask = work[gcol].astype(str).isin(focus)
    work = work[focus_mask].copy()
    if work.empty:
        return None, None

    up_mask   = work[gcol].astype(str).isin(up_set)
    down_mask = work[gcol].astype(str).isin(down_set)
    bg_mask   = work[gcol].astype(str).isin(focus - (up_set | down_set))

    # weighted value per row = mean stoich across this condition's stoich cols
    if st_cols:
        work["_w"] = work[st_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1).fillna(0.0)
    else:
        work["_w"] = 0.0

    def counts_for(mask_set):
        a = int(((mask_set) & (work["_region"]=="5UTR")).sum())
        b = int(((bg_mask)   & (work["_region"]=="5UTR")).sum())
        c = int(((mask_set) & (work["_region"]=="CDS")).sum())
        d = int(((bg_mask)   & (work["_region"]=="CDS")).sum())
        e = int(((mask_set) & (work["_region"]=="3UTR")).sum())
        f = int(((bg_mask)   & (work["_region"]=="3UTR")).sum())

        a_w = float(work.loc[(mask_set) & (work["_region"]=="5UTR"), "_w"].sum())
        b_w = float(work.loc[(bg_mask)  & (work["_region"]=="5UTR"), "_w"].sum())
        c_w = float(work.loc[(mask_set) & (work["_region"]=="CDS"),  "_w"].sum())
        d_w = float(work.loc[(bg_mask)  & (work["_region"]=="CDS"),  "_w"].sum())
        e_w = float(work.loc[(mask_set) & (work["_region"]=="3UTR"), "_w"].sum())
        f_w = float(work.loc[(bg_mask)  & (work["_region"]=="3UTR"),  "_w"].sum())
        return (a,b,c,d,e,f), (a_w,b_w,c_w,d_w,e_w,f_w)

    up_cts, up_w = counts_for(up_mask)
    down_cts, down_w = counts_for(down_mask)
    return (up_cts, up_w), (down_cts, down_w)

def flat_region_counts(triple):
    cts, wcts = triple
    out = {}
    out["5UTR"] = (cts[0], cts[1], wcts[0], wcts[1])
    out["CDS"]  = (cts[2], cts[3], wcts[2], wcts[3])
    out["3UTR"] = (cts[4], cts[5], wcts[4], wcts[5])
    return out

# ----------------------------- per-file (per-mod) -----------------------------
def compute_for_file(tsv_path, rec_by_id, rec_by_name, focus, up_focus, down_focus, score_min, global_perms):
    base = os.path.basename(tsv_path)
    mod = os.path.splitext(base.split("_")[0])[0]

    try:
        df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    except Exception as e:
        warnings.warn(f"Failed to read {tsv_path}: {e}")
        return mod, [], [], [], []

    results_real, results_perm1, results_perm2, results_perm3 = [], [], [], []

    def add_rows(rows_container, set_name, by_reg, cond_label):
        total_a = sum([by_reg[r][0] for r in REGIONS])
        total_b = sum([by_reg[r][1] for r in REGIONS])
        total_aw = sum([by_reg[r][2] for r in REGIONS])
        total_bw = sum([by_reg[r][3] for r in REGIONS])

        for region in REGIONS:
            a, b, a_w, b_w = by_reg[region]
            c = total_a - a
            d = total_b - b
            c_w = max(total_aw - a_w, 0.0)
            d_w = max(total_bw - b_w, 0.0)

            orat, lnOR, se, lo, hi = or_ci_from_counts(a,b,c,d, pseudo=0.5)
            w_or, w_lnOR, w_se, w_lo, w_hi = or_ci_from_counts(a_w, b_w, c_w, d_w, pseudo=0.5)

            rows_container.append({
                "mod": mod, "cond": cond_label, "set": set_name, "region": region,
                "a":a, "b":b, "c":c, "d":d,
                "or": orat, "lnOR": lnOR, "lnOR_se": se, "or_ci_low": lo, "or_ci_high": hi,
                "a_w": a_w, "b_w": b_w, "c_w": c_w, "d_w": d_w,
                "w_or": w_or, "w_lnOR": w_lnOR, "w_lnOR_se": w_se, "w_or_ci_low": w_lo, "w_or_ci_high": w_hi
            })

    for _ck, sample_nums, cond_label, _color in COND_ORDER:
        sums_up, sums_down = summarize_condition(df, sample_nums, focus, up_focus, down_focus, score_min, rec_by_id, rec_by_name)
        if sums_up is not None:
            add_rows(results_real, "up_focus", flat_region_counts(sums_up), cond_label)
        if sums_down is not None:
            add_rows(results_real, "down_focus", flat_region_counts(sums_down), cond_label)

        for pidx, (perm_up, pdown) in enumerate(global_perms, start=1):
            sums_up_p, sums_down_p = summarize_condition(df, sample_nums, focus, perm_up, pdown, score_min, rec_by_id, rec_by_name)
            target = results_perm1 if pidx==1 else results_perm2 if pidx==2 else results_perm3
            if sums_up_p is not None:
                add_rows(target, "up_focus", flat_region_counts(sums_up_p), cond_label)
            if sums_down_p is not None:
                add_rows(target, "down_focus", flat_region_counts(sums_down_p), cond_label)

    return mod, results_real, results_perm1, results_perm2, results_perm3

# ----------------------------- plotting -----------------------------
def save_multi(fig, basepath):
    for ext in (".png",".pdf",".svg"):
        fig.savefig(basepath + ext, dpi=150, bbox_inches="tight")

def plot_bars(df, mod, set_name, scenario_label, outbase, weighted=False):
    regions = REGIONS
    cond_labels = [c[2] for c in COND_ORDER]
    cond_color = {c[2]: c[3] for c in COND_ORDER}

    fig, ax = plt.subplots(figsize=(10,6))
    width = 0.2
    x = np.arange(len(regions))

    val_col   = "w_or"      if weighted else "or"
    ln_col    = "w_lnOR"    if weighted else "lnOR"
    se_col    = "w_lnOR_se" if weighted else "lnOR_se"
    y_label   = "Weighted odds ratio" if weighted else "Odds ratio"

    for i, clab in enumerate(cond_labels):
        sub = df[(df["cond"]==clab) & (df["region"].isin(regions))].copy()
        ys, low_err, high_err = [], [], []
        for r in regions:
            row = sub[sub["region"]==r]
            if row.empty or pd.isna(row[val_col]).all():
                ys.append(np.nan); low_err.append(np.nan); high_err.append(np.nan)
            else:
                orat = float(row[val_col].iloc[0])
                lnOR = float(row[ln_col].iloc[0])
                se   = float(row[se_col].iloc[0])
                lo = math.exp(lnOR - 1.96*se)
                hi = math.exp(lnOR + 1.96*se)
                ys.append(orat)
                low_err.append(max(orat - lo, 0.0))
                high_err.append(max(hi - orat, 0.0))
        pos = x + (i-1.5)*width
        ax.bar(pos, ys, width=width, color=cond_color[clab], label=clab,
               yerr=[low_err, high_err], capsize=3, edgecolor="black", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(regions)
    ax.axhline(1.0, color="#333333", lw=1, ls="--")
    ax.set_ylabel(y_label + " (set vs. focus-BG)")
    ax.set_title(f"{mod}: {set_name} — {scenario_label} ({'WEIGHTED' if weighted else 'UNWEIGHTED'})")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    save_multi(fig, outbase)
    plt.close(fig)

# ----------------------------- IO helpers -----------------------------
def read_focus_table(path, frac_min):
    df = pd.read_csv(path, sep="\t", header=None, names=["gene","frac"])
    df["frac"] = pd.to_numeric(df["frac"], errors="coerce")
    df = df.dropna(subset=["gene","frac"])
    focus = set(df.loc[df["frac"] > frac_min, "gene"].astype(str))
    return focus, df

def read_gene_list(path):
    s=set()
    with open(path,"r") as fh:
        for ln in fh:
            t=ln.strip()
            if t:
                s.add(t.split()[0])
    return s

# ----------------------------- main -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Real + 3 global permutations, per-condition ORs (UP/DOWN) vs focus background; plots for unweighted & weighted; saves PNG/PDF/SVG.")
    ap.add_argument("data_dir")
    ap.add_argument("gtf")
    ap.add_argument("genes_pass_tsv", help="two-column TSV: gene<TAB>fraction")
    ap.add_argument("--fraction-min", type=float, default=0.4)
    ap.add_argument("--up-list", required=True)
    ap.add_argument("--down-list", required=True)
    ap.add_argument("--outdir", default="focus_bg_with_perms")
    ap.add_argument("--score-min", type=float, default=20.0)
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--pattern", default="*_merged_union_matrix_wide__z_gt3.0_annotated.tsv")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    focus, fg = read_focus_table(args.genes_pass_tsv, args.fraction_min)
    up_all = read_gene_list(args.up_list)
    down_all = read_gene_list(args.down_list)
    up_focus = set(g for g in up_all if g in focus)
    down_focus = set(g for g in down_all if g in focus)
    print(f"[INFO] focus genes: {len(focus)}")
    print(f"[INFO] up_focus={len(up_focus)}, down_focus={len(down_focus)}")

    with open(os.path.join(args.outdir,"up_focus_genes.txt"),"w") as fh:
        for g in sorted(up_focus): fh.write(g+"\n")
    with open(os.path.join(args.outdir,"down_focus_genes.txt"),"w") as fh:
        for g in sorted(down_focus): fh.write(g+"\n")

    print("[INFO] parsing GTF ...")
    rec_by_id, rec_by_name = parse_gtf_canonical(args.gtf)
    print(f"[INFO] genes with canonical CDS (by gene_id): {len(rec_by_id)} ; (by gene_name): {len(rec_by_name)}")

    rng = np.random.default_rng(args.seed)
    n_up = len(up_focus); n_down = len(down_focus)
    focus_list = sorted(list(focus))

    def one_perm_pair():
        idx_up = rng.choice(len(focus_list), size=n_up, replace=False)
        perm_up = set(focus_list[i] for i in idx_up)
        remaining = [g for i,g in enumerate(focus_list) if i not in idx_up]
        if n_down > len(remaining):
            perm_down = set(rng.choice(focus_list, size=n_down, replace=False))
        else:
            perm_down = set(rng.choice(remaining, size=n_down, replace=False))
        return perm_up, perm_down

    global_perms = [one_perm_pair() for _ in range(3)]
    for i,(pu,pdown) in enumerate(global_perms, start=1):
        with open(os.path.join(args.outdir, f"perm{i}_up_genes.txt"),"w") as fh:
            for g in sorted(pu): fh.write(g+"\n")
        with open(os.path.join(args.outdir, f"perm{i}_down_genes.txt"),"w") as fh:
            for g in sorted(pdown): fh.write(g+"\n")

    tsvs = sorted(set(
        glob.glob(os.path.join(args.data_dir, "*", args.pattern)) +
        glob.glob(os.path.join(args.data_dir, args.pattern))
    ))
    if not tsvs:
        print("[WARN] no TSVs found; check data_dir/pattern"); return

    for tsv in tsvs:
        print(f"[INFO] read {tsv}")
        mod, real, p1, p2, p3 = compute_for_file(
            tsv, rec_by_id, rec_by_name, focus, up_focus, down_focus, args.score_min, global_perms
        )

        def to_df(rows):
            if not rows:
                return pd.DataFrame(columns=[
                    "mod","cond","set","region","a","b","c","d",
                    "or","lnOR","lnOR_se","or_ci_low","or_ci_high",
                    "a_w","b_w","c_w","d_w","w_or","w_lnOR","w_lnOR_se","w_or_ci_low","w_or_ci_high"
                ])
            return pd.DataFrame(rows).sort_values(["set","region","cond"])

        df_real = to_df(real); df_p1 = to_df(p1); df_p2 = to_df(p2); df_p3 = to_df(p3)
        base = os.path.join(args.outdir, f"{mod}")
        df_real.to_csv(base + "__REAL.tsv", sep="\t", index=False)
        df_p1.to_csv(base + "__PERM1.tsv", sep="\t", index=False)
        df_p2.to_csv(base + "__PERM2.tsv", sep="\t", index=False)
        df_p3.to_csv(base + "__PERM3.tsv", sep="\t", index=False)

        # Plots: UP/DOWN x (REAL, PERM1-3) x (UNWEIGHTED, WEIGHTED)
        for set_name in ["up_focus","down_focus"]:
            for scen_label, dfX, tag in [
                ("REAL", df_real, "REAL"),
                ("PERM1", df_p1, "PERM1"),
                ("PERM2", df_p2, "PERM2"),
                ("PERM3", df_p3, "PERM3"),
            ]:
                sub = dfX[dfX["set"]==set_name].copy()
                if sub.empty:
                    warnings.warn(f"{mod} {set_name} {scen_label}: empty, skipping plots.")
                    continue
                # unweighted
                outbase_u = os.path.join(args.outdir, f"{mod}__{set_name.upper()}__{tag}__UNWEIGHTED")
                plot_bars(sub, mod, set_name.replace("_"," ").upper(), scen_label, outbase_u, weighted=False)
                # weighted
                outbase_w = os.path.join(args.outdir, f"{mod}__{set_name.upper()}__{tag}__WEIGHTED")
                plot_bars(sub, mod, set_name.replace("_"," ").upper(), scen_label, outbase_w, weighted=True)

if __name__ == "__main__":
    main()

