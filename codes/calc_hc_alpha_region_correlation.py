#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt


H_TAXA = [
    "Bacteria.Spirochaetota.Spirochaetia.Spirochaetales.Spirochaetaceae.Treponema",
    "Bacteria.Bacillota.Clostridia.Oscillospirales.Ruminococcaceae.Ruminococcus",
    "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Rikenellaceae.Rikenellaceae RC9 gut group",
    "Bacteria.Bacillota.Clostridia.Oscillospirales.Oscillospiraceae.NK4A214 group",
    "Bacteria.Bacillota.Clostridia.Christensenellales.Christensenellaceae.Christensenellaceae R-7 group",
    "Bacteria.Bacillota.Clostridia.Oscillospirales.Oscillospiraceae.UCG-005",
    "Bacteria.Bacillota.Clostridia.Oscillospirales.Oscillospiraceae.UCG-002",
    "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotellaceae UCG-001",
    "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotellaceae UCG-003",
    "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotellaceae UCG-004",
    "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.p-251-o5.Incertae Sedis",
    "Bacteria.Bacillota.Clostridia.Oscillospirales.[Eubacterium] coprostanoligenes group.Incertae Sedis",
    "Bacteria.Bacillota.Clostridia.Oscillospirales.UCG-010.Incertae Sedis",
]

C_TAXA = [
    "Bacteria.Fusobacteriota.Fusobacteriia.Fusobacteriales.Fusobacteriaceae.Fusobacterium",
    "Bacteria.Bacillota.Clostridia.Peptostreptococcales-Tissierellales.Peptostreptococcaceae.Peptoclostridium",
    "Bacteria.Pseudomonadota.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae.Escherichia-Shigella",
    "Bacteria.Bacillota.Clostridia.Lachnospirales.Lachnospiraceae.Blautia",
]


def load_taxonomic_table(path):
    df = pd.read_csv(path, compression="infer", low_memory=False)

    if "sample" not in df.columns:
        raise ValueError("Column 'sample' was not found.")

    sample_ids = df["sample"].astype(str)

    meta_cols = ["sample"]
    meta_cols += [c for c in df.columns if str(c).startswith("Unnamed:")]

    taxon_cols = [c for c in df.columns if c not in meta_cols]

    mat = df[taxon_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat.index = sample_ids

    if pd.Index(mat.columns).duplicated().any():
        mat = mat.T.groupby(level=0).sum().T

    return mat


def compute_chao1(counts):
    counts = counts.astype(float)
    observed = (counts > 0).sum(axis=1)
    f1 = (counts == 1).sum(axis=1)
    f2 = (counts == 2).sum(axis=1)

    chao1 = observed + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    return chao1


def compute_shannon(counts):
    row_sums = counts.sum(axis=1).replace(0, np.nan)
    rel = counts.div(row_sums, axis=0).fillna(0.0)

    rel_nonzero = rel.replace(0, np.nan)
    shannon = -(rel_nonzero * np.log(rel_nonzero)).sum(axis=1, skipna=True)
    return shannon


def compute_hc_index(counts, eps=5e-5, detection_limit=1e-4):
    row_sums = counts.sum(axis=1).replace(0, np.nan)
    rel = counts.div(row_sums, axis=0).fillna(0.0)

    # detection limit未満は0扱い
    rel = rel.mask(rel < detection_limit, 0.0)

    h_present = [x for x in H_TAXA if x in rel.columns]
    c_present = [x for x in C_TAXA if x in rel.columns]

    if len(h_present) == 0:
        raise ValueError("No H taxa found.")
    if len(c_present) == 0:
        raise ValueError("No C taxa found.")

    h_sum = rel[h_present].sum(axis=1)
    c_sum = rel[c_present].sum(axis=1)

    hc_index = np.log10((h_sum + eps) / (c_sum + eps))

    return hc_index, h_sum, c_sum, h_present, c_present


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]

    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)

    out = np.empty(n)
    out[order] = adj
    return out


def pick_column(df, requested, candidates, label):
    """Return a column name, using the requested name if present, otherwise trying candidates."""
    if requested != "auto" and requested in df.columns:
        return requested

    for c in candidates:
        if c in df.columns:
            print(f"{label} column '{requested}' was not found. Using '{c}' instead.")
            return c

    lower_map = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            picked = lower_map[c.lower()]
            print(f"{label} column '{requested}' was not found. Using '{picked}' instead.")
            return picked

    raise ValueError(
        f"{label} column was not found. Requested='{requested}'. "
        f"Tried candidates={candidates}. Available columns={list(df.columns)[:50]}"
    )


# Helper for normalizing sample IDs for metadata merge
def normalize_sample_id_for_merge(x):
    """Normalize sample IDs for metadata merge.

    Examples:
    - PRJDB10485_DRR243823 -> DRR243823
    - PRJNAxxxx_SRR123456 -> SRR123456
    - DRR304700 -> DRR304700
    """
    s = str(x).strip()
    if "_" in s:
        last = s.split("_")[-1]
        if last.startswith(("SRR", "ERR", "DRR", "SRX", "ERX", "DRX", "SRS", "ERS", "DRS")):
            return last
    return s


# Additional helper functions for adult-human stool filtering and trend line
def normalize_text_for_filter(x) -> str:
    if pd.isna(x):
        return "NA"
    s = str(x).strip()
    if s == "" or s.lower() in ["nan", "none", "na"]:
        return "NA"
    return s


def is_adult_human_subject(subject) -> bool:
    s = normalize_text_for_filter(subject).lower()

    if s == "na":
        return False

    exclude_keywords = [
        "mouse", "mice", "rhesus", "macaque",
        "younger than 2", "infant", "infants", "newborn", "baby", "babies",
        "children", "child", "pediatric", "paediatric", "adolescent", "adolescents",
        "toddlers", "toddler", "youth", "preterm", "very low birth weight",
        "female infant", "young children", "girls",
        "all ages", "newborn to centenarian", "healthy chinese ages 3-100+",
        "indian families"
    ]
    if any(k in s for k in exclude_keywords):
        return False

    include_keywords = [
        "adult", "adults", "elderly", "centenarian", "centenarians",
        "women", "woman", "men", "females", "female", "healthy humans",
        "healthy individuals", "healthy donors", "healthy western individuals",
        "patients", "patiens", "volunteers", "individuals", "donors",
        "long-term acute care", "hepatic cirrhosis", "cirrhotic", "crc",
        "ibd", "ibs", "adenoma", "carcinoma", "cancer", "diabetic",
        "obese", "overweight", "anorexia", "chinese adults", "finnish adult men",
        "african american men"
    ]
    return any(k in s for k in include_keywords)


def add_linear_fit(ax, x, y):
    """Add a least-squares linear trend line to an axis."""
    x = pd.to_numeric(pd.Series(x), errors="coerce")
    y = pd.to_numeric(pd.Series(y), errors="coerce")
    ok = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return
    if x[ok].nunique() < 2:
        return
    slope, intercept = np.polyfit(x[ok], y[ok], 1)
    xs = np.linspace(x[ok].min(), x[ok].max(), 100)
    ys = slope * xs + intercept
    ax.plot(xs, ys, linestyle="--", linewidth=1.5, color="black", label="linear fit")
    ax.legend(frameon=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="taxonomic_table.csv(.gz)")
    parser.add_argument("-m", "--metadata", required=True, help="metadata TSV/CSV containing sample and region")
    parser.add_argument("-p", "--projects", required=True, help="projects.csv containing project-level sample_type and subjects")
    parser.add_argument("-o", "--out_prefix", default="HC_alpha_region")
    parser.add_argument("--sample_col", default="auto", help="Metadata sample ID column. Use 'auto' to detect automatically.")
    parser.add_argument("--region_col", default="region")
    parser.add_argument("--eps", type=float, default=5e-5)
    parser.add_argument("--detection_limit", type=float, default=1e-4)
    parser.add_argument("--min_n", type=int, default=30)
    args = parser.parse_args()

    print("Loading taxonomic table...")
    counts = load_taxonomic_table(args.input)
    print(f"Loaded: {counts.shape[0]} samples x {counts.shape[1]} taxa")

    print("Computing HC-index, Chao1, Shannon...")
    hc_index, h_sum, c_sum, h_present, c_present = compute_hc_index(
        counts,
        eps=args.eps,
        detection_limit=args.detection_limit,
    )

    alpha_df = pd.DataFrame({
        "sample": counts.index.astype(str),
        "HC_index": hc_index.values,
        "H_sum": h_sum.values,
        "C_sum": c_sum.values,
        "Chao1": compute_chao1(counts).values,
        "Shannon": compute_shannon(counts).values,
        "total_counts": counts.sum(axis=1).values,
        "observed_taxa": (counts > 0).sum(axis=1).values,
    })

    meta = pd.read_csv(args.metadata, sep=None, engine="python")
    print("Metadata columns:", list(meta.columns)[:30])

    sample_col = pick_column(
        meta,
        args.sample_col,
        candidates=[
            # taxonomic_table sample IDs in this dataset are often run IDs, so prefer srr/run first.
            "srr", "SRR", "run", "Run",
            "sample", "Sample", "sample_id", "SampleID", "#SampleID", "sample_name",
            "srs", "SRS", "id", "ID"
        ],
        label="Sample ID",
    )

    region_col = pick_column(
        meta,
        args.region_col,
        candidates=[
            "region", "Region", "continent", "Continent", "geo_region",
            "geo_loc_name", "country", "Country", "location", "Location"
        ],
        label="Region",
    )

    meta[sample_col] = meta[sample_col].astype(str)

    # The taxonomic table may contain sample IDs like PRJDB10485_DRR243823,
    # while metadata often stores the run ID only, e.g. DRR243823.
    alpha_df["merge_sample_id"] = alpha_df["sample"].map(normalize_sample_id_for_merge)
    meta["merge_sample_id"] = meta[sample_col].map(normalize_sample_id_for_merge)

    print(f"Using metadata sample column: {sample_col}")
    print(f"Using metadata region column: {region_col}")
    print("First taxonomic table sample IDs:", list(alpha_df["sample"].head(10)))
    print("First normalized taxonomic IDs:", list(alpha_df["merge_sample_id"].head(10)))
    print("First metadata sample IDs:", list(meta[sample_col].head(10)))
    print("First normalized metadata IDs:", list(meta["merge_sample_id"].head(10)))

    meta_cols_for_merge = ["merge_sample_id", region_col]
    for optional_col in ["project", "library_strategy", "library_source", "geo_loc_name", "iso", "srs", "srr"]:
        if optional_col in meta.columns and optional_col not in meta_cols_for_merge:
            meta_cols_for_merge.append(optional_col)

    meta_for_merge = meta[meta_cols_for_merge].drop_duplicates(subset=["merge_sample_id"]).copy()

    merged = alpha_df.merge(
        meta_for_merge,
        on="merge_sample_id",
        how="left",
    )

    n_total = len(merged)
    n_matched = merged[region_col].notna().sum()
    n_unmatched = n_total - n_matched
    print("Merge check:")
    print(f"  Total taxonomic samples: {n_total}")
    print(f"  Matched metadata samples: {n_matched}")
    print(f"  Unmatched samples: {n_unmatched}")

    if n_total > 0 and n_matched / n_total < 0.5:
        print("WARNING: Less than 50% of samples matched metadata.")
        print("         Check whether --sample_col should be 'srr' or 'srs'.")
        print("         The script now normalizes IDs such as PRJDB10485_DRR243823 -> DRR243823.")
        print("         Example: --sample_col srr")

    merged[region_col] = merged[region_col].fillna("unknown")

    print("Region counts after merge:")
    print(merged[region_col].value_counts(dropna=False).head(20))

    print("Loading projects for adult-human stool filtering...")
    projects = pd.read_csv(args.projects, dtype=str)
    required_project_cols = {"project", "sample_type", "subjects"}
    missing_project_cols = required_project_cols - set(projects.columns)
    if missing_project_cols:
        raise ValueError(f"Projects file is missing required columns: {missing_project_cols}")

    if "project" not in merged.columns:
        raise ValueError("Column 'project' was not found after metadata merge. Cannot merge projects.csv.")

    before_project_merge = len(merged)
    merged = merged.merge(projects, on="project", how="left")
    print(f"Rows after project metadata merge: {len(merged)} / {before_project_merge}")

    if "library_strategy" in merged.columns:
        before_amplicon = len(merged)
        merged = merged[merged["library_strategy"].astype(str).str.upper() == "AMPLICON"].copy()
        print(f"After AMPLICON filter: {len(merged)} / {before_amplicon}")

    merged["sample_type_clean"] = merged["sample_type"].apply(normalize_text_for_filter)
    before_stool = len(merged)
    merged = merged[merged["sample_type_clean"].str.lower() == "stool"].copy()
    print(f"After stool filter: {len(merged)} / {before_stool}")

    merged["subjects_clean"] = merged["subjects"].apply(normalize_text_for_filter)
    before_adult = len(merged)
    merged = merged[merged["subjects_clean"].apply(is_adult_human_subject)].copy()
    print(f"After adult-human subject filter: {len(merged)} / {before_adult}")

    # =========================
    # Exclude projects with high proportion of near-zero HC_index
    # =========================
    print("Filtering zero-like HC projects...")

    tmp_proj = merged.dropna(subset=["project", "HC_index"]).copy()

    proj_summary = (
        tmp_proj.assign(is_zero_like=lambda x: x["HC_index"].abs() < 0.01)
        .groupby("project")
        .agg(
            n=("HC_index", "count"),
            zero_like_ratio=("is_zero_like", "mean")
        )
        .reset_index()
    )

    exclude_projects = proj_summary.loc[
        proj_summary["zero_like_ratio"] >= 0.20,
        "project"
    ].tolist()

    print(f"Excluding {len(exclude_projects)} projects (high zero-like HC)")

    before_proj_filter = len(merged)
    merged = merged[~merged["project"].isin(exclude_projects)].copy()
    print(f"After zero-like project filter: {len(merged)} / {before_proj_filter}")

    print("Region counts after adult-human stool filtering:")
    print(merged[region_col].value_counts(dropna=False).head(20))

    # =========================
    # Publication-oriented region filter
    # =========================
    # Remove unknown regions and keep only regions with n >= --min_n.
    before_unknown_filter = len(merged)
    merged = merged[merged[region_col].astype(str).str.lower() != "unknown"].copy()
    print(f"After removing unknown regions: {len(merged)} / {before_unknown_filter}")

    region_counts_for_filter = merged[region_col].value_counts(dropna=False)
    valid_regions = region_counts_for_filter[region_counts_for_filter >= args.min_n].index.tolist()
    excluded_regions = region_counts_for_filter[region_counts_for_filter < args.min_n]

    if len(excluded_regions) > 0:
        print("Excluded regions with n < min_n:")
        for region_name, n_region in excluded_regions.items():
            print(f"  {region_name}: n={n_region}")

    before_region_filter = len(merged)
    merged = merged[merged[region_col].isin(valid_regions)].copy()
    print(f"After min_n region filter (n >= {args.min_n}): {len(merged)} / {before_region_filter}")

    print("Final region counts used for analysis:")
    print(merged[region_col].value_counts(dropna=False).head(20))

    out_alpha = f"{args.out_prefix}_HC_Chao_Shannon_with_region.csv"
    merged.to_csv(out_alpha, index=False)
    print(f"Saved: {out_alpha}")

    results = []

    for region, sub in merged.groupby(region_col):
        sub = sub.dropna(subset=["HC_index", "Chao1", "Shannon"])

        if len(sub) < args.min_n:
            continue

        for metric in ["Chao1", "Shannon"]:
            rho, p = spearmanr(sub["HC_index"], sub[metric], nan_policy="omit")

            results.append({
                "region": region,
                "metric": metric,
                "n": len(sub),
                "spearman_rho": rho,
                "p_value": p,
            })

            fig, ax = plt.subplots(figsize=(5, 4))
            ax.scatter(sub["HC_index"], sub[metric], alpha=0.4)
            add_linear_fit(ax, sub["HC_index"], sub[metric])
            ax.set_xlabel("HC-index")
            ax.set_ylabel(metric)
            ax.set_title(f"{region}: HC-index vs {metric}\nSpearman rho={rho:.3f}, p={p:.2e}")
            plt.tight_layout()

            fig_name = f"{args.out_prefix}_{region}_{metric}_scatter.png"
            fig_name = fig_name.replace("/", "_").replace(" ", "_")
            plt.savefig(fig_name, dpi=300)
            plt.close()

    res_df = pd.DataFrame(results)

    if len(res_df) > 0:
        res_df["p_adj_BH"] = bh_adjust(res_df["p_value"].values)
        res_df = res_df.sort_values(["metric", "p_adj_BH"])

    out_stats = f"{args.out_prefix}_region_correlation_stats.csv"
    res_df.to_csv(out_stats, index=False)
    print(f"Saved: {out_stats}")

    found_df = pd.DataFrame({
        "group": ["H"] * len(h_present) + ["C"] * len(c_present),
        "taxon": h_present + c_present,
    })
    found_df.to_csv(f"{args.out_prefix}_found_taxa.csv", index=False)

    print(f"H taxa found: {len(h_present)} / {len(H_TAXA)}")
    print(f"C taxa found: {len(c_present)} / {len(C_TAXA)}")


if __name__ == "__main__":
    main()