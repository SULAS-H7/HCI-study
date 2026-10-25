#!/usr/bin/env python3

import argparse
import os
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import kruskal

# Dunn検定用
try:
    import scikit_posthocs as sp
except ImportError:
    sp = None


def sanitize_filename(text: str) -> str:
    text = str(text)
    text = re.sub(r"[\\/:*?\"<>|]+", "_", text)
    text = re.sub(r"\s+", "_", text)
    return text


def load_hc(hc_path: str) -> pd.DataFrame:
    df = pd.read_csv(hc_path)
    required = {"SampleID", "HC_index"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"HC file is missing required columns: {missing}")
    return df


def load_meta(meta_path: str) -> pd.DataFrame:
    df = pd.read_csv(meta_path, sep="\t", dtype=str)

    required = {"project", "srr", "library_strategy", "region"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Metadata file is missing required columns: {missing}")

    df["SampleID"] = df["project"].astype(str) + "_" + df["srr"].astype(str)
    return df


def load_projects(projects_path: str) -> pd.DataFrame:
    df = pd.read_csv(projects_path, dtype=str)

    required = {"project", "amplicon", "bead_beating", "kit", "sample_type", "condition"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Projects file is missing required columns: {missing}")

    return df


# ======================== Inserted functions ========================

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
        "mouse", "mice", "mice!!", "rhesus", "macaque",
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


def collapse_condition(condition) -> str:
    s = normalize_text_for_filter(condition)
    sl = s.lower()

    if s == "NA":
        return "NA"

    # Colorectal cancer / adenoma / carcinogenesis related projects
    crc_keywords = [
        "colorectal cancer", "colorectal adenoma", "colorectal carcinogenesis",
        "colorectal cancer surgery", "crc", "adenoma", "carcinoma"
    ]
    if any(k in sl for k in crc_keywords):
        return "CRC/adenoma"

    # IBD / UC related projects
    ibd_keywords = [
        "ibd", "ulcerative colitis", "uc"
    ]
    if any(k in sl for k in ibd_keywords):
        return "IBD/UC"

    # IBS related projects
    if "ibs" in sl:
        return "IBS"

    return "Other"


def summarize_by_region(df: pd.DataFrame) -> pd.DataFrame:
    out = (
        df.groupby("region", dropna=False)["HC_index"]
        .agg(
            n="count",
            mean="mean",
            median="median",
            sd="std",
            min="min",
            q1=lambda x: x.quantile(0.25),
            q3=lambda x: x.quantile(0.75),
            max="max",
        )
        .reset_index()
        .sort_values("median", ascending=False)
    )
    return out


def summarize_by_group(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    tmp = df.dropna(subset=[group_col, "HC_index"]).copy()
    tmp[group_col] = tmp[group_col].astype(str)
    tmp = tmp[tmp[group_col].str.strip() != ""]

    out = (
        tmp.groupby(group_col, dropna=False)["HC_index"]
        .agg(
            n="count",
            mean="mean",
            median="median",
            sd="std",
            min="min",
            q1=lambda x: x.quantile(0.25),
            q3=lambda x: x.quantile(0.75),
            max="max",
        )
        .reset_index()
        .sort_values("median", ascending=False)
    )
    return out


def run_kruskal(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.dropna(subset=["region", "HC_index"]).copy()
    counts = tmp["region"].value_counts()
    valid_regions = counts[counts >= 2].index.tolist()
    tmp = tmp[tmp["region"].isin(valid_regions)]

    if len(valid_regions) < 2:
        return pd.DataFrame([{
            "test": "Kruskal-Wallis",
            "n_regions": len(valid_regions),
            "statistic": np.nan,
            "p_value": np.nan,
            "note": "Not enough groups with n>=2"
        }])

    groups = [tmp.loc[tmp["region"] == r, "HC_index"].values for r in valid_regions]
    stat, p = kruskal(*groups)

    return pd.DataFrame([{
        "test": "Kruskal-Wallis",
        "n_regions": len(valid_regions),
        "statistic": stat,
        "p_value": p,
        "note": ""
    }])


def run_kruskal_by_group(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    tmp = df.dropna(subset=[group_col, "HC_index"]).copy()
    tmp[group_col] = tmp[group_col].astype(str)
    tmp = tmp[tmp[group_col].str.strip() != ""]

    counts = tmp[group_col].value_counts()
    valid_groups = counts[counts >= 2].index.tolist()
    tmp = tmp[tmp[group_col].isin(valid_groups)]

    if len(valid_groups) < 2:
        return pd.DataFrame([{
            "variable": group_col,
            "test": "Kruskal-Wallis",
            "n_groups": len(valid_groups),
            "statistic": np.nan,
            "p_value": np.nan,
            "note": "Not enough groups with n>=2"
        }])

    groups = [tmp.loc[tmp[group_col] == g, "HC_index"].values for g in valid_groups]
    stat, p = kruskal(*groups)

    return pd.DataFrame([{
        "variable": group_col,
        "test": "Kruskal-Wallis",
        "n_groups": len(valid_groups),
        "statistic": stat,
        "p_value": p,
        "note": ""
    }])


def run_dunn(df: pd.DataFrame) -> pd.DataFrame:
    if sp is None:
        return pd.DataFrame([{
            "group1": np.nan,
            "group2": np.nan,
            "p_adjusted": np.nan,
            "note": "scikit-posthocs is not installed. Install with: pip install scikit-posthocs"
        }])

    tmp = df.dropna(subset=["region", "HC_index"]).copy()
    counts = tmp["region"].value_counts()
    valid_regions = counts[counts >= 2].index.tolist()
    tmp = tmp[tmp["region"].isin(valid_regions)]

    if len(valid_regions) < 2:
        return pd.DataFrame([{
            "group1": np.nan,
            "group2": np.nan,
            "p_adjusted": np.nan,
            "note": "Not enough groups with n>=2"
        }])

    dunn_mat = sp.posthoc_dunn(
        tmp,
        val_col="HC_index",
        group_col="region",
        p_adjust="bonferroni"
    )

    long_rows = []
    regions = list(dunn_mat.index)
    for i in range(len(regions)):
        for j in range(i + 1, len(regions)):
            long_rows.append({
                "group1": regions[i],
                "group2": regions[j],
                "p_adjusted": dunn_mat.iloc[i, j]
            })

    out = pd.DataFrame(long_rows).sort_values("p_adjusted")
    return out


def run_dunn_by_group(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if sp is None:
        return pd.DataFrame([{
            "variable": group_col,
            "group1": np.nan,
            "group2": np.nan,
            "p_adjusted": np.nan,
            "note": "scikit-posthocs is not installed. Install with: pip install scikit-posthocs"
        }])

    tmp = df.dropna(subset=[group_col, "HC_index"]).copy()
    tmp[group_col] = tmp[group_col].astype(str)
    tmp = tmp[tmp[group_col].str.strip() != ""]

    counts = tmp[group_col].value_counts()
    valid_groups = counts[counts >= 2].index.tolist()
    tmp = tmp[tmp[group_col].isin(valid_groups)]

    if len(valid_groups) < 2:
        return pd.DataFrame([{
            "variable": group_col,
            "group1": np.nan,
            "group2": np.nan,
            "p_adjusted": np.nan,
            "note": "Not enough groups with n>=2"
        }])

    dunn_mat = sp.posthoc_dunn(
        tmp,
        val_col="HC_index",
        group_col=group_col,
        p_adjust="bonferroni"
    )

    long_rows = []
    groups = list(dunn_mat.index)
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            long_rows.append({
                "variable": group_col,
                "group1": groups[i],
                "group2": groups[j],
                "p_adjusted": dunn_mat.iloc[i, j]
            })

    out = pd.DataFrame(long_rows).sort_values("p_adjusted")
    return out


def make_boxplot(df: pd.DataFrame, kruskal_df: pd.DataFrame, out_png: str):
    tmp = df.dropna(subset=["region", "HC_index"]).copy()

    counts = tmp["region"].value_counts()
    valid_regions = counts[counts >= 1].index.tolist()
    tmp = tmp[tmp["region"].isin(valid_regions)]

    if tmp.empty:
        raise ValueError("No data available for plotting after filtering.")

    summary = (
        tmp.groupby("region")["HC_index"]
        .median()
        .sort_values(ascending=False)
    )
    order = summary.index.tolist()

    plot_data = [tmp.loc[tmp["region"] == r, "HC_index"].values for r in order]

    plt.figure(figsize=(max(10, len(order) * 0.7), 6))
    bp = plt.boxplot(plot_data, labels=order, patch_artist=False, showfliers=False)

    # jitter points
    rng = np.random.default_rng(42)
    for i, r in enumerate(order, start=1):
        y = tmp.loc[tmp["region"] == r, "HC_index"].values
        x = rng.normal(i, 0.06, size=len(y))
        plt.plot(x, y, "o", alpha=0.18, markersize=2)

    pval = kruskal_df.loc[0, "p_value"] if not kruskal_df.empty else np.nan
    if pd.notna(pval):
        title = f"HC index by region (Kruskal-Wallis p = {pval:.3e})"
    else:
        title = "HC index by region"

    plt.title(title)
    plt.ylabel("HC_index")
    plt.xlabel("Region")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# Inserted function: make_region_boxplot_with_hue
def make_region_boxplot_with_hue(
    df: pd.DataFrame,
    hue_col: str,
    kruskal_df: pd.DataFrame,
    out_png: str,
    top_n_hue: int = 8,
):
    tmp = df.dropna(subset=["region", "HC_index"]).copy()

    if hue_col not in tmp.columns:
        print(f"Skip hue plot: {hue_col} was not found")
        return

    tmp[hue_col] = tmp[hue_col].fillna("NA").astype(str)
    tmp[hue_col] = tmp[hue_col].str.strip()
    tmp.loc[tmp[hue_col] == "", hue_col] = "NA"
    tmp.loc[tmp[hue_col].str.lower().isin(["nan", "none", "na"]), hue_col] = "NA"

    if hue_col == "condition_group":
        category_order = ["CRC/adenoma", "IBD/UC", "IBS", "Other", "NA"]
        tmp[f"{hue_col}_plot"] = pd.Categorical(tmp[hue_col], categories=category_order, ordered=True).astype(str)
        hue_order = [x for x in category_order if x in set(tmp[f"{hue_col}_plot"])]
    else:
        # Keep frequent hue categories and collapse rare ones into Other.
        hue_counts = tmp[hue_col].value_counts()
        keep_hues = hue_counts.head(top_n_hue).index.tolist()
        tmp[f"{hue_col}_plot"] = np.where(tmp[hue_col].isin(keep_hues), tmp[hue_col], "Other")

    region_order = (
        tmp.groupby("region")["HC_index"]
        .median()
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    if hue_col == "condition_group":
        # Draw background categories first, then target disease categories on top.
        draw_order = [x for x in ["Other", "NA", "CRC/adenoma", "IBD/UC", "IBS"] if x in hue_order]
    else:
        hue_order = tmp[f"{hue_col}_plot"].value_counts().index.tolist()
        draw_order = hue_order

    plt.figure(figsize=(max(11, len(region_order) * 0.85), 6.5))

    plot_data = [tmp.loc[tmp["region"] == r, "HC_index"].values for r in region_order]
    plt.boxplot(plot_data, labels=region_order, patch_artist=False, showfliers=False)

    # Color settings. For condition_group, keep Other/NA light gray and highlight target groups.
    if hue_col == "condition_group":
        color_map = {
            "CRC/adenoma": "tab:red",
            "IBD/UC": "tab:blue",
            "IBS": "tab:green",
            "Other": "lightgray",
            "NA": "gainsboro",
        }
    else:
        color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
        color_map = {h: color_cycle[i % len(color_cycle)] for i, h in enumerate(hue_order)} if color_cycle else {}

    rng = np.random.default_rng(42)
    region_pos = {r: i for i, r in enumerate(region_order, start=1)}

    for hue in draw_order:
        sub = tmp[tmp[f"{hue_col}_plot"] == hue]
        xs = []
        ys = []
        for r in region_order:
            y = sub.loc[sub["region"] == r, "HC_index"].values
            if len(y) == 0:
                continue
            x = rng.normal(region_pos[r], 0.06, size=len(y))
            xs.extend(x)
            ys.extend(y)
        if len(xs) > 0:
            if color_map:
                if hue_col == "condition_group" and hue in ["Other", "NA"]:
                    plt.scatter(xs, ys, s=8, alpha=0.10, label=hue, color=color_map.get(hue, "lightgray"), zorder=2)
                elif hue_col == "condition_group":
                    plt.scatter(xs, ys, s=18, alpha=0.70, label=hue, color=color_map.get(hue, None), edgecolors="black", linewidths=0.25, zorder=5)
                else:
                    plt.scatter(xs, ys, s=8, alpha=0.25, label=hue, color=color_map[hue], zorder=3)
            else:
                plt.scatter(xs, ys, s=8, alpha=0.25, label=hue, zorder=3)

    pval = kruskal_df.loc[0, "p_value"] if not kruskal_df.empty else np.nan
    if pd.notna(pval):
        title = f"HC index by region colored by {hue_col} (Kruskal-Wallis p = {pval:.3e})"
    else:
        title = f"HC index by region colored by {hue_col}"

    plt.title(title)
    plt.ylabel("HC_index")
    plt.xlabel("Region")
    plt.xticks(rotation=60, ha="right")
    legend_order = hue_order if hue_col == "condition_group" else draw_order
    handles, labels = plt.gca().get_legend_handles_labels()
    handle_map = dict(zip(labels, handles))
    ordered_handles = [handle_map[x] for x in legend_order if x in handle_map]
    ordered_labels = [x for x in legend_order if x in handle_map]
    plt.legend(ordered_handles, ordered_labels, title=hue_col, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


# ================== Inserted function ==================

def make_condition_focused_region_plots(df: pd.DataFrame, outdir: str):
    tmp = df.dropna(subset=["region", "HC_index", "condition_group"]).copy()

    target_conditions = ["CRC/adenoma", "IBD/UC", "IBS"]
    background_conditions = ["Other", "NA"]

    # Keep only regions where at least one target condition exists.
    target_regions = (
        tmp.loc[tmp["condition_group"].isin(target_conditions), "region"]
        .dropna()
        .unique()
        .tolist()
    )

    if len(target_regions) == 0:
        print("No regions with target condition groups. Skipping condition-focused region plots.")
        return

    tmp = tmp[tmp["region"].isin(target_regions)].copy()

    region_order = (
        tmp.groupby("region")["HC_index"]
        .median()
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    color_map = {
        "CRC/adenoma": "tab:red",
        "IBD/UC": "tab:blue",
        "IBS": "tab:green",
        "Other": "lightgray",
        "NA": "gainsboro",
    }

    def draw_one(sub_df: pd.DataFrame, condition_order: list, out_png: str, title_suffix: str):
        if sub_df.empty:
            print(f"No data for {title_suffix}. Skipping.")
            return

        plt.figure(figsize=(max(8, len(region_order) * 1.1), 6.5))

        plot_data = [sub_df.loc[sub_df["region"] == r, "HC_index"].values for r in region_order]
        plt.boxplot(plot_data, labels=region_order, patch_artist=False, showfliers=False)

        rng = np.random.default_rng(42)
        region_pos = {r: i for i, r in enumerate(region_order, start=1)}

        for cond in condition_order:
            sub = sub_df[sub_df["condition_group"] == cond]
            xs = []
            ys = []
            for r in region_order:
                y = sub.loc[sub["region"] == r, "HC_index"].values
                if len(y) == 0:
                    continue
                x = rng.normal(region_pos[r], 0.06, size=len(y))
                xs.extend(x)
                ys.extend(y)

            if len(xs) == 0:
                continue

            if cond in background_conditions:
                plt.scatter(xs, ys, s=9, alpha=0.18, label=cond, color=color_map.get(cond, "lightgray"), zorder=2)
            else:
                plt.scatter(xs, ys, s=20, alpha=0.70, label=cond, color=color_map.get(cond, None), edgecolors="black", linewidths=0.25, zorder=5)

        plt.title(f"HC index by region: {title_suffix}")
        plt.ylabel("HC_index")
        plt.xlabel("Region")
        plt.xticks(rotation=60, ha="right")
        plt.legend(title="condition_group", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
        plt.tight_layout()
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()

    target_df = tmp[tmp["condition_group"].isin(target_conditions)].copy()
    background_df = tmp[tmp["condition_group"].isin(background_conditions)].copy()

    draw_one(
        target_df,
        target_conditions,
        os.path.join(outdir, "HC_index_region_boxplot_condition_focus_target_conditions.png"),
        "target conditions only"
    )

    draw_one(
        background_df,
        background_conditions,
        os.path.join(outdir, "HC_index_region_boxplot_condition_focus_other_NA.png"),
        "Other / NA only in target-condition regions"
    )

    # Save the table used for these focused plots.
    tmp.to_csv(
        os.path.join(outdir, "HC_index_condition_focus_target_regions_table.tsv"),
        sep="\t",
        index=False
    )


def make_boxplot_by_group(df: pd.DataFrame, group_col: str, kruskal_df: pd.DataFrame, out_png: str):
    tmp = df.dropna(subset=[group_col, "HC_index"]).copy()
    tmp[group_col] = tmp[group_col].astype(str)
    tmp = tmp[tmp[group_col].str.strip() != ""]

    if group_col == "kit":
        counts = tmp[group_col].value_counts()
        top_groups = counts.head(15).index.tolist()
        tmp = tmp[tmp[group_col].isin(top_groups)]

    counts = tmp[group_col].value_counts()
    valid_groups = counts[counts >= 1].index.tolist()
    tmp = tmp[tmp[group_col].isin(valid_groups)]

    if tmp.empty:
        raise ValueError(f"No data available for plotting after filtering: {group_col}")

    summary = (
        tmp.groupby(group_col)["HC_index"]
        .median()
        .sort_values(ascending=False)
    )
    order = summary.index.tolist()

    plot_data = [tmp.loc[tmp[group_col] == g, "HC_index"].values for g in order]

    plt.figure(figsize=(max(8, len(order) * 0.7), 6))
    plt.boxplot(plot_data, labels=order, patch_artist=False, showfliers=False)

    rng = np.random.default_rng(42)
    for i, g in enumerate(order, start=1):
        y = tmp.loc[tmp[group_col] == g, "HC_index"].values
        x = rng.normal(i, 0.06, size=len(y))
        plt.plot(x, y, "o", alpha=0.18, markersize=2)

    pval = kruskal_df.loc[0, "p_value"] if not kruskal_df.empty else np.nan
    if pd.notna(pval):
        title = f"HC index by {group_col} (Kruskal-Wallis p = {pval:.3e})"
    else:
        title = f"HC index by {group_col}"

    plt.title(title)
    plt.ylabel("HC_index")
    plt.xlabel(group_col)
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Merge HC index with sample metadata and compare HC index across regions"
    )
    parser.add_argument("-i", "--hc_index", required=True, help="HC_index_compendium.csv")
    parser.add_argument("-m", "--metadata", required=True, help="sample_metadata.tsv")
    parser.add_argument("-p", "--projects", required=True, help="projects.csv")
    parser.add_argument("-o", "--outdir", default="HC_region_results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading HC index...")
    hc = load_hc(args.hc_index)

    print("Loading metadata...")
    meta = load_meta(args.metadata)

    print("Loading projects...")
    projects = load_projects(args.projects)

    print("Merging...")
    merged = hc.merge(meta, on="SampleID", how="inner")
    merged = merged.merge(projects, on="project", how="left")

    merged_out = os.path.join(args.outdir, "HC_index_with_metadata_and_projects.tsv")
    merged.to_csv(merged_out, sep="\t", index=False)

    print(f"Merged rows: {merged.shape[0]}")

    amplicon = merged[merged["library_strategy"] == "AMPLICON"].copy()
    stool = amplicon[amplicon["sample_type"].astype(str).str.lower() == "stool"].copy()

    stool["subjects_clean"] = stool["subjects"].apply(normalize_text_for_filter)
    stool["condition_clean"] = stool["condition"].apply(normalize_text_for_filter)
    stool["condition_group"] = stool["condition_clean"].apply(collapse_condition)

    before_subject_filter = stool.shape[0]
    stool = stool[stool["subjects_clean"].apply(is_adult_human_subject)].copy()
    print(f"Adult human stool rows after subjects filter: {stool.shape[0]} / {before_subject_filter}")

    # numeric safety
    stool["HC_index"] = pd.to_numeric(stool["HC_index"], errors="coerce")

    stool_prefilter_out = os.path.join(args.outdir, "HC_index_with_metadata_projects_AMPLICON_stool_before_zero_project_filter.tsv")
    stool.to_csv(stool_prefilter_out, sep="\t", index=False)

    print(f"AMPLICON stool rows before zero-project filter: {stool.shape[0]}")

    # =========================
    # Exclude projects likely affected by technical bias
    # Definition: projects where >=80% of samples have HC_index near zero.
    # Near zero is defined as abs(HC_index) < 0.2, matching exploratory checks.
    # =========================
    zero_project_summary = (
        stool.dropna(subset=["project", "HC_index"])
        .assign(is_zero_like=lambda x: x["HC_index"].abs() < 0.01)
        .groupby("project", dropna=False)
        .agg(
            n=("HC_index", "count"),
            zero_like_n=("is_zero_like", "sum"),
            zero_like_ratio=("is_zero_like", "mean"),
            median_HC=("HC_index", "median"),
            condition=("condition", lambda x: " / ".join(sorted(set([str(v) for v in x.dropna().unique()])))),
            subjects=("subjects", lambda x: " / ".join(sorted(set([str(v) for v in x.dropna().unique()])))),
            geo_loc_name=("geo_loc_name", lambda x: " / ".join(sorted(set([str(v) for v in x.dropna().unique()])))),
        )
        .reset_index()
        .sort_values(["zero_like_ratio", "n"], ascending=[False, False])
    )

    zero_project_summary_out = os.path.join(args.outdir, "HC_index_zero_like_project_summary.tsv")
    zero_project_summary.to_csv(zero_project_summary_out, sep="\t", index=False)

    excluded_projects = zero_project_summary.loc[
        zero_project_summary["zero_like_ratio"] >= 0.2,
        "project"
    ].tolist()

    print("Excluded project IDs:")
    for p in excluded_projects:
        print(f"  {p}")

    excluded_projects_out = os.path.join(args.outdir, "HC_index_excluded_zero_like_projects.tsv")
    zero_project_summary[zero_project_summary["project"].isin(excluded_projects)].to_csv(
        excluded_projects_out,
        sep="\t",
        index=False
    )

    before_zero_project_filter = stool.shape[0]
    stool = stool[~stool["project"].isin(excluded_projects)].copy()
    print(f"Excluded zero-like projects: {len(excluded_projects)}")
    print(f"AMPLICON stool rows after zero-project filter: {stool.shape[0]} / {before_zero_project_filter}")

    # =========================
    # Publication-oriented region filter
    # =========================
    # Exclude regions with low interpretability or very small sample size.
    exclude_regions = ["unknown", "Oceania"]
    before_region_filter = stool.shape[0]
    stool = stool[~stool["region"].astype(str).str.lower().isin([x.lower() for x in exclude_regions])].copy()
    print(f"Excluded regions {exclude_regions}: {stool.shape[0]} / {before_region_filter}")
    print("Final region counts used for output:")
    print(stool["region"].value_counts(dropna=False))

    stool_out = os.path.join(args.outdir, "HC_index_with_metadata_projects_AMPLICON_stool.tsv")
    stool.to_csv(stool_out, sep="\t", index=False)

    # =========================
    # Low / High HC group analysis
    # =========================
    low = stool[stool["HC_index"] < -5].copy()
    high = stool[stool["HC_index"] > 5].copy()

    print(f"Low HC samples (< -5): {low.shape[0]}")
    print(f"High HC samples (> 5): {high.shape[0]}")

    def save_group_summary(df, label, outdir):
        # region distribution (count + fraction)
        region = df["region"].value_counts().reset_index()
        region.columns = ["region", "count"]
        total = region["count"].sum()
        region["fraction"] = region["count"] / total
        region.to_csv(os.path.join(outdir, f"{label}_region_distribution.tsv"), sep="\t", index=False)

        # project distribution (with fraction)
        project = df["project"].value_counts().reset_index()
        project.columns = ["project", "count"]
        total = project["count"].sum()
        project["fraction"] = project["count"] / total
        project.to_csv(os.path.join(outdir, f"{label}_project_distribution.tsv"), sep="\t", index=False)

        # technical factors (count + fraction)
        for col in ["amplicon", "kit", "bead_beating"]:
            tmp = df[col].value_counts().reset_index()
            tmp.columns = [col, "count"]
            total = tmp["count"].sum()
            tmp["fraction"] = tmp["count"] / total
            tmp.to_csv(os.path.join(outdir, f"{label}_{col}_distribution.tsv"), sep="\t", index=False)

        # geo location (count + fraction)
        if "geo_loc_name" in df.columns:
            geo = df["geo_loc_name"].value_counts().reset_index()
            geo.columns = ["geo_loc_name", "count"]
            total = geo["count"].sum()
            geo["fraction"] = geo["count"] / total
            geo.to_csv(os.path.join(outdir, f"{label}_geo_distribution.tsv"), sep="\t", index=False)

    save_group_summary(low, "lowHC", args.outdir)
    save_group_summary(high, "highHC", args.outdir)

    # region summary and tests
    summary_df = summarize_by_region(stool)
    summary_out = os.path.join(args.outdir, "HC_index_region_summary_stool.tsv")
    summary_df.to_csv(summary_out, sep="\t", index=False)

    kw_df = run_kruskal(stool)
    kw_out = os.path.join(args.outdir, "HC_index_region_kruskal_stool.tsv")
    kw_df.to_csv(kw_out, sep="\t", index=False)

    dunn_df = run_dunn(stool)
    dunn_out = os.path.join(args.outdir, "HC_index_region_dunn_stool.tsv")
    dunn_df.to_csv(dunn_out, sep="\t", index=False)

    plot_out = os.path.join(args.outdir, "HC_index_region_boxplot_stool.png")
    make_boxplot(stool, kw_df, plot_out)

    # region boxplots colored by project-level annotations
    for hue_col in ["subjects_clean", "condition_group"]:
        hue_plot_out = os.path.join(args.outdir, f"HC_index_region_boxplot_stool_colored_by_{sanitize_filename(hue_col)}.png")
        make_region_boxplot_with_hue(stool, hue_col, kw_df, hue_plot_out, top_n_hue=8)

    # Focused plots: only regions containing CRC/adenoma, IBD/UC, or IBS projects.
    # Target condition groups and Other/NA are plotted separately for easier comparison.
    make_condition_focused_region_plots(stool, args.outdir)

    # technical factor summaries and tests
    for var in ["amplicon", "bead_beating", "kit"]:
        var_summary = summarize_by_group(stool, var)
        var_summary_out = os.path.join(args.outdir, f"HC_index_{sanitize_filename(var)}_summary_stool.tsv")
        var_summary.to_csv(var_summary_out, sep="\t", index=False)

        var_kw = run_kruskal_by_group(stool, var)
        var_kw_out = os.path.join(args.outdir, f"HC_index_{sanitize_filename(var)}_kruskal_stool.tsv")
        var_kw.to_csv(var_kw_out, sep="\t", index=False)

        var_dunn = run_dunn_by_group(stool, var)
        var_dunn_out = os.path.join(args.outdir, f"HC_index_{sanitize_filename(var)}_dunn_stool.tsv")
        var_dunn.to_csv(var_dunn_out, sep="\t", index=False)

        var_plot_out = os.path.join(args.outdir, f"HC_index_{sanitize_filename(var)}_boxplot_stool.png")
        make_boxplot_by_group(stool, var, var_kw, var_plot_out)

    print("Done.")
    print(f"Saved merged table: {merged_out}")
    print(f"Saved AMPLICON stool table before zero-project filter: {stool_prefilter_out}")
    print(f"Saved zero-like project summary: {zero_project_summary_out}")
    print(f"Saved excluded zero-like projects: {excluded_projects_out}")
    print(f"Saved AMPLICON stool table after zero-project filter: {stool_out}")
    print(f"Saved region summary: {summary_out}")
    print(f"Saved region Kruskal-Wallis result: {kw_out}")
    print(f"Saved region Dunn test result: {dunn_out}")
    print(f"Saved region boxplot: {plot_out}")
    print("Saved region boxplots colored by adult-human subjects and collapsed condition groups.")
    print("Saved condition-focused region plots for target condition regions.")
    print("Saved amplicon / bead_beating / kit summaries, Kruskal-Wallis, Dunn, and boxplots.")


if __name__ == "__main__":
    main()