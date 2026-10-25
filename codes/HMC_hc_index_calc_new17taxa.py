#!/usr/bin/env python3

import argparse
from typing import List, Tuple

import numpy as np
import pandas as pd

# Compendium-specific taxon names in original dot format
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


def load_compendium_table(path: str) -> pd.DataFrame:
    """
    Read Human Microbiome Compendium taxonomic_table.csv(.gz)
    Expected format:
      col1 = row number (unnamed)
      col2 = sample
      col3.. = taxa
    Return:
      sample x taxon dataframe
    """
    df = pd.read_csv(path, compression="infer", low_memory=False)

    if "sample" not in df.columns:
        raise ValueError("Column 'sample' was not found. Please check the input file.")

    # Keep sample column + all taxon columns
    sample_ids = df["sample"].astype(str)

    meta_cols = []
    for col in df.columns:
        if col == "sample":
            continue
        # typical first unnamed index column
        if str(col).startswith("Unnamed:"):
            meta_cols.append(col)

    taxon_cols = [c for c in df.columns if c not in ["sample"] + meta_cols]

    if len(taxon_cols) == 0:
        raise ValueError("No taxon columns found.")

    mat = df[taxon_cols].copy()
    mat.index = sample_ids

    # numeric conversion
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Keep original compendium dot-format taxon names as-is.
    # If duplicated taxa somehow exist, merge them by summing columns.
    if pd.Index(mat.columns).duplicated().any():
        mat = mat.T.groupby(level=0).sum().T

    return mat


def compute_hc_index(
    mat: pd.DataFrame,
    h_taxa: List[str],
    c_taxa: List[str],
    eps: float = 5e-5,
) -> Tuple[pd.DataFrame, List[str], List[str], pd.Series, pd.Series]:
    h_present = [x for x in h_taxa if x in mat.columns]
    c_present = [x for x in c_taxa if x in mat.columns]

    if len(h_present) == 0:
        raise ValueError(f"No H-taxa found in input table. Example columns: {list(mat.columns[:10])}")
    if len(c_present) == 0:
        raise ValueError(f"No C-taxa found in input table. Example columns: {list(mat.columns[:10])}")

    # convert to relative abundance (per-sample normalization)
    row_sums = mat.sum(axis=1)
    mat = mat.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)

    h_sum = mat[h_present].sum(axis=1)
    c_sum = mat[c_present].sum(axis=1)
    hc_index = np.log10((h_sum + eps) / (c_sum + eps))

    h_contrib = mat[h_present].div(h_sum.replace(0, np.nan), axis=0)
    c_contrib = mat[c_present].div(c_sum.replace(0, np.nan), axis=0)

    h_importance = h_contrib.mean(axis=0, skipna=True)
    c_importance = c_contrib.mean(axis=0, skipna=True)

    out = pd.DataFrame({
        "SampleID": mat.index,
        "H_sum": h_sum.values,
        "C_sum": c_sum.values,
        "HC_index": hc_index.values,
    })

    h_dominant = mat[h_present].idxmax(axis=1)
    c_dominant = mat[c_present].idxmax(axis=1)

    out["H_dominant_taxa"] = np.where(h_sum.to_numpy() > 0, h_dominant.to_numpy(), np.nan)
    out["C_dominant_taxa"] = np.where(c_sum.to_numpy() > 0, c_dominant.to_numpy(), np.nan)

    return out, h_present, c_present, h_importance, c_importance


def save_found_taxa(output_csv: str, h_present: List[str], c_present: List[str]) -> None:
    found_df = pd.DataFrame({
        "group": ["H"] * len(h_present) + ["C"] * len(c_present),
        "taxon": h_present + c_present,
    })
    found_out = output_csv.replace(".csv", "_found_taxa.csv")
    found_df.to_csv(found_out, index=False)
    print(f"Saved found taxa list: {found_out}")


def save_missing_taxa(output_csv: str, h_present: List[str], c_present: List[str]) -> None:
    missing_h = [x for x in H_TAXA if x not in h_present]
    missing_c = [x for x in C_TAXA if x not in c_present]

    miss_df = pd.DataFrame({
        "group": ["H"] * len(missing_h) + ["C"] * len(missing_c),
        "taxon": missing_h + missing_c,
    })
    miss_out = output_csv.replace(".csv", "_missing_taxa.csv")
    miss_df.to_csv(miss_out, index=False)
    print(f"Saved missing taxa list: {miss_out}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute HC-index from Human Microbiome Compendium taxonomic_table.csv(.gz)"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input taxonomic_table.csv or taxonomic_table.csv.gz"
    )
    parser.add_argument(
        "--output", "-o",
        default="HC_index_compendium.csv",
        help="Output CSV file"
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=1e-4,
        help="Pseudo-count added to H_sum and C_sum before log10 ratio"
    )

    args = parser.parse_args()

    print("Loading compendium taxonomic table...")
    mat = load_compendium_table(args.input)
    print(f"Loaded matrix: {mat.shape[0]} samples x {mat.shape[1]} taxa")

    out_df, h_present, c_present, h_importance, c_importance = compute_hc_index(
        mat=mat,
        h_taxa=H_TAXA,
        c_taxa=C_TAXA,
        eps=args.eps,
    )

    out_df.to_csv(args.output, index=False)

    contrib_df = pd.DataFrame({
        "taxon": list(h_importance.index) + list(c_importance.index),
        "contribution": list(h_importance.values) + list(c_importance.values),
        "group": ["H"] * len(h_importance) + ["C"] * len(c_importance)
    }).sort_values(["group", "contribution"], ascending=[True, False])

    contrib_out = args.output.replace(".csv", "_contribution.csv")
    contrib_df.to_csv(contrib_out, index=False)
    print(f"Saved contribution summary: {contrib_out}")

    save_found_taxa(args.output, h_present, c_present)
    save_missing_taxa(args.output, h_present, c_present)

    print(f"Input samples: {mat.shape[0]}")
    print(f"Input taxa: {mat.shape[1]}")
    print(f"H taxa found: {len(h_present)} / {len(H_TAXA)}")
    print(f"C taxa found: {len(c_present)} / {len(C_TAXA)}")
    print(f"Samples with H_sum > 0: {(out_df['H_sum'] > 0).sum()} / {len(out_df)}")
    print(f"Samples with C_sum > 0: {(out_df['C_sum'] > 0).sum()} / {len(out_df)}")
    print(f"Non-NA H_dominant_taxa: {out_df['H_dominant_taxa'].notna().sum()} / {len(out_df)}")
    print(f"Non-NA C_dominant_taxa: {out_df['C_dominant_taxa'].notna().sum()} / {len(out_df)}")

    print(f"\nSaved: {args.output}")
    print(out_df.head())


if __name__ == "__main__":
    main()