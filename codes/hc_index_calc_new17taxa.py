#!/usr/bin/env python3

import argparse
from typing import List, Tuple

import pandas as pd
import numpy as np

# C-tax IndValue 0.85, H-tax IndValue 0.9
H_TAXA = [
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__UCG-010;g__UCG-010",
"d__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Spirochaetales;f__Spirochaetaceae;g__Treponema",
"d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotellaceae_UCG-004",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Ruminococcus",
"d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Rikenellaceae_RC9_gut_group",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__NK4A214_group",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Christensenellales;f__Christensenellaceae;g__Christensenellaceae_R-7_group",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__[Eubacterium]_coprostanoligenes_group;g__[Eubacterium]_coprostanoligenes_group",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__UCG-005",
"d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotellaceae_UCG-003",
"d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__p-251-o5;g__p-251-o5",
"d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotellaceae_UCG-001",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__UCG-002",
]

C_TAXA = [
"d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Peptostreptococcaceae;g__Peptoclostridium",
"d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia",
]


def load_table(path: str) -> pd.DataFrame:
    """
    Read table.csv in the format:
      first column = taxon
      remaining columns = samples
    Return sample x taxon dataframe.
    """
    df = pd.read_csv(path, dtype=str)
    if df.shape[1] < 2:
        raise ValueError("Input table must have at least 2 columns: taxon + samples")

    taxon_col = df.columns[0]
    df = df.set_index(taxon_col)

    # numeric conversion
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # transpose: sample x taxon
    mat = df.T
    mat.index = mat.index.astype(str)
    return mat


def compute_hc_index(
    mat: pd.DataFrame,
    h_taxa: List[str],
    c_taxa: List[str],
    eps: float = 5e-5,
    detection_limit: float = 1e-4,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    h_present = [x for x in h_taxa if x in mat.columns]
    c_present = [x for x in c_taxa if x in mat.columns]

    if len(h_present) == 0:
        raise ValueError("No H-taxa found in input table")
    if len(c_present) == 0:
        raise ValueError("No C-taxa found in input table")

    # =========================
    # Detection limit handling
    # Values below detection_limit are treated as zero.
    # Default rationale: at 10,000 reads, 1 read corresponds to 1e-4.
    # =========================
    if detection_limit is not None and detection_limit > 0:
        mat = mat.where(mat >= detection_limit, 0.0)

    h_sum = mat[h_present].sum(axis=1)
    c_sum = mat[c_present].sum(axis=1)

    hc_index = np.log10((h_sum + eps) / (c_sum + eps))

    # =========================
    # Contribution calculation (relative contribution within H and C)
    # =========================
    h_contrib = mat[h_present].div(h_sum.replace(0, np.nan), axis=0)
    c_contrib = mat[c_present].div(c_sum.replace(0, np.nan), axis=0)

    # mean contribution per taxa (global importance)
    h_importance = h_contrib.mean(axis=0, skipna=True)
    c_importance = c_contrib.mean(axis=0, skipna=True)

    out = pd.DataFrame({
        "SampleID": mat.index,
        "H_sum": h_sum.values,
        "C_sum": c_sum.values,
        "HC_index": hc_index.values,
    })

    # Determine dominant taxa based on raw abundance within H and C groups
    # If the group sum is zero, set dominant taxa to NaN
    h_dominant = mat[h_present].idxmax(axis=1)
    c_dominant = mat[c_present].idxmax(axis=1)

    out["H_dominant_taxa"] = np.where(h_sum.to_numpy() > 0, h_dominant.to_numpy(), np.nan)
    out["C_dominant_taxa"] = np.where(c_sum.to_numpy() > 0, c_dominant.to_numpy(), np.nan)

    return out, h_present, c_present, h_importance, c_importance


def main():
    parser = argparse.ArgumentParser(
        description="Compute HC-index from genus abundance table"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input table.csv (first column = taxon, remaining columns = samples)"
    )
    parser.add_argument(
        "--output", "-o",
        default="HC_index_11taxa.csv",
        help="Output CSV file"
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=5e-5,
        help="Pseudo-count added to H_sum and C_sum before log10 ratio"
    )
    parser.add_argument(
        "--detection-limit",
        type=float,
        default=1e-4,
        help="Values below this relative abundance are treated as zero; use 0 to disable"
    )

    args = parser.parse_args()

    mat = load_table(args.input)
    out_df, h_present, c_present, h_importance, c_importance = compute_hc_index(
        mat=mat,
        h_taxa=H_TAXA,
        c_taxa=C_TAXA,
        eps=args.eps,
        detection_limit=args.detection_limit,
    )

    out_df.to_csv(args.output, index=False)

    # =========================
    # Save contribution summary
    # =========================
    contrib_df = pd.DataFrame({
        "taxon": list(h_importance.index) + list(c_importance.index),
        "contribution": list(h_importance.values) + list(c_importance.values),
        "group": ["H"] * len(h_importance) + ["C"] * len(c_importance)
    }).sort_values("contribution", ascending=False)

    contrib_out = args.output.replace(".csv", "_contribution.csv")
    contrib_df.to_csv(contrib_out, index=False)

    print(f"Saved contribution summary: {contrib_out}")

    print(f"Input samples: {mat.shape[0]}")
    print(f"H taxa found: {len(h_present)} / {len(H_TAXA)}")
    print(f"C taxa found: {len(c_present)} / {len(C_TAXA)}")
    print(f"Epsilon: {args.eps}")
    print(f"Detection limit: {args.detection_limit}")
    print(f"Samples with H_sum > 0: {(out_df['H_sum'] > 0).sum()} / {len(out_df)}")
    print(f"Samples with C_sum > 0: {(out_df['C_sum'] > 0).sum()} / {len(out_df)}")
    print(f"Non-NA H_dominant_taxa: {out_df['H_dominant_taxa'].notna().sum()} / {len(out_df)}")
    print(f"Non-NA C_dominant_taxa: {out_df['C_dominant_taxa'].notna().sum()} / {len(out_df)}")

    missing_h = [x for x in H_TAXA if x not in h_present]
    missing_c = [x for x in C_TAXA if x not in c_present]

    if missing_h:
        print("\nH taxa not found:")
        for x in missing_h:
            print(f"  - {x}")

    if missing_c:
        print("\nC taxa not found:")
        for x in missing_c:
            print(f"  - {x}")

    print(f"\nSaved: {args.output}")
    print(out_df.head())


if __name__ == "__main__":
    main()