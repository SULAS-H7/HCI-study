import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# =========================
# Input data from 5 cohorts
# =========================
cohorts = {
    "Original": ("/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Original1803_HC_taxa_evaluation/original_table.csv", "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Original1803_HC_taxa_evaluation/original_alpha.csv"),
    "Kyotango": ("/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Kyotango/Kyotango_table.csv", "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Kyotango/Kyotango_alpha.csv"),
    "AGP": ("/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/AGP_over3000/AGP_table.csv", "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/AGP_over3000/AGP_alpha.csv"),
    "Hadza": ("/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Hadza/hadza_table.csv", "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Hadza/hadza_alpha.csv"),
    "FD_DB": ("/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/FD_DB/SearchBestCombinationTaxa/FD_table_filtered_by_alpha.csv", "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/FD_DB/SearchBestCombinationTaxa/FD_alpha.csv"),
    "Chinese": ("/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/PRJNA385551_china_cohort/Level6.csv", "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/PRJNA385551_china_cohort/alpha.csv")
}

# =========================
# Calculation of HC-index
# =========================
# C-tax IndValue 0.85, H-tax IndValue 0.9
H_taxa = [
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

C_taxa = [
"d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Peptostreptococcales-Tissierellales;f__Peptostreptococcaceae;g__Peptoclostridium",
"d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella",
"d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia",
]

eps = 5e-5

results = []

def detect_context_column(alpha_df):
    candidates = ["context", "Context", "source_alpha", "Source", "dataset", "Dataset"]
    for col in candidates:
        if col in alpha_df.columns:
            return col
    return None

# =========================
# Processing each cohort
# =========================
for name, (table_file, alpha_file) in cohorts.items():
    print(f"\nProcessing {name}")

    table = pd.read_csv(table_file)
    alpha = pd.read_csv(alpha_file)

    genus_col = table.columns[0]
    mat = table.set_index(genus_col).T
    mat.index.name = "SampleID"

    mat.index = mat.index.astype(str)
    alpha["SampleID"] = alpha["SampleID"].astype(str)

    common = mat.index.intersection(alpha["SampleID"])
    mat = mat.loc[common].apply(pd.to_numeric, errors="coerce").fillna(0)
    alpha = alpha.set_index("SampleID").loc[common]

    # Use only detected taxa
    H_present = [x for x in H_taxa if x in mat.columns]
    C_present = [x for x in C_taxa if x in mat.columns]

    if len(H_present) == 0 or len(C_present) == 0:
        print(f"Skip {name} (no taxa found)")
        continue

    context_col = detect_context_column(alpha.reset_index())

    if context_col is None:
        subsets = [(name, alpha.index)]
    else:
        subsets = []
        alpha_reset = alpha.reset_index()
        for ctx_value, sub_df in alpha_reset.groupby(context_col):
            subset_name = f"{name}_{ctx_value}"
            subsets.append((subset_name, sub_df["SampleID"].astype(str).tolist()))
        # 全体集計も残す
        subsets.insert(0, (name, alpha.index.tolist()))

    for subset_name, subset_ids in subsets:
        subset_ids = [x for x in subset_ids if x in mat.index and x in alpha.index]
        if len(subset_ids) < 3:
            print(f"Skip {subset_name} (n < 3)")
            continue

        mat_sub = mat.loc[subset_ids]
        alpha_sub = alpha.loc[subset_ids]

        # === Exclude abnormal samples（e.g. more than 70% occupied by single taxa） ===
        max_abundance = mat_sub.max(axis=1)
        mask_good = max_abundance <= 0.7

        # フィルタ適用
        mat_sub = mat_sub.loc[mask_good]
        alpha_sub = alpha_sub.loc[mask_good]

        # === detection limit（treat as 0 when less than 1e-4） ===
        mat_sub = mat_sub.where(mat_sub >= 1e-4, 0)

        # check number of samples
        if mat_sub.shape[0] < 3:
            print(f"Skip {subset_name} after filtering (n < 3)")
            continue

        H_sum = mat_sub[H_present].sum(axis=1).to_numpy()
        C_sum = mat_sub[C_present].sum(axis=1).to_numpy()

        HC = np.log10((H_sum + eps) / (C_sum + eps))

        shannon = alpha_sub["shannon"].to_numpy()
        chao = alpha_sub["Chao"].to_numpy()

        rho_s, p_s = spearmanr(HC, shannon)
        rho_c, p_c = spearmanr(HC, chao)

        # correlation figures（each cohort / each context in case AGP）
        plt.figure(figsize=(5, 4))
        plt.scatter(HC, shannon, s=12, alpha=0.7)
        plt.xlabel("HC index")
        plt.ylabel("Shannon")
        plt.title(f"{subset_name}: HC-index vs Shannon\nrho={rho_s:.3f}, p={p_s:.3g}")
        plt.tight_layout()
        plt.savefig(f"{subset_name}_HC_vs_Shannon.png", dpi=200)
        plt.close()

        plt.figure(figsize=(5, 4))
        plt.scatter(HC, chao, s=12, alpha=0.7)
        plt.xlabel("HC index")
        plt.ylabel("Chao")
        plt.title(f"{subset_name}: HC-index vs Chao\nrho={rho_c:.3f}, p={p_c:.3g}")
        plt.tight_layout()
        plt.savefig(f"{subset_name}_HC_vs_Chao.png", dpi=200)
        plt.close()

        results.append({
            "cohort": subset_name,
            "base_cohort": name,
            "context": subset_name.replace(f"{name}_", "") if subset_name != name else "all",
            "rho_shannon": rho_s,
            "p_shannon": p_s,
            "rho_chao": rho_c,
            "p_chao": p_c,
            "n": len(HC)
        })

# =========================
# Making DataFrame
# =========================
res_df = pd.DataFrame(results)

if res_df.empty:
    raise ValueError("No valid cohort/context subsets were analyzed. Check context column names and sample overlap.")

# p adjustement（BH）
res_df["p_shannon_adj"] = multipletests(res_df["p_shannon"], method="fdr_bh")[1]
res_df["p_chao_adj"] = multipletests(res_df["p_chao"], method="fdr_bh")[1]

print("\n=== Results ===")
print(res_df)

res_df.to_csv("HC_index_validation_results.csv", index=False)

# =========================
# Forest plot
# =========================
def plot_forest(df, value_col, title, outfile):
    y = np.arange(len(df))

    plt.figure(figsize=(6, 4))
    plt.errorbar(df[value_col], y, xerr=0, fmt='o')
    plt.yticks(y, df["cohort"])
    plt.axvline(0, linestyle='--')
    plt.xlabel("Spearman rho")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()

plot_forest(res_df, "rho_shannon", "HC-index vs Shannon", "forest_shannon.png")
plot_forest(res_df, "rho_chao", "HC-index vs Chao", "forest_chao.png")

print("\nForest plots and per-cohort scatter plots saved.")