library(dplyr)
library(ggplot2)
library(readr)

# =========================
# ファイル指定
# =========================
files <- c(
  "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Original1803_HC_taxa_evaluation/O_Analysis_new17/hc_index17.csv",
  "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Kyotango/O_Analysis_new17/hc_index17.csv",
  "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/AGP_over3000/O_Analysis_new17/hc_index17.csv",
  "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Hadza/hc_index17.csv",
  "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/FD_DB/O_Analysis_new17/hc_index17.csv",
  "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/PRJNA385551_china_cohort/hc_index_new17.csv"
)

cohort_names <- c(
  "Original",
  "Kyotango",
  "AGP",
  "Hadza",
  "Flora_Scan",
  "Chinese"
)

# =========================
# 読み込み
# =========================
df_list <- lapply(seq_along(files), function(i) {
  dat <- read_csv(files[i], show_col_types = FALSE,
                  col_types = cols(.default = "c"))  # ← ここ重要
  
  dat$Cohort <- cohort_names[i]
  
  # HC_indexだけ数値に戻す
  dat$HC_index <- as.numeric(dat$HC_index)
  
  dat
})

df_all <- bind_rows(df_list)

# 念のため
df_all <- df_all %>%
  filter(!is.na(HC_index))

# =========================
# ヒストグラム
# =========================
p <- ggplot(df_all, aes(x = HC_index)) +
  geom_histogram(bins = 40, color = "black", fill = "grey80") +
  facet_wrap(~ Cohort, ncol = 2, scales = "fixed") +
  theme_bw() +
  labs(
    title = "Distribution of HC-index across cohorts",
    x = "HC-index",
    y = "Count"
  )

ggsave("HC_index_histograms_4cohorts.png", p, width = 8, height = 6, dpi = 300)

p <- ggplot(df_all, aes(x = HC_index)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40,
                 color = "black",
                 fill = "grey80") +
  facet_wrap(~ Cohort, ncol = 2, scales = "fixed") +
  theme_bw() +
  labs(
    title = "Distribution of HC-index across cohorts",
    x = "HC-index",
    y = "Density"
  )

ggsave("HC_index_histograms_4cohorts_density.png", p, width = 8, height = 6, dpi = 300)



library(dplyr)
library(ggplot2)
library(readr)

# =========================
# ファイル指定
# =========================
# plot_hc_by_region.py で出力された stool + metadata + projects 結合ファイルを指定 
input_file <- "/Users/ryoinoue/Library/CloudStorage/Dropbox/DataTransfer/2025_VariousAnimals2/H_Ctaxa_validation/CohortData/Human_Microbiome_Compendium/HC_region_results6/HC_index_with_metadata_projects_AMPLICON_stool.tsv"

# =========================
# 読み込み
# =========================
df <- read_tsv(input_file, show_col_types = FALSE, col_types = cols(.default = "c"))

# HC_indexを数値化、region欠損を除外
df <- df %>%
  mutate(
    HC_index = as.numeric(HC_index),
    region = ifelse(is.na(region) | region == "" | region == "NA", "unknown", region)
  ) %>%
  filter(!is.na(HC_index), !is.na(region))

# regionの表示順：中央値が高い順
region_order <- df %>%
  group_by(region) %>%
  summarise(
    n = n(),
    median_HC = median(HC_index, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_HC)) %>%
  pull(region)

df <- df %>%
  mutate(region = factor(region, levels = region_order))

# regionごとのサンプル数も保存
region_summary <- df %>%
  group_by(region) %>%
  summarise(
    n = n(),
    mean_HC = mean(HC_index, na.rm = TRUE),
    median_HC = median(HC_index, na.rm = TRUE),
    sd_HC = sd(HC_index, na.rm = TRUE),
    min_HC = min(HC_index, na.rm = TRUE),
    max_HC = max(HC_index, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(region_summary, "HC_index_region_histogram_summary.tsv")

# =========================
# ヒストグラム：Count
# =========================
p_count <- ggplot(df, aes(x = HC_index)) +
  geom_histogram(bins = 60, color = "black", fill = "grey80") +
  facet_wrap(~ region, ncol = 3, scales = "fixed") +
  theme_bw() +
  labs(
    title = "Distribution of HC-index by region",
    x = "HC-index",
    y = "Count"
  ) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave("HC_index_histograms_by_region_count.png", p_count, width = 10, height = 8, dpi = 300)

# =========================
# ヒストグラム：Density
# =========================
p_density <- ggplot(df, aes(x = HC_index)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60,
                 color = "black",
                 fill = "grey80") +
  facet_wrap(~ region, ncol = 3, scales = "fixed") +
  theme_bw() +
  labs(
    title = "Distribution of HC-index by region",
    x = "HC-index",
    y = "Density"
  ) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave("HC_index_histograms_by_region_density.png", p_density, width = 10, height = 8, dpi = 300)

# =========================
# 密度曲線：region重ね描き
# =========================
p_overlay <- ggplot(df, aes(x = HC_index, color = region)) +
  geom_density(linewidth = 0.8, na.rm = TRUE) +
  theme_bw() +
  labs(
    title = "Density curves of HC-index by region",
    x = "HC-index",
    y = "Density",
    color = "Region"
  )

ggsave("HC_index_density_overlay_by_region.png", p_overlay, width = 9, height = 6, dpi = 300)


df_adult <- df %>%
  filter(grepl("adult|adults|elderly|healthy adults", subjects, ignore.case = TRUE))


ggplot(df_adult, aes(x = HC_index)) +
  geom_histogram() +
  facet_wrap(~ region)
