library(dplyr)

read_qiime_tsv <- function(path) {
  df <- read.delim(path, sep = "\t", skip = 1, check.names = FALSE)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  t(as.matrix(df))
}


mat90  <- read_qiime_tsv("diet90_L6_rel.tsv")
mat125 <- read_qiime_tsv("diet125_L6_rel.tsv")

dim(mat90)
dim(mat125)

head(rownames(mat90))
head(rownames(mat125))

head(colnames(mat90))
head(colnames(mat125))

common_samples_90_125 <- intersect(rownames(mat90), rownames(mat125))

length(common_samples_90_125)
head(common_samples_90_125)

only90 <- setdiff(rownames(mat90), rownames(mat125))

mat90_only <- mat90[only90, , drop = FALSE]

dim(mat90_only)

all_taxa <- union(colnames(mat125), colnames(mat90_only))

fill_taxa <- function(m, taxa) {
  miss <- setdiff(taxa, colnames(m))
  if (length(miss) > 0) {
    add <- matrix(
      0,
      nrow = nrow(m),
      ncol = length(miss),
      dimnames = list(rownames(m), miss)
    )
    m <- cbind(m, add)
  }
  m[, taxa, drop = FALSE]
}

mat125_f <- fill_taxa(mat125, all_taxa)
mat90_f  <- fill_taxa(mat90_only, all_taxa)

dim(mat125_f)
dim(mat90_f)

mat_diet_merged <- rbind(mat125_f, mat90_f)

dim(mat_diet_merged)

meta_diet <- read.delim("meta_candidates_diet.tsv", sep = "\t", check.names = FALSE,
                        stringsAsFactors = FALSE)

nrow(meta_diet)
head(meta_diet[, 1:min(6, ncol(meta_diet))])

# 追加メタデータ: 年齢・性別を AGP 元メタデータから取得
meta_agp <- read.delim("../10317_20260402-160702.txt", sep = "\t", check.names = FALSE,
                       stringsAsFactors = FALSE)

meta_agp_small <- meta_agp %>%
  transmute(
    sample_name,
    age_use = dplyr::coalesce(na_if(age_years, ""), na_if(host_age, ""), na_if(age_cat, "")),
    sex_use = dplyr::coalesce(na_if(sex, ""), na_if(gender_v2, "")),
    age_cat,
    bmi_use = bmi_cat
  )

# meta_candidates_diet.tsv の #SampleID はすでに base ID なので、そのまま結合する
meta_diet <- meta_diet %>%
  left_join(meta_agp_small, by = c("#SampleID" = "sample_name"))

sum(is.na(meta_diet$age_use))
sum(is.na(meta_diet$sex_use))
head(meta_diet[, c("#SampleID", "age_use", "sex_use", "age_cat", "bmi_use")])

common_samples <- intersect(rownames(mat_diet_merged), meta_diet$`#SampleID`)

length(common_samples)

mat_sub  <- mat_diet_merged[common_samples, , drop = FALSE]
meta_sub <- meta_diet[match(common_samples, meta_diet$`#SampleID`), , drop = FALSE]

all(rownames(mat_sub) == meta_sub$`#SampleID`)

# === 1属のabundanceが80%以上のサンプルを除外 ===
max_abundance <- apply(mat_sub, 1, max)
keep_samples <- max_abundance <= 0.8

# フィルタ適用
mat_sub  <- mat_sub[keep_samples, , drop = FALSE]
meta_sub <- meta_sub[keep_samples, , drop = FALSE]

cat("Removed samples (>=80% single genus):", sum(!keep_samples), "\n")
cat("Remaining samples:", nrow(mat_sub), "\n")


df_final <- cbind(meta_sub, mat_sub)

dim(df_final)
head(df_final[, 1:min(10, ncol(df_final))])

# 1. alpha diversity を読む
shannon125 <- read.delim("shannon_export125/alpha-diversity.tsv", sep = "\t", check.names = FALSE,
                         stringsAsFactors = FALSE)
chao125    <- read.delim("chao1_export125/alpha-diversity.tsv", sep = "\t", check.names = FALSE,
                         stringsAsFactors = FALSE)

shannon90  <- read.delim("shannon_export90/alpha-diversity.tsv", sep = "\t", check.names = FALSE,
                         stringsAsFactors = FALSE)
chao90     <- read.delim("chao1_export90/alpha-diversity.tsv", sep = "\t", check.names = FALSE,
                         stringsAsFactors = FALSE)

colnames(shannon125) <- c("#SampleID", "Shannon")
colnames(chao125)    <- c("#SampleID", "Chao1")
colnames(shannon90)  <- c("#SampleID", "Shannon")
colnames(chao90)     <- c("#SampleID", "Chao1")

# 2. 125 と 90 をそれぞれまとめる
alpha125 <- shannon125 %>%
  left_join(chao125, by = "#SampleID") %>%
  mutate(source_alpha = "125nt")

alpha90 <- shannon90 %>%
  left_join(chao90, by = "#SampleID") %>%
  mutate(source_alpha = "90nt")

# 3. 125 を優先して統合
only90_ids <- setdiff(alpha90$`#SampleID`, alpha125$`#SampleID`)
alpha90_only <- alpha90 %>%
  filter(`#SampleID` %in% only90_ids)

alpha_merged <- bind_rows(alpha125, alpha90_only)

# 4. 確認
nrow(alpha125)
nrow(alpha90)
nrow(alpha_merged)

table(alpha_merged$source_alpha)

df_final <- df_final %>%
  left_join(alpha_merged, by = "#SampleID")

summary(df_final$Shannon)
summary(df_final$Chao1)

sum(is.na(df_final$Shannon))
sum(is.na(df_final$Chao1))

table(df_final$source_alpha, useNA = "ifany")

meta_with_alpha <- meta_diet %>%
  left_join(alpha_merged, by = "#SampleID")

sum(rownames(mat_diet_merged) %in% alpha_merged$`#SampleID`)
sum(alpha_merged$`#SampleID` %in% rownames(mat_diet_merged))
setdiff(rownames(mat_diet_merged), alpha_merged$`#SampleID`)[1:10]

length(setdiff(rownames(mat_diet_merged), alpha_merged$`#SampleID`))

dim(df_final)

summary(df_final$Shannon)
summary(df_final$Chao1)

sum(is.na(df_final$Shannon))
sum(is.na(df_final$Chao1))

table(df_final$source_alpha, useNA = "ifany")

saveRDS(df_final, "AGP_diet_merged_90_125_L6.rds")
write.table(df_final,
            file = "AGP_diet_merged_90_125_L6.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
write.csv(df_final,
          file = "AGP_diet_merged_90_125_L6.csv",
          row.names = FALSE)
