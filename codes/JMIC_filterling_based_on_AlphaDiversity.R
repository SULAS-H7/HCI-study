library(data.table)

# 入力ファイル
infile  <- "L6_rel.tsv"
outfile <- "FD_table.csv"

# 1行目のコメントを飛ばして読む
df <- fread(
  infile,
  sep = "\t",
  header = TRUE,
  skip = 1,
  quote = "",
  data.table = FALSE,
  check.names = FALSE
)

# 出力
write.csv(df, file = outfile, row.names = FALSE, quote = FALSE)

cat("Saved:", outfile, "\n")
cat("Dimensions:", nrow(df), "rows x", ncol(df), "cols\n")


library(data.table)

# ----------------------------
# 入力ファイル
# ----------------------------
table_file   <- "FD_table.csv"
shannon_file <- "shannon.tsv"
chao_file    <- "chao.tsv"

# 出力
out_file <- "FD_table_filtered_by_alpha.csv"

# ----------------------------
# table 読み込み
# 1列目は taxa 名、2列目以降が SampleID
# ----------------------------
tab <- fread(
  table_file,
  data.table = FALSE,
  check.names = FALSE
)

# ----------------------------
# alpha diversity 側の SampleID を読む
# QIIME2 export の alpha-diversity.tsv は
# 1列目が sample-id のことが多い
# ----------------------------
shannon <- fread(shannon_file, data.table = FALSE, check.names = FALSE)
chao    <- fread(chao_file, data.table = FALSE, check.names = FALSE)

# 先頭列名を確認
print(colnames(shannon))
print(colnames(chao))

# 1列目を SampleID とみなす
ids_shannon <- shannon[[1]]
ids_chao    <- chao[[1]]

# 両方に残っているサンプルだけ使う
keep_ids <- intersect(ids_shannon, ids_chao)

cat("Samples in Shannon:", length(ids_shannon), "\n")
cat("Samples in Chao   :", length(ids_chao), "\n")
cat("Samples in both   :", length(keep_ids), "\n")

# ----------------------------
# table 側から残す列を選ぶ
# 1列目は taxa 列なので必ず残す
# ----------------------------
table_sample_ids <- colnames(tab)[-1]

keep_ids_in_table <- intersect(table_sample_ids, keep_ids)

cat("Samples found in table and alpha files:", length(keep_ids_in_table), "\n")

tab_filtered <- tab[, c(colnames(tab)[1], keep_ids_in_table), drop = FALSE]

# ----------------------------
# 保存
# ----------------------------
write.csv(tab_filtered, out_file, row.names = FALSE, quote = FALSE)

cat("Saved:", out_file, "\n")
cat("Dimensions:", nrow(tab_filtered), "rows x", ncol(tab_filtered), "cols\n")
