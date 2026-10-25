

library(dplyr)
library(FSA)
library(ggpubr)
library(gridExtra)
library(grid)
library(nnet)

# =========================
# 出力設定
# =========================
# 探索的解析なので、図示・抽出の基準は BH 調整後 p < 0.1 とする
# 論文化・報告書では FDR < 0.1 と明記する
FDR_CUTOFF <- 0.10
RHO_CUTOFF <- 0.20
MAX_PLOTS <- 50

out_dirs <- c("categorical", "boxplot", "correlation", "ageband", "adjusted")
invisible(lapply(out_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

safe_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("<", "lt", x)
  x <- gsub(">", "gt", x)
  x <- gsub("[^A-Za-z0-9_\u3040-\u30ff\u4e00-\u9fff.-]+", "_", x)
  x
}

# Helper to wrap long labels for plot readability
wrap_label <- function(x, width = 45) {
  x <- as.character(x)
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

make_ageband_dirs <- function(response_val, band_val, sex_val) {
  response_safe <- safe_filename(response_val)
  band_safe <- safe_filename(band_val)
  sex_safe <- safe_filename(sex_val)
  base_dir <- file.path("ageband", response_safe, band_safe, sex_safe)
  dirs <- file.path(base_dir, c("boxplot", "categorical", "correlation"))
  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))
  list(
    base = base_dir,
    boxplot = file.path(base_dir, "boxplot"),
    categorical = file.path(base_dir, "categorical"),
    correlation = file.path(base_dir, "correlation")
  )
}

# =========================
# Dunn検定結果を ggpubr 用の注釈データに変換
# =========================
make_dunn_annotations <- function(dunn_res, y_values, padj_col = "P.adj") {
  if (is.null(dunn_res) || nrow(dunn_res) == 0) return(NULL)
  if (!(padj_col %in% names(dunn_res))) return(NULL)
  
  ann <- dunn_res %>%
    dplyr::filter(!is.na(.data[[padj_col]]), .data[[padj_col]] < FDR_CUTOFF)
  
  if (nrow(ann) == 0) return(NULL)
  
  comps <- strsplit(as.character(ann$Comparison), " - ")
  ann$group1 <- vapply(comps, `[`, character(1), 1)
  ann$group2 <- vapply(comps, `[`, character(1), 2)
  
  y_max <- max(y_values, na.rm = TRUE)
  y_min <- min(y_values, na.rm = TRUE)
  y_span <- y_max - y_min
  if (!is.finite(y_span) || y_span == 0) y_span <- max(abs(y_max), 1)
  
  ann$y.position <- y_max + seq_len(nrow(ann)) * (0.08 * y_span)
  ann$p.adj <- ann[[padj_col]]
  ann$p.adj.signif <- dplyr::case_when(
    ann$p.adj < 0.001 ~ "***",
    ann$p.adj < 0.01 ~ "**",
    ann$p.adj < 0.05 ~ "*",
    ann$p.adj < FDR_CUTOFF ~ "†",
    TRUE ~ "ns"
  )
  
  ann
}

# =========================
# 標準化残差テーブルを作成（|stdres| > 2 のみ）
# =========================
make_residual_table <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  
  tab <- res_df %>%
    dplyr::filter(!is.na(stdres), abs(stdres) > 2) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        stdres > 2 ~ "over",
        stdres < -2 ~ "under",
        TRUE ~ "ns"
      ),
      stdres = round(stdres, 2)
    ) %>%
    dplyr::select(group, category, stdres, direction)
  
  if (nrow(tab) == 0) {
    tab <- data.frame(
      group = "-",
      category = "No residuals |stdres| > 2",
      stdres = "",
      direction = ""
    )
  }
  
  tab
}


df <- read.csv("HC_index_and_metadata.csv", stringsAsFactors = FALSE)

# =========================
# Dataset-wide preprocessing
# =========================
invalid_response_values <- c("not provided", "not collected", "not applicable")
char_cols_all <- names(df)[sapply(df, is.character)]

if (length(char_cols_all) > 0) {
  invalid_mask <- apply(
    df[, char_cols_all, drop = FALSE],
    1,
    function(x) any(tolower(trimws(x)) %in% invalid_response_values, na.rm = TRUE)
  )
  message(paste0("Removing rows with invalid questionnaire responses: ", sum(invalid_mask, na.rm = TRUE)))
  df <- df[!invalid_mask, , drop = FALSE]
}

age_filter_candidates <- c("Age", "age", "age_years")
age_filter_col <- age_filter_candidates[age_filter_candidates %in% names(df)][1]

if (!is.na(age_filter_col)) {
  age_num_tmp <- suppressWarnings(as.numeric(trimws(as.character(df[[age_filter_col]]))))
  df[[age_filter_col]] <- age_num_tmp
  n_before_age <- nrow(df)
  df <- df[is.na(age_num_tmp) | age_num_tmp > 20, , drop = FALSE]
  message(paste0("Age filter applied using ", age_filter_col, ": removed ", n_before_age - nrow(df), " rows (Age <= 20)."))
}

response_candidates_pre <- c("HC_index", "Chao1", "Chao", "shannon", "Shannon")
for (rv_tmp in response_candidates_pre[response_candidates_pre %in% names(df)]) {
  df[[rv_tmp]] <- suppressWarnings(as.numeric(trimws(as.character(df[[rv_tmp]]))))
}

# =========================
# BMI外れ値の除外（BMI < 15 または BMI > 40 を除外）
# 既存の BMI 列を使用
# =========================
if ("BMI" %in% names(df)) {
  df <- df %>%
    mutate(BMI = suppressWarnings(as.numeric(BMI)))
  
  message("BMI summary before filtering:")
  print(summary(df$BMI))
  
  n_before_bmi <- nrow(df)
  df <- df %>%
    filter(is.na(BMI) | (BMI >= 15 & BMI <= 40))
  n_after_bmi <- nrow(df)
  
  message(paste0("BMI filter applied: removed ", n_before_bmi - n_after_bmi,
                 " rows (BMI < 15 or BMI > 40)."))
} else {
  message("BMI column not found. Skipping BMI filter.")
}

# =========================
# Outcome-side exclusion list
# =========================
# Age/Sex are used as covariates in adjusted models, so they should not be tested as outcomes.
# Height/Weight are also excluded from main outcome-side screening because they are structural body-size variables.
outcome_exclude_vars <- c(
  "Age", "age", "age_years",
  "Sex", "SEX", "Gender", "gender", "Sex_clean",
  "Height", "Hight", "Weight"
)

numeric_cols <- names(df)[sapply(df, is.numeric)]

# 目的変数・共変量・構造的変数は outcome 側から除外する
numeric_cols <- setdiff(
  numeric_cols,
  c("HC_index", "Chao1", "Chao", "shannon", "Shannon", outcome_exclude_vars)
)

# =========================
# 目的変数の設定（HC_index に加えて Chao, Shannon にも対応）
# =========================
response_candidates <- c("HC_index", "Chao1", "Chao", "shannon", "Shannon")
response_vars <- response_candidates[response_candidates %in% names(df)]

if (length(response_vars) == 0) {
  stop("No response variables (HC_index / Chao / Shannon) found in dataset")
}

# =========================
# 各指標ごとの四分位グループ作成
# =========================
for (rv in response_vars) {
  grp_name <- paste0(rv, "_group")
  df[[grp_name]] <- ntile(df[[rv]], 4) %>%
    factor(labels = c("Q1", "Q2", "Q3", "Q4"))
}

# =========================
# 結果保存
# =========================
results <- list()

# =========================
# カテゴリ列と図示対象変数を先に定義
# Age Category別解析でも使うので、前半で定義しておく
# =========================
cat_cols <- names(df)[sapply(df, is.character)]
cat_cols <- cat_cols[!is.na(cat_cols)]
cat_cols <- setdiff(cat_cols, c("SampleID", outcome_exclude_vars))

# カテゴリ変数の図は、手入力ではなく有意差のあるものだけ自動抽出して出力する
vars <- character(0)

# 解析部分
for (rv in response_vars) {
  grp <- paste0(rv, "_group")
  
  for (col in numeric_cols) {
    if (col %in% outcome_exclude_vars) next
    tmp <- df %>%
      select(all_of(grp), all_of(col)) %>%
      filter(!is.na(.data[[grp]]), !is.na(.data[[col]]))
    
    if (nrow(tmp) < 20) next
    
    kw <- kruskal.test(tmp[[col]] ~ tmp[[grp]])
    
    dunn <- tryCatch({
      dunnTest(tmp[[col]] ~ tmp[[grp]], method = "bh")$res
    }, error = function(e) NULL)
    
    results[[paste(rv, col, sep = "__")]] <- list(
      response = rv,
      variable = col,
      kruskal_p = kw$p.value,
      dunn = dunn
    )
  }
}

# =========================
# Age Category 別解析（サンプル数>=50のみ）
# =========================
# # 想定される列名（存在するものを自動検出）
# age_candidates <- c("age_cat", "age_category", "age_group", "age")
# age_col <- age_candidates[age_candidates %in% names(df)][1]
# 
# if (!is.na(age_col)) {
#   message(paste0("Using age column: ", age_col))
# 
#   # サンプル数が50以上のカテゴリのみ採用
#   age_counts <- df %>%
#     filter(!is.na(.data[[age_col]])) %>%
#     count(.data[[age_col]], name = "n")
# 
#   valid_ages <- age_counts %>%
#     filter(n >= 50) %>%
#     pull(1)
# 
#   message("Valid age categories (n>=50):")
#   print(valid_ages)
# 
#   # 各Age Categoryごとに同じ解析を実行
#   for (age_val in valid_ages) {
#     message(paste0("\n=== Age Category: ", age_val, " ==="))
# 
#     df_age <- df %>%
#       filter(.data[[age_col]] == age_val)
# 
#     # ---- 数値データ解析 ----
#     results_age <- list()
# 
#     for (col in numeric_cols) {
#       tmp <- df_age %>%
#         select(HC_group, all_of(col)) %>%
#         filter(!is.na(HC_group), !is.na(.data[[col]]))
# 
#       if (nrow(tmp) < 20) next
# 
#       kw <- kruskal.test(tmp[[col]] ~ tmp$HC_group)
# 
#       dunn <- tryCatch({
#         dunnTest(tmp[[col]] ~ tmp$HC_group, method = "bh")$res
#       }, error = function(e) NULL)
# 
#       results_age[[col]] <- list(
#         kruskal_p = kw$p.value,
#         dunn = dunn
#       )
#     }
# 
#     if (length(results_age) > 0) {
#       summary_age <- data.frame(
#         variable = names(results_age),
#         p_value = sapply(results_age, function(x) x$kruskal_p)
#       )
# 
#       summary_age$p_adj <- p.adjust(summary_age$p_value, method = "BH")
#       summary_age <- summary_age %>% arrange(p_adj)
# 
#       out_num <- paste0("HC_vs_nutrition_", age_val, ".csv")
#       write.csv(summary_age, out_num, row.names = FALSE)
#     }
# 
#     # ---- 数値データの可視化（Age Category別Boxplot） ----
#     if (exists("summary_age") && nrow(summary_age) > 0) {
#       top_vars_age <- summary_age$variable
#       if (length(top_vars_age) > 50) top_vars_age <- top_vars_age[1:50]
#       top_vars_age <- top_vars_age[!is.na(top_vars_age)]
#       top_vars_age <- top_vars_age[top_vars_age %in% names(df_age)]
# 
#       for (v in top_vars_age) {
#         if (is.na(v) || !(v %in% names(df_age))) next
# 
#         plot_df_num <- df_age %>%
#           filter(!is.na(HC_group), !is.na(.data[[v]]))
# 
#         if (nrow(plot_df_num) < 5) next
# 
#         p_num <- ggplot(plot_df_num, aes(x = HC_group, y = .data[[v]], fill = HC_group)) +
#           geom_boxplot(outlier.shape = NA) +
#           geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
#           theme_bw() +
#           labs(title = paste0(v, " (", age_val, ")"))
# 
#         ggsave(paste0("box_", v, "_", age_val, ".png"), p_num, width = 5, height = 4)
#       }
#     }
# 
#     # ---- カテゴリデータ解析 ----
#     cat_results_age <- list()
# 
#     for (col in cat_cols) {
#       if (is.na(col) || !(col %in% names(df_age))) next
# 
#       tmp <- df_age %>%
#         select(HC_group, all_of(col)) %>%
#         filter(!is.na(HC_group), !is.na(.data[[col]]))
# 
#       tbl <- table(tmp[[col]], tmp$HC_group)
#       if (min(dim(tbl)) < 2) next
# 
#       chisq <- chisq.test(tbl)
#       cat_results_age[[col]] <- chisq$p.value
#     }
# 
#     if (length(cat_results_age) > 0) {
#       cat_df_age <- data.frame(
#         variable = names(cat_results_age),
#         p_value = unlist(cat_results_age)
#       )
# 
#       cat_df_age$p_adj <- p.adjust(cat_df_age$p_value, method = "BH")
#       cat_df_age <- cat_df_age %>% arrange(p_adj)
# 
#       out_cat <- paste0("HC_vs_categorical_", age_val, ".csv")
#       write.csv(cat_df_age, out_cat, row.names = FALSE)
#     }
# 
#     # ---- 可視化（指定vars） ----
#     for (v in vars) {
#       if (is.na(v) || !(v %in% names(df_age))) next
# 
#       plot_df <- df_age %>%
#         filter(!is.na(HC_group), !is.na(.data[[v]])) %>%
#         group_by(HC_group, .data[[v]]) %>%
#         summarise(n = n(), .groups = "drop") %>%
#         group_by(HC_group) %>%
#         mutate(prop = n / sum(n))
# 
#       p <- ggplot(plot_df, aes(x = HC_group, y = prop, fill = .data[[v]])) +
#         geom_col(position = "fill") +
#         theme_bw() +
#         labs(
#           title = paste0(v, " (", age_val, ")"),
#           x = "HC-index quartile",
#           y = "Proportion"
#         )
# 
#       ggsave(paste0("stacked_", v, "_", age_val, ".png"), p, width = 5, height = 4)
#     }
#   }
# } else {
#   message("No age category column found. Skipping age-stratified analysis.")
# }

# =========================
# 年齢3群（<40, 40-60, >60）別解析
# =========================
# 数値の年齢列がある場合のみ実行
age_numeric_candidates <- c("age_years", "age", "Age")

# =========================
# Sex列の自動検出と正規化
# =========================
sex_candidates <- c("Sex", "SEX", "Gender", "gender")
sex_col <- sex_candidates[sex_candidates %in% names(df)][1]

if (!is.na(sex_col)) {
  df[["Sex_clean"]] <- case_when(
    df[[sex_col]] %in% c("Male", "M", "1") ~ "Male",
    df[[sex_col]] %in% c("Female", "F", "2") ~ "Female",
    TRUE ~ NA_character_
  )
} else {
  df[["Sex_clean"]] <- NA_character_
}
age_num_col <- age_numeric_candidates[age_numeric_candidates %in% names(df)][1]

if (!is.na(age_num_col) && is.numeric(df[[age_num_col]])) {
  message(paste0("Using numeric age column for custom age bands: ", age_num_col))
  
  df_ageband <- df %>%
    mutate(
      age_band3 = case_when(
        .data[[age_num_col]] < 40 ~ "<40",
        .data[[age_num_col]] <= 60 ~ "40-60",
        .data[[age_num_col]] > 60 ~ ">60",
        TRUE ~ NA_character_
      )
    )
  
  ageband_counts <- df_ageband %>%
    filter(!is.na(age_band3)) %>%
    count(age_band3, name = "n")
  
  message("Custom age bands:")
  print(ageband_counts)
  
  # 年齢3群 × 性別の件数を表示
  age_sex_counts <- df_ageband %>%
    filter(!is.na(age_band3), !is.na(Sex_clean)) %>%
    count(age_band3, Sex_clean, name = "n") %>%
    arrange(age_band3, Sex_clean)
  
  message("Counts by custom age band and sex:")
  print(age_sex_counts)
  
  age_sex_table <- with(
    df_ageband %>% filter(!is.na(age_band3), !is.na(Sex_clean)),
    table(age_band3, Sex_clean)
  )
  message("Contingency table: age band x sex")
  print(age_sex_table)
  
  valid_bands <- ageband_counts %>%
    filter(n >= 50) %>%
    pull(age_band3)
  
  message("Valid custom age bands (n>=50):")
  print(valid_bands)
  
  for (band_val in valid_bands) {
    message(paste0("\n=== Custom Age Band: ", band_val, " ==="))
    
    sex_levels <- unique(df$Sex_clean)
    sex_levels <- sex_levels[!is.na(sex_levels)]
    sex_levels <- c("All", sex_levels)
    sex_levels <- unique(sex_levels)
    
    for (sex_val in sex_levels) {
      message(paste0("--- Sex: ", sex_val))
      
      if (sex_val == "All") {
        df_band <- df_ageband %>% filter(age_band3 == band_val)
      } else {
        df_band <- df_ageband %>% filter(age_band3 == band_val, Sex_clean == sex_val)
      }
      
      if (nrow(df_band) < 30) next
      
      for (rv in response_vars) {
        grp <- paste0(rv, "_group")
        age_out_dirs <- make_ageband_dirs(rv, band_val, sex_val)
        
        # ---- 数値データ解析 ----
        results_band <- list()
        
        for (col in numeric_cols) {
          if (col %in% outcome_exclude_vars) next
          tmp <- df_band %>%
            select(all_of(grp), all_of(col)) %>%
            filter(!is.na(.data[[grp]]), !is.na(.data[[col]]))
          
          if (nrow(tmp) < 20) next
          
          kw <- kruskal.test(tmp[[col]] ~ tmp[[grp]])
          
          dunn <- tryCatch({
            dunnTest(tmp[[col]] ~ tmp[[grp]], method = "bh")$res
          }, error = function(e) NULL)
          
          results_band[[col]] <- list(
            response = rv,
            variable = col,
            kruskal_p = kw$p.value,
            dunn = dunn
          )
        }
        
        if (length(results_band) > 0) {
          summary_band <- data.frame(
            response = rv,
            variable = names(results_band),
            p_value = sapply(results_band, function(x) x$kruskal_p)
          )
          
          summary_band$p_adj <- p.adjust(summary_band$p_value, method = "BH")
          summary_band$significant_FDR <- summary_band$p_adj < FDR_CUTOFF
          summary_band <- summary_band %>% arrange(p_adj)
          
          write.csv(summary_band, file.path(age_out_dirs$boxplot, "numeric_quartile_summary.csv"), row.names = FALSE)
        }
        
        # ---- 数値データの可視化（年齢3群別Boxplot） ----
        if (exists("summary_band") && nrow(summary_band) > 0) {
          top_vars_band <- summary_band$variable[summary_band$p_adj < FDR_CUTOFF]
          if (length(top_vars_band) > 50) top_vars_band <- top_vars_band[1:50]
          top_vars_band <- top_vars_band[!is.na(top_vars_band)]
          top_vars_band <- top_vars_band[top_vars_band %in% names(df_band)]
          
          for (v in top_vars_band) {
            if (is.na(v) || !(v %in% names(df_band))) next
            
            plot_df_num <- df_band %>%
              filter(!is.na(.data[[grp]]), !is.na(.data[[v]]))
            
            if (nrow(plot_df_num) < 5) next
            
            dunn_ann_band <- make_dunn_annotations(results_band[[v]]$dunn, plot_df_num[[v]])
            
            p_num <- ggplot(plot_df_num, aes(x = .data[[grp]], y = .data[[v]], fill = .data[[grp]])) +
              geom_boxplot(outlier.shape = NA) +
              geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
              stat_compare_means(method = "kruskal.test") +
              theme_bw() +
              labs(
                title = paste0(v, " vs ", rv, " (", band_val, ", ", sex_val, ")"),
                subtitle = paste0("Kruskal-Wallis BH adj p=", signif(summary_band$p_adj[summary_band$variable == v], 3), "; plotted if FDR < ", FDR_CUTOFF),
                x = paste0(rv, " quartile"),
                y = v
              )
            
            if (!is.null(dunn_ann_band)) {
              p_num <- p_num +
                stat_pvalue_manual(
                  dunn_ann_band,
                  label = "p.adj.signif",
                  xmin = "group1",
                  xmax = "group2",
                  y.position = "y.position",
                  tip.length = 0.01,
                  inherit.aes = FALSE
                )
            }
            
            ggsave(file.path(age_out_dirs$boxplot, paste0("box_", safe_filename(v), ".png")), p_num, width = 5, height = 4)
          }
        }
        
        # ---- カテゴリデータ解析 ----
        cat_results_band <- list()
        
        for (col in cat_cols) {
          if (is.na(col) || col %in% outcome_exclude_vars || !(col %in% names(df_band))) next
          
          tmp <- df_band %>%
            select(all_of(grp), all_of(col)) %>%
            filter(!is.na(.data[[grp]]), !is.na(.data[[col]]))
          
          tbl <- table(tmp[[col]], tmp[[grp]])
          if (min(dim(tbl)) < 2) next
          
          chisq <- suppressWarnings(chisq.test(tbl))
          cat_results_band[[col]] <- chisq$p.value
        }
        
        sig_vars_band_cat <- character(0)
        if (length(cat_results_band) > 0) {
          cat_df_band <- data.frame(
            response = rv,
            variable = names(cat_results_band),
            p_value = unlist(cat_results_band)
          )
          
          cat_df_band$p_adj <- p.adjust(cat_df_band$p_value, method = "BH")
          cat_df_band$significant_FDR <- cat_df_band$p_adj < FDR_CUTOFF
          cat_df_band <- cat_df_band %>% arrange(p_adj)
          
          write.csv(cat_df_band, file.path(age_out_dirs$categorical, "categorical_summary.csv"), row.names = FALSE)
          sig_vars_band_cat <- cat_df_band$variable[cat_df_band$p_adj < FDR_CUTOFF]
        }
        
        # ---- カテゴリデータの可視化 ----
        for (v in sig_vars_band_cat) {
          if (is.na(v) || !(v %in% names(df_band))) next
          
          plot_df <- df_band %>%
            filter(!is.na(.data[[grp]]), !is.na(.data[[v]])) %>%
            group_by(.data[[grp]], .data[[v]]) %>%
            summarise(n = n(), .groups = "drop") %>%
            group_by(.data[[grp]]) %>%
            mutate(prop = n / sum(n))
          
          p_age_cat <- if (exists("cat_df_band") && v %in% cat_df_band$variable) {
            cat_df_band$p_adj[match(v, cat_df_band$variable)]
          } else {
            NA_real_
          }
          p <- ggplot(plot_df, aes(x = .data[[grp]], y = prop, fill = .data[[v]])) +
            geom_col(position = "fill") +
            theme_bw() +
            labs(
              title = paste0(v, " vs ", rv, " (", band_val, ", ", sex_val, ")"),
              subtitle = paste0("Chi-square BH adj p=", signif(p_age_cat, 3), "; plotted if FDR < ", FDR_CUTOFF),
              x = paste0(rv, " quartile"),
              y = "Proportion"
            )
          
          res_df <- tryCatch({
            tmp <- df_band %>%
              filter(!is.na(.data[[grp]]), !is.na(.data[[v]]))
            tbl <- table(tmp[[v]], tmp[[grp]])
            if (min(dim(tbl)) < 2) return(NULL)
            chisq <- suppressWarnings(chisq.test(tbl))
            res_df <- as.data.frame(as.table(tbl))
            names(res_df) <- c("category", "group", "n")
            res_df$stdres <- as.vector(chisq$stdres)
            res_df
          }, error = function(e) NULL)
          res_tab <- make_residual_table(res_df)
          tab_grob <- tableGrob(res_tab, rows = NULL)
          combo_grob <- arrangeGrob(p, tab_grob, ncol = 2, widths = c(2.4, 1.6))
          
          ggsave(
            file.path(age_out_dirs$categorical, paste0("stacked_", safe_filename(v), ".png")),
            combo_grob,
            width = 15,
            height = 4.5
          )
        }
        
        # ---- Spearman相関解析（年齢3群別） ----
        cor_results_band <- list()
        
        for (col in numeric_cols) {
          if (col %in% outcome_exclude_vars) next
          tmp <- df_band %>%
            select(all_of(rv), all_of(col)) %>%
            filter(!is.na(.data[[rv]]), !is.na(.data[[col]]))
          
          if (nrow(tmp) < 20) next
          
          cor_test <- suppressWarnings(cor.test(tmp[[rv]], tmp[[col]], method = "spearman", exact = FALSE))
          
          cor_results_band[[col]] <- list(
            response = rv,
            variable = col,
            rho = unname(cor_test$estimate),
            p_value = cor_test$p.value
          )
        }
        
        if (length(cor_results_band) > 0) {
          cor_df_band <- data.frame(
            response = rv,
            variable = names(cor_results_band),
            rho = vapply(cor_results_band, function(x) as.numeric(x$rho), numeric(1)),
            p_value = vapply(cor_results_band, function(x) as.numeric(x$p_value), numeric(1)),
            stringsAsFactors = FALSE
          )
          
          cor_df_band <- cor_df_band %>%
            filter(!is.na(rho), !is.na(p_value), is.finite(rho), is.finite(p_value))
          
          if (nrow(cor_df_band) == 0) next
          
          cor_df_band$p_adj <- p.adjust(cor_df_band$p_value, method = "BH")
          cor_df_band$abs_rho <- abs(cor_df_band$rho)
          cor_df_band$significant_FDR <- cor_df_band$p_adj < FDR_CUTOFF
          cor_df_band$large_enough_effect <- cor_df_band$abs_rho >= RHO_CUTOFF
          cor_df_band <- cor_df_band %>% arrange(desc(abs_rho), p_adj)
          
          write.csv(cor_df_band, file.path(age_out_dirs$correlation, "numeric_correlation.csv"), row.names = FALSE)
          
          sig_vars_band <- cor_df_band$variable[cor_df_band$p_adj < FDR_CUTOFF & cor_df_band$abs_rho >= RHO_CUTOFF]
          
          for (v in sig_vars_band) {
            if (!(v %in% names(df_band))) next
            
            idx <- which(cor_df_band$variable == v)[1]
            plot_df_cor <- df_band %>%
              filter(!is.na(.data[[rv]]), !is.na(.data[[v]]))
            
            if (nrow(plot_df_cor) < 20) next
            
            p_cor <- ggplot(plot_df_cor, aes(x = .data[[rv]], y = .data[[v]])) +
              geom_point(alpha = 0.4) +
              geom_smooth(method = "lm") +
              theme_bw() +
              labs(
                title = paste0(v, " vs ", rv, " (", band_val, ", ", sex_val, ")"),
                subtitle = paste0(
                  "Spearman rho=", round(cor_df_band$rho[idx], 3),
                  ", |rho|=", round(cor_df_band$abs_rho[idx], 3),
                  ", BH adj p=", signif(cor_df_band$p_adj[idx], 3),
                  "; plotted if FDR < ", FDR_CUTOFF,
                  " and |rho| >= ", RHO_CUTOFF
                ),
                x = rv,
                y = v
              )
            
            ggsave(
              file.path(age_out_dirs$correlation, paste0("correlation_", safe_filename(v), ".png")),
              p_cor,
              width = 5,
              height = 4
            )
          }
        }
      }
    }  # end sex loop
  }  # end sex loop
} else {
  message("No numeric age column found for custom age-band analysis. Skipping custom age-band analysis.")
}

# =========================
# p値まとめ
# =========================
summary_df <- data.frame(
  response = sapply(results, function(x) x$response),
  variable = sapply(results, function(x) x$variable),
  p_value = sapply(results, function(x) x$kruskal_p)
)

summary_df$p_adj <- p.adjust(summary_df$p_value, method = "BH")
summary_df$significant_FDR <- summary_df$p_adj < FDR_CUTOFF

summary_df <- summary_df %>%
  arrange(p_adj)

write.csv(summary_df, file.path("boxplot", "HC_vs_nutrition_summary.csv"), row.names = FALSE)

head(summary_df, 20)


library(ggplot2)

top_vars <- summary_df$variable[summary_df$p_adj < FDR_CUTOFF]
if (length(top_vars) > 50) top_vars <- top_vars[1:50]
top_vars <- top_vars[!is.na(top_vars)]
top_vars <- top_vars[top_vars %in% names(df)]

plot_idx <- which(summary_df$p_adj < FDR_CUTOFF)
if (length(plot_idx) > MAX_PLOTS) plot_idx <- plot_idx[seq_len(MAX_PLOTS)]

for (i in plot_idx) {
  v <- summary_df$variable[i]
  rv <- summary_df$response[i]
  grp <- paste0(rv, "_group")
  
  if (!(v %in% names(df))) next
  
  plot_df <- df %>% filter(!is.na(.data[[grp]]), !is.na(.data[[v]]))
  
  dunn_ann <- make_dunn_annotations(results[[paste(rv, v, sep = "__")]]$dunn, plot_df[[v]])
  
  p <- ggplot(plot_df, aes(x = .data[[grp]], y = .data[[v]], fill = .data[[grp]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    stat_compare_means(method = "kruskal.test") +
    theme_bw() +
    labs(
      title = paste0(v, " vs ", rv),
      subtitle = paste0("Kruskal-Wallis BH adj p=", signif(summary_df$p_adj[i], 3), "; plotted if FDR < ", FDR_CUTOFF),
      x = paste0(rv, " quartile"),
      y = v
    )
  
  if (!is.null(dunn_ann)) {
    p <- p +
      stat_pvalue_manual(
        dunn_ann,
        label = "p.adj.signif",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        inherit.aes = FALSE
      )
  }
  
  out_dir <- file.path("boxplot", safe_filename(rv))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, paste0("box_", safe_filename(v), ".png")), p, width = 5, height = 4)
}


# cat_cols は上で定義済み

cat_results <- list()

for (rv in response_vars) {
  grp <- paste0(rv, "_group")
  
  for (col in cat_cols) {
    if (is.na(col) || col %in% outcome_exclude_vars || !(col %in% names(df))) next
    
    tmp <- df %>%
      select(all_of(grp), all_of(col)) %>%
      filter(!is.na(.data[[grp]]), !is.na(.data[[col]]))
    
    tbl <- table(tmp[[col]], tmp[[grp]])
    
    if (min(dim(tbl)) < 2) next
    
    chisq <- suppressWarnings(chisq.test(tbl))
    
    cat_results[[paste(rv, col, sep = "__")]] <- chisq$p.value
  }
}

cat_df <- if (length(cat_results) > 0) {
  data.frame(
    variable = sapply(strsplit(names(cat_results), "__"), `[`, 2),
    response = sapply(strsplit(names(cat_results), "__"), `[`, 1),
    p_value = unlist(cat_results)
  )
} else {
  data.frame(variable = character(0), response = character(0), p_value = numeric(0))
}

if (nrow(cat_df) > 0) {
  cat_df$p_adj <- p.adjust(cat_df$p_value, method = "BH")
  cat_df <- cat_df %>% arrange(p_adj)
  cat_df$significant_FDR <- cat_df$p_adj < FDR_CUTOFF
  sig_cat_vars <- cat_df$variable[cat_df$p_adj < FDR_CUTOFF]
} else {
  cat_df$p_adj <- numeric(0)
  sig_cat_vars <- character(0)
}

write.csv(cat_df, file.path("categorical", "HC_vs_categorical_summary.csv"), row.names = FALSE)

head(cat_df, 20)

# vars は上で定義済み
cat_plot_idx <- which(cat_df$p_adj < FDR_CUTOFF)
if (length(cat_plot_idx) > MAX_PLOTS) cat_plot_idx <- cat_plot_idx[seq_len(MAX_PLOTS)]

for (i in cat_plot_idx) {
  v <- cat_df$variable[i]
  rv <- cat_df$response[i]
  grp <- paste0(rv, "_group")
  
  if (is.na(v) || !(v %in% names(df))) next
  
  plot_df <- df %>%
    filter(!is.na(.data[[grp]]), !is.na(.data[[v]])) %>%
    group_by(.data[[grp]], .data[[v]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[grp]]) %>%
    mutate(prop = n / sum(n))
  
  p_adj_cat <- cat_df$p_adj[i]
  p <- ggplot(plot_df,
              aes(x = .data[[grp]], y = prop, fill = .data[[v]])) +
    geom_col(position = "fill") +
    theme_bw() +
    labs(
      title = paste0(v, " vs ", rv),
      subtitle = paste0("Chi-square BH adj p=", signif(p_adj_cat, 3), "; plotted if FDR < ", FDR_CUTOFF),
      x = paste0(rv, " quartile"),
      y = "Proportion"
    )
  
  res_df <- tryCatch({
    tmp <- df %>%
      filter(!is.na(.data[[grp]]), !is.na(.data[[v]]))
    tbl <- table(tmp[[v]], tmp[[grp]])
    if (min(dim(tbl)) < 2) return(NULL)
    chisq <- suppressWarnings(chisq.test(tbl))
    res_df <- as.data.frame(as.table(tbl))
    names(res_df) <- c("category", "group", "n")
    res_df$stdres <- as.vector(chisq$stdres)
    res_df
  }, error = function(e) NULL)
  res_tab <- make_residual_table(res_df)
  tab_grob <- tableGrob(res_tab, rows = NULL)
  combo_grob <- arrangeGrob(p, tab_grob, ncol = 2, widths = c(2.4, 1.6))
  
  out_dir <- file.path("categorical", safe_filename(rv))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, paste0("stacked_", safe_filename(v), ".png")), combo_grob, width = 15, height = 4.5)
}

# =========================
# Spearman相関解析（HC-index等 vs 数値データ）
# =========================
cor_results <- list()

for (rv in response_vars) {
  for (col in numeric_cols) {
    if (col %in% outcome_exclude_vars) next
    tmp <- df %>%
      select(all_of(rv), all_of(col))
    tmp <- tmp %>% filter(!is.na(.data[[rv]]), !is.na(.data[[col]]))
    if (nrow(tmp) < 20) next
    cor_test <- suppressWarnings(cor.test(tmp[[rv]], tmp[[col]], method = "spearman", exact = FALSE))
    cor_results[[paste(rv, col, sep="__")]] <- list(
      response = rv,
      variable = col,
      rho = unname(cor_test$estimate),
      p_value = cor_test$p.value
    )
  }
}

cor_df <- if (length(cor_results) > 0) {
  data.frame(
    response = vapply(cor_results, function(x) as.character(x$response), character(1)),
    variable = vapply(cor_results, function(x) as.character(x$variable), character(1)),
    rho = vapply(cor_results, function(x) as.numeric(x$rho), numeric(1)),
    p_value = vapply(cor_results, function(x) as.numeric(x$p_value), numeric(1)),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(response = character(0), variable = character(0), rho = numeric(0), p_value = numeric(0))
}

cor_df <- cor_df %>%
  filter(!is.na(rho), !is.na(p_value), is.finite(rho), is.finite(p_value))

if (nrow(cor_df) == 0) {
  message("No valid overall Spearman correlation results.")
} else {
  cor_df$p_adj <- p.adjust(cor_df$p_value, method = "BH")
  cor_df$abs_rho <- abs(cor_df$rho)
  cor_df$significant_FDR <- cor_df$p_adj < FDR_CUTOFF
  cor_df$large_enough_effect <- cor_df$abs_rho >= RHO_CUTOFF
  cor_df <- cor_df %>% arrange(desc(abs_rho), p_adj)
  
  write.csv(cor_df, file.path("correlation", "HC_vs_numeric_correlation.csv"), row.names = FALSE)
  
  sig_vars <- cor_df$variable[cor_df$p_adj < FDR_CUTOFF & cor_df$abs_rho >= RHO_CUTOFF]
  sig_responses <- cor_df$response[cor_df$p_adj < FDR_CUTOFF & cor_df$abs_rho >= RHO_CUTOFF]
  
  for (i in which(cor_df$p_adj < FDR_CUTOFF & cor_df$abs_rho >= RHO_CUTOFF)) {
    v <- cor_df$variable[i]
    rv <- cor_df$response[i]
    if (!(v %in% names(df))) next
    
    plot_df <- df %>%
      filter(!is.na(.data[[rv]]), !is.na(.data[[v]]))
    
    p <- ggplot(plot_df, aes(x = .data[[rv]], y = .data[[v]])) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "lm") +
      theme_bw() +
      labs(
        title = paste0(v, " vs ", rv),
        subtitle = paste0(
          "Spearman rho=", round(cor_df$rho[i], 3),
          ", |rho|=", round(cor_df$abs_rho[i], 3),
          ", BH adj p=", signif(cor_df$p_adj[i], 3),
          "; plotted if FDR < ", FDR_CUTOFF,
          " and |rho| >= ", RHO_CUTOFF
        )
      )
    
    out_dir <- file.path("correlation", safe_filename(rv))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, paste0("correlation_", safe_filename(v), ".png")), p, width = 5, height = 4)
  }
}


# =========================
# Age / Sex adjusted analysis
# =========================
message("\n=== Running adjusted models (Age + Sex) ===")

# Ensure Age numeric
if ("Age" %in% names(df)) {
  df$Age <- suppressWarnings(as.numeric(df$Age))
}

# Keep only rows with covariates
adj_df <- df %>%
  filter(!is.na(Age), !is.na(Sex_clean))

adj_num_results <- list()
adj_num_plot_data <- list()

for (rv in response_vars) {
  if (!(rv %in% names(adj_df))) next
  
  for (col in numeric_cols) {
    if (col %in% outcome_exclude_vars) next
    if (!(col %in% names(adj_df))) next
    
    tmp <- adj_df %>%
      select(all_of(rv), Age, Sex_clean, all_of(col)) %>%
      filter(!is.na(.data[[col]]), !is.na(.data[[rv]]))
    
    if (nrow(tmp) < 50) next
    if (length(unique(tmp[[col]])) < 5) next
    
    # --- Use clean variable names for modeling ---
    model_df <- tmp
    model_df$outcome_value <- model_df[[col]]
    model_df$response_value <- model_df[[rv]]
    
    fit <- tryCatch(
      lm(outcome_value ~ response_value + Age + Sex_clean, data = model_df),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    summ <- summary(fit)
    if (!is.null(summ$r.squared) && is.finite(summ$r.squared) && summ$r.squared > 0.99) next
    coef_tab <- summ$coefficients
    
    if (!("response_value" %in% rownames(coef_tab))) next
    
    beta_hc <- coef_tab["response_value", "Estimate"]
    se_hc <- coef_tab["response_value", "Std. Error"]
    
    result_key <- paste(rv, col, sep = "__")
    
    adj_num_results[[result_key]] <- data.frame(
      response = rv,
      variable = col,
      beta = beta_hc,
      se = se_hc,
      ci_low = beta_hc - 1.96 * se_hc,
      ci_high = beta_hc + 1.96 * se_hc,
      p_value = coef_tab["response_value", "Pr(>|t|)"]
    )
    
    adj_num_plot_data[[result_key]] <- list(
      response = rv,
      variable = col,
      tmp = model_df,
      coef_tab = coef_tab
    )
  }
}

if (length(adj_num_results) > 0) {
  adj_num_df <- bind_rows(adj_num_results)
  adj_num_df <- adj_num_df %>%
    group_by(response) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(response, p_adj)
  
  # ---- Plot only significant numeric results ----
  sig_keys_num <- paste(adj_num_df$response[adj_num_df$p_adj < FDR_CUTOFF], adj_num_df$variable[adj_num_df$p_adj < FDR_CUTOFF], sep = "__")
  
  for (key in sig_keys_num) {
    if (!(key %in% names(adj_num_plot_data))) next
    
    rv <- adj_num_plot_data[[key]]$response
    v <- adj_num_plot_data[[key]]$variable
    tmp <- adj_num_plot_data[[key]]$tmp
    coef_tab <- adj_num_plot_data[[key]]$coef_tab
    
    resid_fit <- tryCatch(
      lm(outcome_value ~ Age + Sex_clean, data = tmp),
      error = function(e) NULL
    )
    
    if (is.null(resid_fit)) next
    
    tmp$adjusted_residual <- resid(resid_fit)
    
    p_adj_num <- ggplot(tmp, aes(x = response_value, y = adjusted_residual)) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "lm") +
      theme_bw() +
      labs(
        title = wrap_label(paste0(v, " vs ", rv, " adjusted for Age/Sex"), 60),
        subtitle = paste0(
          "Adjusted beta=", signif(coef_tab["response_value", "Estimate"], 3),
          ", p=", signif(coef_tab["response_value", "Pr(>|t|)"], 3)
        ),
        x = rv,
        y = paste0(v, " residual after Age/Sex adjustment")
      )
    
    ggsave(
      file.path("adjusted", paste0("numeric_adjusted_", safe_filename(rv), "_", safe_filename(v), ".png")),
      p_adj_num,
      width = 5,
      height = 4
    )
  }
  
  write.csv(adj_num_df,
            file.path("adjusted", "numeric_adjusted_results.csv"),
            row.names = FALSE)
  
  # ---- Numeric adjusted forest plot ----
  forest_num_df <- adj_num_df %>%
    filter(p_adj < FDR_CUTOFF) %>%
    mutate(label = paste0(response, ": ", variable)) %>%
    arrange(response, beta) %>%
    mutate(label = factor(label, levels = label))
  
  if (nrow(forest_num_df) > 0) {
    p_forest_num <- ggplot(forest_num_df, aes(x = beta, y = label)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
      geom_point(size = 2) +
      theme_bw() +
      labs(
        title = "Age/Sex-adjusted numeric associations for HC-index / alpha diversity",
        subtitle = paste0("Linear model beta for response variable; plotted if BH adj p < ", FDR_CUTOFF),
        x = "Adjusted beta per 1-unit increase in response variable",
        y = NULL
      )
    
    ggsave(
      file.path("adjusted", "forest_numeric_adjusted_beta.png"),
      p_forest_num,
      width = 7,
      height = max(4, 0.35 * nrow(forest_num_df) + 1)
    )
  }
}

adj_cat_results <- list()
adj_cat_plot_data <- list()
adj_cat_plot_diagnostics <- list()

for (rv in response_vars) {
  if (!(rv %in% names(adj_df))) next
  
  for (col in cat_cols) {
    if (col %in% outcome_exclude_vars) next
    if (!(col %in% names(adj_df))) next
    
    tmp <- adj_df %>%
      select(all_of(rv), Age, Sex_clean, all_of(col)) %>%
      filter(!is.na(.data[[col]]), !is.na(.data[[rv]]))
    
    # Convert to clean categorical values
    tmp[[col]] <- trimws(as.character(tmp[[col]]))
    tmp[[col]][tmp[[col]] %in% c("", "NA", "N/A", "na", "n/a", "Unknown", "unknown", "None", "none", "Not available", "not available", "Not provided", "not provided")] <- NA
    tmp <- tmp %>% filter(!is.na(.data[[col]]))
    
    if (nrow(tmp) < 50) next
    
    # If a variable contains Yes/No plus other sparse or unknown-like categories,
    # keep only Yes/No so it can be modeled as a binary outcome.
    vals_lower <- tolower(unique(tmp[[col]]))
    if (all(c("yes", "no") %in% vals_lower) && length(vals_lower) > 2) {
      tmp <- tmp %>% filter(tolower(.data[[col]]) %in% c("yes", "no"))
    }
    
    tmp[[col]] <- as.factor(tmp[[col]])
    tmp[[col]] <- droplevels(tmp[[col]])
    
    if (nlevels(tmp[[col]]) < 2) next
    
    # Binary → logistic
    if (nlevels(tmp[[col]]) == 2) {
      model_df <- tmp
      model_df$outcome_bin <- droplevels(as.factor(model_df[[col]]))
      model_df$response_value <- model_df[[rv]]
      
      fit <- tryCatch(
        glm(outcome_bin ~ response_value + Age + Sex_clean,
            data = model_df,
            family = binomial),
        error = function(e) NULL
      )
      
      if (is.null(fit) || !isTRUE(fit$converged)) next
      
      coef_tab <- summary(fit)$coefficients
      
      if (!("response_value" %in% rownames(coef_tab))) next
      
      beta_hc <- coef_tab["response_value", "Estimate"]
      se_hc <- coef_tab["response_value", "Std. Error"]
      
      # Likelihood ratio test for binary model: full model vs covariates-only model
      fit_reduced_bin <- tryCatch(
        glm(outcome_bin ~ Age + Sex_clean,
            data = model_df,
            family = binomial),
        error = function(e) NULL
      )
      
      lrt_stat_bin <- NA_real_
      lrt_df_bin <- NA_real_
      lrt_p_bin <- NA_real_
      
      if (!is.null(fit_reduced_bin)) {
        ll_full_bin <- as.numeric(logLik(fit))
        ll_reduced_bin <- as.numeric(logLik(fit_reduced_bin))
        df_full_bin <- attr(logLik(fit), "df")
        df_reduced_bin <- attr(logLik(fit_reduced_bin), "df")
        lrt_stat_bin <- 2 * (ll_full_bin - ll_reduced_bin)
        lrt_df_bin <- df_full_bin - df_reduced_bin
        lrt_p_bin <- if (is.finite(lrt_stat_bin) && lrt_df_bin > 0) {
          pchisq(lrt_stat_bin, df = lrt_df_bin, lower.tail = FALSE)
        } else {
          NA_real_
        }
      }
      
      adj_cat_results[[paste(rv, col, sep = "__")]] <- data.frame(
        response = rv,
        variable = col,
        type = "binary",
        beta = beta_hc,
        se = se_hc,
        ci_low = beta_hc - 1.96 * se_hc,
        ci_high = beta_hc + 1.96 * se_hc,
        odds_ratio = exp(beta_hc),
        or_ci_low = exp(beta_hc - 1.96 * se_hc),
        or_ci_high = exp(beta_hc + 1.96 * se_hc),
        lrt_stat = lrt_stat_bin,
        lrt_df = lrt_df_bin,
        lrt_p_value = lrt_p_bin,
        p_value = coef_tab["response_value", "Pr(>|z|)"]
      )
      
      result_key <- paste(rv, col, sep = "__")
      
      adj_cat_plot_data[[result_key]] <- list(
        response = rv,
        variable = col,
        tmp = model_df,
        fit = fit,
        coef_tab = coef_tab
      )
      
      adj_cat_plot_diagnostics[[result_key]] <- data.frame(
        response = rv,
        variable = col,
        type = "binary",
        n = nrow(model_df),
        n_levels = nlevels(model_df$outcome_bin),
        model_converged = isTRUE(fit$converged),
        stored_for_plot = TRUE,
        reason = "binary model stored; plot depends on BH-adjusted p-value"
      )
      
    } else {
      # Multi-category → multinomial logistic regression
      # This is preferable to numeric encoding because nominal categories do not have a natural order.
      model_df <- tmp
      model_df$outcome_multinom <- droplevels(as.factor(model_df[[col]]))
      model_df$response_value <- model_df[[rv]]
      
      fit_full <- tryCatch(
        nnet::multinom(outcome_multinom ~ response_value + Age + Sex_clean,
                       data = model_df,
                       trace = FALSE),
        error = function(e) NULL
      )
      
      fit_reduced <- tryCatch(
        nnet::multinom(outcome_multinom ~ Age + Sex_clean,
                       data = model_df,
                       trace = FALSE),
        error = function(e) NULL
      )
      
      if (is.null(fit_full) || is.null(fit_reduced)) next
      
      ll_full <- as.numeric(logLik(fit_full))
      ll_reduced <- as.numeric(logLik(fit_reduced))
      df_full <- attr(logLik(fit_full), "df")
      df_reduced <- attr(logLik(fit_reduced), "df")
      lrt_stat <- 2 * (ll_full - ll_reduced)
      lrt_df <- df_full - df_reduced
      lrt_p <- if (is.finite(lrt_stat) && lrt_df > 0) pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE) else NA_real_
      
      if (is.na(lrt_p)) next
      
      result_key <- paste(rv, col, sep = "__")
      
      adj_cat_results[[result_key]] <- data.frame(
        response = rv,
        variable = col,
        type = "multinomial",
        beta = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        odds_ratio = NA_real_,
        or_ci_low = NA_real_,
        or_ci_high = NA_real_,
        lrt_stat = lrt_stat,
        lrt_df = lrt_df,
        lrt_p_value = lrt_p,
        p_value = lrt_p
      )
      
      adj_cat_plot_data[[result_key]] <- list(
        response = rv,
        variable = col,
        tmp = model_df,
        fit = fit_full,
        coef_tab = NULL
      )
      
      adj_cat_plot_diagnostics[[result_key]] <- data.frame(
        response = rv,
        variable = col,
        type = "multinomial",
        n = nrow(model_df),
        n_levels = nlevels(model_df$outcome_multinom),
        model_converged = NA,
        stored_for_plot = TRUE,
        reason = "multi-category variable modeled by multinomial logistic regression; plot depends on BH-adjusted p-value"
      )
    }
  }
}

if (length(adj_cat_results) > 0) {
  adj_cat_df <- bind_rows(adj_cat_results)
  adj_cat_df <- adj_cat_df %>%
    group_by(response) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(response, p_adj)
  
  write.csv(adj_cat_df,
            file.path("adjusted", "categorical_adjusted_results.csv"),
            row.names = FALSE)
  
  # =========================
  # LRT-based comparison among HC-index, Chao, and Shannon for categorical outcomes
  # =========================
  cat_lrt_df <- adj_cat_df %>%
    filter(!is.na(lrt_stat), is.finite(lrt_stat)) %>%
    select(response, variable, type, lrt_stat, lrt_df, lrt_p_value, p_adj)
  
  if (nrow(cat_lrt_df) > 0) {
    write.csv(
      cat_lrt_df,
      file.path("adjusted", "categorical_LRT_statistics.csv"),
      row.names = FALSE
    )
    
    # Long-format distribution plot of LRT statistics
    p_lrt_box <- ggplot(cat_lrt_df, aes(x = response, y = lrt_stat)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
      theme_bw() +
      labs(
        title = "LRT statistics for categorical outcomes",
        subtitle = "Higher LRT indicates stronger explanatory contribution beyond Age/Sex",
        x = "Response variable",
        y = "Likelihood ratio statistic"
      )
    
    ggsave(
      file.path("adjusted", "categorical_LRT_distribution_boxplot.png"),
      p_lrt_box,
      width = 7,
      height = 4.5
    )
    
    p_lrt_density <- ggplot(cat_lrt_df, aes(x = lrt_stat, color = response)) +
      geom_density(linewidth = 1, na.rm = TRUE) +
      theme_bw() +
      labs(
        title = "Distribution of LRT statistics for categorical outcomes",
        subtitle = "Higher values indicate stronger association after Age/Sex adjustment",
        x = "Likelihood ratio statistic",
        y = "Density",
        color = "Response"
      )
    
    ggsave(
      file.path("adjusted", "categorical_LRT_distribution_density.png"),
      p_lrt_density,
      width = 7,
      height = 4.5
    )
    
    # Summary table by response
    cat_lrt_summary <- cat_lrt_df %>%
      group_by(response) %>%
      summarise(
        n_tests = n(),
        mean_LRT = mean(lrt_stat, na.rm = TRUE),
        median_LRT = median(lrt_stat, na.rm = TRUE),
        max_LRT = max(lrt_stat, na.rm = TRUE),
        n_FDR_sig = sum(!is.na(p_adj) & p_adj < FDR_CUTOFF),
        .groups = "drop"
      ) %>%
      arrange(desc(median_LRT))
    
    write.csv(
      cat_lrt_summary,
      file.path("adjusted", "categorical_LRT_summary_by_response.csv"),
      row.names = FALSE
    )
  }
  
  if (length(adj_cat_plot_diagnostics) > 0) {
    adj_cat_diag_df <- bind_rows(adj_cat_plot_diagnostics)
    adj_cat_diag_df <- adj_cat_diag_df %>%
      left_join(adj_cat_df %>% select(response, variable, type, p_value, p_adj),
                by = c("response", "variable", "type")) %>%
      mutate(
        passes_FDR = !is.na(p_adj) & p_adj < FDR_CUTOFF,
        will_plot = stored_for_plot & passes_FDR
      )
    
    write.csv(
      adj_cat_diag_df,
      file.path("adjusted", "categorical_adjusted_plot_diagnostics.csv"),
      row.names = FALSE
    )
  }
  
  # ---- Plot only significant binary categorical results ----
  sig_keys_cat <- adj_cat_df %>%
    filter(type == "binary", p_adj < FDR_CUTOFF) %>%
    mutate(key = paste(response, variable, sep = "__")) %>%
    pull(key)
  
  message(paste0("Adjusted categorical plots to generate: ", length(sig_keys_cat)))
  
  for (key in sig_keys_cat) {
    if (!(key %in% names(adj_cat_plot_data))) next
    
    rv <- adj_cat_plot_data[[key]]$response
    v <- adj_cat_plot_data[[key]]$variable
    tmp <- adj_cat_plot_data[[key]]$tmp
    fit <- adj_cat_plot_data[[key]]$fit
    coef_tab <- adj_cat_plot_data[[key]]$coef_tab
    
    age_median <- 50
    sex_levels_for_plot <- unique(tmp$Sex_clean)
    sex_levels_for_plot <- sex_levels_for_plot[!is.na(sex_levels_for_plot)]
    
    # Marginalized prediction: Age is fixed at the median, Sex is averaged over.
    # This treats Sex_clean as a confounder rather than displaying sex-specific curves.
    pred_base <- data.frame(
      response_value = seq(min(tmp$response_value, na.rm = TRUE), max(tmp$response_value, na.rm = TRUE), length.out = 100),
      Age = age_median
    )
    
    pred_all <- do.call(
      rbind,
      lapply(sex_levels_for_plot, function(sx) {
        z <- pred_base
        z$Sex_clean <- sx
        z
      })
    )
    
    pred_all$pred_prob <- tryCatch(
      predict(fit, newdata = pred_all, type = "response"),
      error = function(e) NA_real_
    )
    
    pred_grid <- pred_all %>%
      group_by(response_value) %>%
      summarise(pred_prob = mean(pred_prob, na.rm = TRUE), .groups = "drop")
    
    if (any(!is.na(pred_grid$pred_prob))) {
      p_adj_cat <- ggplot(pred_grid, aes(x = response_value, y = pred_prob)) +
        geom_line(linewidth = 1) +
        theme_bw() +
        labs(
          title = wrap_label(paste0(v, " marginal adjusted probability vs ", rv), 60),
          subtitle = paste0(
            "Age fixed at median=", signif(age_median, 3),
            "; Sex averaged; beta=", signif(coef_tab["response_value", "Estimate"], 3),
            ", p=", signif(coef_tab["response_value", "Pr(>|z|)"], 3)
          ),
          x = rv,
          y = paste0("Marginal adjusted probability of ", levels(tmp$outcome_bin)[2])
        )
      
      ggsave(
        file.path("adjusted", paste0("categorical_marginal_adjusted_", safe_filename(rv), "_", safe_filename(v), ".png")),
        p_adj_cat,
        width = 7,
        height = 4.5
      )
    }
  }
  
  # ---- Plot significant multinomial categorical results ----
  sig_keys_multinom <- adj_cat_df %>%
    filter(type == "multinomial", p_adj < FDR_CUTOFF) %>%
    mutate(key = paste(response, variable, sep = "__")) %>%
    pull(key)
  
  message(paste0("Adjusted multinomial categorical plots to generate: ", length(sig_keys_multinom)))
  
  for (key in sig_keys_multinom) {
    if (!(key %in% names(adj_cat_plot_data))) next
    
    rv <- adj_cat_plot_data[[key]]$response
    v <- adj_cat_plot_data[[key]]$variable
    model_df <- adj_cat_plot_data[[key]]$tmp
    fit <- adj_cat_plot_data[[key]]$fit
    
    age_median <- 50
    sex_levels_for_plot <- unique(model_df$Sex_clean)
    sex_levels_for_plot <- sex_levels_for_plot[!is.na(sex_levels_for_plot)]
    
    # Marginalized prediction: Age is fixed at the median, Sex is averaged over.
    # This treats Sex_clean as a confounder rather than displaying sex-specific facets.
    pred_grid <- expand.grid(
      response_value = seq(min(model_df$response_value, na.rm = TRUE), max(model_df$response_value, na.rm = TRUE), length.out = 100),
      Age = age_median,
      Sex_clean = sex_levels_for_plot
    )
    
    probs <- tryCatch(
      predict(fit, newdata = pred_grid, type = "probs"),
      error = function(e) NULL
    )
    
    if (is.null(probs)) next
    
    if (is.vector(probs)) {
      probs <- cbind(probs)
      colnames(probs) <- levels(model_df$outcome_multinom)[1]
    }
    
    prob_df <- cbind(pred_grid, as.data.frame(probs, check.names = FALSE))
    prob_cols <- setdiff(colnames(prob_df), c("response_value", "Age", "Sex_clean"))
    
    prob_long <- do.call(
      rbind,
      lapply(prob_cols, function(cat_name) {
        data.frame(
          response_value = prob_df$response_value,
          Age = prob_df$Age,
          Sex_clean = prob_df$Sex_clean,
          category = cat_name,
          pred_prob = prob_df[[cat_name]],
          stringsAsFactors = FALSE
        )
      })
    ) %>%
      group_by(response_value, category) %>%
      summarise(pred_prob = mean(pred_prob, na.rm = TRUE), .groups = "drop")
    
    p_multinom <- ggplot(prob_long, aes(x = response_value, y = pred_prob, color = category)) +
      geom_line(linewidth = 1) +
      theme_bw() +
      scale_color_discrete(labels = function(x) wrap_label(x, 32)) +
      guides(color = guide_legend(ncol = 1)) +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 13),
        plot.subtitle = element_text(size = 10)
      ) +
      labs(
        title = wrap_label(paste0(v, " marginal adjusted category probabilities vs ", rv), 65),
        subtitle = paste0("Age fixed at median=", signif(age_median, 3), "; Sex averaged; multinomial model; plotted if BH adj p < ", FDR_CUTOFF),
        x = rv,
        y = "Marginal adjusted predicted probability",
        color = wrap_label(v, 35)
      )
    
    ggsave(
      file.path("adjusted", paste0("categorical_multinomial_marginal_adjusted_", safe_filename(rv), "_", safe_filename(v), ".png")),
      p_multinom,
      width = 10,
      height = 5.5
    )
  }
  
  # ---- Binary categorical adjusted forest plot ----
  forest_cat_df <- adj_cat_df %>%
    filter(type == "binary", p_adj < FDR_CUTOFF, !is.na(odds_ratio)) %>%
    mutate(label = paste0(response, ": ", variable)) %>%
    arrange(response, odds_ratio) %>%
    mutate(label = factor(label, levels = label))
  
  if (nrow(forest_cat_df) > 0) {
    p_forest_cat <- ggplot(forest_cat_df, aes(x = odds_ratio, y = label)) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_errorbarh(aes(xmin = or_ci_low, xmax = or_ci_high), height = 0.2) +
      geom_point(size = 2) +
      scale_x_log10() +
      theme_bw() +
      theme(plot.title = element_text(size = 13), plot.subtitle = element_text(size = 10)) +
      labs(
        title = wrap_label("Age/Sex-adjusted binary categorical associations for HC-index / alpha diversity", 70),
        subtitle = paste0("Logistic model odds ratio for response variable; plotted if BH adj p < ", FDR_CUTOFF),
        x = "Adjusted odds ratio per 1-unit increase in response variable (log scale)",
        y = NULL
      )
    
    ggsave(
      file.path("adjusted", "forest_categorical_adjusted_OR.png"),
      p_forest_cat,
      width = 7,
      height = max(4, 0.35 * nrow(forest_cat_df) + 1)
    )
  }
}

# =========================
# Age / Sex adjusted correlation:
# HC-index vs alpha diversity
# =========================
message("\n=== Running adjusted HC-index vs alpha diversity correlations (Age + Sex) ===")

alpha_response_vars <- setdiff(response_vars, "HC_index")
adj_alpha_results <- list()

if ("HC_index" %in% names(adj_df) && length(alpha_response_vars) > 0) {
  for (alpha_var in alpha_response_vars) {
    if (!(alpha_var %in% names(adj_df))) next
    
    tmp <- adj_df %>%
      select(HC_index, all_of(alpha_var), Age, Sex_clean) %>%
      filter(!is.na(HC_index), !is.na(.data[[alpha_var]]), !is.na(Age), !is.na(Sex_clean))
    
    if (nrow(tmp) < 50) next
    if (length(unique(tmp$HC_index)) < 5) next
    if (length(unique(tmp[[alpha_var]])) < 5) next
    
    # Age/Sexの影響をそれぞれから除去した残差同士の相関を取る
    fit_hc <- tryCatch(
      lm(HC_index ~ Age + Sex_clean, data = tmp),
      error = function(e) NULL
    )
    
    fit_alpha <- tryCatch(
      lm(tmp[[alpha_var]] ~ Age + Sex_clean, data = tmp),
      error = function(e) NULL
    )
    
    if (is.null(fit_hc) || is.null(fit_alpha)) next
    
    tmp$HC_index_adj_resid <- resid(fit_hc)
    tmp$alpha_adj_resid <- resid(fit_alpha)
    
    cor_adj <- tryCatch(
      cor.test(tmp$HC_index_adj_resid, tmp$alpha_adj_resid, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(cor_adj)) next
    
    # 比較用：未補正Spearman相関
    cor_raw <- tryCatch(
      cor.test(tmp$HC_index, tmp[[alpha_var]], method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    
    raw_rho <- if (!is.null(cor_raw)) unname(cor_raw$estimate) else NA_real_
    raw_p <- if (!is.null(cor_raw)) cor_raw$p.value else NA_real_
    
    adj_alpha_results[[alpha_var]] <- data.frame(
      alpha_variable = alpha_var,
      n = nrow(tmp),
      raw_spearman_rho = raw_rho,
      raw_p_value = raw_p,
      adjusted_spearman_rho = unname(cor_adj$estimate),
      adjusted_p_value = cor_adj$p.value
    )
    
    # Adjusted residual plot
    p_adj_alpha <- ggplot(tmp, aes(x = HC_index_adj_resid, y = alpha_adj_resid)) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "lm") +
      theme_bw() +
      labs(
        title = paste0("HC-index vs ", alpha_var, " adjusted for Age/Sex"),
        subtitle = paste0(
          "Adjusted Spearman rho=", round(unname(cor_adj$estimate), 3),
          ", p=", signif(cor_adj$p.value, 3)
        ),
        x = "HC-index residual after Age/Sex adjustment",
        y = paste0(alpha_var, " residual after Age/Sex adjustment")
      )
    
    ggsave(
      file.path("adjusted", paste0("adjusted_correlation_HC_index_vs_", safe_filename(alpha_var), ".png")),
      p_adj_alpha,
      width = 5,
      height = 4
    )
  }
}

if (length(adj_alpha_results) > 0) {
  adj_alpha_df <- bind_rows(adj_alpha_results)
  adj_alpha_df$adjusted_p_adj <- p.adjust(adj_alpha_df$adjusted_p_value, method = "BH")
  adj_alpha_df$raw_p_adj <- p.adjust(adj_alpha_df$raw_p_value, method = "BH")
  adj_alpha_df <- adj_alpha_df %>% arrange(adjusted_p_adj)
  
  write.csv(
    adj_alpha_df,
    file.path("adjusted", "HC_index_vs_alpha_adjusted_correlation.csv"),
    row.names = FALSE
  )
  
  # Forest-style summary plot for adjusted HC-index vs alpha correlations
  forest_alpha_df <- adj_alpha_df %>%
    mutate(alpha_variable = factor(alpha_variable, levels = alpha_variable))
  
  p_forest_alpha <- ggplot(forest_alpha_df, aes(x = adjusted_spearman_rho, y = alpha_variable)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      title = "Age/Sex-adjusted correlation between HC-index and alpha diversity",
      subtitle = "Spearman correlation between residuals after Age/Sex adjustment",
      x = "Adjusted Spearman rho",
      y = NULL
    )
  
  ggsave(
    file.path("adjusted", "forest_HC_index_vs_alpha_adjusted_correlation.png"),
    p_forest_alpha,
    width = 6,
    height = max(3, 0.4 * nrow(forest_alpha_df) + 1)
  )
}

message("Adjusted analysis completed.")

# =========================
# Model comparison: HC-index vs alpha diversity (ΔR²)
# =========================
message("\n=== Running model comparison (HC vs Chao/Shannon) ===")

model_comp_results <- list()

if ("HC_index" %in% names(adj_df)) {
  
  alpha_vars <- setdiff(response_vars, "HC_index")
  
  for (alpha_var in alpha_vars) {
    
    if (!(alpha_var %in% names(adj_df))) next
    
    for (outcome in numeric_cols) {
      
      if (outcome %in% outcome_exclude_vars) next
      if (!(outcome %in% names(adj_df))) next
      
      tmp <- adj_df %>%
        select(all_of(outcome), HC_index, all_of(alpha_var), Age, Sex_clean) %>%
        filter(
          !is.na(.data[[outcome]]),
          !is.na(HC_index),
          !is.na(.data[[alpha_var]]),
          !is.na(Age),
          !is.na(Sex_clean)
        )
      
      if (nrow(tmp) < 50) next
      
      # --- 単独モデル ---
      fit_hc <- tryCatch(
        lm(tmp[[outcome]] ~ HC_index + Age + Sex_clean, data = tmp),
        error = function(e) NULL
      )
      
      fit_alpha <- tryCatch(
        lm(tmp[[outcome]] ~ tmp[[alpha_var]] + Age + Sex_clean, data = tmp),
        error = function(e) NULL
      )
      
      # --- 同時モデル ---
      fit_combined <- tryCatch(
        lm(tmp[[outcome]] ~ HC_index + tmp[[alpha_var]] + Age + Sex_clean, data = tmp),
        error = function(e) NULL
      )
      
      if (is.null(fit_hc) || is.null(fit_alpha) || is.null(fit_combined)) next
      
      r2_hc <- summary(fit_hc)$adj.r.squared
      r2_alpha <- summary(fit_alpha)$adj.r.squared
      r2_combined <- summary(fit_combined)$adj.r.squared
      
      delta_r2 <- r2_combined - r2_alpha
      
      model_comp_results[[paste(outcome, alpha_var, sep="__")]] <- data.frame(
        outcome = outcome,
        alpha_metric = alpha_var,
        R2_HC = r2_hc,
        R2_alpha = r2_alpha,
        R2_combined = r2_combined,
        delta_R2_HC_over_alpha = delta_r2
      )
    }
  }
}

if (length(model_comp_results) > 0) {
  
  model_comp_df <- bind_rows(model_comp_results) %>%
    arrange(desc(delta_R2_HC_over_alpha))
  
  write.csv(
    model_comp_df,
    file.path("adjusted", "model_comparison_deltaR2.csv"),
    row.names = FALSE
  )
  
  # --- 可視化（ΔR² forest風） ---
  model_comp_plot_df <- model_comp_df %>%
    mutate(label = paste0(outcome, " | ", alpha_metric)) %>%
    arrange(delta_R2_HC_over_alpha) %>%
    mutate(label = factor(label, levels = label))
  
  p_delta <- ggplot(model_comp_plot_df,
                    aes(x = delta_R2_HC_over_alpha, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      title = "ΔR²: Added value of HC-index over alpha diversity",
      subtitle = "Positive values indicate HC-index adds explanatory power beyond alpha diversity",
      x = "Δ adjusted R² (HC_index added to alpha model)",
      y = NULL
    )
  
  ggsave(
    file.path("adjusted", "deltaR2_forest.png"),
    p_delta,
    width = 7,
    height = max(4, 0.35 * nrow(model_comp_plot_df) + 1)
  )
}

message("Model comparison completed.")