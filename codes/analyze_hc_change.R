#!/usr/bin/env Rscript
# Reanalysis of Inoue et al. (2025), doi:10.3390/microorganisms13092068.
# Run with Rscript scripts/analyze_hc_change.R from the repository root.
# Also supports invocation by absolute script path from another working directory.
suppressPackageStartupMessages(library(ggplot2))
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo_root <- if (length(script_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
} else getwd()
d <- read.csv(file.path(repo_root, "data", "hc_index_microorganisms_2025_2068.csv"))
stopifnot(all(c("participant_id", "week", "fiber_group", "hc_index") %in% names(d)),
          !anyDuplicated(d[c("participant_id", "week")]),
          all(d$week %in% c(0, 2, 4)), all(d$fiber_group %in% c("Low", "High")),
          all(is.finite(d$hc_index)))
stopifnot(all(vapply(split(d$fiber_group, d$participant_id), function(g) length(unique(g)) == 1L, logical(1))))
baseline <- d[d$week == 0, c("participant_id", "hc_index")]
names(baseline)[2] <- "baseline_HC"
stopifnot(all(d$participant_id %in% baseline$participant_id))
d <- merge(d, baseline, by = "participant_id", sort = FALSE)
d$Subject <- d$participant_id
d$Period <- paste0(d$week, "w")
d$HC_index <- d$hc_index
d$group <- factor(d$fiber_group, levels = c("Low", "High"))
# Independent ANCOVA at weeks 2 and 4, with a common baseline slope.
# Predictions use the pooled baseline mean among participants observed at that week.
# Within-group change subtracts this mean from the predicted follow-up HC-index.
# Standard errors treat the reference baseline mean as fixed.
# Two-sided t-tests are unadjusted for multiple comparisons; intervals are 95% CIs.
rows <- lapply(c("Final", "Mid"), function(time_label) {
  period <- if (time_label == "Final") "4w" else "2w"
  x <- d[d$Period == period & complete.cases(d[c("HC_index", "baseline_HC", "group")]), ]
  stopifnot(all(table(x$group) > 0), !anyDuplicated(x$Subject))
  fit <- lm(HC_index ~ baseline_HC + group, data = x)
  baseline_mean <- mean(x$baseline_HC)
  newdata <- data.frame(baseline_HC = baseline_mean,
                        group = factor(c("Low", "High"), levels = levels(x$group)))
  design <- model.matrix(delete.response(terms(fit)), newdata)
  contrasts <- rbind(design, design[2, ] - design[1, ])
  estimate <- as.vector(contrasts %*% coef(fit)) - c(baseline_mean, baseline_mean, 0)
  se <- sqrt(diag(contrasts %*% vcov(fit) %*% t(contrasts)))
  df <- df.residual(fit)
  margin <- qt(0.975, df) * se
  data.frame(time_label = time_label, group = c("Low", "High", "High_vs_Low"),
             n = c(sum(x$group == "Low"), sum(x$group == "High"), nrow(x)),
             mean_baseline = baseline_mean, adjusted_change = estimate, se = se,
             df = df, p_value = 2 * pt(abs(estimate / se), df, lower.tail = FALSE),
             conf_low = estimate - margin, conf_high = estimate + margin,
             row.names = NULL)
})
result <- do.call(rbind, rows)

out_dir <- file.path(repo_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(out_dir, "adjusted_change_and_group_difference.csv"), row.names = FALSE)
plot_dir <- out_dir
format_p <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p))
}
plot_df <- result[result$group %in% c("Low", "High"), ]
plot_df$time_label <- factor(plot_df$time_label, levels = c("Mid", "Final"),
                             labels = c("2wk", "4wk"))
plot_df$group <- factor(plot_df$group, levels = c("Low", "High"))
plot_df$p_label <- paste0("p=", format_p(plot_df$p_value))
plot_df$label_y <- plot_df$conf_high + 0.055

p <- ggplot(plot_df, aes(x = time_label, y = adjusted_change, color = group)) +
  geom_hline(yintercept = 0, color = "grey72", linewidth = 0.4) +
  geom_pointrange(
    aes(ymin = conf_low, ymax = conf_high),
    position = position_dodge(width = 0.45),
    linewidth = 0.65,
    size = 0.65
  ) +
  geom_text(
    aes(y = label_y, label = p_label),
    position = position_dodge(width = 0.45),
    size = 3.1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(Low = "#2f6f73", High = "#6f9f3a"), labels = c(Low = "Low complex fiber", High = "High complex fiber")) +
  labs(
    x = NULL,
    y = "Baseline-adjusted change in HC-index",
    color = "Fiber group",
    title = "Inoue et al. (2025): baseline-adjusted HC-index change",
    subtitle = "Within-group estimates from follow-up HC-index ~ baseline HC-index + High/Low"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(plot_dir, "baseline_adjusted_hc_change.png"),
  p,
  width = 6.8,
  height = 4.8,
  dpi = 300
)
ggsave(
  file.path(plot_dir, "baseline_adjusted_hc_change.pdf"),
  p,
  width = 6.8,
  height = 4.8,
  device = cairo_pdf
)



# Participant-level change correlations, pooled across fiber groups.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
dat <- read_csv(file.path(repo_root, "data", "hc_index_microorganisms_2025_2068.csv"),
                show_col_types = FALSE) %>% rename(HC_index = hc_index)
ratio_specs <- tibble(
  ratio_metric = c("Propionate_over_Acetate", "nButyrate_over_Acetate", "Propionate_over_nButyrate"),
  numerator = c("Propionate", "nButyrate", "Propionate"),
  denominator = c("Acetate", "Acetate", "nButyrate")
)
missing_ratio_cols <- setdiff(unique(c(ratio_specs$numerator, ratio_specs$denominator)), names(dat))
if (length(missing_ratio_cols) > 0) {
  stop("Missing columns required for ratio metrics: ", paste(missing_ratio_cols, collapse = ", "))
}
for (i in seq_len(nrow(ratio_specs))) {
  num <- ratio_specs$numerator[i]
  den <- ratio_specs$denominator[i]
  out <- ratio_specs$ratio_metric[i]
  dat[[out]] <- ifelse(is.na(dat[[den]]) | dat[[den]] == 0, NA_real_, dat[[num]] / dat[[den]])
}

analysis_df <- dat %>% mutate(Period = factor(paste0(week, "w"), levels = c("0w", "2w", "4w")))
numeric_metrics <- names(dat)[vapply(dat, is.numeric, logical(1))]
metadata_metrics <- setdiff(numeric_metrics, c("HC_index", "week"))
ratio_metrics <- ratio_specs$ratio_metric
comparisons <- tibble(
  period_ref = c("0w", "0w", "2w"),
  period_test = c("2w", "4w", "4w"),
  comparison = c("2w_vs_0w", "4w_vs_0w", "4w_vs_2w")
)

make_delta_table <- function(df, metrics, p_ref, p_test) {
  wide <- df %>%
    filter(Period %in% c(p_ref, p_test)) %>%
    select(participant_id, fiber_group, Period, all_of(c("HC_index", metrics))) %>%
    distinct() %>%
    pivot_wider(names_from = Period, values_from = all_of(c("HC_index", metrics)), names_sep = "__")

  hc_ref <- paste0("HC_index__", p_ref)
  hc_test <- paste0("HC_index__", p_test)
  wide <- wide %>%
    filter(!is.na(.data[[hc_ref]]), !is.na(.data[[hc_test]])) %>%
    mutate(delta_HC_index = .data[[hc_test]] - .data[[hc_ref]])

  bind_rows(lapply(metrics, function(metric) {
    ref_col <- paste0(metric, "__", p_ref)
    test_col <- paste0(metric, "__", p_test)
    if (!all(c(ref_col, test_col) %in% names(wide))) {
      return(NULL)
    }
    wide %>%
      transmute(
        participant_id,
        fiber_group,
        metric = metric,
        delta_HC_index,
        delta_metadata = .data[[test_col]] - .data[[ref_col]]
      )
  }))
}

delta_long <- bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
  make_delta_table(
    analysis_df,
    metadata_metrics,
    comparisons$period_ref[i],
    comparisons$period_test[i]
  ) %>%
    mutate(
      comparison = comparisons$comparison[i],
      period_ref = comparisons$period_ref[i],
      period_test = comparisons$period_test[i],
      .before = 1
    )
}))
write_csv(delta_long, file.path(out_dir, "hc_delta_metadata_delta_long.csv"))

cor_one <- function(df, method) {
  ok <- complete.cases(df$delta_HC_index, df$delta_metadata)
  x <- df$delta_HC_index[ok]
  y <- df$delta_metadata[ok]
  if (length(x) < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(n = length(x), estimate = NA_real_, p_value = NA_real_))
  }
  ct <- cor.test(x, y, method = method, exact = FALSE)
  tibble(n = length(x), estimate = unname(ct$estimate), p_value = ct$p.value)
}

cor_results <- delta_long %>%
  group_by(comparison, period_ref, period_test, metric) %>%
  group_modify(~ bind_rows(
    cor_one(.x, "spearman") %>% mutate(method = "spearman"),
    cor_one(.x, "pearson") %>% mutate(method = "pearson")
  )) %>%
  ungroup() %>%
  group_by(method, comparison) %>%
  mutate(p_fdr_by_comparison = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(method) %>%
  mutate(p_fdr_all_comparisons = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(method, p_value)
write_csv(cor_results, file.path(out_dir, "hc_delta_metadata_delta_correlations.csv"))

top_pearson <- cor_results %>%
  filter(method == "pearson", !is.na(p_value)) %>%
  group_by(comparison) %>%
  slice_min(order_by = p_value, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    metric = factor(metric, levels = rev(unique(metric[order(p_value)]))),
    comparison = factor(comparison, levels = comparisons$comparison)
  )
write_csv(top_pearson, file.path(out_dir, "hc_delta_metadata_delta_top_pearson_heatmap_data.csv"))

p_heat_pearson <- ggplot(top_pearson, aes(x = comparison, y = metric, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("r %.2f\np %.3g", estimate, p_value)), size = 2.8, lineheight = 0.9) +
  scale_fill_gradient2(low = "#2f6f73", mid = "white", high = "#b35a3a", midpoint = 0, name = "Pearson r") +
  labs(
    x = NULL,
    y = NULL,
    title = "Inoue et al. (2025): HC-index and metadata changes"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(plot_dir, "hc_delta_metadata_delta_top_pearson_heatmap.png"), p_heat_pearson, width = 7.5, height = 7.5, dpi = 300)
ggsave(file.path(plot_dir, "hc_delta_metadata_delta_top_pearson_heatmap.pdf"), p_heat_pearson, width = 7.5, height = 7.5)


writeLines(c(paste("R", getRversion()),
             vapply(c("ggplot2", "dplyr", "tidyr", "readr"),
                    function(pkg) paste(pkg, packageVersion(pkg)), character(1))),
           file.path(out_dir, "software_versions.txt"))
message("Analysis complete: adjusted changes and change-correlation heatmap saved in results/.")
