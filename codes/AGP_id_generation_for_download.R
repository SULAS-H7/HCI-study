library(data.table)
library(dplyr)

df <- fread(
  "10317_20260402-160702.txt",
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = TRUE,
  data.table = FALSE
)

meta_raw <- df %>%
  transmute(
    sample_name = as.character(sample_name),
    diet_type,
    red_meat_frequency,
    vegetable_frequency,
    meat_eggs_frequency,
    fruit_frequency,
    IBD = ibd_diagnosed_as_type
  )

meta_clean <- meta_raw %>%
  mutate(across(everything(), ~ na_if(., ""))) %>%
  mutate(across(everything(), ~ na_if(., "Not provided"))) %>%
  mutate(across(everything(), ~ na_if(., "not provided"))) %>%
  mutate(across(everything(), ~ na_if(., "Unspecified"))) %>%
  mutate(across(everything(), ~ na_if(., "unspecified"))) %>%
  distinct(sample_name, .keep_all = TRUE) %>%
  mutate(
    IBD_flag = case_when(
      IBD %in% c(
        "Ulcerative Colitis",
        "Colonic Crohn's disease",
        "Ileal Crohn's disease",
        "Ileal and Colonic Crohn's disease"
      ) ~ "IBD",
      IBD == "control sample" ~ "Control",
      TRUE ~ NA_character_
    )
  )

meta_ibd_all <- meta_clean %>%
  filter(!is.na(IBD_flag))

table(meta_ibd_all$IBD_flag, useNA = "ifany")

set.seed(123)

meta_ibd_cases <- meta_ibd_all %>%
  filter(IBD_flag == "IBD")

meta_ibd_controls <- meta_ibd_all %>%
  filter(IBD_flag == "Control") %>%
  sample_n(min(2000, n()))

meta_ibd_sub <- bind_rows(meta_ibd_cases, meta_ibd_controls) %>%
  distinct(sample_name, .keep_all = TRUE)

table(meta_ibd_sub$IBD_flag)

meta_diet <- meta_clean %>%
  filter(
    !is.na(sample_name),
    !is.na(diet_type),
    !is.na(red_meat_frequency),
    !is.na(vegetable_frequency),
    !is.na(meat_eggs_frequency),
    !is.na(fruit_frequency)
  )

nrow(meta_diet)
table(meta_diet$diet_type)
table(meta_diet$IBD_flag, useNA = "ifany")

set.seed(123)

meta_diet_sub <- meta_diet %>%
  group_by(diet_type) %>%
  sample_n(min(1000, n())) %>%
  ungroup()

nrow(meta_diet_sub)
table(meta_diet_sub$diet_type)
table(meta_diet_sub$IBD_flag, useNA = "ifany")

meta_export_diet <- meta_diet_sub
colnames(meta_export_diet)[colnames(meta_export_diet) == "sample_name"] <- "#SampleID"

write.table(
  meta_export_diet,
  file = "meta_candidates_diet.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  meta_export_diet$`#SampleID`,
  file = "sample_ids_candidates_diet.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


meta_export_ibd <- meta_ibd_sub
colnames(meta_export_ibd)[colnames(meta_export_ibd) == "sample_name"] <- "#SampleID"

write.table(
  meta_export_ibd,
  file = "meta_candidates_ibd.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  meta_export_ibd$`#SampleID`,
  file = "sample_ids_candidates_ibd.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
