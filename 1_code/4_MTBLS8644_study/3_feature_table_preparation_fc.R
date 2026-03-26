library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(tidymass)

load("data_cleaning/POS/object_pos2.rda")
load("data_cleaning/NEG/object_neg2.rda")

# ---- helper: compute log2FC for each feature ----
compute_log2fc <- function(expr_mat, case_ids, control_ids) {
  case_mat <- expr_mat[, case_ids, drop = FALSE]
  ctrl_mat <- expr_mat[, control_ids, drop = FALSE]
  
  mu_case <- rowMeans(case_mat, na.rm = TRUE)
  mu_ctrl <- rowMeans(ctrl_mat, na.rm = TRUE)
  
  mu_ctrl[mu_ctrl == 0 | is.na(mu_ctrl)] <- NA_real_
  
  log2(mu_case / mu_ctrl)
}

# =========================
# POS
# =========================
control_sample_id <-
  object_pos2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Control") %>%
  dplyr::pull(sample_id)

case_sample_id <-
  object_pos2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Case") %>%
  dplyr::pull(sample_id)

expression_data_pos_raw <- object_pos2@expression_data

control_sample_id_pos <- intersect(control_sample_id, colnames(expression_data_pos_raw))
case_sample_id_pos    <- intersect(case_sample_id,    colnames(expression_data_pos_raw))

log2fc_pos <- compute_log2fc(expression_data_pos_raw, case_sample_id_pos, control_sample_id_pos)
mean_intensity_pos <- rowMeans(expression_data_pos_raw, na.rm = TRUE)

feature_table_pos <-
  object_pos2@variable_info %>%
  mutate(condition = log2fc_pos,
         polarity = "positive",
         mean_intensity = mean_intensity_pos)

# =========================
# NEG
# =========================
control_sample_id <-
  object_neg2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Control") %>%
  dplyr::pull(sample_id)

case_sample_id <-
  object_neg2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "Case") %>%
  dplyr::pull(sample_id)

expression_data_neg_raw <- object_neg2@expression_data

control_sample_id_neg <- intersect(control_sample_id, colnames(expression_data_neg_raw))
case_sample_id_neg    <- intersect(case_sample_id,    colnames(expression_data_neg_raw))

log2fc_neg <- compute_log2fc(expression_data_neg_raw, case_sample_id_neg, control_sample_id_neg)
mean_intensity_neg <- rowMeans(expression_data_neg_raw, na.rm = TRUE)

feature_table_neg <-
  object_neg2@variable_info %>%
  mutate(condition = log2fc_neg,
         polarity = "negative",
         mean_intensity = mean_intensity_neg)

# =========================
# Combine + save
# =========================
feature_table <- bind_rows(feature_table_pos, feature_table_neg)
rownames(feature_table) <- 1:nrow(feature_table)

save(feature_table, file = "~/Projects/fmsea_manuscript/3_data_analysis/4_MTBLS265/feature_table.rda")