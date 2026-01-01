library(r4projects)
setwd(get_project_wd())
rm(list = ls())

load("2_data/3_Xia_Case_Study_v2/data_cleaning/POS/object_pos2.rda")
load("2_data/3_Xia_Case_Study_v2/data_cleaning/NEG/object_neg2.rda")

# ---- helper: compute SNR for each feature ----
compute_snr <- function(expr_mat, case_ids, control_ids) {
  # expr_mat: features x samples (rows = features, cols = samples)
  case_mat <- expr_mat[, case_ids, drop = FALSE]
  ctrl_mat <- expr_mat[, control_ids, drop = FALSE]
  
  mu_case <- rowMeans(case_mat, na.rm = TRUE)
  mu_ctrl <- rowMeans(ctrl_mat, na.rm = TRUE)
  
  sd_case <- apply(case_mat, 1, sd, na.rm = TRUE)
  sd_ctrl <- apply(ctrl_mat, 1, sd, na.rm = TRUE)
  
  denom <- sd_case + sd_ctrl
  denom[denom == 0 | is.na(denom)] <- NA_real_
  
  (mu_case - mu_ctrl) / denom
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

# object_pos <-
#   mutate_fc(object = object_pos2,
#             control_sample_id = control_sample_id,
#             case_sample_id = case_sample_id,
#             mean_median = "mean")

# raw + log2(x+1)
expression_data_pos_raw <- object_pos2@expression_data
expression_data_pos <- log2(as.matrix(expression_data_pos_raw) + 1)

# robust sample matching
control_sample_id_pos <- intersect(control_sample_id, colnames(expression_data_pos))
case_sample_id_pos    <- intersect(case_sample_id,    colnames(expression_data_pos))

snr_pos <- compute_snr(expression_data_pos, case_sample_id_pos, control_sample_id_pos)

mean_intensity_pos <- rowMeans(expression_data_pos_raw, na.rm = TRUE)

feature_table_pos <- 
  object_pos2@variable_info %>% 
  mutate(condition = snr_pos,
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

# object_neg <-
#   mutate_fc(object = object_neg2,
#             control_sample_id = control_sample_id,
#             case_sample_id = case_sample_id,
#             mean_median = "mean")

expression_data_neg_raw <- object_neg2@expression_data
expression_data_neg <- log2(as.matrix(expression_data_neg_raw) + 1)

control_sample_id_neg <- intersect(control_sample_id, colnames(expression_data_neg))
case_sample_id_neg    <- intersect(case_sample_id,    colnames(expression_data_neg))

snr_neg <- compute_snr(expression_data_neg, case_sample_id_neg, control_sample_id_neg)

mean_intensity_neg <- rowMeans(expression_data_neg_raw, na.rm = TRUE)

feature_table_neg <- 
  object_neg2@variable_info %>% 
  mutate(condition = snr_neg,
         polarity = "negative",
         mean_intensity = mean_intensity_neg)

# =========================
# Combine + save
# =========================
feature_table <- bind_rows(feature_table_pos, feature_table_neg)
rownames(feature_table) <- 1:nrow(feature_table)

save(feature_table, file = "3_data_analysis/3_Xia_case_study/feature_table_snr.rda")
