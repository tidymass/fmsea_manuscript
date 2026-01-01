library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)

load("3_data_analysis/1_MSUD/POS/object")
object_pos <- object

load("3_data_analysis/1_MSUD/NEG/object")
object_neg <- object

outlier_samples =
  object_pos %>%
  `+`(1) %>% 
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples


outlier_samples =
  object_neg %>%
  `+`(1) %>% 
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

object_pos <- 
  impute_mv(object = object_pos, method = "knn")

object_pos <- 
  normalize_data(object_pos, method = "median")

object_pos2 <- 
  integrate_data(object_pos, method = "subject_median")

object_neg <- 
  impute_mv(object = object_neg, method = "knn")

object_neg <- 
  normalize_data(object_neg, method = "median")

object_neg2 <- 
  integrate_data(object_neg, method = "subject_median")



control_sample_id =
  object_pos2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Control") %>%
  pull(sample_id)

case_sample_id =
  object_pos2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Case") %>%
  pull(sample_id)

object_pos <-
  mutate_fc(object = object_pos2, 
            control_sample_id = control_sample_id, 
            case_sample_id = case_sample_id, 
            mean_median = "mean")

expression_data_pos <- object_pos2@expression_data

mean_intensity <-
  expression_data_pos %>%
  rowMeans(na.rm = TRUE)

feature_table_pos <- 
  object_pos@variable_info %>% 
  mutate(condition = log2(fc),
         polarity = "positive",
         mean_intensity = mean_intensity) %>% 
  dplyr::select(-fc)

control_sample_id =
  object_neg2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Control") %>%
  pull(sample_id)

case_sample_id =
  object_neg2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Case") %>%
  pull(sample_id)

object_neg <-
  mutate_fc(object = object_neg2, 
            control_sample_id = control_sample_id, 
            case_sample_id = case_sample_id, 
            mean_median = "mean")

expression_data_neg <- object_neg2@expression_data

mean_intensity <-
  expression_data_neg %>%
  rowMeans(na.rm = TRUE)

feature_table_neg <- 
  object_neg@variable_info %>% 
  mutate(condition = log2(fc),
         polarity = "negative",
         mean_intensity = mean_intensity) %>% 
  dplyr::select(-fc)

feature_table <- 
  bind_rows(feature_table_pos, feature_table_neg)

rownames(feature_table) <- 1:nrow(feature_table)

save(feature_table, file = "3_data_analysis/2_feature_table/feature_table2.rda")
