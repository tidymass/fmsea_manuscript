library(r4projects)
setwd(get_project_wd())
rm(list = ls())

load("2_data/3_Xia_Case_Study_v2/data_cleaning/POS/object_pos2.rda")
load("2_data/3_Xia_Case_Study_v2/data_cleaning/NEG/object_neg2.rda")

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

save(feature_table, file = "3_data_analysis/3_Xia_case_study/feature_table.rda")
