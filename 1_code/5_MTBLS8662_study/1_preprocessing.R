library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)
library(tidyverse)

setwd("/Volumes/Western/projects/fmsea_cases/MTBLS8662")

process_data(
  path = "Result/POS",
  polarity = "positive",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 4,
  output_tic = FALSE,
  output_bpc = FALSE,
  output_rt_correction_plot = FALSE,
  min_fraction = 0.5
)

process_data(
  path = "Result/NEG",
  polarity = "negative",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 4,
  output_tic = FALSE,
  output_bpc = FALSE,
  output_rt_correction_plot = FALSE,
  min_fraction = 0.5
)

load("Result/POS/Result/object")
object_pos <- object
object_pos

load("Result/NEG/Result/object")
object_neg <- object
object_neg

massdataset::show_sample_missing_values(object = object_pos,
                                        color_by = "group",
                                        order_by = "injection.order",
                                        percentage = TRUE) +
  theme(axis.text.x = element_text(size = 2)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

outlier_samples =
  object_pos %>%
  `+`(1) %>% 
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

get_mv_number(object_pos)

object_pos <- 
  impute_mv(object = object_pos, method = "knn")

get_mv_number(object_pos)

get_mv_number(object_neg)

object_neg <- 
  impute_mv(object = object_neg, method = "knn")

get_mv_number(object_neg)

object_pos <- 
  normalize_data(object_pos, method = "median")

object_pos2 <- 
  integrate_data(object_pos, method = "subject_median")

object_neg <- 
  normalize_data(object_neg, method = "median")

object_neg2 <- 
  integrate_data(object_neg, method = "subject_median")

save(object_pos2, file = "data_cleaning/POS/object_pos2.rda")
save(object_neg2, file = "data_cleaning/NEG/object_neg2.rda")

# load("2_data/3_Xia_Case_Study_v2/data_cleaning/POS/object_pos2.rda")
# load("2_data/3_Xia_Case_Study_v2/data_cleaning/NEG/object_neg2.rda")
# 
# control_sample_id =
#   object_pos2 %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   filter(group == "Control") %>%
#   pull(sample_id)
# 
# case_sample_id =
#   object_pos2 %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   filter(group == "Case") %>%
#   pull(sample_id)
# 
# object_pos <-
#   mutate_fc(object = object_pos2, 
#             control_sample_id = control_sample_id, 
#             case_sample_id = case_sample_id, 
#             mean_median = "mean")
# 
# object_pos <- mutate_p_value(
#   object = object_pos,
#   control_sample_id = control_sample_id,
#   case_sample_id = case_sample_id,
#   method = "wilcox.test",           
#   p_adjust_method = "BH"        # Benjamini-Hochberg FDR 校正
# )
# 
# control_sample_id =
#   object_neg2 %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   filter(group == "Control") %>%
#   pull(sample_id)
# 
# case_sample_id =
#   object_neg2 %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   filter(group == "Case") %>%
#   pull(sample_id)
# 
# object_neg <-
#   mutate_fc(object = object_neg2,
#             control_sample_id = control_sample_id,
#             case_sample_id = case_sample_id,
#             mean_median = "mean")
# 
# object_neg <- mutate_p_value(
#   object = object_neg,
#   control_sample_id = control_sample_id,
#   case_sample_id = case_sample_id,
#   method = "wilcox.test",           
#   p_adjust_method = "BH"        # Benjamini-Hochberg FDR 校正
# )
#   
# feature_table_pos <- object_pos@variable_info
# feature_table_neg <- object_neg@variable_info
# feature_table <- rbind(feature_table_pos, feature_table_neg)
# 
# rownames(feature_table) <- NULL
# 
# feature_table2 <- 
#   feature_table %>%
#   dplyr::mutate(
#     condition = log2(fc),
#     polarity = ifelse(str_detect(variable_id, "POS"), "positive", "negative")) %>% 
#   dplyr::select(variable_id, mz, rt, condition, polarity)






# metabolite_name <- "M582T46_POS"
# metabolite_index <- which(rownames(object_pos@expression_data) == metabolite_name)
# 
# # 获取 case 样本在 expression_data 中的列索引
# case_index <- match(case_sample_id, colnames(object_pos@expression_data))
# 
# # 将该代谢物在 case 样本中的值乘以 10
# object_pos@expression_data[metabolite_index, case_index] <- 
#   object_pos@expression_data[metabolite_index, case_index] * 0.1















