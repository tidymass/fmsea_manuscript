library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)
library(tidyverse)

process_data(
  path = "2_data/8_mice_diet/hilic/POS",
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
  path = "2_data/8_mice_diet/hilic/NEG",
  polarity = "negative",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 4,
  output_tic = FALSE,
  output_bpc = FALSE,
  output_rt_correction_plot = FALSE,
  min_fraction = 0.5
)

load("2_data/8_mice_diet/hilic/POS/Result/object")
object_pos <- object
object_pos

load("2_data/8_mice_diet/hilic/NEG/Result/object")
object_neg <- object
object_neg

object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(group)

control_id =
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Control") %>%
  pull(sample_id)

case_id =
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Case") %>%
  pull(sample_id)

object_pos =
  object_pos %>%
  mutate_variable_na_freq(according_to_samples = control_id) |> 
  mutate_variable_na_freq(according_to_samples = case_id)


object_pos <-
  object_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.5 | na_freq.1 < 0.5)

##
object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(group)

control_id =
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Control") %>%
  pull(sample_id)

case_id =
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Case") %>%
  pull(sample_id)

object_neg =
  object_neg %>%
  mutate_variable_na_freq(according_to_samples = control_id) |> 
  mutate_variable_na_freq(according_to_samples = case_id)



object_neg <-
  object_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.5 | na_freq.1 < 0.5)

outlier_samples =
  object_pos %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples =
  object_neg %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

object_pos <-
  impute_mv(object = object_pos, method = "knn")

object_neg <-
  impute_mv(object = object_neg, method = "knn")

object_pos <-
  normalize_data(object_pos, method = "median")

object_pos2 <-
  integrate_data(object_pos, method = "subject_median")

object_neg <-
  normalize_data(object_neg, method = "median")

object_neg2 <-
  integrate_data(object_neg, method = "subject_median")

save(object_pos2, file = "2_data/8_mice_diet/data_cleaning/POS/object_pos2.rda")
save(object_neg2, file = "2_data/8_mice_diet/data_cleaning/NEG/object_neg2.rda")
