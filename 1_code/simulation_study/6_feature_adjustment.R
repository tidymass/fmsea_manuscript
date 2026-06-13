library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(tidyverse)

load("3_data_analysis/3_Xia_case_study/ranking_table_kegg.rda")
load("2_data/3_Xia_Case_Study_v2/feature_ids/ids_kegg.rda")

ranking_table <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% ids, 
                                 ranking_weight + 3, 
                                 ranking_weight))

load("3_data_analysis/3_Xia_case_study/human_pathway_database.rda")
load("3_data_analysis/3_Xia_case_study/original_score_annotation_kegg.rda")

results <- perform_fmsea_analysis(
  pathway_database,
  annotation_table = original_score_annotation,         
  ranking_table,
  n_cores = 20,
  min_compounds = 5,
  max_compounds = 300,
  id_col = "KEGG_ID",
  n_perm = 1000,
  seed = 123,
  fdr_threshold = 0.05,
  max_iter = 10,
  verbose = TRUE
)

