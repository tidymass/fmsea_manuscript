library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)

setwd("~/Projects/fmsea_manuscript")

load("3_data_analysis/5_MTBLS8662/feature_table.rda")
load("2_data/2_KEGG_MS1/kegg_compound_ms1.rda")

annotation_table_final <-
  annotate_feature_table(
    feature_table = feature_table,
    column = "rp",
    metabolite_database = kegg_compound_ms1, # select ms1 database
    database_type = "KEGG", # select database type :KEGG or HMDB
    ms1_match_ppm = 15,
    mfc_rt_tol = 10,
    isotope_number = 3
  )


calculate_redundancy(annotation_table = annotation_table_final)

annotation_table_final2 =
  featuremsea::remove_redundancy(
    annotation_table = annotation_table_final)

calculate_redundancy(annotation_table = annotation_table_final2)

results <- process_annotation_table(
  annotation_table_final2 = annotation_table_final2,
  database_type = "KEGG" # KEGG or HMDB here
)

ranking_table <- results$ranking_table
original_score_annotation <- results$original_score_annotation

save(ranking_table, file = "3_data_analysis/5_MTBLS8662/ranking_table_kegg_snr.rda")
save(original_score_annotation, 
     file = "3_data_analysis/5_MTBLS8662/original_score_annotation_kegg.rda", 
     compress = "xz", 
     compression_level = 9)

# load("2_data/3_Xia_Case_Study_v2/fmsea_input/human_pathway_database.rda")
# load("2_data/3_Xia_Case_Study_v2/fmsea_input/original_score_annotation.rda")
# load("2_data/3_Xia_Case_Study_v2/fmsea_input/ranking_table.rda")
# 
# pathway_dataset <- human_pathway_all_clean %>% 
#   dplyr::rename(KEGG_ID = KEGG_ID_clean)

annotation_table <- original_score_annotation

# for test
pathway_database <- metpath::kegg_hsa_pathway

results <- perform_fmsea_analysis(
  pathway_database,
  annotation_table,         
  ranking_table,
  threads = 6,
  min.compounds.num = 5,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  perm.num = 1000,
  seed = 123,
  fdr.thr = 0.05,
  max.iter.num = 1,
  verbose = TRUE
)




