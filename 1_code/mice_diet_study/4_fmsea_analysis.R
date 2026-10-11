library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)

load("3_data_analysis/8_mice_diet/feature_table.rda")
load("2_data/2_KEGG_MS1/kegg_compound_ms1.rda")

annotation_table_final <-
  annotate_feature_table(
    feature_table = feature_table,
    column = "hilic",
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

save(ranking_table, file = "3_data_analysis/8_mice_diet/ranking_table_kegg_snr_pre.rda")
save(original_score_annotation,
     file = "3_data_analysis/8_mice_diet/original_score_annotation_kegg_pre.rda",
     compress = "xz",
     compression_level = 9)

load("3_data_analysis/8_mice_diet/ranking_table_kegg_snr_pre.rda")
load("3_data_analysis/8_mice_diet/original_score_annotation_kegg_pre.rda")

annotation_table <- original_score_annotation

load("3_data_analysis/8_mice_diet/mmu_pathway.rda")
pathway_database <- mmu_pathway

results_original <- perform_fmsea_analysis(
  pathway_database,
  annotation_table,         
  ranking_table,
  threads = 6,
  min.compounds.num = 15,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  perm.num = 10000,
  seed = 123,
  fdr.thr = 0.05,
  max.iter.num = 1,
  verbose = TRUE
)

results_iterations <- perform_fmsea_analysis(
  pathway_database,
  annotation_table,         
  ranking_table,
  threads = 6,
  min.compounds.num = 15,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  perm.num = 10000,
  seed = 123,
  fdr.thr = 0.05,
  # max.iter.num = 1,
  verbose = TRUE
)

save(results_original, file = "3_data_analysis/8_mice_diet/results_original.rda")
save(results_iterations, file = "3_data_analysis/8_mice_diet/results_iterations.rda")
