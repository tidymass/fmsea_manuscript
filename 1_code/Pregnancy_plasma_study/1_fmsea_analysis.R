library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)

load("3_data_analysis/7_pregnancy_plasma/feature_table.rda")
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

annotation_table <- original_score_annotation

# for test
pathway_database <- metpath::kegg_hsa_pathway

results <- perform_fmsea_analysis(
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
  max.iter.num = 3,
  verbose = TRUE
)

save(results, file = "3_data_analysis/7_pregnancy_plasma/results_iteration3.rda ")
