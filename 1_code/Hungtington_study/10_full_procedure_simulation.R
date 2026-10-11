library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)

load("3_data_analysis/5_MTBLS749/ranking_table_kegg_snr_pre.rda")
load("3_data_analysis/5_MTBLS749/original_score_annotation_kegg_pre.rda")

annotation_table <- original_score_annotation

# for test
pathway_database <- metpath::kegg_hsa_pathway

results_iterations <- perform_fmsea_analysis(
  pathway_database,
  annotation_table,         
  ranking_table,
  threads = 50,
  min.compounds.num = 15,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  perm.num = 10000,
  seed = 123,
  fdr.thr = 0.1,
  verbose = TRUE
)

permutation_results <- perform_fmsea_iterative_permutation(
  pathway_database,
  annotation_table,
  ranking_table,
  outer.perm.num = 1000,      
  inner.perm.num = 10000,    
  threads = 50,
  min.compounds.num = 15,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  seed = 123,
  fdr.thr = 0.1,
  verbose = TRUE
)

save(permutation_results, file = "~/tidymass/fmsea_manuscript/3_data_analysis/permutation_results.rda")

results <- test$results

kegg_pathway <- metpath::kegg_hsa_pathway
kegg_database <- convert_kegg2fmsea(kegg_pathway)

results2 <- results %>% 
  left_join(kegg_database[,c("KEGG_ID", "pathway_id")], by = "pathway_id")

significant_modules <- results2 %>% 
  dplyr::select(pathway_id, pathway_name, KEGG_ID, ES, NES_iter_aware, p_value_iter_aware, FDR_iter_aware) %>% 
  dplyr::rename(NES = NES_iter_aware,
                p_value = p_value_iter_aware,
                FDR = FDR_iter_aware) %>% 
  dplyr::filter(FDR < 0.1)

write.csv(significant_modules, file = "3_data_analysis/5_MTBLS749/test_permutation.csv")
