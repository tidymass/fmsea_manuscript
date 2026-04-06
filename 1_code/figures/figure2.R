library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)
library(featuremsea)
library(tidyverse)

load("3_data_analysis/3_Xia_case_study/ranking_table_kegg_snr.rda")
load("3_data_analysis/3_Xia_case_study/original_score_annotation_kegg.rda")
load("2_data/3_Xia_Case_Study_v2/pathway_ratio/pathway_database.rda")
load("2_data/3_Xia_Case_Study_v2/feature_ids/hsa04978.rda")

kegg_database <- pathway_database

set.seed(123) 

ranking_table_new <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% ids, 
                                 ranking_weight + rnorm(n(), mean = 1, sd = 0.25), 
                                 ranking_weight))

# load("3_data_analysis/3_Xia_case_study/human_pathway_database.rda")
pathway_database <- metpath::kegg_hsa_pathway

# pathway_database <- imetpd_pathway

results <- perform_fmsea_analysis(
  pathway_database,
  annotation_table,         
  ranking_table = ranking_table_new,
  threads = 6,
  min.compounds.num = 5,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  perm.num = 10000,
  seed = 123,
  fdr.thr = 0.05,
  max.iter.num = 1,
  verbose = TRUE
)

# kegg_database <- metpath::kegg_hsa_pathway
# kegg_database <- convert_kegg2fmsea(kegg_database)

kegg_pathway <- data.frame(
  pathway_id = names(results@res_list),
  p_value = sapply(results@res_list, function(x) x$p_value[1]),
  n_leading_edge = sapply(results@res_list, function(x) length(unique(x$leading_edge$KEGG_ID))),
  row.names = NULL
)

kegg_pathway$fdr <- p.adjust(kegg_pathway$p_value, method = "BH")

kegg_pathway <- kegg_pathway %>%
  left_join(kegg_database %>% dplyr::select(pathway_id, pathway_name, KEGG_ID, match_ratio), by = "pathway_id")

kegg_pathway <- kegg_pathway %>% 
  dplyr::mutate(metabolite_all_number = str_count(KEGG_ID, "\\{\\}") + 1)

# kegg_pathway <- kegg_pathway %>% 
#   dplyr::mutate(mapped_percentage = n_leading_edge / metabolite_all_number)

kegg_pathway <- 
  kegg_pathway %>% 
  dplyr::rename(mapped_number = n_leading_edge,
                p_value_adjust = fdr,
                mapped_percentage = match_ratio,
                all_number = metabolite_all_number)

load("2_data/3_Xia_Case_Study_v2/result_kegg/result_kegg.rda")
result_kegg@result <- kegg_pathway

p <- 
enrich_scatter_plot(
  object     = result_kegg,
  x_axis     = "mapped_percentage",
  point_size = "all_number"
)
p

ggsave(p,
       filename = "4_manuscript/Figures/Figure_2/hsa04978.pdf",
       width = 5,
       height = 4)



load("3_data_analysis/3_Xia_case_study/ranking_table_kegg_snr.rda")
load("3_data_analysis/3_Xia_case_study/original_score_annotation_kegg.rda")
load("2_data/3_Xia_Case_Study_v2/pathway_ratio/pathway_database.rda")

kegg_database <- pathway_database
ranking_table_new <- ranking_table

pathway_database <- metpath::kegg_hsa_pathway


# pathway_database <- imetpd_pathway

results <- perform_fmsea_analysis(
  pathway_database,
  annotation_table,         
  ranking_table = ranking_table_new,
  threads = 10,
  min.compounds.num = 5,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  perm.num = 10000,
  seed = 123,
  fdr.thr = 0.05,
  max.iter.num = 1,
  verbose = TRUE
)

kegg_pathway <- data.frame(
  pathway_id = names(results@res_list),
  p_value = sapply(results@res_list, function(x) x$p_value[1]),
  n_leading_edge = sapply(results@res_list, function(x) length(unique(x$leading_edge$KEGG_ID))),
  row.names = NULL
)

kegg_pathway$fdr <- p.adjust(kegg_pathway$p_value, method = "BH")

kegg_pathway <- kegg_pathway %>%
  left_join(kegg_database %>% dplyr::select(pathway_id, pathway_name, KEGG_ID, match_ratio), by = "pathway_id")

kegg_pathway <- kegg_pathway %>% 
  dplyr::mutate(metabolite_all_number = str_count(KEGG_ID, "\\{\\}") + 1)

# kegg_pathway <- kegg_pathway %>% 
#   dplyr::mutate(mapped_percentage = n_leading_edge / metabolite_all_number)

kegg_pathway <- 
  kegg_pathway %>% 
  dplyr::rename(mapped_number = n_leading_edge,
                p_value_adjust = fdr,
                mapped_percentage = match_ratio,
                all_number = metabolite_all_number)

load("2_data/3_Xia_Case_Study_v2/result_kegg/result_kegg.rda")
result_kegg@result <- kegg_pathway

p <- 
  enrich_scatter_plot(
    object     = result_kegg,
    x_axis     = "mapped_percentage",
    point_size = "all_number"
  )
p

ggsave(p,
       filename = "4_manuscript/Figures/Figure_2/original.pdf",
       width = 5,
       height = 4)



