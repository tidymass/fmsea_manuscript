library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(dplyr)

load("3_data_analysis/5_MTBLS749/results_iterations.rda")

significant_modules <- results_iterations@significant_modules

# ---- leading-edge features of each significant module ----
leading_edge_features <-
  lapply(significant_modules$pathway_id, function(id) {
    results_iterations@res_list[[id]]$leading_edge %>%
      mutate(pathway_id = id)
  }) %>%
  bind_rows() %>%
  left_join(significant_modules %>% dplyr::select(pathway_id, pathway_name),
            by = "pathway_id") %>%
  dplyr::select(pathway_id, pathway_name, everything())
rownames(leading_edge_features) <- NULL

# ---- number of leading-edge features per module ----
leading_edge_summary <-
  leading_edge_features %>%
  group_by(pathway_id) %>%
  summarise(
    n_leading_edge_feature = n_distinct(variable_id),
    n_leading_edge_compound = n_distinct(KEGG_ID),
    n_pos = n_distinct(variable_id[grepl("_POS$", variable_id)]),
    n_neg = n_distinct(variable_id[grepl("_NEG$", variable_id)]),
    .groups = "drop"
  )

leading_edge_summary <-
  significant_modules %>%
  dplyr::select(pathway_id, pathway_name, NES, FDR) %>%
  left_join(leading_edge_summary, by = "pathway_id") %>%
  arrange(desc(n_leading_edge_feature))

print(as.data.frame(leading_edge_summary))

# feature list of a single module, e.g.
# leading_edge_features %>% dplyr::filter(pathway_id == "hsa00020")
