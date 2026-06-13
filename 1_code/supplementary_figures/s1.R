library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)
library(featuremsea)
library(tidyverse)
library(ggplot2)
library(ggrepel)


# Selection of target metabolite sets
load("2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_pos2.rda")
load("2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_neg2.rda")

annotation_pos <- object_pos2@annotation_table
annotation_neg <- object_neg2@annotation_table

annotation_original <- rbind(annotation_pos, annotation_neg)

annotation_table <- 
  rbind(annotation_pos, annotation_neg) %>% 
  dplyr::select(variable_id, Compound.name, KEGG.ID) %>% 
  dplyr::filter(!is.na(KEGG.ID)) %>% 
  distinct()

KEGG_unique <- unique(annotation_table$KEGG.ID)

# obtain kegg pathway database

kegg_database <- metpath::kegg_hsa_pathway
pathway_database <- convert_kegg2fmsea(kegg_database)


# calculate the matched percentage for all kegg hsa pathways
pathway_database <- pathway_database %>%
  mutate(
    
    KEGG_list = str_split(KEGG_ID, "\\{\\}"),
    
    total_count = map_int(KEGG_list, length),
    
    matched_count = map_int(KEGG_list, ~ sum(.x %in% KEGG_unique)),
    
    matched_percentage = ifelse(total_count > 0, matched_count / total_count, 0)
  ) 

# find corresponding feature ids for target pathways
metabolite_ids <- pathway_database %>%
  dplyr::filter(pathway_id == "hsa05230") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  .[[1]] %>%
  .[. %in% KEGG_unique]

hsa05230_feature_ids <- annotation_table %>%
  dplyr::filter(KEGG.ID %in% metabolite_ids) %>%
  dplyr::select(variable_id, KEGG.ID) %>%
  distinct() %>%
  pull(variable_id) %>%
  unique()

metabolite_ids <- pathway_database %>%
  dplyr::filter(pathway_id == "hsa04978") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  .[[1]] %>%
  .[. %in% KEGG_unique]

hsa04978_feature_ids <- annotation_table %>%
  dplyr::filter(KEGG.ID %in% metabolite_ids) %>%
  dplyr::select(variable_id, KEGG.ID) %>%
  distinct() %>%
  pull(variable_id) %>%
  unique()

metabolite_ids <- pathway_database %>%
  dplyr::filter(pathway_id == "hsa00970") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  .[[1]] %>%
  .[. %in% KEGG_unique]

hsa00970_feature_ids <- annotation_table %>%
  dplyr::filter(KEGG.ID %in% metabolite_ids) %>%
  dplyr::select(variable_id, KEGG.ID) %>%
  distinct() %>%
  pull(variable_id) %>%
  unique()


metabolite_ids <- pathway_database %>%
  dplyr::filter(pathway_id == "hsa04974") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  .[[1]] %>%
  .[. %in% KEGG_unique]

hsa04974_feature_ids <- annotation_table %>%
  dplyr::filter(KEGG.ID %in% metabolite_ids) %>%
  dplyr::select(variable_id, KEGG.ID) %>%
  distinct() %>%
  pull(variable_id) %>%
  unique()
# perform featureMSEA and visualization

load("3_data_analysis/3_Xia_case_study/ranking_table_kegg_snr.rda")
load("3_data_analysis/3_Xia_case_study/original_score_annotation_kegg.rda")

kegg_database <- pathway_database


# for hsa04974
set.seed(123) 
ranking_table_new <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% hsa04974_feature_ids, 
                                 ranking_weight + rnorm(n(), mean = 1, sd = 0.25), 
                                 ranking_weight))

load("3_data_analysis/3_Xia_case_study/original_score_annotation_kegg.rda")

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

significant_modules <- results@significant_modules

all_modules <- do.call(rbind, lapply(names(results@res_list), function(id) {
  data.frame(
    id      = id,
    p_value = results@res_list[[id]]$p_value,
    NES     = results@res_list[[id]]$NES,
    stringsAsFactors = FALSE
  )
}))

all_modules$FDR <- p.adjust(all_modules$p_value, method = "BH")

all_modules$neg_log10_FDR <- -log10(all_modules$FDR)

all_modules <- all_modules %>% 
  left_join(significant_modules %>% dplyr::select(pathway_id, pathway_name), by = c("id" = "pathway_id"))


sig <- subset(all_modules, FDR < 0.05)

p <- 
ggplot(all_modules, aes(x = NES, y = neg_log10_FDR)) +
  geom_point(aes(color = FDR < 0.05), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D7534E"),
                     labels = c("FDR ≥ 0.05", "FDR < 0.05"),
                     name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  geom_text_repel(data = sig,
                  aes(label = pathway_name),
                  size = 3,
                  max.overlaps = 30,
                  box.padding = 0.4,
                  segment.color = "grey60",
                  segment.size = 0.3) +
  labs(x = "NES", y = expression(-log[10](FDR))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.line = element_line(color = "black"))

ggsave(p, filename = "4_manuscript/Supplementary_figures/4_scatter_plots/4.pdf", width = 6, height = 5)

# for hsa00970
set.seed(123) 
ranking_table_new <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% hsa00970_feature_ids, 
                                 ranking_weight + rnorm(n(), mean = 1, sd = 0.25), 
                                 ranking_weight))

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

significant_modules <- results@significant_modules

all_modules <- do.call(rbind, lapply(names(results@res_list), function(id) {
  data.frame(
    id      = id,
    p_value = results@res_list[[id]]$p_value,
    NES     = results@res_list[[id]]$NES,
    stringsAsFactors = FALSE
  )
}))

all_modules$FDR <- p.adjust(all_modules$p_value, method = "BH")

all_modules$neg_log10_FDR <- -log10(all_modules$FDR)

all_modules <- all_modules %>% 
  left_join(significant_modules %>% dplyr::select(pathway_id, pathway_name), by = c("id" = "pathway_id"))


sig <- subset(all_modules, FDR < 0.05)

p <- 
  ggplot(all_modules, aes(x = NES, y = neg_log10_FDR)) +
  geom_point(aes(color = FDR < 0.05), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D7534E"),
                     labels = c("FDR ≥ 0.05", "FDR < 0.05"),
                     name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  geom_text_repel(data = sig,
                  aes(label = pathway_name),
                  size = 3,
                  max.overlaps = 30,
                  box.padding = 0.4,
                  segment.color = "grey60",
                  segment.size = 0.3) +
  labs(x = "NES", y = expression(-log[10](FDR))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.line = element_line(color = "black"))

ggsave(p, filename = "4_manuscript/Supplementary_figures/4_scatter_plots/1.pdf", width = 6, height = 5)

# for hsa05230
set.seed(123) 
ranking_table_new <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% hsa05230_feature_ids, 
                                 ranking_weight + rnorm(n(), mean = 1, sd = 0.25), 
                                 ranking_weight))

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

significant_modules <- results@significant_modules

all_modules <- do.call(rbind, lapply(names(results@res_list), function(id) {
  data.frame(
    id      = id,
    p_value = results@res_list[[id]]$p_value,
    NES     = results@res_list[[id]]$NES,
    stringsAsFactors = FALSE
  )
}))

all_modules$FDR <- p.adjust(all_modules$p_value, method = "BH")

all_modules$neg_log10_FDR <- -log10(all_modules$FDR)

all_modules <- all_modules %>% 
  left_join(significant_modules %>% dplyr::select(pathway_id, pathway_name), by = c("id" = "pathway_id"))


sig <- subset(all_modules, FDR < 0.05)

p <- 
  ggplot(all_modules, aes(x = NES, y = neg_log10_FDR)) +
  geom_point(aes(color = FDR < 0.05), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D7534E"),
                     labels = c("FDR ≥ 0.05", "FDR < 0.05"),
                     name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  geom_text_repel(data = sig,
                  aes(label = pathway_name),
                  size = 3,
                  max.overlaps = 30,
                  box.padding = 0.4,
                  segment.color = "grey60",
                  segment.size = 0.3) +
  labs(x = "NES", y = expression(-log[10](FDR))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.line = element_line(color = "black"))

ggsave(p, filename = "4_manuscript/Supplementary_figures/4_scatter_plots/2.pdf", width = 6, height = 5)

# for hsa04978
set.seed(123) 
ranking_table_new <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% hsa04978_feature_ids, 
                                 ranking_weight + rnorm(n(), mean = 1, sd = 0.25), 
                                 ranking_weight))

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

significant_modules <- results@significant_modules

all_modules <- do.call(rbind, lapply(names(results@res_list), function(id) {
  data.frame(
    id      = id,
    p_value = results@res_list[[id]]$p_value,
    NES     = results@res_list[[id]]$NES,
    stringsAsFactors = FALSE
  )
}))

all_modules$FDR <- p.adjust(all_modules$p_value, method = "BH")

all_modules$neg_log10_FDR <- -log10(all_modules$FDR)

all_modules <- all_modules %>% 
  left_join(significant_modules %>% dplyr::select(pathway_id, pathway_name), by = c("id" = "pathway_id"))


sig <- subset(all_modules, FDR < 0.05)

p <- 
  ggplot(all_modules, aes(x = NES, y = neg_log10_FDR)) +
  geom_point(aes(color = FDR < 0.05), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D7534E"),
                     labels = c("FDR ≥ 0.05", "FDR < 0.05"),
                     name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  geom_text_repel(data = sig,
                  aes(label = pathway_name),
                  size = 3,
                  max.overlaps = 30,
                  box.padding = 0.4,
                  segment.color = "grey60",
                  segment.size = 0.3) +
  labs(x = "NES", y = expression(-log[10](FDR))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.line = element_line(color = "black"))

ggsave(p, filename = "4_manuscript/Supplementary_figures/4_scatter_plots/3.pdf", width = 6, height = 5)
