library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(tidyverse)
library(tidymass)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(patchwork)
library(purrr)

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

original_results <- perform_fmsea_analysis(
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

iteration_results <- perform_fmsea_analysis(
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

save(original_results, file = "3_data_analysis/7_pregnancy_plasma/original_results.rda")
save(iteration_results, file = "3_data_analysis/7_pregnancy_plasma/iteration_results.rda")

# visualization
results <- original_results
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

all_modules <- all_modules %>% 
  dplyr::rename(pathway_id = id)


all_modules$log10_FDR <- -log10(all_modules$FDR)
fdr_threshold <- -log10(0.05)

# Merge pathway name
all_modules <- merge(all_modules, significant_modules[, c("pathway_id", "pathway_name")],
                     by = "pathway_id", all.x = TRUE)



p <- 
  ggplot(all_modules, aes(x = NES, y = log10_FDR)) +
  geom_point(aes(fill = log10_FDR > fdr_threshold), 
             alpha = 0.7, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60"), guide = "none") +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_hline(yintercept = fdr_threshold,
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_text_repel(
    data = subset(all_modules, log10_FDR > fdr_threshold),
    aes(label = pathway_name),
    size = 3,
    color = "black",
    max.overlaps = 20
  ) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(
    x = "NES",
    y = expression(-log[10](FDR)),
    title = "Pathway Enrichment Analysis"
  ) +
  theme_bw()

p

ggsave(
  filename = "4_manuscript/Supplementary_figures/pregnancy_plasma/pathways_original.pdf",
  plot     = p,
  width    = 10,
  height   = 6)


original_ids <- original_results@significant_modules %>% pull(pathway_id) %>% unique()
iterations_ids <- iteration_results@significant_modules %>% pull(pathway_id) %>% unique()
different_ids <- setdiff(iterations_ids, original_ids)

results <- iteration_results
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

all_modules <- all_modules %>% 
  dplyr::rename(pathway_id = id)


all_modules$log10_FDR <- -log10(all_modules$FDR)
fdr_threshold <- -log10(0.05)

# Merge pathway name
all_modules <- merge(all_modules, significant_modules[, c("pathway_id", "pathway_name")],
                     by = "pathway_id", all.x = TRUE)

p <- 
  ggplot(all_modules, aes(x = NES, y = log10_FDR)) +
  geom_point(aes(fill = case_when(
    pathway_id %in% different_ids ~ "diff",
    log10_FDR > fdr_threshold ~ "significant",
    TRUE ~ "normal"
  )), 
  alpha = 0.7, size = 3, shape = 21, color = "black") +
  scale_fill_manual(
    values = c(
      "diff"        = "#27AE60",
      "significant" = "red",
      "normal"      = "grey60"
    ),
    guide = "none"
  ) +
  geom_vline(xintercept = 0,
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_hline(yintercept = fdr_threshold,
             linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_text_repel(
    data = subset(all_modules, log10_FDR > fdr_threshold),
    aes(label = pathway_name),
    size = 3,
    color = "black",
    max.overlaps = 20
  ) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(
    x = "NES",
    y = expression(-log[10](FDR)),
    title = "Pathway Enrichment Analysis"
  ) +
  theme_bw()

p

ggsave(
  filename = "4_manuscript/Supplementary_figures/pregnancy_plasma/pathways_iteration.pdf",
  plot     = p,
  width    = 10,
  height   = 6)


# module visualization
load("3_data_analysis/7_pregnancy_plasma/original_results.rda")
load("3_data_analysis/7_pregnancy_plasma/iteration_results.rda")

plasma_res <- analyze_matrix_relevance(
  results = iteration_results,
  sample_source = "plasma",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  provider = "openai",
  model = "gpt-5.5"
)

plasma_significant_modules <- plasma_res@significant_modules
save(plasma_significant_modules, file = "3_data_analysis/7_pregnancy_plasma/plasma_significant_modules.rda", compress = "xz")

res <- analyze_topic_relevance(iteration_results,
                               research_topic = "Pregnancy",
                               provider = "openai",
                               api_key = Sys.getenv("OPENAI_API_KEY"))
research_relevance_results <- res@significant_modules
save(research_relevance_results, file = "3_data_analysis/7_pregnancy_plasma/research_relevance_results.rda", compress = "xz")


target_id <- "hsa00140"

p <- plot_fmsea_plot(fmsea_obj = iteration_results, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/pregnancy_plasma/hsa00140.pdf",
       width = 5,
       height = 5)



target_id <- "hsa04913"

p <- plot_fmsea_plot(fmsea_obj = iteration_results, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/pregnancy_plasma/hsa04913.pdf",
       width = 5,
       height = 5)

target_id <- "hsa00830"

p <- plot_fmsea_plot(fmsea_obj = iteration_results, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/pregnancy_plasma/hsa00830.pdf",
       width = 5,
       height = 5)

target_id <- "hsa00120"

p <- plot_fmsea_plot(fmsea_obj = iteration_results, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/pregnancy_plasma/hsa00120.pdf",
       width = 5,
       height = 5)


