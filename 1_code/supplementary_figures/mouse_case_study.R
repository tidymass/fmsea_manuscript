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


# visualization
results <- results_original
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
  filename = "4_manuscript/Supplementary_figures/HFD_mice/pathways_original.pdf",
  plot     = p,
  width    = 10,
  height   = 6)


original_ids <- results_original@significant_modules %>% pull(pathway_id) %>% unique()
iterations_ids <- results_iterations@significant_modules %>% pull(pathway_id) %>% unique()
different_ids <- setdiff(iterations_ids, original_ids)

results <- results_iterations
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
  filename = "4_manuscript/Supplementary_figures/HFD_mice/pathways_iteration.pdf",
  plot     = p,
  width    = 10,
  height   = 6)

plasma_res <- analyze_matrix_relevance(
  results = results_iterations,
  sample_source = "plasma",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  provider = "openai",
  model = "gpt-5.5"
)

plasma_significant_modules <- plasma_res@significant_modules
save(plasma_significant_modules, file = "3_data_analysis/8_mice_diet/plasma_significant_modules.rda", compress = "xz")

res <- analyze_topic_relevance(results_iterations,
                               research_topic = "high-fat diet",
                               provider = "openai",
                               api_key = Sys.getenv("OPENAI_API_KEY"))
research_relevance_results <- res@significant_modules
save(research_relevance_results, file = "3_data_analysis/8_mice_diet/research_relevance_results.rda", compress = "xz")


target_id <- "mmu00350"

p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/HFD_mice/mmu00350.pdf",
       width = 5,
       height = 5)



target_id <- "mmu00640"

p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/HFD_mice/mmu00640.pdf",
       width = 5,
       height = 5)

target_id <- "mmu00564"

p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/HFD_mice/mmu00564.pdf",
       width = 5,
       height = 5)

target_id <- "mmu00650"

p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/HFD_mice/mmu00650.pdf",
       width = 5,
       height = 5)

