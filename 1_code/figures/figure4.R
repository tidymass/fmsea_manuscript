library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(tidymass)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(patchwork)
library(purrr)

load("2_data/6_Huntington/data_cleaning/POS/object_pos2.rda")
load("2_data/6_Huntington/data_cleaning/NEG/object_neg2.rda")

# =========================
# Helper functions
# =========================

compute_log2fc <- function(expr_mat, case_ids, control_ids) {
  rowMeans(expr_mat[, case_ids, drop = FALSE], na.rm = TRUE) -
    rowMeans(expr_mat[, control_ids, drop = FALSE], na.rm = TRUE)
}

compute_pvalue <- function(expr_mat, case_ids, control_ids) {
  apply(expr_mat, 1, function(x) {
    cv <- x[case_ids];    cv <- cv[!is.na(cv)]
    tv <- x[control_ids]; tv <- tv[!is.na(tv)]
    if (length(cv) < 2 || length(tv) < 2) return(NA_real_)
    tryCatch(t.test(cv, tv)$p.value, error = function(e) NA_real_)
  })
}

# =========================
# POS
# =========================

case_ids_pos <- object_pos2 %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::filter(group == "Case") %>% dplyr::pull(sample_id)

ctrl_ids_pos <- object_pos2 %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::filter(group == "Control") %>% dplyr::pull(sample_id)

expr_pos <- log2(as.matrix(object_pos2@expression_data) + 1)
case_ids_pos <- intersect(case_ids_pos, colnames(expr_pos))
ctrl_ids_pos <- intersect(ctrl_ids_pos, colnames(expr_pos))

feature_table_pos <- object_pos2@variable_info %>%
  mutate(log2fc = compute_log2fc(expr_pos, case_ids_pos, ctrl_ids_pos),
         pvalue = compute_pvalue(expr_pos, case_ids_pos, ctrl_ids_pos),
         polarity = "positive")

# =========================
# NEG
# =========================

case_ids_neg <- object_neg2 %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::filter(group == "Case") %>% dplyr::pull(sample_id)

ctrl_ids_neg <- object_neg2 %>%
  activate_mass_dataset("sample_info") %>%
  dplyr::filter(group == "Control") %>% dplyr::pull(sample_id)

expr_neg <- log2(as.matrix(object_neg2@expression_data) + 1)
case_ids_neg <- intersect(case_ids_neg, colnames(expr_neg))
ctrl_ids_neg <- intersect(ctrl_ids_neg, colnames(expr_neg))

feature_table_neg <- object_neg2@variable_info %>%
  mutate(log2fc = compute_log2fc(expr_neg, case_ids_neg, ctrl_ids_neg),
         pvalue = compute_pvalue(expr_neg, case_ids_neg, ctrl_ids_neg),
         polarity = "negative")

# =========================
# Combine + FDR + Plot
# =========================

feature_table <-
  bind_rows(feature_table_pos, feature_table_neg) %>%
  mutate(
    FDR           = p.adjust(pvalue, method = "BH"),
    neg_log10_FDR = -log10(FDR),
    significance  = case_when(
      FDR < 0.05 & log2fc > 0 ~ "Up",
      FDR < 0.05 & log2fc < 0 ~ "Down",
      TRUE                     ~ "NS"
    )
  )

n_up   <- sum(feature_table$significance == "Up",   na.rm = TRUE)
n_down <- sum(feature_table$significance == "Down", na.rm = TRUE)
n_total <- nrow(feature_table)

up_ids <- feature_table %>% filter(significance == "Up") %>% pull(variable_id) %>% unique()
down_ids <- feature_table %>% filter(significance == "Down") %>% pull(variable_id) %>% unique()

p <- 
ggplot(feature_table, aes(x = log2fc, y = neg_log10_FDR, color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Up: ", n_up, " (", round(n_up/n_total*100, 1), "%)"),
           hjust = 1.1, vjust = 2, color = "#E64B35", size = 3.5) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("Down: ", n_down, " (", round(n_down/n_total*100, 1), "%)"),
           hjust = -0.1, vjust = 2, color = "#2166AC", size = 3.5) +
  scale_color_manual(values = c(Up = "#E64B35", Down = "#2166AC", NS = "grey70"), name = NULL) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = expression(log[2]~FC), y = expression(-log[10](FDR))) +
  theme_bw()

ggsave(filename = "4_manuscript/Figures/Figure_4/volcano_plot.pdf",
       plot     = p,
       width    = 5,
       height   = 5)

load("2_data/6_Huntington/metabolite_annotation/object_pos2_mona_hmdb_mp_massbank.rda")
load("2_data/6_Huntington/metabolite_annotation/object_neg2_mona_hmdb_mp_massbank.rda")

annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table

annotation_table <- rbind(annotation_table_pos, annotation_table_neg)

ms2_ids <- annotation_table %>% 
  pull(variable_id) %>% unique()

length(intersect(up_ids, ms2_ids))
length(intersect(down_ids, ms2_ids))

total_up <- length(up_ids)
total_down <- length(down_ids)

up_ms2 <- length(intersect(up_ids, ms2_ids))
down_ms2 <- length(intersect(down_ids, ms2_ids))

make_venn <- function(n_total, n_sub, fill_sub, label_sub) {
  r_big   <- 1
  r_small <- r_big * sqrt(n_sub / n_total)
  x_small <- -r_big + r_small  # left inscribed
  
  df <- data.frame(
    x     = c(0, x_small),
    y     = c(0, 0),
    r     = c(r_big, r_small),
    fill  = c("grey70", fill_sub)
  )
  
  ggplot() +
    geom_circle(data = df,
                aes(x0 = x, y0 = y, r = r, fill = fill),
                alpha = 0.7, color = "white", linewidth = 0.5) +
    scale_fill_identity() +
    coord_fixed() +
    theme_void()
}

p_up   <- make_venn(total_up,   up_ms2,   "#E64B35", "Up-MS2")
p_down <- make_venn(total_down, down_ms2, "#2166AC", "Down-MS2")

ggsave("4_manuscript/Figures/Figure_4/venn_up.pdf", plot = p_up, width = 8, height = 4)
ggsave("4_manuscript/Figures/Figure_4/venn_down.pdf", plot = p_down, width = 8, height = 4)

# original pathways
original_pathways <- read_csv("2_data/6_Huntington/paper_pathways.csv")

# Sort pathways by FDR ascending (lower FDR = higher position)
original_pathways$Pathway <- factor(
  original_pathways$Pathway,
  levels = original_pathways$Pathway[order(original_pathways$FDR, decreasing = TRUE)]
)

p <- 
ggplot(original_pathways, aes(x = -log10(FDR), y = Pathway, size = Impact)) +
  geom_point(color = "#1565C0", alpha = 0.8) +
  scale_size_continuous(
    name = "Impact",
    range = c(5, 10)
  ) +
  labs(
    title = "Pathway Impact Analysis",
    x = expression(-log[10](FDR)),
    y = "Pathway"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "4_manuscript/Figures/Figure_4/paper_pathways.pdf",
  plot     = p,
  width    = 10,
  height   = 10
)

# load("2_data/6_Huntington/metabolite_annotation/object_pos2.rda")
# load("2_data/6_Huntington/metabolite_annotation/object_neg2.rda")
# 
# # pos
# sample_info_pos <- object_pos2@sample_info
# expression_data_pos <- object_pos2@expression_data
# 
# sample_info_pos <- sample_info_pos[match(colnames(expression_data_pos), sample_info_pos$sample_id), ]
# 
# case_samples    <- sample_info_pos$sample_id[sample_info_pos$group == "Case"]
# control_samples <- sample_info_pos$sample_id[sample_info_pos$group == "Control"]
# 
# results <- apply(expression_data_pos, 1, function(x) {
#   case_vals    <- x[case_samples]
#   control_vals <- x[control_samples]
#   
#   log2fc <- log2(mean(case_vals, na.rm = TRUE) / mean(control_vals, na.rm = TRUE))
#   p_val  <- tryCatch(
#     t.test(case_vals, control_vals)$p.value,
#     error = function(e) NA_real_
#   )
#   
#   c(log2FoldChange = log2fc, p_value = p_val)
# })
# 
# result_pos <- data.frame(
#   feature_id     = rownames(expression_data_pos),
#   log2FoldChange = results["log2FoldChange", ],
#   p_value        = results["p_value", ],
#   row.names      = NULL
# )
# 
# # neg
# sample_info_neg <- object_neg2@sample_info
# expression_data_neg <- object_neg2@expression_data
# sample_info_neg <- sample_info_neg[match(colnames(expression_data_neg), sample_info_neg$sample_id), ]
# case_samples    <- sample_info_neg$sample_id[sample_info_neg$group == "Case"]
# control_samples <- sample_info_neg$sample_id[sample_info_neg$group == "Control"]
# results <- apply(expression_data_neg, 1, function(x) {
#   case_vals    <- x[case_samples]
#   control_vals <- x[control_samples]
#   
#   log2fc <- log2(mean(case_vals, na.rm = TRUE) / mean(control_vals, na.rm = TRUE))
#   p_val  <- tryCatch(
#     t.test(case_vals, control_vals)$p.value,
#     error = function(e) NA_real_
#   )
#   
#   c(log2FoldChange = log2fc, p_value = p_val)
# })
# 
# result_neg <- data.frame(
#   feature_id     = rownames(expression_data_neg),
#   log2FoldChange = results["log2FoldChange", ],
#   p_value        = results["p_value", ],
#   row.names      = NULL
# )
# 
# results <- rbind(result_pos, result_neg)
# 
# library(ggplot2)
# 
# results$p_adj <- p.adjust(results$p_value, method = "BH")
# 
# results$significance <- "Not Significant"
# results$significance[results$log2FoldChange >  1 & results$p_adj < 0.05] <- "Up"
# results$significance[results$log2FoldChange < -1 & results$p_adj < 0.05] <- "Down"
# 
# # Calculate proportions
# total    <- nrow(results)
# n_up     <- sum(results$significance == "Up")
# n_down   <- sum(results$significance == "Down")
# pct_up   <- round(n_up   / total * 100, 2)
# pct_down <- round(n_down / total * 100, 2)
# 
# label_up   <- paste0("Up: ",   n_up,   " (", pct_up,   "%)")
# label_down <- paste0("Down: ", n_down, " (", pct_down, "%)")
# 
# p <- 
#   ggplot(results, aes(x = log2FoldChange, y = -log10(p_adj), fill = significance)) +
#   geom_point(alpha = 0.6, size = 3, shape = 21, color = "black", stroke = 0.3) +
#   scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "Not Significant" = "grey60")) +
#   geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
#   annotate("text", x =  Inf, y = Inf, label = label_up,
#            hjust = 1.1, vjust = 2, color = "#E64B35", size = 4, fontface = "bold") +
#   annotate("text", x = -Inf, y = Inf, label = label_down,
#            hjust = -0.1, vjust = 2, color = "#4DBBD5", size = 4, fontface = "bold") +
#   labs(
#     title = "Volcano Plot",
#     x = expression(log[2]~"Fold Change"),
#     y = expression(-log[10]~"(FDR)"),
#     fill = ""
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     legend.position = "top"
#   )

load("3_data_analysis/5_MTBLS749/feature_table.rda")
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

# save(ranking_table, file = "3_data_analysis/5_MTBLS749/ranking_table_kegg_snr_pre.rda")
# save(original_score_annotation,
#      file = "3_data_analysis/5_MTBLS749/original_score_annotation_kegg_pre.rda",
#      compress = "xz",
#      compression_level = 9)

load("3_data_analysis/5_MTBLS749/ranking_table_kegg_snr_pre.rda")
load("3_data_analysis/5_MTBLS749/original_score_annotation_kegg_pre.rda")

annotation_table <- original_score_annotation

# for test
pathway_database <- metpath::kegg_hsa_pathway

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
  fdr.thr = 0.1,
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
  fdr.thr = 0.1,
  # max.iter.num = 1,
  verbose = TRUE
)

save(results_original, file = "3_data_analysis/5_MTBLS749/results_original.rda", compress = "xz")
save(results_iterations, file = "3_data_analysis/5_MTBLS749/results_iterations.rda", compress = "xz")

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
fdr_threshold <- -log10(0.1)

# Merge pathway name
all_modules <- merge(all_modules, significant_modules[, c("pathway_id", "pathway_name")],
                     by = "pathway_id", all.x = TRUE)


# p <- 
#   ggplot(all_modules, aes(x = NES, y = log10_FDR)) +
#   # Layer 1: points not in csf_significant_modules, or non-significant, fixed size
#   geom_point(
#     data = ~ subset(., !(pathway_id %in% csf_significant_modules$pathway_id & log10_FDR > fdr_threshold)),
#     aes(fill = log10_FDR > fdr_threshold),
#     alpha = 0.7, size = 2, shape = 21, color = "black"
#   ) +
#   # Layer 2: significant points with confidence score, size mapped to confidence score
#   geom_point(
#     data = ~ merge(
#       subset(., pathway_id %in% csf_significant_modules$pathway_id & log10_FDR > fdr_threshold),
#       csf_significant_modules[, c("pathway_id", "matrix_confidence_score")],
#       by = "pathway_id"
#     ),
#     aes(fill = log10_FDR > fdr_threshold,
#         size = matrix_confidence_score),
#     alpha = 0.7, shape = 21, color = "black"
#   ) +
#   scale_fill_manual(
#     values = c("TRUE" = "red", "FALSE" = "grey60"),
#     guide = "none"
#   ) +
#   scale_size_continuous(
#     name   = "Confidence Score",
#     breaks = c(0, 25, 50, 75, 100),
#     range  = c(1, 6)
#   ) +
#   geom_vline(xintercept = 0,
#              linetype = "dashed", color = "grey40", linewidth = 0.8) +
#   geom_hline(yintercept = fdr_threshold,
#              linetype = "dashed", color = "grey40", linewidth = 0.8) +
#   geom_text_repel(
#     data = subset(all_modules, log10_FDR > fdr_threshold),
#     aes(label = pathway_name),
#     size = 3,
#     color = "black",
#     max.overlaps = 20
#   ) +
#   scale_x_continuous(limits = c(0, NA)) +
#   labs(
#     x     = "NES",
#     y     = expression(-log[10](FDR)),
#     title = "Pathway Enrichment Analysis"
#   ) +
#   theme_bw()
# 
# ggsave(
#   filename = "4_manuscript/Figures/Figure_4/pathways_iterations.pdf",
#   plot     = p,
#   width    = 10,
#   height   = 6
# )



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
  filename = "4_manuscript/Figures/Figure_4/pathways_original.pdf",
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
fdr_threshold <- -log10(0.1)

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
  filename = "4_manuscript/Figures/Figure_4/pathways_iterations.pdf",
  plot     = p,
  width    = 10,
  height   = 6)



# urine_res <- analyze_matrix_relevance(
#   results = results_iterations,
#   sample_source = "urine",
#   api_key = Sys.getenv("OPENAI_API_KEY"),
#   provider = "openai",
#   model = "gpt-5.5"
# )
# 
# plasma_res <- analyze_matrix_relevance(
#   results = results_iterations,
#   sample_source = "plasma",
#   api_key = Sys.getenv("OPENAI_API_KEY"),
#   provider = "openai",
#   model = "gpt-5.5"
# )

# csf_res <- analyze_matrix_relevance(
#   results = results_iterations,
#   sample_source = "cerebrospinal fluid",
#   api_key = Sys.getenv("OPENAI_API_KEY"),
#   provider = "openai",
#   model = "gpt-5.5"
# )
# 
# csf_significant_modules <- csf_res@significant_modules
# 
# p <- 
#   ggplot(all_modules, aes(x = NES, y = log10_FDR)) +
#   geom_point(
#     data = ~ subset(., !(pathway_id %in% csf_significant_modules$pathway_id & log10_FDR > fdr_threshold)),
#     aes(fill = case_when(
#       pathway_id %in% different_ids ~ "diff",
#       log10_FDR > fdr_threshold    ~ "significant",
#       TRUE                          ~ "normal"
#     )),
#     alpha = 0.7, size = 2, shape = 21, color = "black"
#   ) +
#   geom_point(
#     data = ~ merge(subset(., pathway_id %in% csf_significant_modules$pathway_id & log10_FDR > fdr_threshold),
#                    csf_significant_modules[, c("pathway_id", "matrix_confidence_score")],
#                    by = "pathway_id"),
#     aes(fill = case_when(
#       pathway_id %in% different_ids ~ "diff",
#       TRUE                          ~ "significant"
#     ),
#     size = matrix_confidence_score),
#     alpha = 0.7, shape = 21, color = "black"
#   ) +
#   scale_fill_manual(
#     values = c(
#       "diff"        = "#27AE60",
#       "significant" = "red",
#       "normal"      = "grey60"
#     ),
#     guide = "none"
#   ) +
#   scale_size_continuous(
#     name   = "Confidence Score",
#     breaks = c(0, 25, 50, 75, 100),
#     range  = c(1, 6)
#   ) +
#   geom_vline(xintercept = 0,
#              linetype = "dashed", color = "grey40", linewidth = 0.8) +
#   geom_hline(yintercept = fdr_threshold,
#              linetype = "dashed", color = "grey40", linewidth = 0.8) +
#   geom_text_repel(
#     data = subset(all_modules, log10_FDR > fdr_threshold),
#     aes(label = pathway_name),
#     size = 3,
#     color = "black",
#     max.overlaps = 20
#   ) +
#   scale_x_continuous(limits = c(0, NA)) +
#   labs(
#     x     = "NES",
#     y     = expression(-log[10](FDR)),
#     title = "Pathway Enrichment Analysis"
#   ) +
#   theme_bw()
# urine_modules <- urine_res@significant_modules
# plasma_modules <- plasma_res@significant_modules
# csf_modules <- csf_res@significant_modules
# 
# # boxplot
# library(ggplot2)
# library(dplyr)
# 
# # Merge data from plasma, urine, and CSF modules
# plasma_data <- plasma_modules %>%
#   dplyr::select(pathway_name, matrix_confidence_score) %>%
#   dplyr::mutate(source = "Plasma", x_num = 1)
# 
# urine_data <- urine_modules %>%
#   dplyr::select(pathway_name, matrix_confidence_score) %>%
#   mutate(source = "Urine", x_num = 2)
# 
# csf_data <- csf_modules %>%
#   dplyr::select(pathway_name, matrix_confidence_score) %>%
#   mutate(source = "CSF", x_num = 3)
# 
# combined_data <- bind_rows(plasma_data, urine_data, csf_data) %>%
#   mutate(source = factor(source, levels = c("Plasma", "Urine", "CSF")))
# 
# # Plot
# p <- ggplot() +
#   # Boxplot layer
#   geom_boxplot(
#     data = combined_data,
#     aes(x = x_num, y = matrix_confidence_score,
#         fill = source, group = x_num),
#     width = 0.4,
#     outlier.shape = NA,
#     alpha = 0.6
#   ) +
#   # Jittered points for each group
#   geom_jitter(
#     data = combined_data,
#     aes(x = x_num, y = matrix_confidence_score,
#         color = source),
#     width = 0.08,
#     size = 2,
#     alpha = 0.9
#   ) +
#   # Set colors
#   scale_fill_manual(
#     values = c("Plasma" = "#C44E52", "Urine" = "#E6B800", "CSF" = "#2E86C1")
#   ) +
#   scale_color_manual(
#     values = c("Plasma" = "#C0392B", "Urine" = "#D4AC0D", "CSF" = "#1A5276")
#   ) +
#   # Map numeric x to named axis labels
#   scale_x_continuous(
#     breaks = c(1, 2, 3),
#     labels = c("Plasma", "Urine", "CSF"),
#     limits = c(0.5, 3.5)
#   ) +
#   coord_cartesian(clip = "off") +
#   labs(
#     x = NULL,
#     y = "Matrix Confidence Score",
#     title = "Pathway Matrix Confidence Score: Plasma vs Urine vs CSF"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     panel.grid.minor = element_blank(),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size = 13),
#     plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
#   )
# 
# 

# module visualization
load("3_data_analysis/5_MTBLS749/results_original.rda")
load("3_data_analysis/5_MTBLS749/results_iterations.rda")

target_id <- "hsa00350"

csf_res <- analyze_matrix_relevance(
  results = results_iterations,
  sample_source = "cerebrospinal fluid",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  provider = "openai",
  model = "gpt-5.5"
)

csf_significant_modules <- csf_res@significant_modules
save(csf_significant_modules, file = "3_data_analysis/5_MTBLS749/csf_significant_modules_iterations.rda", compress = "xz")

res <- analyze_topic_relevance(results_iterations,
                               research_topic = "Huntington's disease",
                               provider = "openai",
                               api_key = Sys.getenv("OPENAI_API_KEY"))
research_relevance_results <- res@significant_modules
save(research_relevance_results, file = "3_data_analysis/5_MTBLS749/research_relevance_results.rda", compress = "xz")



p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Figures/Figure_4/hsa00350.pdf",
       width = 5,
       height = 5)



target_id <- "hsa00350"

p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Figures/Figure_4/hsa00350.pdf",
       width = 5,
       height = 5)

target_id <- "hsa04080"

p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Figures/Figure_4/hsa04080.pdf",
       width = 5,
       height = 5)

target_id <- "hsa00190"

p <- plot_fmsea_plot(fmsea_obj = results_iterations, 
                     pathway_id = target_id, 
                     title = NULL)

ggsave(p,
       filename = "4_manuscript/Figures/Figure_4/hsa00190.pdf",
       width = 5,
       height = 5)
