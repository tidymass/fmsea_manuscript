library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)
library(featuremsea)
library(tidyverse)
library(ggrepel)

load("2_data/3_Xia_Case_Study_v2/data_cleaning/POS/object_pos2.rda")
load("2_data/3_Xia_Case_Study_v2/data_cleaning/NEG/object_neg2.rda")

# pos
sample_info_pos <- object_pos2@sample_info
expression_data_pos <- object_pos2@expression_data

sample_info_pos <- sample_info_pos[match(colnames(expression_data_pos), sample_info_pos$sample_id), ]

case_samples    <- sample_info_pos$sample_id[sample_info_pos$group == "Case"]
control_samples <- sample_info_pos$sample_id[sample_info_pos$group == "Control"]

results <- apply(expression_data_pos, 1, function(x) {
  case_vals    <- x[case_samples]
  control_vals <- x[control_samples]
  
  log2fc <- log2(mean(case_vals, na.rm = TRUE) / mean(control_vals, na.rm = TRUE))
  p_val  <- tryCatch(
    t.test(case_vals, control_vals)$p.value,
    error = function(e) NA_real_
  )
  
  c(log2FoldChange = log2fc, p_value = p_val)
})

result_pos <- data.frame(
  feature_id     = rownames(expression_data_pos),
  log2FoldChange = results["log2FoldChange", ],
  p_value        = results["p_value", ],
  row.names      = NULL
)

# neg
sample_info_neg <- object_neg2@sample_info
expression_data_neg <- object_neg2@expression_data
sample_info_neg <- sample_info_neg[match(colnames(expression_data_neg), sample_info_neg$sample_id), ]
case_samples    <- sample_info_neg$sample_id[sample_info_neg$group == "Case"]
control_samples <- sample_info_neg$sample_id[sample_info_neg$group == "Control"]
results <- apply(expression_data_neg, 1, function(x) {
  case_vals    <- x[case_samples]
  control_vals <- x[control_samples]
  
  log2fc <- log2(mean(case_vals, na.rm = TRUE) / mean(control_vals, na.rm = TRUE))
  p_val  <- tryCatch(
    t.test(case_vals, control_vals)$p.value,
    error = function(e) NA_real_
  )
  
  c(log2FoldChange = log2fc, p_value = p_val)
})

result_neg <- data.frame(
  feature_id     = rownames(expression_data_neg),
  log2FoldChange = results["log2FoldChange", ],
  p_value        = results["p_value", ],
  row.names      = NULL
)

results <- rbind(result_pos, result_neg)

results$p_adj <- p.adjust(results$p_value, method = "BH")

results$significance <- "Not Significant"
results$significance[results$log2FoldChange >  1 & results$p_adj < 0.05] <- "Up"
results$significance[results$log2FoldChange < -1 & results$p_adj < 0.05] <- "Down"

# Calculate proportions
total    <- nrow(results)
n_up     <- sum(results$significance == "Up")
n_down   <- sum(results$significance == "Down")
pct_up   <- round(n_up   / total * 100, 2)
pct_down <- round(n_down / total * 100, 2)

label_up   <- paste0("Up: ",   n_up,   " (", pct_up,   "%)")
label_down <- paste0("Down: ", n_down, " (", pct_down, "%)")

p <- 
  ggplot(results, aes(x = log2FoldChange, y = -log10(p_adj), fill = significance)) +
  geom_point(alpha = 0.6, size = 3, shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "Not Significant" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("text", x =  Inf, y = Inf, label = label_up,
           hjust = 1.1, vjust = 2, color = "#E64B35", size = 4, fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = label_down,
           hjust = -0.1, vjust = 2, color = "#4DBBD5", size = 4, fontface = "bold") +
  labs(
    title = "Volcano Plot",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"(FDR)"),
    fill = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

ggsave(
  filename = "4_manuscript/Figures/Figure_3/volcano.pdf",
  plot     = p,
  width    = 12,
  height   = 6)





# p <- 
# ggplot(results, aes(x = log2FoldChange, y = -log10(p_adj), color = significance)) +
#   geom_point(alpha = 0.6, size = 1.5) +
#   scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "Not Significant" = "grey60")) +
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
#     color = ""
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     legend.position = "top"
#   )


# Use original healthy cohort to perform featureMSEA
load("3_data_analysis/3_Xia_case_study/feature_table_snr.rda")
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

save(ranking_table, file = "3_data_analysis/3_Xia_case_study/ranking_table_kegg_snr.rda")
save(annotation_table, 
     file = "3_data_analysis/3_Xia_case_study/original_score_annotation_kegg.rda", 
     compress = "xz", 
     compression_level = 9)

pathway_database <- metpath::kegg_hsa_pathway


results <- perform_fmsea_analysis(
  pathway_database,
  annotation_table,         
  ranking_table = ranking_table,
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

p <- 
ggplot(all_modules, aes(x = id, y = neg_log10_FDR)) +
  geom_point(shape = 21, fill = "gray60", color = "black", size = 3) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "black", 
             linewidth = 0.8) +
  labs(x = "ID", 
       y = expression(-log[10](FDR)),
       title = "FDR Significance Plot") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "4_manuscript/Figures/Figure_3/scatter_pathways.pdf",
  plot     = p,
  width    = 12,
  height   = 6)



###
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

p <-
ggplot(pathway_database, aes(x = matched_percentage, y = total_count)) +
  geom_point(shape = 21, fill = "grey", color = "black", size = 2, alpha = 0.8) +
  geom_point(
    data = subset(pathway_database, matched_percentage > 0.25 & total_count > 15),
    shape = 21, fill = "#E74C3C", color = "black", size = 2.5
  ) +
  geom_text_repel(
    data = subset(pathway_database, matched_percentage > 0.25 & total_count > 15),
    aes(label = pathway_name),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "grey70", linewidth = 0.8) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "grey70", linewidth = 0.8) +
  labs(
    x = "Matched Percentage",
    y = "Total Number",
    title = "Pathway Distribution"
  ) +
  theme_bw()

ggsave(
  filename = "4_manuscript/Figures/Figure_3/kegg_pathway_distribution.pdf",
  plot     = p,
  width    = 8,
  height   = 6)

# highlight perturbed features
plot_fmsea_plot <- function(fmsea_obj,
                            pathway_id,
                            title = NULL,
                            highlight_feature_ids = NULL) {
  
  # --- 0. Configuration: Colors ---
  col_pos  <- "#E64B35"  # Red
  col_neg  <- "#3498DB"  # Blue
  col_line <- "#00A087"  # Teal
  
  # --- 1. Data Extraction ---
  if (!pathway_id %in% names(fmsea_obj@res_list)) {
    stop(paste("Pathway ID", pathway_id, "not found."))
  }
  
  # 1.1 Get Steps Data (X-axis source)
  step_data <- fmsea_obj@res_list[[pathway_id]]$steps
  
  # Preprocess P5/P6 columns: Replace NA with 0 immediately
  if("a_ij" %in% colnames(step_data)) step_data$a_ij[is.na(step_data$a_ij)] <- 0
  if("pair_weight" %in% colnames(step_data)) step_data$pair_weight[is.na(step_data$pair_weight)] <- 0
  
  # --- 2. Link with Ranking Table (for Direction) ---
  if (nrow(fmsea_obj@ranking_table) == 0) stop("ranking_table is empty.")
  
  # Reconstruct global order based on absolute weight
  rank_ref <- fmsea_obj@ranking_table %>%
    dplyr::arrange(dplyr::desc(abs(ranking_weight))) %>%
    dplyr::mutate(
      order_pos = dplyr::row_number(),
      real_ranking_weight = ranking_weight,
      norm_weight = abs(ranking_weight) / max(abs(ranking_weight))
    ) %>%
    dplyr::select(order_pos, real_ranking_weight, norm_weight)
  
  # Merge to get signed weights for colors
  plot_data <- step_data %>%
    dplyr::left_join(rank_ref, by = "order_pos")
  
  # --- 3. Sub-plot Data Prep ---
  
  # P2: Hits
  hits_data <- plot_data %>%
    dplyr::filter(!is.na(variable_id_hit) & variable_id_hit != "") %>%
    dplyr::mutate(
      hit_color = dplyr::case_when(
        !is.null(highlight_feature_ids) & variable_id_hit %in% highlight_feature_ids ~ "highlight",
        TRUE ~ "default"
      )
    )
  
  # P4: Global Weight Background
  min_pos <- min(plot_data$order_pos)
  max_pos <- max(plot_data$order_pos)
  weight_dist_data <- rank_ref %>%
    dplyr::filter(order_pos >= min_pos & order_pos <= max_pos)
  
  # --- 4. Stats & Table ---
  stats_df <- fmsea_obj@significant_modules
  id_col <- if("MFM_id" %in% colnames(stats_df)) "MFM_id" else "pathway_id"
  p_stats <- stats_df[stats_df[[id_col]] == pathway_id, ]
  
  # Helper to safe extract
  safe_val <- function(x, digits = 3, sci = FALSE) {
    if (length(x) == 0 || is.na(x)) return("-")
    if (is.numeric(x)) {
      if (sci) format(x, scientific = TRUE, digits = digits) else round(x, digits)
    } else as.character(x)
  }
  
  # Metadata
  p_name <- if (nrow(p_stats) > 0) (if ("pathway_name" %in% names(p_stats)) p_stats$pathway_name else p_stats$MFM_name) else pathway_id
  le_pos <- plot_data$order_pos[which.max(abs(plot_data$ES_cum))]
  
  # Table Data
  table_df <- data.frame(
    Metric = c("ID", "Name", "LE Pos", "P-val", "FDR", "ES", "NES"),
    Value = c(
      pathway_id,
      ifelse(nchar(p_name) > 25, paste0(substr(p_name, 1, 22), "..."), p_name),
      safe_val(le_pos, 0),
      safe_val(p_stats$p_value, sci = TRUE),
      safe_val(p_stats$FDR, sci = TRUE),
      safe_val(p_stats$ES, 4),
      safe_val(p_stats$NES, 4)
    ), stringsAsFactors = FALSE
  )
  
  # Table Grob
  tbl_grob <- gridExtra::tableGrob(
    table_df, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core = list(
        fg_params = list(hjust = 0, x = 0.1, fontsize = 5),
        bg_params = list(fill = c("white", "grey95"), col = "grey60", lwd = 0.2)
      ),
      colhead = list(
        fg_params = list(fontsize = 5, fontface = "bold"),
        bg_params = list(fill = "grey85", col = "grey40")
      )
    )
  )
  tbl_grob <- gtable::gtable_add_grob(
    tbl_grob,
    grid::rectGrob(gp = grid::gpar(fill = NA, col = "grey30")),
    t = 1, l = 1, b = nrow(tbl_grob), r = ncol(tbl_grob), z = 0
  )
  
  # --- 5. Visualization Panels ---
  x_lim <- c(min_pos, max_pos)
  
  common_x_scale <- ggplot2::scale_x_continuous(
    limits = x_lim,
    expand = c(0, 0)
  )
  
  theme_strip <- ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x  = ggplot2::element_blank(),
      axis.text.x   = ggplot2::element_blank(),
      axis.ticks.x  = ggplot2::element_blank(),
      panel.grid    = ggplot2::element_blank(),
      plot.margin   = ggplot2::margin(1, 5, 1, 1),
      axis.title.y  = ggplot2::element_text(size = 9, angle = 90, vjust = 0.5),
      axis.text.y   = ggplot2::element_text(size = 7)
    )
  
  # P1: ES Curve
  y_rng <- range(plot_data$ES_cum)
  tbl_y <- if (abs(max(y_rng)) >= abs(min(y_rng))) {
    c(max(y_rng) * 0.5, max(y_rng) * 0.95)
  } else {
    c(min(y_rng) * 0.1, min(y_rng) * 0.5)
  }
  
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = order_pos, y = ES_cum)) +
    ggplot2::geom_vline(xintercept = le_pos, linetype = "dashed", color = col_pos, alpha = 0.8) +
    ggplot2::geom_line(color = col_line, linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    ggplot2::annotation_custom(
      tbl_grob,
      xmin = x_lim[2] * 0.6, xmax = x_lim[2] * 1.02,
      ymin = tbl_y[1],        ymax = tbl_y[2]
    ) +
    ggplot2::labs(title = if (is.null(title)) p_name else title, y = "ES") +
    common_x_scale +
    theme_strip +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", size = 12),
      axis.title.y = ggplot2::element_text(size = 11)
    )
  
  # P2: Hits (grey by default, red for highlight_feature_ids)
  p2 <- ggplot2::ggplot(hits_data, ggplot2::aes(x = order_pos, y = 1)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = 0, ymax = 1, color = hit_color)
    ) +
    ggplot2::scale_color_manual(
      values = c("highlight" = col_pos, "default" = "grey70"),
      guide  = "none"
    ) +
    ggplot2::labs(y = "Hit") +
    common_x_scale +
    theme_strip +
    ggplot2::theme(
      axis.line.x  = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )
  
  # P3: Step
  p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = order_pos, y = step)) +
    ggplot2::geom_area(fill = "grey85", alpha = 0.5) +
    ggplot2::geom_line(color = "grey40", linewidth = 0.4) +
    ggplot2::labs(y = "Step") +
    common_x_scale +
    theme_strip
  
  # P4: Ranking Weight
  p4 <- ggplot2::ggplot(weight_dist_data, ggplot2::aes(x = order_pos, y = abs(real_ranking_weight))) +
    ggplot2::geom_line(color = "#C0392B", linewidth = 0.6) +
    ggplot2::labs(y = "Weight") +
    common_x_scale +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    theme_strip +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA))
  
  # P5: a_ij
  p5 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = order_pos, y = a_ij)) +
    ggplot2::geom_line(color = "#6C3483", linewidth = 0.6) +
    ggplot2::labs(y = "a_ij") +
    common_x_scale +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    theme_strip +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA))
  
  # P6: Pair Weight
  p6 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = order_pos, y = pair_weight)) +
    ggplot2::geom_line(color = "#BA4A00", linewidth = 0.6) +
    ggplot2::labs(x = "Ordered Features", y = "Pair W.") +
    common_x_scale +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    theme_strip +
    ggplot2::theme(
      plot.margin  = ggplot2::margin(1, 5, 5, 1),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.text.x  = ggplot2::element_text(size = 8),
      axis.ticks.x = ggplot2::element_line(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA)
    )
  
  # --- 6. Assembly ---
  combined_plot <- (p1 / p2 / p3 / p4 / p5 / p6) +
    patchwork::plot_layout(heights = c(3, 0.4, 1.2, 0.6, 0.6, 0.6))
  
  return(combined_plot)
}

###
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

# for hsa00970
set.seed(123) 
ranking_table_new <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% hsa00970_feature_ids, 
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

results@significant_modules

target_id <- "hsa00970"

p <- plot_fmsea_plot(
  fmsea_obj            = results,
  pathway_id           = "hsa00970",
  title                = NULL,
  highlight_feature_ids = hsa00970_feature_ids
)

ggsave(p,
       filename = "4_manuscript/Figures/Figure_3/hsa00970_enrichment_result_highlight.pdf",
       width = 6,
       height = 4)


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

results@significant_modules

target_id <- "hsa04978"

p <- plot_fmsea_plot(
  fmsea_obj            = results,
  pathway_id           = "hsa04978",
  title                = NULL,
  highlight_feature_ids = hsa04978_feature_ids
)

ggsave(p,
       filename = "4_manuscript/Figures/Figure_3/hsa04978_enrichment_result_highlight.pdf",
       width = 6,
       height = 4)


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

results@significant_modules

target_id <- "hsa05230"

p <- plot_fmsea_plot(
  fmsea_obj            = results,
  pathway_id           = "hsa05230",
  title                = NULL,
  highlight_feature_ids = hsa05230_feature_ids
)

ggsave(p,
       filename = "4_manuscript/Figures/Figure_3/hsa05230_enrichment_result_highlight.pdf",
       width = 6,
       height = 4)


# for hsa04974
set.seed(123) 
ranking_table_new <- ranking_table %>%
  mutate(ranking_weight = ifelse(variable_id %in% hsa04974_feature_ids, 
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

results@significant_modules

results@significant_modules[["pathway_id"]]

hsa04974_ids <- kegg_database %>% 
  dplyr::filter(pathway_id == "hsa04974") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  unlist() %>%
  unique()

hsa00970_ids <- kegg_database %>% 
  dplyr::filter(pathway_id == "hsa00970") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  unlist() %>%
  unique()

hsa04978_ids <- kegg_database %>% 
  dplyr::filter(pathway_id == "hsa04978") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  unlist() %>%
  unique()


length(intersect(hsa04974_ids, hsa00970_ids))
length(intersect(hsa04974_ids, hsa04978_ids))


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


## Heatmap

# Sample data (replace with actual data)
ranking_table <- data.frame(
  variable_id = paste0("V", 1:30),
  ranking_weight = runif(30, 0, 1)  # range set to 0~2
)

# Sort
ranking_table <- ranking_table[order(ranking_table$ranking_weight, decreasing = TRUE), ]
ranking_table$order <- seq_len(nrow(ranking_table))

p <-
ggplot(ranking_table, aes(x = 1, y = reorder(variable_id, ranking_weight), fill = ranking_weight)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",        # close to 0 → white
    high = "#08306b",         # close to 2 → dark blue
    limits = c(0, 2),     # fixed color scale range
    name = "Ranking\nWeight"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin  = margin(10, 10, 10, 10)
  ) +
  labs(title = "Variable Ranking Heatmap")

ggsave("4_manuscript/Figures/Figure_2/heatmap1.pdf", plot = p, width = 2, height = 6)



set.seed(42)

# Base data
ranking_table <- data.frame(
  variable_id = paste0("V", 1:30),
  ranking_weight = runif(30, 0, 1)
)

# Randomly select 10 points and add normal perturbation
selected_idx <- sample(1:30, 10)
ranking_table$ranking_weight[selected_idx] <-
  ranking_table$ranking_weight[selected_idx] + rnorm(10, mean = 0.7, sd = 0.1)

# Clip to range 0~2
ranking_table$ranking_weight <- pmin(pmax(ranking_table$ranking_weight, 0), 2)

# Sort
ranking_table <- ranking_table[order(ranking_table$ranking_weight, decreasing = TRUE), ]
ranking_table$order <- seq_len(nrow(ranking_table))

# Heatmap
p2 <-
ggplot(ranking_table, aes(x = 1, y = reorder(variable_id, ranking_weight), fill = ranking_weight)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",
    high = "#08306b",
    limits = c(0, 2),
    name = "Ranking\nWeight"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin  = margin(10, 10, 10, 10)
  ) +
  labs(title = "Variable Ranking Heatmap")

ggsave("4_manuscript/Figures/Figure_2/heatmap2.pdf", plot = p2, width = 2, height = 6)


### normal distribution line

# Normal distribution data
x <- seq(-4, 4, length.out = 1000)
y <- dnorm(x, mean = 0, sd = 1)
df <- data.frame(x = x, y = y)

p3 <- 
ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#4faed9", linewidth = 1.2) +
  geom_area(fill = "white", alpha = 0.5) +
  theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.title   = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = "Normal Distribution",
    x = "x",
    y = "Density"
  )

ggsave("4_manuscript/Figures/Figure_2/normal_distribution.pdf", plot = p3, width = 4, height = 4)
