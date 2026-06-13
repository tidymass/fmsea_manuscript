library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)
library(featuremsea)
library(tidyverse)
library(eulerr)

# Selection of target metabolite sets
load("2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_pos2.rda")
load("2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_neg2.rda")

load("3_data_analysis/3_Xia_case_study/ranking_table_kegg_snr.rda")

annotation_pos <- object_pos2@annotation_table
annotation_neg <- object_neg2@annotation_table

annotation_original <- rbind(annotation_pos, annotation_neg)

annotation_table <- 
  rbind(annotation_pos, annotation_neg) %>% 
  dplyr::select(variable_id, Compound.name, KEGG.ID) %>% 
  dplyr::filter(!is.na(KEGG.ID)) %>% 
  distinct()

kegg_database <- metpath::kegg_hsa_pathway
pathway_database <- convert_kegg2fmsea(kegg_database)

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


metabolite_ids <- pathway_database %>%
  dplyr::filter(pathway_id == "hsa04974") %>%
  pull(KEGG_ID) %>%
  str_split("\\{\\}") %>%
  unlist %>%
  unique()

hsa04974_feature_ids <- annotation_table %>%
  dplyr::filter(KEGG.ID %in% metabolite_ids) %>%
  dplyr::select(variable_id, KEGG.ID) %>%
  distinct() %>%
  pull(variable_id) %>%
  unique()

hsa04974_metabolite_ids <- annotation_table %>%
  dplyr::filter(KEGG.ID %in% metabolite_ids) %>%
  dplyr::select(variable_id, KEGG.ID) %>%
  distinct() %>%
  pull(KEGG.ID) %>%
  unique()

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

results@significant_modules

# hsa04974
target_id <- "hsa04974"

p <- plot_fmsea_plot(
  fmsea_obj            = results,
  pathway_id           = "hsa04974",
  title                = NULL,
  highlight_feature_ids = hsa04974_feature_ids
)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/hsa04974/hsa04974_enrichment_result_highlight.pdf",
       width = 6,
       height = 4)

# hsa00970 
target_id <- "hsa00970"

p <- plot_fmsea_plot(
  fmsea_obj            = results,
  pathway_id           = "hsa00970",
  title                = NULL,
  highlight_feature_ids = hsa04974_feature_ids
)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/hsa04974/hsa00970_enrichment_result_highlight.pdf",
       width = 6,
       height = 4)

# hsa04978
target_id <- "hsa04978"

p <- plot_fmsea_plot(
  fmsea_obj            = results,
  pathway_id           = "hsa04978",
  title                = NULL,
  highlight_feature_ids = hsa04974_feature_ids
)

ggsave(p,
       filename = "4_manuscript/Supplementary_figures/hsa04974/hsa04978_enrichment_result_highlight.pdf",
       width = 6,
       height = 4)



kegg_database <- convert_kegg2fmsea(kegg_database)




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

# 1. Color Palette Definition (Consistent with previous colors)
# ==============================================================================
color_04974 <- "#4E79A7"  # Muted Blue
color_00970 <- "#F28E2B"  # Muted Orange
color_04978 <- "#E15759"  # Muted Red

# ==============================================================================
# 2. Venn/Euler Diagram 1: hsa04974 vs hsa00970 (Area Proportional)
# ==============================================================================
list1 <- list(hsa04974 = hsa04974_ids, hsa00970 = hsa00970_ids)

# Fit the Euler model to calculate exact proportional areas
fit1 <- euler(list1)

# Customize the plot with academic publication style
p1 <- plot(fit1,
           quantities = list(type = "counts", cex = 1.2, fontfamily = "sans"), # Show compound counts inside
           fills = list(fill = c(color_04974, color_00970), alpha = 0.6),      # Apply designated colors with transparency
           edges = list(col = "white", lex = 2),                               # Clean white borders
           labels = list(cex = 1.2, fontface = "bold", fontfamily = "sans"),   # Category text style
           main = list(label = "Proportional Overlap: hsa04974 vs hsa00970", 
                       cex = 1.3, fontface = "bold", fontfamily = "sans"))


pdf("4_manuscript/Supplementary_figures/hsa04974/venn1.pdf", width = 6, height = 6)
print(p1)
dev.off()

# 3. Venn/Euler Diagram 2: hsa04974 vs hsa04978 (Area Proportional)
# ==============================================================================
list2 <- list(hsa04974 = hsa04974_ids, hsa04978 = hsa04978_ids)

# Fit the Euler model
fit2 <- euler(list2)

# Customize the plot
p2 <- plot(fit2,
           quantities = list(type = "counts", cex = 1.2, fontfamily = "sans"),
           fills = list(fill = c(color_04974, color_04978), alpha = 0.6),      # hsa04974 is still blue, hsa04978 is red
           edges = list(col = "white", lex = 2),
           labels = list(cex = 1.2, fontface = "bold", fontfamily = "sans"),
           main = list(label = "Proportional Overlap: hsa04974 vs hsa04978", 
                       cex = 1.3, fontface = "bold", fontfamily = "sans"))

pdf("4_manuscript/Supplementary_figures/hsa04974/venn2.pdf", width = 6, height = 6)
print(p2)
dev.off()

length(intersect(hsa04974_metabolite_ids, hsa00970_ids))
length(intersect(hsa04974_metabolite_ids, hsa04978_ids))



# Input data
plot_df <- data.frame(
  Pathway = c("Protein digestion and absorption",
              "Mineral absorption",
              "Aminoacyl-tRNA biosynthesis"),
  Total_metabolites = c(47, 29, 52),
  Overlap_with_adjusted_A = c(19, 9, 15)
)

# Calculate metabolites not overlapping with the 19 adjusted metabolites in pathway A
plot_df <- plot_df %>%
  mutate(
    Not_overlap = Total_metabolites - Overlap_with_adjusted_A
  )

# Convert to long format
plot_long <- plot_df %>%
  pivot_longer(
    cols = c(Not_overlap, Overlap_with_adjusted_A),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c("Not_overlap", "Overlap_with_adjusted_A"),
      labels = c("Other metabolites", "Metabolites overlapping with the 19 adjusted metabolites in pathway A")
    ),
    Pathway = factor(
      Pathway,
      levels = c("Protein digestion and absorption",
                 "Mineral absorption",
                 "Aminoacyl-tRNA biosynthesis")
    )
  )

# Plot
p <- ggplot(plot_long, aes(x = Pathway, y = Count, fill = Category)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "Other metabolites" = "#BFC7C9",
      "Metabolites overlapping with the 19 adjusted metabolites in pathway A" = "#D55E00"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of metabolites",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10)
  )

p

# Save figure
ggsave(
  filename = "4_manuscript/Supplementary_figures/hsa04974/overlap.pdf",
  plot = p,
  width = 7,
  height = 4.5,
  dpi = 300
)


