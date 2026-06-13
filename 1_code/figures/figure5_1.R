library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(tidyverse)
library(ggpubr)

load("3_data_analysis/5_MTBLS749/original_score_annotation_kegg_pre.rda")
load("3_data_analysis/5_MTBLS749/ranking_table_kegg_snr_pre.rda")
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
  fdr.thr = 0.1,
  max.iter.num = 3,
  verbose = TRUE
)

save(results, file = "3_data_analysis/5_MTBLS749/results_iteration1.rda")
save(results, file = "3_data_analysis/5_MTBLS749/results_iteration2.rda")
save(results, file = "3_data_analysis/5_MTBLS749/results_iteration3.rda")


load("2_data/6_Huntington/metabolite_annotation/object_pos2_mona_hmdb_mp_massbank.rda")
load("2_data/6_Huntington/metabolite_annotation/object_neg2_mona_hmdb_mp_massbank.rda")
load("2_data/2_KEGG_MS1/kegg_compound_ms1.rda")


annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table

annotation_table_filtering <- rbind(annotation_table_pos, annotation_table_neg)

metabolite_database <- kegg_compound_ms1@spectra.info %>% 
  dplyr::rename(MS2_KEGG.ID = KEGG.ID) %>% 
  dplyr::select(MS2_KEGG.ID, CAS.ID) %>% 
  dplyr::filter(!is.na(MS2_KEGG.ID) & !is.na(CAS.ID)) %>% 
  distinct()

metabolite_database$CAS.ID <- gsub("\\s.*", "", metabolite_database$CAS.ID)

annotation_table_filtering <- 
  annotation_table_filtering %>% 
  left_join(metabolite_database, by = "CAS.ID")

annotation_table_filtering$KEGG.ID <- sub("\\{\\}.*$", "", annotation_table_filtering$KEGG.ID)

annotation_table_filtering <- annotation_table_filtering %>%
  mutate(
    KEGG.ID = if_else(
      !is.na(KEGG.ID) &
        !is.na(MS2_KEGG.ID) &
        KEGG.ID != MS2_KEGG.ID,
      MS2_KEGG.ID,
      KEGG.ID
    )
  )

annotation_table_filtering <- annotation_table_filtering %>% 
  filter(!is.na(KEGG.ID)) %>% 
  dplyr::select(variable_id, KEGG.ID) %>% 
  distinct()

feature_metabolite_count <- results@feature_metabolite_count

variable_ids <- unique(feature_metabolite_count$variable_id)

annotation_table_filtering2 <- annotation_table_filtering %>% 
  filter(variable_id %in% variable_ids)

weighted_annotation_table <- results@annotation_table_weighting

long_table <- weighted_annotation_table %>%
  rownames_to_column("feature_id") %>%
  pivot_longer(
    cols = -feature_id,
    names_to = "KEGG.ID",
    values_to = "Score"
  ) %>%
  dplyr::filter(Score != 0)

# long_table_ranked <- long_table %>%
#   group_by(feature_id) %>%
#   mutate(Rank = dense_rank(desc(Score))) %>%
#   ungroup()

long_table_ranked <- long_table %>%
  group_by(feature_id) %>%
  mutate(Rank = min_rank(desc(Score))) %>%
  ungroup()
# 
# long_table_ranked <- long_table %>% 
#   group_by(feature_id) %>% 
#   mutate(Rank = rank(-Score, ties.method = "max")) %>% 
#   ungroup()

# Join to annotation_table_filtering2, keep matching rows only
annotation_table_filtering3 <- annotation_table_filtering2 %>%
  left_join(long_table_ranked, by = c("variable_id" = "feature_id", "KEGG.ID")) %>%
  group_by(variable_id) %>%
  slice_min(Rank, n = 1, with_ties = FALSE) %>%
  ungroup()

length(unique(annotation_table_filtering3$variable_id))

annotation_table_filtering3 %>% 
  filter(Rank == 3) %>% 
  pull(variable_id) %>% 
  unique() %>% 
  length()

save(annotation_table_filtering3, file = "3_data_analysis/5_MTBLS749/annotation_table_filtering3.rda")
save(long_table_ranked, file = "3_data_analysis/5_MTBLS749/long_table_ranked.rda")

# visualization
top3_accuracy <- read.csv("3_data_analysis/5_MTBLS749/top3_accuracy.csv", 
               check.names = FALSE,
               row.names = 1)

df <- top3_accuracy %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("iteration") %>%
  mutate(
    across(c(Features_total, Top_1, Top_2, Top_3), as.numeric)
  )

df_plot <- df %>%
  mutate(
    iteration = str_replace(iteration, "Iteration", "Iteration "),
    `Top 1` = Top_1 / Features_total * 100,
    `Top 2` = (Top_1 + Top_2) / Features_total * 100,
    `Top 3` = (Top_1 + Top_2 + Top_3) / Features_total * 100
  ) %>%
  filter(iteration == "Iteration 3") %>%
  dplyr::select(`Top 1`, `Top 2`, `Top 3`) %>%
  pivot_longer(
    cols = everything(),
    names_to = "metric",
    values_to = "accuracy"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Top 1", "Top 2", "Top 3")
    )
  )

p <- ggplot(df_plot, aes(x = metric, y = accuracy, fill = metric)) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = sprintf("%.1f%%", accuracy)),
    vjust = -0.5,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Top 1" = "#4C78A8",
      "Top 2" = "#4C78A8",
      "Top 3" = "#4C78A8"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Accuracy (%)",
    title = "Top-3 accuracy"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    legend.position = "none"
  )

p

ggsave("4_manuscript/Figures/Figure_5/top_3_accuracy.pdf", plot = p, width = 14, height = 6)



# top_count_df <- long_table %>%
#   dplyr::filter(feature_id %in% variable_ids) %>%
#   dplyr::group_by(feature_id) %>%
#   dplyr::mutate(
#     Rank = dense_rank(desc(Score))
#   ) %>%
#   dplyr::filter(Rank %in% c(1, 2, 3)) %>%
#   dplyr::count(feature_id, Rank, name = "num") %>%
#   dplyr::mutate(
#     Rank = paste0("top", Rank, "_num")
#   ) %>%
#   tidyr::pivot_wider(
#     names_from = Rank,
#     values_from = num,
#     values_fill = 0
#   ) %>%
#   dplyr::right_join(
#     tibble::tibble(feature_id = variable_ids),
#     by = "feature_id"
#   ) %>%
#   dplyr::mutate(
#     top1_num = tidyr::replace_na(top1_num, 0),
#     top2_num = tidyr::replace_na(top2_num, 0),
#     top3_num = tidyr::replace_na(top3_num, 0)
#   ) %>%
#   dplyr::select(feature_id, top1_num, top2_num, top3_num)

# top_count_df <- long_table %>%
#   dplyr::filter(feature_id %in% variable_ids) %>%
#   dplyr::group_by(feature_id) %>%
#   dplyr::mutate(
#     Rank = rank(-Score, ties.method = "min")
#   ) %>%
#   dplyr::filter(Rank %in% c(1, 2, 3)) %>%
#   dplyr::count(feature_id, Rank, name = "num") %>%
#   dplyr::mutate(
#     Rank = paste0("top", Rank, "_num")
#   ) %>%
#   tidyr::pivot_wider(
#     names_from = Rank,
#     values_from = num,
#     values_fill = 0
#   ) %>%
#   dplyr::right_join(
#     tibble::tibble(feature_id = variable_ids),
#     by = "feature_id"
#   ) %>%
#   dplyr::mutate(
#     top1_num = tidyr::replace_na(top1_num, 0),
#     top2_num = tidyr::replace_na(top2_num, 0),
#     top3_num = tidyr::replace_na(top3_num, 0)
#   ) %>%
#   dplyr::select(feature_id, top1_num, top2_num, top3_num)

# top_count_df <- long_table %>%
#   dplyr::filter(feature_id %in% variable_ids) %>%
#   dplyr::group_by(feature_id) %>%
#   dplyr::mutate(
#     Rank = rank(-Score, ties.method = "min")
#   ) %>%
#   dplyr::summarise(
#     top1_num = sum(Rank <= 1, na.rm = TRUE),
#     top2_num = sum(Rank <= 2, na.rm = TRUE),
#     top3_num = sum(Rank <= 3, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   dplyr::right_join(
#     tibble::tibble(feature_id = variable_ids),
#     by = "feature_id"
#   ) %>%
#   dplyr::mutate(
#     top1_num = tidyr::replace_na(top1_num, 0),
#     top2_num = tidyr::replace_na(top2_num, 0),
#     top3_num = tidyr::replace_na(top3_num, 0)
#   ) %>%
#   dplyr::select(feature_id, top1_num, top2_num, top3_num)
# 
# plot_df <- top_count_df %>%
#   tidyr::pivot_longer(
#     cols = c(top1_num, top2_num, top3_num),
#     names_to = "Top",
#     values_to = "Num"
#   ) %>%
#   dplyr::mutate(
#     Top = dplyr::case_when(
#       Top == "top1_num" ~ "Top 1",
#       Top == "top2_num" ~ "Top 2",
#       Top == "top3_num" ~ "Top 3"
#     ),
#     Top = factor(Top, levels = c("Top 1", "Top 2", "Top 3"))
#   )
# 
# # summary_df <- plot_df %>%
# #   group_by(Top) %>%
# #   summarise(
# #     mean_num = mean(Num, na.rm = TRUE),
# #     median_num = median(Num, na.rm = TRUE),
# #     .groups = "drop"
# #   )
# 
# summary_df <- plot_df %>%
#   dplyr::group_by(Top) %>%
#   dplyr::summarise(
#     mean_num = mean(Num, na.rm = TRUE),
#     median_num = median(Num, na.rm = TRUE),
#     sd_num = sd(Num, na.rm = TRUE),
#     se_num = sd_num / sqrt(dplyr::n()),
#     .groups = "drop"
#   ) %>%
#   dplyr::mutate(
#     dplyr::across(
#       .cols = c(mean_num, median_num, sd_num, se_num),
#       .fns = ~ round(.x, 2)
#     )
#   )
# 
# p <- ggplot(plot_df, aes(x = Top, y = Num, fill = Top)) +
#   geom_boxplot(
#     width = 0.45,
#     outlier.shape = NA,
#     alpha = 0.75,
#     color = "grey25",
#     linewidth = 0.6
#   ) +
#   stat_summary(
#     fun = mean,
#     geom = "point",
#     shape = 23,
#     size = 3.2,
#     fill = "white",
#     color = "black"
#   ) +
#   geom_text(
#     data = summary_df,
#     aes(
#       x = Top,
#       y = 14.5,
#       label = paste0("Mean = ", round(mean_num, 2))
#     ),
#     inherit.aes = FALSE,
#     size = 4
#   ) +
#   coord_cartesian(ylim = c(0, 15)) +
#   scale_y_continuous(
#     breaks = seq(0, 15, 3)
#   ) +
#   scale_fill_manual(
#     values = c(
#       "Top 1" = "#4C72B0",
#       "Top 2" = "#55A868",
#       "Top 3" = "#C44E52"
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Number of annotations",
#     title = "Top-ranked annotation count per feature"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     axis.text = element_text(color = "black"),
#     axis.title.y = element_text(color = "black")
#   )
# p


calculate_top_count <- function(input_long_table, group_name, variable_ids) {
  
  input_long_table %>%
    dplyr::filter(feature_id %in% variable_ids) %>%
    dplyr::group_by(feature_id) %>%
    dplyr::mutate(
      Rank = rank(-Score, ties.method = "min")
    ) %>%
    dplyr::summarise(
      top1_num = sum(Rank <= 1, na.rm = TRUE),
      top2_num = sum(Rank <= 2, na.rm = TRUE),
      top3_num = sum(Rank <= 3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::right_join(
      tibble::tibble(feature_id = variable_ids),
      by = "feature_id"
    ) %>%
    dplyr::mutate(
      top1_num = tidyr::replace_na(top1_num, 0),
      top2_num = tidyr::replace_na(top2_num, 0),
      top3_num = tidyr::replace_na(top3_num, 0),
      Group = group_name
    )
}

weighted_annotation_table <- results@annotation_table_weighting

long_table <- weighted_annotation_table %>%
  rownames_to_column("feature_id") %>%
  pivot_longer(
    cols = -feature_id,
    names_to = "KEGG.ID",
    values_to = "Score"
  ) %>%
  dplyr::filter(Score != 0)

original_long_table <- original_score_annotation %>%
  rownames_to_column("feature_id") %>%
  pivot_longer(
    cols = -feature_id,
    names_to = "KEGG.ID",
    values_to = "Score"
  ) %>%
  dplyr::filter(Score != 0)

top_count_new <- calculate_top_count(
  input_long_table = long_table,
  group_name = "featureMSEA",
  variable_ids = variable_ids
)

top_count_original <- calculate_top_count(
  input_long_table = original_long_table,
  group_name = "Original",
  variable_ids = variable_ids
)

plot_df <- dplyr::bind_rows(
  top_count_new,
  top_count_original
) %>%
  tidyr::pivot_longer(
    cols = c(top1_num, top2_num, top3_num),
    names_to = "Top",
    values_to = "Num"
  ) %>%
  dplyr::mutate(
    Top = dplyr::case_when(
      Top == "top1_num" ~ "Top 1",
      Top == "top2_num" ~ "Top 2",
      Top == "top3_num" ~ "Top 3"
    ),
    Top = factor(Top, levels = c("Top 1", "Top 2", "Top 3")),
    Group = factor(Group, levels = c("Original", "featureMSEA"))
  )

# p <- ggplot(plot_df, aes(x = Top, y = Num, fill = Group)) +
#   geom_boxplot(
#     position = position_dodge(width = 0.8),
#     width = 0.6,
#     outlier.shape = NA,
#     alpha = 0.75
#   ) +
#   stat_compare_means(
#     aes(group = Group),
#     method = "wilcox.test",
#     label = "p.signif",
#     hide.ns = FALSE
#   ) +
#   scale_y_continuous(
#     limits = c(0, 15),
#     breaks = seq(0, 15, 3)
#   ) +
#   scale_fill_manual(
#     values = c(
#       "Original" = "#4E79A7",
#       "featureMSEA" = "#F28E2B"
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Number of annotations",
#     fill = NULL
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     axis.text.x = element_text(size = 13, color = "black"),
#     axis.text.y = element_text(size = 12, color = "black"),
#     axis.title.y = element_text(size = 14, color = "black"),
#     legend.position = "top",
#     legend.text = element_text(size = 12)
#   )

p <- ggplot(plot_df, aes(x = Top, y = Num, fill = Group)) +
  # 1. Draw half-violin plot (Slightly reduced nudge to narrow the gap)
  geom_half_violin(
    position = position_dodge(width = 0.8),
    side = "r",           # Violin bulge facing right
    nudge = 0.07,         # Fine-tuned from 0.1 to 0.07 to bring it closer
    alpha = 0.4,          # Soft transparency
    color = "#4A4A4A",    # Light dark grey border
    trim = TRUE,          # Restrict to actual data range
    scale = "width"       # Keep maximum width consistent across groups
  ) +
  # 2. Nested narrow boxplot (Slightly widened to balance the layout)
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.18,         # Adjusted from 0.15 to 0.18 for better visual density
    outlier.shape = NA,   # Hide outliers since they are already shown by the violin
    alpha = 0.7,          # Deeper color for visual hierarchy
    color = "#333333",
    size = 0.6
  ) +
  # 3. Journal-level precise statistical significance test
  stat_compare_means(
    aes(group = Group),
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
      symbols = c("***", "**", "*", "ns")
    ),
    bracket.size = 0.4,   # Explicitly draw the comparison brackets
    color = "#333333",
    label.y = 17.5        # Unified height for neat alignment
  ) +
  # 4. Optimize y-axis to provide enough breathing room for top brackets
  scale_y_continuous(
    limits = c(0, 20),    # Raise upper limit to 20 to avoid compression
    breaks = seq(0, 18, 3),
    expand = expansion(mult = c(0, 0.02)) # No padding at the bottom, slight padding at the top
  ) +
  # 5. Academic Morandi low-saturation color palette (Nature style)
  scale_fill_manual(
    values = c(
      "Original" = "#6A8CAF",     # Muted blue
      "featureMSEA" = "#E09F67"   # Sunset orange
    )
  ) +
  # 6. Refined labels and typography theme
  labs(
    x = NULL,
    y = "Number of annotations",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    # Axis text and title optimization
    axis.text.x = element_text(size = 13, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, color = "black", margin = margin(r = 12)),
    # Minimalist and elegant legend design
    legend.position = "top", 
    legend.text = element_text(size = 12, face = "bold"),
    legend.key.spacing.x = unit(0.6, "cm"),
    # Adjust plot margin to prevent elements from being cut off
    plot.margin = margin(15, 20, 15, 15, "pt")
  )

ggsave(
  filename = "4_manuscript/Figures/Figure_5/top3_redundancy_2.pdf",
  plot = p,
  width = 6,
  height = 4
)
