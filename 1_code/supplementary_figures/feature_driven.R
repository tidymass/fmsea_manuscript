library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(dplyr)
library(tidyr)
library(ggplot2)

load("3_data_analysis/5_MTBLS749/results_iterations.rda")

significant_modules <- results_iterations@significant_modules

# ---- per-module feature counts and top-k contributions ----
# contribution of a leading-edge feature = its pair_weight, which is proportional
# to the increment it adds to the running enrichment score
module_summary <-
  lapply(significant_modules$pathway_id, function(id) {
    res <- results_iterations@res_list[[id]]
    leading_edge <- res$leading_edge
    hits <- res$steps %>% dplyr::filter(type == "pair_slot")
    w <- sort(leading_edge$pair_weight, decreasing = TRUE) / sum(leading_edge$pair_weight)
    data.frame(
      pathway_id = id,
      n_hit = n_distinct(hits$variable_id_hit),
      n_leading_edge = n_distinct(leading_edge$variable_id),
      top1 = w[1],
      top2 = sum(w[1:2]),
      top3 = sum(w[1:3])
    )
  }) %>%
  bind_rows() %>%
  left_join(significant_modules %>% dplyr::select(pathway_id, pathway_name),
            by = "pathway_id")

theme_feature_driven <-
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# ---- A: number of hit features and of leading-edge features ----
plot_feature_number <- function(n, y_label, fill, color) {
  ggplot(data.frame(n = n), aes(x = "", y = n)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6, fill = fill) +
    geom_jitter(position = position_jitter(width = 0.12, height = 0, seed = 123),
                size = 2, alpha = 0.8, color = color) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = y_label) +
    theme_feature_driven
}

p_hit <- plot_feature_number(module_summary$n_hit, "Number of hit features",
                             fill = "#9ECAE1", color = "#4292C6")
p_hit

ggsave("4_manuscript/Supplementary_figures/feature_driven/hit_feature_number.pdf",
       p_hit, width = 2.2, height = 4.5)

p_leading_edge <- plot_feature_number(module_summary$n_leading_edge,
                                      "Number of leading-edge features",
                                      fill = "#3182BD", color = "#08519C")
p_leading_edge

ggsave("4_manuscript/Supplementary_figures/feature_driven/leading_edge_feature_number.pdf",
       p_leading_edge, width = 2.2, height = 4.5)

# ---- B: cumulative contribution of the top 1, 2 and 3 leading-edge features ----
contribution_data <-
  module_summary %>%
  dplyr::select(pathway_name, top1, top2, top3) %>%
  pivot_longer(-pathway_name, names_to = "top", values_to = "proportion") %>%
  mutate(top = factor(top,
                      levels = c("top1", "top2", "top3"),
                      labels = c("Top 1", "Top 2", "Top 3")))

p_contribution <-
  ggplot(contribution_data, aes(x = top, y = proportion, fill = top)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = top),
              position = position_jitter(width = 0.12, height = 0, seed = 123),
              size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Top 1" = "#08519C",
                               "Top 2" = "#3182BD",
                               "Top 3" = "#9ECAE1")) +
  scale_color_manual(values = c("Top 1" = "#08306B",
                                "Top 2" = "#08519C",
                                "Top 3" = "#4292C6")) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cumulative contribution") +
  theme_feature_driven

p_contribution

ggsave("4_manuscript/Supplementary_figures/feature_driven/top_feature_contribution.pdf",
       p_contribution, width = 4, height = 4.5)
