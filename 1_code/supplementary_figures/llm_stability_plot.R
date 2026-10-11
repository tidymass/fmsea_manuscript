library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

load("3_data_analysis/9_llm_stability/llm_stability_results.rda")

out_dir <- "4_manuscript/Supplementary_figures/llm_stability"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_levels <- c("GPT-5.5", "Qwen3-32B", "Qwen3.6-35B-A3B")
label_colors <- c("Plasma" = "#C44E52", "Urine" = "#E6B800")

theme_llm_stability <-
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 11)
  )

# metabolite sets ordered by their mean score over all models, labels and runs
pathway_levels <- llm_scores %>%
  group_by(pathway_name) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_score) %>%
  pull(pathway_name)

# ---- A: distribution of the 10 scores of each metabolite set ----
score_counts <- llm_scores %>%
  filter(!is.na(score)) %>%
  count(model, sample_source, pathway_name, score) %>%
  mutate(
    model = factor(model, levels = model_levels),
    sample_source = factor(str_to_title(sample_source), levels = names(label_colors)),
    pathway_name = factor(pathway_name, levels = pathway_levels)
  )

p_scores <- ggplot(score_counts,
                   aes(x = score, y = pathway_name, size = n, color = sample_source)) +
  geom_point(position = position_dodge(width = 0.6, reverse = TRUE), alpha = 0.85) +
  scale_size_area(max_size = 5, breaks = c(1, 5, 10), name = "Number of runs") +
  scale_color_manual(values = label_colors, name = "Matrix label") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(-8, 108)) +
  facet_wrap(~ model, nrow = 1) +
  labs(x = "Matrix confidence score", y = NULL) +
  theme_llm_stability +
  theme(panel.grid.major.y = element_line(color = "grey92"))

p_scores

ggsave(file.path(out_dir, "score_distribution.pdf"),
       p_scores, width = 12, height = 4.5)

# ---- B: pairwise weighted kappa between runs, within and between models ----
comparison_levels <- c(
  paste("Within", model_levels),
  "GPT-5.5 vs Qwen3-32B", "GPT-5.5 vs Qwen3.6-35B-A3B", "Qwen3-32B vs Qwen3.6-35B-A3B"
)

# comparisons listed from top to bottom: within-model first, then between models
kappa_data <- pairwise_kappa %>%
  mutate(
    comparison = factor(comparison, levels = rev(comparison_levels)),
    type = if_else(str_detect(comparison, "^Within"), "Within model", "Between models")
  )

p_kappa <- ggplot(kappa_data, aes(x = weighted_kappa, y = comparison)) +
  geom_boxplot(aes(fill = type), width = 0.55, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = type),
              position = position_jitter(width = 0, height = 0.15, seed = 123),
              size = 1, alpha = 0.6) +
  scale_fill_manual(values = c("Within model" = "#9ECAE1", "Between models" = "#D9D9D9")) +
  scale_color_manual(values = c("Within model" = "#4292C6", "Between models" = "#737373")) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Weighted kappa between two runs", y = NULL) +
  theme_llm_stability +
  theme(panel.grid.major.y = element_blank(), legend.position = "none")

p_kappa

ggsave(file.path(out_dir, "pairwise_kappa.pdf"),
       p_kappa, width = 7, height = 4)

# ---- C: plasma vs urine comparison of the 9 metabolite sets in each run ----
direction_colors <- c("Plasma higher" = "#C44E52", "Equal" = "#BDBDBD", "Urine higher" = "#E6B800")

direction_data <- label_diff_by_run %>%
  dplyr::select(model, run, `Plasma higher` = n_plasma_higher, Equal = n_equal,
                `Urine higher` = n_urine_higher) %>%
  pivot_longer(-c(model, run), names_to = "direction", values_to = "n_sets") %>%
  mutate(
    model = factor(model, levels = model_levels),
    direction = factor(direction, levels = rev(names(direction_colors)))
  )

p_direction <- ggplot(direction_data, aes(x = factor(run), y = n_sets, fill = direction)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = direction_colors, breaks = names(direction_colors), name = NULL) +
  scale_y_continuous(breaks = seq(0, 9, 3), expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ model, nrow = 1) +
  labs(x = "Run", y = "Number of metabolite sets") +
  theme_llm_stability +
  theme(panel.grid.major.x = element_blank())

p_direction

ggsave(file.path(out_dir, "plasma_vs_urine_by_run.pdf"),
       p_direction, width = 12, height = 3.5)
