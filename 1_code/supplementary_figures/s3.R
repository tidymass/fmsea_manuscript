library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)
library(tidyverse)
library(ggforce)

# =========================
# Load data
# =========================

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
# Combine + FDR
# =========================

feature_table <-
  bind_rows(feature_table_pos, feature_table_neg) %>%
  mutate(
    FDR          = p.adjust(pvalue, method = "BH"),
    significance = case_when(
      FDR < 0.05 & log2fc > 0 ~ "Up",
      FDR < 0.05 & log2fc < 0 ~ "Down",
      TRUE                     ~ "NS"
    )
  )

# =========================
# Up-regulated features Venn
# =========================

fdr_thresholds <- c(0.1, 0.05, 0.01)

total      <- nrow(feature_table)
up_counts  <- sapply(fdr_thresholds, function(thresh) {
  sum(feature_table$FDR < thresh & feature_table$log2fc > 0)
})

all_counts <- c(total, up_counts)
all_labels <- c("All features", paste0("FDR < ", fdr_thresholds))

R     <- 1
radii <- R * sqrt(all_counts / max(all_counts))

circle_df <- data.frame(
  label = factor(all_labels, levels = all_labels),
  count = all_counts,
  r     = radii,
  x0    = 0,
  y0    = -R + radii
)

fills_up <- c("All features" = "#D0D0D0",
              "FDR < 0.1"    = "#FDAE6B",
              "FDR < 0.05"   = "#E6550D",
              "FDR < 0.01"   = "#A63603")

circle_df <- arrange(circle_df, desc(r))

legend_df <- data.frame(
  label = paste0("FDR < ", fdr_thresholds),
  fill  = c("#FDAE6B", "#E6550D", "#A63603"),
  x     = 1.4,
  y     = c(0.4, 0.15, -0.1),
  xend  = 1.7,
  yend  = c(0.55, 0.3, 0.05)
)

p_up <- ggplot() +
  geom_circle(
    data = circle_df,
    aes(x0 = x0, y0 = y0, r = r, fill = label),
    alpha = 0.75, color = "white", linewidth = 0.7
  ) +
  geom_rect(
    data = legend_df,
    aes(xmin = x, xmax = xend, ymin = y, ymax = yend, fill = label),
    color = NA
  ) +
  geom_text(
    data = legend_df,
    aes(x = xend + 0.05, y = (y + yend) / 2, label = label),
    hjust = 0, vjust = 0.5, size = 3.5
  ) +
  scale_fill_manual(values = fills_up) +
  coord_fixed(xlim = c(-1.2, 2.4), ylim = c(-1.2, 1.2)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(
    title    = "Significant Up-regulated Features",
    subtitle = "Circle area ∝ feature count"
  )

ggsave("4_manuscript/Supplementary_figures/features_venn/venn_features_up.pdf",
       plot = p_up, width = 5, height = 5)

# =========================
# Down-regulated features Venn
# =========================

down_counts <- sapply(fdr_thresholds, function(thresh) {
  sum(feature_table$FDR < thresh & feature_table$log2fc < 0)
})

all_counts <- c(total, down_counts)

circle_df <- data.frame(
  label = factor(all_labels, levels = all_labels),
  count = all_counts,
  r     = R * sqrt(all_counts / max(all_counts)),
  x0    = 0,
  y0    = -R + R * sqrt(all_counts / max(all_counts))
)

fills_down <- c("All features" = "#D0D0D0",
                "FDR < 0.1"    = "#9ECAE1",
                "FDR < 0.05"   = "#3182BD",
                "FDR < 0.01"   = "#08519C")

circle_df <- arrange(circle_df, desc(r))

legend_df <- data.frame(
  label = paste0("FDR < ", fdr_thresholds),
  fill  = c("#9ECAE1", "#3182BD", "#08519C"),
  x     = 1.4,
  y     = c(0.4, 0.15, -0.1),
  xend  = 1.7,
  yend  = c(0.55, 0.3, 0.05)
)

p_down <- ggplot() +
  geom_circle(
    data = circle_df,
    aes(x0 = x0, y0 = y0, r = r, fill = label),
    alpha = 0.75, color = "white", linewidth = 0.7
  ) +
  geom_rect(
    data = legend_df,
    aes(xmin = x, xmax = xend, ymin = y, ymax = yend, fill = label),
    color = NA
  ) +
  geom_text(
    data = legend_df,
    aes(x = xend + 0.05, y = (y + yend) / 2, label = label),
    hjust = 0, vjust = 0.5, size = 3.5
  ) +
  scale_fill_manual(values = fills_down) +
  coord_fixed(xlim = c(-1.2, 2.4), ylim = c(-1.2, 1.2)) +
  theme_void() +
  theme(legend.position = "none") +
  labs(
    title    = "Significant Down-regulated Features",
    subtitle = "Circle area ∝ feature count"
  )

ggsave("4_manuscript/Supplementary_figures/features_venn/venn_features_down.pdf",
       plot = p_down, width = 5, height = 5)
