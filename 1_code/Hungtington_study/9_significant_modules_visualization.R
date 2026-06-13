library(ggplot2)
library(dplyr)
library(forcats)

load("3_data_analysis/5_MTBLS749/results.rda")

significant_modules <- results@significant_modules

p_bubble <- significant_modules %>%
  mutate(
    FDR = ifelse(FDR == 0, 1e-10, FDR),
    neg_log10_FDR = -log10(FDR),
    direction = ifelse(NES > 0, "Activated", "Suppressed")
  ) %>%
  arrange(NES) %>%
  mutate(pathway_name = fct_inorder(pathway_name)) %>%
  ggplot(aes(
    x = NES,
    y = pathway_name,
    size = neg_log10_FDR,
    color = direction
  )) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_color_manual(
    values = c("Activated" = "#D74E3C", "Suppressed" = "#2474CA")
  ) +
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(3, 10)
  ) +
  scale_x_continuous(limits = c(1.0, 1.5)) +
  guides(color = "none") +
  labs(
    title = "Enriched Metabolic Modules",
    x = "NES",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor  = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

p_bubble

ggsave("3_data_analysis/5_MTBLS749/featureMSEA_bubble.pdf", p_bubble, width = 10, height = 6)
