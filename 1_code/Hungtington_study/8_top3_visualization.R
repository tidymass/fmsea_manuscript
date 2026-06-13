library(ggplot2)
library(tidyr)
library(dplyr)

df <- data.frame(
  iteration = c("Iteration 1", "Iteration 2", "Iteration 3"),
  Top1      = c(9/22, 9/24, 11/26) * 100,
  Top2      = c((9+0)/22, (9+3)/24, (11+2)/26) * 100,
  Top3      = c((9+0+1)/22, (9+3+1)/24, (11+2+1)/26) * 100
)

df_long <- df %>%
  pivot_longer(cols = c(Top1, Top2, Top3),
               names_to = "metric", values_to = "acc") %>%
  mutate(metric = factor(metric, levels = c("Top1", "Top2", "Top3")))

p <- 
ggplot(df_long, aes(x = metric, y = acc, color = iteration, group = iteration)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.1f%%", acc)),
            vjust = -1.2, size = 5, show.legend = FALSE) +
  scale_x_discrete(labels = c("Top-1", "Top-2", "Top-3")) +
  scale_y_continuous(limits = c(30, 60), labels = function(x) paste0(x, "%")) +
  scale_color_manual(values = c("Iteration 1" = "#185fa5",
                                "Iteration 2" = "#e07b39",
                                "Iteration 3" = "#3a9e6f")) +
  labs(x = NULL, y = "Accuracy (%)", color = NULL,
       title = "Top-3 Annotation Accuracy across Iterations") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top")

ggsave("3_data_analysis/5_MTBLS749/top_k_accuracy.pdf", plot = p, width = 8, height = 6)
