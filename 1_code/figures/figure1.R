library(pheatmap)
library(ggplotify)
library(ggplot2)

set.seed(42)
mat <- matrix(rnorm(300), nrow = 4, ncol = 3)

p <- as.ggplot(
  pheatmap(mat,
           color = colorRampPalette(c("#1A6FAF", "#FAFAF5", "#C0392B"))(100),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = FALSE,
           legend = FALSE,
           border_color = NA,
           cellwidth = 10,
           cellheight = 10,
           silent = TRUE)
)

ggsave("4_manuscript/Figures/Figure_1/heatmap3.pdf", plot = p,
       width = 30, height = 300, units = "mm",
       device = "pdf")

# enrichment plot

library(ggplot2)

data <- data.frame(
  label = c("F1-M1", "F1-M3", "F2", "F3-M4", "F4"),
  value = c(0.6, 0.3, -0.5, 0.1, -0.5)
)

data$color <- ifelse(data$value >= 0, "#ED1B23", "#1A75BC")
data$label <- factor(data$label, levels = rev(data$label))

n <- nrow(data)

line_x <- c(0, cumsum(data$value))

line_data <- data.frame(
  x = line_x,
  y = c(n + 0.6, n:1)
)

x_min <- min(min(line_x), 0) - 0.1
x_max <- max(max(abs(data$value)), max(line_x)) + 0.1

p <- 
ggplot() +
  geom_col(data = data,
           aes(y = label, x = abs(value), fill = color),
           width = 0.5) +
  geom_path(data = line_data,
            aes(x = x, y = y),
            color = "#FF8C00", linewidth = 0.8) +
  scale_fill_identity() +
  scale_x_continuous(labels = NULL, breaks = NULL,
                     limits = c(x_min, x_max),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.9))) +
  labs(x = "Step", y = "Feature-metabolite pair") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave("4_manuscript/Figures/Figure_1/enrichment_plot.pdf", plot = p, width = 6, height = 4)
