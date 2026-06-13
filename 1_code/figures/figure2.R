library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(ggplot2)
library(dplyr)
library(ggrepel)
#
# setwd("~/Projects/fmsea_manuscript")

load("3_data_analysis/7_pregnancy_plasma/feature_table.rda")
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
  perm.num = 1000,
  seed = 123,
  fdr.thr = 0.05,
  max.iter.num = 1,
  verbose = TRUE
)


urine_res <- analyze_matrix_relevance(
  results = results,
  sample_source = "urine",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  provider = "openai",
  model = "gpt-5.5"
)

plasma_res <- analyze_matrix_relevance(
  results = results,
  sample_source = "plasma",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  provider = "openai",
  model = "gpt-5.5"
)

save(urine_res, file = "2_data/7_pregnancy_plasma/urine_res.rda")
save(plasma_res, file = "2_data/7_pregnancy_plasma/plasma_res.rda")

urine_modules <- urine_res@significant_modules
plasma_modules <- plasma_res@significant_modules


# boxplot

# Merge data from plasma and urine modules
plasma_data <- plasma_modules %>%
  select(pathway_name, matrix_confidence_score) %>%
  mutate(source = "Plasma", x_num = 1)

urine_data <- urine_modules %>%
  select(pathway_name, matrix_confidence_score) %>%
  mutate(source = "Urine", x_num = 2)

combined_data <- bind_rows(plasma_data, urine_data)

# Prepare wide-format data for paired segments
paired_data <- plasma_modules %>%
  select(pathway_name, plasma_score = matrix_confidence_score) %>%
  left_join(
    urine_modules %>% select(pathway_name, urine_score = matrix_confidence_score),
    by = "pathway_name"
  )

# Assign horizontal x offsets for tied Plasma points (same score → spread along x)
paired_data <- paired_data %>%
  arrange(plasma_score) %>%
  group_by(plasma_score) %>%
  mutate(
    plasma_rank = row_number(),
    plasma_n    = n(),
    plasma_x    = 1 + (plasma_rank - (plasma_n + 1) / 2) * 0.08
  ) %>%
  ungroup()

# Assign horizontal x offsets for tied Urine points (same score → spread along x)
paired_data <- paired_data %>%
  arrange(urine_score) %>%
  group_by(urine_score) %>%
  mutate(
    urine_rank = row_number(),
    urine_n    = n(),
    urine_x    = 2 + (urine_rank - (urine_n + 1) / 2) * 0.08
  ) %>%
  ungroup()

# Plot
p <- ggplot() +
  # Boxplot layer using numeric x positions
  geom_boxplot(
    data = combined_data,
    aes(x = x_num, y = matrix_confidence_score,
        fill = source, group = x_num),
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.6
  ) +
  # Paired connecting lines
  geom_segment(
    data = paired_data,
    aes(
      x = plasma_x, xend = urine_x,
      y = plasma_score, yend = urine_score
    ),
    color = "gray50",
    alpha = 0.5,
    linewidth = 0.4
  ) +
  # Horizontally spread points for Plasma
  geom_point(
    data = paired_data,
    aes(x = plasma_x, y = plasma_score),
    size = 2,
    color = "#C0392B",
    alpha = 0.9
  ) +
  # Horizontally spread points for Urine
  geom_point(
    data = paired_data,
    aes(x = urine_x, y = urine_score),
    size = 2,
    color = "#D4AC0D",
    alpha = 0.9
  ) +
  # Pathway labels to the left of each Plasma point
  geom_text_repel(
    data = paired_data,
    aes(
      x = plasma_x,
      y = plasma_score,
      label = pathway_name
    ),
    hjust         = 1,
    direction     = "y",
    nudge_x       = -0.15,
    segment.size  = 0.3,
    segment.color = "gray60",
    size          = 2.5,
    color         = "gray30",
    max.overlaps  = Inf
  ) +
  # Set colors for Plasma (red) and Urine (yellow)
  scale_fill_manual(
    values = c("Plasma" = "#C44E52", "Urine" = "#E6B800")
  ) +
  # Map numeric x to named axis labels
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Plasma", "Urine"),
    limits = c(0.3, 2.5)
  ) +
  # Prevent label clipping outside plot area
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Matrix Confidence Score",
    title = "Pathway Matrix Confidence Score: Plasma vs Urine"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 13),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 150)
  )

p

ggsave(
  filename = "4_manuscript/Figures/Figure_2/plasma_urine_boxplot.pdf",
  plot     = p,
  width    = 8,
  height   = 6)
