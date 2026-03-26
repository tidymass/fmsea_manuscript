library(featuremsea)
library(dplyr)
library(ggplot2)

load("demo_data/7_fmsea_results/results_kegg.rda")

urine_res <- analyze_matrix_relevance(
  results = results,
  sample_source = "urine",
  api_key
)

plasma_res <- analyze_matrix_relevance(
  results = results,
  sample_source = "plasma",
  api_key
)

combined <- bind_rows(urine_res, plasma_res)

ggplot(combined, aes(x = pathway_name, y = matrix_relevance_score, fill = matrix_source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
  scale_fill_manual(values = c("urine" = "#FACE68", "plasma" = "#FA6868")) +
  labs(x = "Pathway", y = "Matrix Relevance Score", fill = "Sample Source") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
