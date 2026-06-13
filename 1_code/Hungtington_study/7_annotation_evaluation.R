library(featuremsea)
library(tidyverse)

load("2_data/6_Huntington/metabolite_annotation/object_pos2.rda")
load("2_data/6_Huntington/metabolite_annotation/object_neg2.rda")

annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table

annotation_table_filtering <- rbind(annotation_table_pos, annotation_table_neg) %>% 
  filter(!is.na(KEGG.ID)) %>% 
  dplyr::select(variable_id, KEGG.ID) %>% 
  distinct()

# annotation_table_pos_filtering <- annotation_table_pos %>% 
#   filter(!is.na(KEGG.ID)) %>% 
#   dplyr::select(variable_id, KEGG.ID) %>% 
#   distinct()

# load("3_data_analysis/5_MTBLS749/results.rda")

load("3_data_analysis/5_MTBLS749/results_original.rda")

results <- results_original
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
#   mutate(Rank = rank(-Score, ties.method = "max")) %>%
#   ungroup()

long_table_ranked <- long_table %>%
  group_by(feature_id) %>%
  mutate(Rank = dense_rank(desc(Score))) %>%
  ungroup()

# Join to annotation_table_filtering2, keep matching rows only
annotation_table_filtering3 <- annotation_table_filtering2 %>%
  left_join(long_table_ranked, by = c("variable_id" = "feature_id", "KEGG.ID")) %>%
  group_by(variable_id) %>%
  slice_min(Rank, n = 1, with_ties = FALSE) %>%
  ungroup()

length(unique(annotation_table_filtering3$variable_id))

annotation_table_filtering3 %>% 
  filter(Rank == 1) %>% 
  pull(variable_id) %>% 
  unique() %>% 
  length()

# visualization

df <- data.frame(
  k    = c(1, 2, 3),
  topk = c("Top-1", "Top-2", "Top-3"),
  acc  = c(10/22, 15/22, 17/22) * 100
)

ggplot(df, aes(x = topk, y = acc)) +
  geom_col(aes(fill = topk), width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/22)", acc, c(10, 15, 17))),
            vjust = -0.4, size = 6) +
  scale_fill_manual(values = c("#b5d4f4", "#378add", "#185fa5")) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Accuracy", title = "Top-3 Annotation Accuracy") +
  theme_minimal(base_size = 16)


# length(unique(annotation_table_filtering2$variable_id))
# 
# annotation_table_filtering2_rank1 <- annotation_table_filtering2 %>% 
#   filter(Rank == 1)
# 
# length(unique(annotation_table_filtering2_rank1$variable_id))
# 
# feature_metabolite_count_filtering_merged <- feature_metabolite_count_filtering %>% 
#   left_join(annotation_table_filtering, by = "variable_id")
# 
# length(unique(feature_metabolite_count_filtering_merged$variable_id))
# 
# feature_metabolite_count_filtering_merged_selected <- feature_metabolite_count_filtering_merged %>% 
#   dplyr::filter(KEGG_ID == KEGG.ID)
