library(tidymass)
library(tidyverse)

load("2_data/6_Huntington/metabolite_annotation/object_pos2.rda")
load("2_data/6_Huntington/metabolite_annotation/object_neg2.rda")

object_pos2 <- 
  object_pos2 %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2)

object_neg2 <- 
  object_neg2 %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2)

expression_data_pos <- object_pos2@expression_data
colnames(expression_data_pos) <- gsub("_POS", "", colnames(expression_data_pos))
expression_data_pos <- expression_data_pos[, 1:(ncol(expression_data_pos) - 15)]


expression_data_neg <- object_neg2@expression_data
colnames(expression_data_neg) <- gsub("_NEG", "", colnames(expression_data_neg))
expression_data_neg <- expression_data_neg[, 1:(ncol(expression_data_neg) - 12)]

expression_data <- rbind(expression_data_pos, expression_data_neg)

# object <- 
#   merge_mass_dataset(x = object_pos2, 
#                      y = object_neg2, 
#                      sample_direction = "inner",
#                      variable_direction = "full", 
#                      sample_by = "sample_id", 
#                      variable_by = c("variable_id", "mz", "rt"))

annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table
annotation_table <- rbind(annotation_table_pos, annotation_table_neg)

annotation_talbe_filtering <- annotation_table %>% dplyr::select(variable_id, KEGG.ID) 
annotation_talbe_filtering2 <- annotation_talbe_filtering %>% dplyr::filter(!is.na(KEGG.ID)) %>% distinct()

sample_info <- object_pos2@sample_info

nrow(sample_info %>% dplyr::filter(group == "Case"))
nrow(sample_info %>% dplyr::filter(group == "Control"))

group_row <- data.frame(matrix(c(rep("Case", 42), rep("Control", 40)), nrow = 1))
colnames(group_row) <- colnames(expression_data)
rownames(group_row) <- "Group"


# Match by row name and variable_id
expression_data$variable_id <- rownames(expression_data)

# Merge: one-to-many will generate new rows automatically
expression_data_kegg <- merge(expression_data, 
                              annotation_talbe_filtering2, 
                              by = "variable_id", 
                              all.x = TRUE)

expression_data <- expression_data_kegg %>% dplyr::filter(!is.na(KEGG.ID))
expression_data <- expression_data[!duplicated(expression_data$KEGG.ID), ]
rownames(expression_data) <- expression_data$KEGG.ID
expression_data$variable_id <- NULL
expression_data$KEGG.ID <- NULL

expression_data <- rbind(group_row, expression_data)

write.csv(expression_data, "3_data_analysis/5_MTBLS749/expression_data_kegg.csv", row.names = TRUE)

