library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

load("2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_pos2.rda")
load("2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_neg2.rda")

annotation_pos <- object_pos2@annotation_table
annotation_neg <- object_neg2@annotation_table

annotation_table <- 
  rbind(annotation_pos, annotation_neg) %>% 
  dplyr::select(variable_id, Compound.name, KEGG.ID) %>% 
  dplyr::filter(!is.na(KEGG.ID)) %>% 
  distinct()

KEGG_unique <- unique(annotation_table$KEGG.ID)

load("3_data_analysis/3_Xia_case_study/human_pathway_database.rda")

# # Calculate match ratio of HMDB_ID per row against HMDB_unique
# pathway_database <- pathway_database %>%
#   mutate(
#     # Split HMDB_ID in each row (using {} as separator)
#     hmdb_list = str_split(HMDB_ID, "\\{\\}"),
#
#     # Count total HMDB IDs
#     total_count = map_int(hmdb_list, length),
#
#     # Count matches against HMDB_unique
#     matched_count = map_int(hmdb_list, ~ sum(.x %in% HMDB_unique)),
#
#     # Calculate match ratio
#     match_ratio = ifelse(total_count > 0, matched_count / total_count, 0)
#   )

# Calculate match ratio of KEGG_ID per row against KEGG_unique
pathway_database <- pathway_database %>%
  mutate(

    KEGG_list = str_split(KEGG_ID, "\\{\\}"),

    # Count total KEGG IDs
    total_count = map_int(KEGG_list, length),

    # Count matches against KEGG_unique
    matched_count = map_int(KEGG_list, ~ sum(.x %in% KEGG_unique)),

    # Calculate match ratio
    match_ratio = ifelse(total_count > 0, matched_count / total_count, 0)
  )

test <- pathway_database %>% 
  dplyr::filter(MFM_id == "MFM_347")

kegg_ids <- str_split(test$KEGG_ID, "\\{\\}")[[1]]
Match_ids <- kegg_ids[kegg_ids %in% KEGG_unique]

annotation_table_filtering <- annotation_table %>% 
  dplyr::filter(KEGG.ID %in% Match_ids)

ids <- unique(annotation_table_filtering$variable_id)
save(ids, file = "2_data/3_Xia_Case_Study_v2/feature_ids/ids_kegg.rda")







