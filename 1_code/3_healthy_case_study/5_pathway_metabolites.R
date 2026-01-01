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

library(dplyr)
library(stringr)

# # 计算每行 HMDB_ID 与 HMDB_unique 的匹配占比
# pathway_database <- pathway_database %>%
#   mutate(
#     # 分割每行的 HMDB_ID（使用 {} 作为分隔符）
#     hmdb_list = str_split(HMDB_ID, "\\{\\}"),
#     
#     # 计算总的 HMDB ID 数量
#     total_count = map_int(hmdb_list, length),
#     
#     # 计算与 HMDB_unique 匹配的数量
#     matched_count = map_int(hmdb_list, ~ sum(.x %in% HMDB_unique)),
#     
#     # 计算匹配占比
#     match_ratio = ifelse(total_count > 0, matched_count / total_count, 0)
#   ) 

# 计算每行 KEGG_ID 与 KEGG_unique 的匹配占比
pathway_database <- pathway_database %>%
  mutate(

    KEGG_list = str_split(KEGG_ID, "\\{\\}"),
    
    # 计算总的 KEGG ID 数量
    total_count = map_int(KEGG_list, length),
    
    # 计算与 KEGG_unique 匹配的数量
    matched_count = map_int(KEGG_list, ~ sum(.x %in% KEGG_unique)),
    
    # 计算匹配占比
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







