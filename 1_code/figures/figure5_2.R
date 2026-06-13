library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(dplyr)
library(httr)
library(jsonlite)

load("3_data_analysis/5_MTBLS749/annotation_table_filtering3.rda")
load("3_data_analysis/5_MTBLS749/long_table_ranked.rda")
load("2_data/1_HMDB_MS1/hmdb_metabolite.rda")

load("2_data/6_Huntington/metabolite_annotation/object_pos2_mona_hmdb_mp_massbank.rda")
load("2_data/6_Huntington/metabolite_annotation/object_neg2_mona_hmdb_mp_massbank.rda")



# annotation_table_not_top3 <- annotation_table_filtering3 %>% 
#   filter(Rank != 1 | is.na(Rank)) %>% 
#   dplyr::rename(MS2_KEGG.ID = KEGG.ID)

annotation_table_not_top3 <- annotation_table_filtering3 %>% 
  filter(!(Rank %in% c(1, 2, 3))) %>% 
  dplyr::rename(MS2_KEGG.ID = KEGG.ID)

ids <- unique(annotation_table_not_top3$variable_id)
# ms1_results <- long_table_ranked %>% 
#   filter(feature_id %in% ids) %>% 
#   dplyr::filter(Rank == 1) %>% 
#   dplyr::rename(MS1_KEGG.ID = KEGG.ID)

top3_ms1_results <- long_table_ranked %>%
  filter(feature_id %in% ids) %>%
  dplyr::filter(Rank == 1 | Rank == 2 | Rank == 3) %>%
  dplyr::rename(MS1_KEGG.ID = KEGG.ID)

# annotation_table_not_top1 <- annotation_table_not_top1 %>% 
#   left_join(
#     ms1_results %>% dplyr::select(feature_id, MS1_KEGG.ID),
#     by = c("variable_id" = "feature_id")
#   )

hmdb_kegg <- hmdb_metabolite %>% 
  dplyr::select(accession, kegg_id) %>% 
  dplyr::rename(HMDB.ID = accession, KEGG.ID = kegg_id) %>%
  dplyr::filter(!is.na(HMDB.ID) & !is.na(KEGG.ID)) %>% 
  distinct()

# annotation_table_not_top1 <- annotation_table_not_top1 %>% 
#   left_join(
#     hmdb_kegg,
#     by = c("MS2_KEGG.ID" = "KEGG.ID")
#   ) %>% 
#   dplyr::rename(MS2_HMDB.ID = HMDB.ID) %>% 
#   left_join(
#     hmdb_kegg,
#     by = c("MS1_KEGG.ID" = "KEGG.ID")
#   ) %>% 
#   dplyr::rename(MS1_HMDB.ID = HMDB.ID)

annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table

annotation_table_filtering <- rbind(annotation_table_pos, annotation_table_neg)

hmdb_kegg <- annotation_table_filtering %>% 
  dplyr::select(HMDB.ID, KEGG.ID) %>% 
  dplyr::filter(!is.na(HMDB.ID) & !is.na(KEGG.ID)) %>% 
  distinct()

annotation_table_not_top3 <- annotation_table_not_top3 %>%
  left_join(
    hmdb_kegg,
    by = c("MS2_KEGG.ID" = "KEGG.ID")
  )

write.csv(annotation_table_not_top1, "3_data_analysis/5_MTBLS749/annotation_table.csv", row.names = FALSE)
annotation_table_not_top1_hmdb <- read.csv("3_data_analysis/5_MTBLS749/annotation_table_manually_checked.csv")
  
# # Match by MS2_HMDB.ID
# annotation_ms2 <- annotation_table_not_top1_hmdb %>%
#   left_join(
#     hmdb_metabolite %>% dplyr::select(accession, kingdom, super_class, class, sub_class, description),
#     by = c("MS2_HMDB.ID" = "accession")
#   ) %>%
#   dplyr::rename(
#     ms2_kingdom     = kingdom,
#     ms2_super_class = super_class,
#     ms2_class       = class,
#     ms2_sub_class   = sub_class,
#     ms2_description = description
#   )

annotation_table_not_top3_combined <- annotation_table_not_top3 %>% 
  left_join(
    hmdb_metabolite %>% dplyr::select(accession, kingdom, super_class, class, sub_class, description),
    by = c("HMDB.ID" = "accession")
  )

kegg_hmdb_map <- 
  metorigindb_database_final %>% 
  dplyr::select(KEGG_ID, HMDB_ID) %>%
  mutate(
    KEGG_ID = as.character(KEGG_ID),
    HMDB_ID = as.character(HMDB_ID)
  ) %>%
  filter(!is.na(KEGG_ID), !is.na(HMDB_ID)) %>%     # keep valid IDs
  distinct(KEGG_ID, HMDB_ID) %>%                   # remove duplicates
  group_by(KEGG_ID) %>%
  summarise(HMDB_ID = dplyr::first(HMDB_ID), .groups = "drop")

top3_ms1_results <- top3_ms1_results %>% 
  left_join(
    kegg_hmdb_map,
    by = c("MS1_KEGG.ID" = "KEGG_ID"))

write.csv(top3_ms1_results, "3_data_analysis/5_MTBLS749/top3_ms1_results.csv", row.names = FALSE)

top3_ms1_results_combined <- top3_ms1_results %>% 
  left_join(
    hmdb_metabolite %>% dplyr::select(accession, kingdom, super_class, class, sub_class, description),
    by = c("HMDB_ID" = "accession")
  )




annotation_combined <- annotation_table_not_top1_hmdb %>%
  # First join: match MS2
  left_join(
    hmdb_metabolite %>% dplyr::select(accession, kingdom, super_class, class, sub_class, description),
    by = c("MS2_HMDB.ID" = "accession")
  ) %>%
  dplyr::rename(
    ms2_kingdom     = kingdom,
    ms2_super_class = super_class,
    ms2_class       = class,
    ms2_sub_class   = sub_class,
    ms2_description = description
  ) %>%
  # Second join: match MS1
  left_join(
    hmdb_metabolite %>% dplyr::select(accession, kingdom, super_class, class, sub_class, description),
    by = c("MS1_HMDB.ID" = "accession")
  ) %>%
  dplyr::rename(
    ms1_kingdom     = kingdom,
    ms1_super_class = super_class,
    ms1_class       = class,
    ms1_sub_class   = sub_class,
    ms1_description = description
  ) %>%
  # Compute match columns
  dplyr::mutate(
    kingdom_match     = dplyr::coalesce(ms2_kingdom     == ms1_kingdom,     FALSE),
    super_class_match = dplyr::coalesce(ms2_super_class == ms1_super_class, FALSE),
    class_match       = dplyr::coalesce(ms2_class       == ms1_class,       FALSE),
    sub_class_match   = dplyr::coalesce(ms2_sub_class   == ms1_sub_class,   FALSE)
  ) %>%
  dplyr::relocate(variable_id, .before = everything())

###
# Function to retrieve compound description via PubChem CID
get_pubchem_description <- function(cid) {
  if (is.na(cid) || cid == "" || is.null(cid)) return(NA_character_)
  
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
    cid,
    "/description/JSON"
  )
  
  tryCatch({
    res <- GET(url, timeout(10))
    if (status_code(res) != 200) return(NA_character_)
    
    content <- fromJSON(content(res, "text", encoding = "UTF-8"))
    
    # Extract description text (take first non-empty entry)
    descs <- content$InformationList$Information$Description
    descs <- descs[!is.na(descs) & nchar(trimws(descs)) > 0]
    
    if (length(descs) == 0) return(NA_character_)
    return(descs[1])
    
  }, error = function(e) {
    message("Error for CID ", cid, ": ", e$message)
    return(NA_character_)
  })
}

# Batch query with cache (avoid repeated requests for the same CID)
get_descriptions_cached <- function(cid_vector) {
  unique_cids <- unique(na.omit(cid_vector))
  cache <- setNames(
    sapply(unique_cids, function(cid) {
      Sys.sleep(0.2)  # avoid rate limiting
      get_pubchem_description(cid)
    }),
    as.character(unique_cids)
  )

  # Map back to the original vector
  result <- cache[as.character(cid_vector)]
  result[is.na(cid_vector)] <- NA
  unname(result)
}

# Main pipeline
annotation_combined2 <- annotation_combined %>%
  mutate(
    MS2_PUBCHEM_Description = get_descriptions_cached(MS2_PUBCHEM.ID),
    MS1_PUBCHEM_Description = get_descriptions_cached(MS1_PUBCHEM.ID)
  )

save(annotation_combined2, file = "3_data_analysis/5_MTBLS749/annotation_combined2.rda")








