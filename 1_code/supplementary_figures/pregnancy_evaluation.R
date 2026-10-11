library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)

load("3_data_analysis/7_pregnancy_plasma/iteration_results.rda")

load("3_data_analysis/7_pregnancy_plasma/metabolite_annotation/object_pos2_mona_hmdb_mp_massbank.rda")
load("3_data_analysis/7_pregnancy_plasma/metabolite_annotation/object_neg2_mona_hmdb_mp_massbank.rda")
load("2_data/2_KEGG_MS1/kegg_compound_ms1.rda")

results <- iteration_results

annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table

annotation_table_filtering <- rbind(annotation_table_pos, annotation_table_neg)

metabolite_database <- kegg_compound_ms1@spectra.info %>% 
  dplyr::rename(MS2_KEGG.ID = KEGG.ID) %>% 
  dplyr::select(MS2_KEGG.ID, CAS.ID) %>% 
  dplyr::filter(!is.na(MS2_KEGG.ID) & !is.na(CAS.ID)) %>% 
  distinct()

metabolite_database$CAS.ID <- gsub("\\s.*", "", metabolite_database$CAS.ID)

annotation_table_filtering <- 
  annotation_table_filtering %>% 
  left_join(metabolite_database, by = "CAS.ID")

annotation_table_filtering$KEGG.ID <- sub("\\{\\}.*$", "", annotation_table_filtering$KEGG.ID)

annotation_table_filtering <- annotation_table_filtering %>%
  mutate(
    KEGG.ID = if_else(
      !is.na(KEGG.ID) &
        !is.na(MS2_KEGG.ID) &
        KEGG.ID != MS2_KEGG.ID,
      MS2_KEGG.ID,
      KEGG.ID
    )
  )

annotation_table_full <- annotation_table_filtering

annotation_table_filtering <- annotation_table_filtering %>% 
  filter(!is.na(KEGG.ID)) %>% 
  dplyr::select(variable_id, KEGG.ID) %>% 
  distinct()

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
#   mutate(Rank = dense_rank(desc(Score))) %>%
#   ungroup()

long_table_ranked <- long_table %>%
  group_by(feature_id) %>%
  mutate(Rank = min_rank(desc(Score))) %>%
  ungroup()
# 
# long_table_ranked <- long_table %>% 
#   group_by(feature_id) %>% 
#   mutate(Rank = rank(-Score, ties.method = "max")) %>% 
#   ungroup()

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

load("2_data/3_Xia_Case_Study_v2/ms2_original/mona_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/hmdb_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/mpsnyder_rplc_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/massbank_ms2.rda")  

top1_ids <- annotation_table_filtering3 %>% 
  filter(Rank == 1)

top1_annotations <- top1_ids %>%
  left_join(annotation_table_full, by = c("variable_id", "KEGG.ID")) %>%
  filter(!is.na(Lab.ID), !is.na(CE), !is.na(SS), SS > 0) %>%
  group_by(variable_id) %>%
  slice_max(SS, n = 1, with_ties = FALSE) %>%
  ungroup()

no_ms2 <- setdiff(top1_ids$variable_id, top1_annotations$variable_id)
length(no_ms2)

unique(top1_annotations$Database)

db_list <- list(
  "MoNA_2022-04-27"              = mona_ms2,
  "Michael_Snyder_RPLC_20220424" = mpsnyder_rplc_ms2,
  "HMDB_2022-04-11"              = hmdb_ms2,
  "MassBank_2022-04-27"                = massbank_ms2
)

setdiff(unique(top1_annotations$Database), names(db_list))

ms2_data_list <- list(
  POS = object_pos2@ms2_data[[1]],
  NEG = object_neg2@ms2_data[[1]]
)

plot_ms2_match <- function(ann_row) {
  polarity      <- if (endsWith(ann_row$variable_id, "POS")) "POS" else "NEG"
  polarity_slot <- if (polarity == "POS") "Spectra.positive" else "Spectra.negative"
  ms2_data      <- ms2_data_list[[polarity]]
  
  # 实验谱：优先 ms2_spectrum_id，退回 variable_id
  idx <- match(ann_row$ms2_spectrum_id, ms2_data@ms2_spectrum_id)
  if (is.na(idx)) idx <- which(ms2_data@variable_id == ann_row$variable_id)[1]
  if (is.na(idx)) stop("Experimental spectrum not found")
  spec <- as.data.frame(ms2_data@ms2_spectra[[idx]])
  
  # 标准谱
  db <- db_list[[ann_row$Database]]
  if (is.null(db)) stop("Unknown database: ", ann_row$Database)
  spec_std <- db@spectra.data[[polarity_slot]][[ann_row$Lab.ID]][[ann_row$CE]]
  if (is.null(spec_std)) stop("Standard spectrum not found: ", ann_row$Lab.ID, " / ", ann_row$CE)
  
  masstools::ms2_plot(spec,
                      as.data.frame(spec_std),
                      spectrum1_name = ann_row$variable_id,
                      spectrum2_name = ann_row$Compound.name)
}

out_dir <- "4_manuscript/Supplementary_figures/pregnancy_plasma/MS2_comparison"

plot_log <- purrr::map_dfr(seq_len(nrow(top1_annotations)), function(i) {
  ann_row <- top1_annotations[i, ]
  status <- tryCatch({
    p <- plot_ms2_match(ann_row)
    ggplot2::ggsave(file.path(out_dir, paste0(ann_row$variable_id, ".pdf")),
                    plot = p, width = 8, height = 6)
    "OK"
  }, error = function(e) conditionMessage(e))
  
  tibble::tibble(variable_id   = ann_row$variable_id,
                 KEGG.ID       = ann_row$KEGG.ID,
                 Compound.name = ann_row$Compound.name,
                 Database      = ann_row$Database,
                 Lab.ID        = ann_row$Lab.ID,
                 CE            = ann_row$CE,
                 SS            = ann_row$SS,
                 status        = status)
})
