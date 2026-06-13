library(r4projects)
setwd(get_project_wd())

rm(list = ls())

load("3_data_analysis/5_MTBLS749/results_iteration3.rda")
load("3_data_analysis/5_MTBLS749/annotation_table_filtering3.rda")
load("3_data_analysis/5_MTBLS749/long_table_ranked.rda")

load("2_data/6_Huntington/metabolite_annotation/object_pos2_mona_hmdb_mp_massbank.rda")
load("2_data/6_Huntington/metabolite_annotation/object_neg2_mona_hmdb_mp_massbank.rda")

load("2_data/3_Xia_Case_Study_v2/ms2_original/mpsnyder_rplc_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/mona_ms2.rda")

load("2_data/3_Xia_Case_Study_v2/ms2_original/mpsnyder_rplc_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/mona_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/hmdb_ms2.rda")

annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table

annotation_table_filtering <- rbind(annotation_table_pos, annotation_table_neg)

top1_annotations <- annotation_table_filtering3 %>% 
  filter(Rank == 1)

top1_annotations <- top1_annotations %>% 
  left_join(annotation_table_filtering, by = c("variable_id", "KEGG.ID"))

ms2_data_pos <- object_pos2@ms2_data[[1]]
ms2_data_neg <- object_neg2@ms2_data[[1]]

top1_annotations_filtering <- top1_annotations %>% 
  filter(!variable_id %in% c("M191T47_1_NEG", "M182T51_POS", "M268T228_POS", "M166T94_POS", "M305T53_POS"))

top1_annotations_filtering2 <- top1_annotations_filtering %>%
  group_by(variable_id) %>%
  slice_max(SS, n = 1, with_ties = FALSE) %>%
  ungroup()

# for M129T206_NEG
which(ms2_data_neg@variable_id == "M129T206_NEG")

spec_1 <- ms2_data_neg@ms2_spectra[[48]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.negative"]][["RPLC_738"]][["NCE25"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M129T206_NEG",
                         spectrum2_name = "4-methyl-2-oxo-pentanoic acid")

ggsave("4_manuscript/Supplementary_figures/MS2/M129T206_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M150T44_POS
which(ms2_data_pos@variable_id == "M150T44_POS")

spec_1 <- ms2_data_pos@ms2_spectra[[97]] %>% as.data.frame()

spec_1_standard <- mona_ms2@spectra.data[["Spectra.positive"]][["EMBL-MCF_spec93736"]][["35.0 eV"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M150T44_POS",
                         spectrum2_name = "L-Methionine")

ggsave("4_manuscript/Supplementary_figures/MS2/M150T44_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)

plot_ms2_match <- function(variable_id, save = FALSE, out_dir = "4_manuscript/Supplementary_figures/MS2") {
  
  ms2_data <- if (endsWith(variable_id, "POS")) ms2_data_pos else ms2_data_neg
  
  # Define polarity_slot BEFORE switch
  polarity_slot <- if (endsWith(variable_id, "POS")) "Spectra.positive" else "Spectra.negative"
  
  idx  <- which(ms2_data@variable_id == variable_id)
  spec <- ms2_data@ms2_spectra[[idx]] %>% as.data.frame()
  
  ann              <- top1_annotations_filtering2 %>% filter(variable_id == !!variable_id)
  compound_name    <- ann$Compound.name
  database         <- ann$Database
  spec_id          <- ann$Lab.ID
  collision_energy <- ann$CE
  
  spec_standard <- switch(database,
                          "MoNA_2022-04-27"              = mona_ms2@spectra.data[[polarity_slot]][[spec_id]][[collision_energy]],
                          "Michael_Snyder_RPLC_20220424" = mpsnyder_rplc_ms2@spectra.data[[polarity_slot]][[spec_id]][[collision_energy]],
                          "HMDB_2022-04-11"              = hmdb_ms2@spectra.data[[polarity_slot]][[spec_id]][[collision_energy]],
                          stop("Unknown database: ", database)
  )
  
  p <- masstools::ms2_plot(spec,
                           spec_standard,
                           spectrum1_name = variable_id,
                           spectrum2_name = compound_name)
  
  if (save) {
    path <- file.path(out_dir, paste0(variable_id, ".pdf"))
    ggsave(path, plot = p, width = 8, height = 6, dpi = 300)
    message("Saved to: ", path)
  }
  
  invisible(p)
}

# for M164T93_1_NEG
p <- plot_ms2_match("M164T93_1_NEG")
ggsave("4_manuscript/Supplementary_figures/MS2/M164T93_1_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M192T217_POS
p <- plot_ms2_match("M192T217_POS")
ggsave("4_manuscript/Supplementary_figures/MS2/M192T217_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M203T165_NEG
p <- plot_ms2_match("M203T165_NEG")
ggsave("4_manuscript/Supplementary_figures/MS2/M203T165_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M205T170_POS
p <- plot_ms2_match("M205T170_POS")
ggsave("4_manuscript/Supplementary_figures/MS2/M205T170_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M267T68_NEG
p <- plot_ms2_match("M267T68_NEG")
ggsave("4_manuscript/Supplementary_figures/MS2/M267T68_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M283T93_NEG
p <- plot_ms2_match("M283T93_NEG")
ggsave("4_manuscript/Supplementary_figures/MS2/M283T93_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M285T54_NEG
p <- plot_ms2_match("M285T54_NEG")
ggsave("4_manuscript/Supplementary_figures/MS2/M285T54_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M303T54_NEG
p <- plot_ms2_match("M303T54_NEG")
ggsave("4_manuscript/Supplementary_figures/MS2/M303T54_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M89T44_2_NEG
p <- plot_ms2_match("M89T44_2_NEG")
ggsave("4_manuscript/Supplementary_figures/MS2/M89T44_2_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M75T47_NEG

which(ms2_data_neg@variable_id == "M75T47_NEG")

spec_1 <- ms2_data_neg@ms2_spectra[[6]] %>% as.data.frame()

spec_1_standard <- mona_ms2@spectra.data[["Spectra.negative"]][["FiehnHILIC001234"]][["HCD (NCE 20-30-40%)"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M75T47_NEG",
                         spectrum2_name = "Glycolic acid")

ggsave("4_manuscript/Supplementary_figures/MS2/M75T47_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)
