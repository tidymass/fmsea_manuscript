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

annotation_table_pos <- object_pos2@annotation_table
annotation_table_neg <- object_neg2@annotation_table

annotation_table_filtering <- rbind(annotation_table_pos, annotation_table_neg)

top1_annotations <- annotation_table_filtering3 %>% 
  filter(Rank == 1)

top1_annotations <- top1_annotations %>% 
  left_join(annotation_table_filtering, by = c("variable_id", "KEGG.ID"))

ms2_data_pos <- object_pos2@ms2_data[[1]]
ms2_data_neg <- object_neg2@ms2_data[[1]]

# for M191T47_1_NEG
which(ms2_data_neg@variable_id == "M191T47_1_NEG")

spec_1 <- ms2_data_neg@ms2_spectra[[152]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.negative"]][["RPLC_128"]][["NCE25"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M191T47_1_NEG",
                         spectrum2_name = "Citric acid")

ggsave("4_manuscript/Figures/Figure_5/M191T47_1_NEG.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M182T51_POS
which(ms2_data_pos@variable_id == "M182T51_POS")

spec_1 <- ms2_data_pos@ms2_spectra[[193]] %>% as.data.frame()

spec_1_standard <- mona_ms2@spectra.data[["Spectra.positive"]][["EMBL-MCF_spec353490"]][["30.0 eV"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M182T51_POS",
                         spectrum2_name = "L-Tyrosine")

ggsave("4_manuscript/Figures/Figure_5/M182T51_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M268T228_POS
which(ms2_data_pos@variable_id == "M268T228_POS")

spec_1 <- ms2_data_pos@ms2_spectra[[441]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_283"]][["NCE25"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M268T228_POS",
                         spectrum2_name = "Adenosine")

ggsave("4_manuscript/Figures/Figure_5/M268T228_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)

# for M305T53_POS
which(ms2_data_pos@variable_id == "M305T53_POS")

spec_1 <- ms2_data_pos@ms2_spectra[[523]] %>% as.data.frame()

spec_1_standard <- mona_ms2@spectra.data[["Spectra.positive"]][["FiehnHILIC000550"]][["HCD (NCE 20-30-40%)"]]

p <- masstools::ms2_plot(spec_1, 
                         spec_1_standard,
                         spectrum1_name = "M305T53_POS",
                         spectrum2_name = "N-Acetylaspartylglutamic acid")

ggsave("4_manuscript/Figures/Figure_5/M305T53_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)


