library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)

load("3_data_analysis/7_pregnancy_plasma/feature_table.rda")

load("2_data/6_Huntington/data_cleaning/POS/object_pos2.rda")
load("2_data/6_Huntington/data_cleaning/NEG/object_neg2.rda")

variable_info_pos <- feature_table |> 
  dplyr::filter(polarity == "positive")

variable_info_neg <- feature_table |> 
  dplyr::filter(polarity == "negative")

object_pos2@variable_info <- variable_info_pos
object_neg2@variable_info <- variable_info_neg

object_pos2 <-
  mutate_ms2(
    object = object_pos2,
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/7_pregnancy_plasma/denmark_project/rp/POS/MS2/mgf"
  )

object_neg2 <-
  mutate_ms2(
    object = object_neg2,
    column = "rp",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/7_pregnancy_plasma/denmark_project/rp/NEG/MS2/mgf"
  )

load("~/Tidymass/ms2_original/mpsnyder_rplc_ms2.rda")
load("~/Tidymass/ms2_original/hmdb_ms2.rda")
load("~/Tidymass/ms2_original/massbank_ms2.rda")
load("~/Tidymass/ms2_original/mona_ms2.rda")

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "rp",
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "rp",
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "rp",
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "rp",
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

# neg
object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

save(object_pos2, file = "3_data_analysis/7_pregnancy_plasma/metabolite_annotation/object_pos2_mona_hmdb_mp_massbank.rda")
save(object_neg2, file = "3_data_analysis/7_pregnancy_plasma/metabolite_annotation/object_neg2_mona_hmdb_mp_massbank.rda")
