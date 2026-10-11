library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)

load("2_data/8_mice_diet/data_cleaning/POS/object_pos2.rda")
load("2_data/8_mice_diet/data_cleaning/NEG/object_neg2.rda")

load("~/Tidymass/ms2_original/mpsnyder_hilic_ms2.rda")
load("~/Tidymass/ms2_original/hmdb_ms2.rda")
load("~/Tidymass/ms2_original/massbank_ms2.rda")
load("~/Tidymass/ms2_original/mona_ms2.rda")

object_pos2 <-
  mutate_ms2(
    object = object_pos2,
    column = "hilic",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/8_mice_diet/MS2/POS"
  )

object_neg2 <-
  mutate_ms2(
    object = object_neg2,
    column = "hilic",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/8_mice_diet/MS2/NEG"
  )

# pos
object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "hilic",
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "hilic",
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "hilic",
    database = mpsnyder_hilic_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "hilic",
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
    column = "hilic",
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "hilic",
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "hilic",
    database = mpsnyder_hilic_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "hilic",
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

save(object_pos2, file = "3_data_analysis/8_mice_diet/metabolite_annotation/object_pos2_mona_hmdb_mp_massbank.rda")
save(object_neg2, file = "3_data_analysis/8_mice_diet/metabolite_annotation/object_neg2_mona_hmdb_mp_massbank.rda")
