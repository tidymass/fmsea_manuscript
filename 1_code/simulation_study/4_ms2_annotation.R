library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidymass)

load("2_data/3_Xia_Case_Study_v2/data_cleaning/POS/object_pos2.rda")
load("2_data/3_Xia_Case_Study_v2/data_cleaning/NEG/object_neg2.rda")

object_pos2 <-
  mutate_ms2(
    object = object_pos2,
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/3_Xia_Case_Study_v2/MS2/mgf/POS"
  )

object_neg2 <-
  mutate_ms2(
    object = object_neg2,
    column = "rp",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/3_Xia_Case_Study_v2/MS2/mgf/NEG"
  )

load("2_data/3_Xia_Case_Study_v2/ms2_original/gnps_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/hmdb_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/massbank_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/metlin_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/mona_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/mpsnyder_rplc_ms2.rda")
load("2_data/3_Xia_Case_Study_v2/ms2_original/nist_ms2.rda")


object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    database = gnps_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    database = metlin_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

 
object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    database = nist_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

###
object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    database = gnps_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "negative",
    column = "rp",
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "negative",
    column = "rp",
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "negative",
    column = "rp",
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    database = metlin_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "negative",
    column = "rp",
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "negative",
    column = "rp",
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "negative",
    column = "rp",
    adduct.table = NULL
  )

object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    database = nist_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "negative",
    column = "rp",
    adduct.table = NULL
  )

save(object_pos2, file = "2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_pos2.rda")
save(object_neg2, file = "2_data/3_Xia_Case_Study_v2/metabolite_annotation/object_neg2.rda")
