load("2_data/6_Huntington/data_cleaning/POS/object_pos2.rda")
load("2_data/6_Huntington/data_cleaning/NEG/object_neg2.rda")

library(tidymass)

object_pos2 <-
  mutate_ms2(
    object = object_pos2,
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/6_Huntington/MS2/POS"
  )


object_neg2 <-
  mutate_ms2(
    object = object_neg2,
    column = "rp",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/6_Huntington/MS2/NEG"
  )



object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "rp",
    database = gnps_ms2,
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
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

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




# neg
object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = gnps_ms2,
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
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

# object_neg2 <-
#   annotate_metabolites(
#     object = object_neg2,
#     ms1.match.ppm = 15,
#     polarity = "negative",
#     column = "rp",
#     database = metlin_ms2,
#     based_on = c("ms1", "ms2"),
#     adduct.table = NULL
#   )

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

# object_neg2 <-
#   annotate_metabolites(
#     object = object_neg2,
#     ms1.match.ppm = 15,
#     polarity = "negative",
#     column = "rp",
#     database = nist_ms2,
#     based_on = c("ms1", "ms2"),
#     adduct.table = NULL
#   )

save(object_pos2, file = "2_data/6_Huntington/metabolite_annotation/object_pos2.rda")
save(object_neg2, file = "2_data/6_Huntington/metabolite_annotation/object_neg2.rda")



