library(r4projects)
setwd(get_project_wd())
rm(list = ls())

source("1_code/annotate_isotope.R")
source("1_code/annotate_feature_table.R")
source("1_code/utils.R")
source("1_code/functions_data_checking.R")
source("1_code/utilities.R")
source("1_code/functions_mfc.R")
source("1_code/functions_redundancy_new.R")

library(tidyverse)
library(tidymass)
library(metid)

####load the correlation with GA
load("3_data_analysis/2_feature_table/feature_table2.rda")
# load("annotation_table_ms2.rda")
load("2_data/1_HMDB_MS1/hmdb_compound_ms1.rda")

annotation_table_final <-
  annotate_feature_table(
    feature_table = feature_table,
    column = "rp",
    metabolite_database = hmdb_compound_ms1,
    ms1_match_ppm = 15,
    mfc_rt_tol = 10,
    isotope_number = 3
  )

###remove redundant annotation according to compound
library(plyr)

calculate_redundancy(annotation_table = annotation_table_final)

annotation_table_final2 =
  remove_redundancy_new(
    annotation_table = annotation_table_final, 
    lab_id_high_cutoff = 0)

calculate_redundancy(annotation_table = annotation_table_final2)


save(annotation_table_final2, file = "/Users/yijiang/Desktop/annotation_table_final.rda")
