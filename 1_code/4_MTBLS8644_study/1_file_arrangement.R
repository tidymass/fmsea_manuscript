library(r4projects)
setwd(get_project_wd())
rm(list = ls())

setwd("/Volumes/Western/projects/fmsea_cases/MTBLS265")

files <- list.files("FILES", pattern = "\\.mzML$", full.names = TRUE)

# 筛选包含NEG和youth的文件 -> Control
neg_youth_files <- files[grepl("NEG", files) & grepl("youth", files)]

# 筛选包含NEG和elder的文件 -> Case
neg_elder_files <- files[grepl("NEG", files) & grepl("elder", files)]

# 复制文件到对应目录
file.copy(neg_youth_files, "Result/NEG/Control/")
file.copy(neg_elder_files, "Result/NEG/Case/")

pos_youth_files <- files[grepl("POS", files) & grepl("youth", files)]
pos_elder_files <- files[grepl("POS", files) & grepl("elder", files)]

file.copy(pos_youth_files, "Result/POS/Control/")
file.copy(pos_elder_files, "Result/POS/Case/")



# library(fs)
# library(tidyverse)
# library(readr)
# sample_info <- read_tsv("s_MTBLS265.txt") 
# 
# healthy <- sample_info %>%
#   filter(`Factor Value[Cohort]` == "Healthy") %>%
#   pull(`Source Name`)
# 
# early_stage <- sample_info %>%
#   filter(`Factor Value[Cohort]` == "early stage (pre-diabetes)") %>%
#   pull(`Source Name`)
# 
# healthy <- gsub("BX", "Bx", healthy)
# early_stage <- gsub("BX", "Bx", early_stage)
# 
# source_dir <- "FILES/DERIVED_FILES/NEG"
# all_files <- list.files(source_dir, full.names = TRUE)
# 
# for (f in all_files) {
#   filename <- basename(f)
#   
#   if (filename %in% healthy) {
#     file.copy(f, file.path("Result/NEG/Control", filename))
#     cat("Control:", filename, "\n")
#     
#   } else if (filename %in% early_stage) {
#     file.copy(f, file.path("Result/NEG/Case", filename))
#     cat("Case:", filename, "\n")
#     
#   } else {
#     cat("未匹配:", filename, "\n")
#   }
# }
# 
# source_dir <- "FILES/DERIVED_FILES/POS"
# all_files <- list.files(source_dir, full.names = TRUE)
# 
# for (f in all_files) {
#   filename <- basename(f)
#   
#   if (filename %in% healthy) {
#     file.copy(f, file.path("Result/POS/Control", filename))
#     cat("Control:", filename, "\n")
#     
#   } else if (filename %in% early_stage) {
#     file.copy(f, file.path("Result/POS/Case", filename))
#     cat("Case:", filename, "\n")
#     
#   } else {
#     cat("未匹配:", filename, "\n")
#   }
# }
