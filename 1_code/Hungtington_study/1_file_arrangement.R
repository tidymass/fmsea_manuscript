library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)
library(readr)
df <- read_delim("~/cases/MTBLS749_Huntington/s_MTBLS749.txt")

manifest_neg_ids <- df %>% dplyr::filter(`Factor Value[Phenotype]` == "manifestHD" & str_detect(`Sample Name`, "NEG")) %>% pull(`Sample Name`)
manifest_pos_ids <- df %>% dplyr::filter(`Factor Value[Phenotype]` == "manifestHD" & str_detect(`Sample Name`, "POS")) %>% pull(`Sample Name`)

control_neg_ids <- df %>% dplyr::filter(`Factor Value[Phenotype]` == "control" & str_detect(`Sample Name`, "NEG")) %>% pull(`Sample Name`)
control_pos_ids <- df %>% dplyr::filter(`Factor Value[Phenotype]` == "control" & str_detect(`Sample Name`, "POS")) %>% pull(`Sample Name`)


# Destination folder
dest_dir <- "~/cases/MTBLS749_Huntington/Result_manifest_control/POS/Case"

# Source folder
source_dir <- "~/cases/MTBLS749_Huntington/FILES"

# manifest_pos_ids is the filtered list of Sample Names
# e.g., c("1_Rep1_POS", "2_Rep1_POS", ...)

# Build full file paths (assuming .mzML extension)
files_to_move <- file.path(source_dir, paste0(manifest_pos_ids, ".mzML"))

# Check which files actually exist
existing_files <- files_to_move[file.exists(files_to_move)]
missing_files  <- files_to_move[!file.exists(files_to_move)]



# Move files
if (length(existing_files) > 0) {
  result <- file.rename(
    from = existing_files,
    to   = file.path(dest_dir, basename(existing_files))
  )
  cat("Successfully moved:", sum(result), "files\n")
  cat("Failed:", sum(!result), "files\n")
}
