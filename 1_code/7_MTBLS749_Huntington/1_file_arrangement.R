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


# 目标文件夹
dest_dir <- "~/cases/MTBLS749_Huntington/Result_manifest_control/POS/Case"

# 源文件夹
source_dir <- "~/cases/MTBLS749_Huntington/FILES"

# manifest_pos_ids 已经是过滤好的 Sample Name 列表
# 例如: c("1_Rep1_POS", "2_Rep1_POS", ...)

# 构建完整文件路径（假设文件扩展名为 .mzML）
files_to_move <- file.path(source_dir, paste0(manifest_pos_ids, ".mzML"))

# 检查哪些文件实际存在
existing_files <- files_to_move[file.exists(files_to_move)]
missing_files  <- files_to_move[!file.exists(files_to_move)]



# 移动文件
if (length(existing_files) > 0) {
  result <- file.rename(
    from = existing_files,
    to   = file.path(dest_dir, basename(existing_files))
  )
  cat("成功移动:", sum(result), "个文件\n")
  cat("失败:", sum(!result), "个文件\n")
}
