library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(massdatabase)
library(metpath)

download_kegg_pathway(path = "2_data/8_mice_diet/kegg_mmu", sleep = 1, organism = "mmu")

data <- read_kegg_pathway(path = "2_data/8_mice_diet/kegg_mmu")

mmu_pathway <- convert_kegg2metpath(data = data, path = "2_data/8_mice_diet/kegg_mmu", threads = 5)

save(mmu_pathway, file = "3_data_analysis/8_mice_diet/mmu_pathway.rda", compress = 'xz')