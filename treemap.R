# Treemap
library(treemap)
library(dplyr)
library(openxlsx)

setwd('~/Desktop/Bulk_BCRrep_Pipeline/clone_expansion_dataMerge')

file_list <- list.files(path = './', pattern = '\\_clone-pass.tsv$', full.names = F)
# file_name <- 'A1_L1_m.process_clone-pass.tsv'
for (file_name in file_list) {
  file_pre <- sapply(strsplit(file_name, "_"), `[`, 1)
  data <- read.table(file_name, sep = '\t', header = T)
  plot <- data %>%
    select(vgene, clone_id) %>%
    group_by(vgene, clone_id) %>%
    mutate(n_seqs = n()) %>%
    mutate(note = paste0(clone_id, '(', n_seqs, ')')) %>%
    distinct()
  
  pdf(paste0('treemap/', file_pre, '.pdf'))
  treemap(plot, index = c('vgene', 'note'), vSize = 'n_seqs', title = file_pre, fontsize.labels = c(20, 10))
  dev.off()
}


