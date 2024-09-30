library(dplyr)

setwd('E:/HJY_BCRrep_0628/analysis/HJY_BCRrep_QC/productive/withUMI')

file_list <- list.files(path = './', pattern = '\\.process.tsv$', full.names = F)
for (file_name in file_list){
  file_pre <- sub("\\.tsv$", "", file_name)
  data <- read.table(file_name, sep = '\t', header = T)
  
  barcode_table <- read.table(paste0('./barcode_id/', file_pre, '_clone-pass.tsv'), sep = '\t', header = T)
  barcode_table$barcode_id <- barcode_table$clone_id
  barcode_table <- barcode_table[, c('sequence_id', 'barcode_id')]
  
  cdr3_table <- read.table(paste0('./cdr3_id/', file_pre, '_clone-pass.tsv'), sep = '\t', header = T)
  cdr3_table$cdr3_id <- cdr3_table$clone_id
  cdr3_table <- cdr3_table[, c('sequence_id', 'cdr3_id')]
  
  
  data <- merge(data, barcode_table, by = 'sequence_id', all.x = T)
  data <- merge(data, cdr3_table, by = 'sequence_id', all.x = T)
  
  data <- data %>%
    arrange(barcode_id, cdr3_id)
  
  test <- data %>%
    distinct(barcode_id, cdr3_id, .keep_all = T)
  write.table(test, paste0('no_PCRrepeat/', file_name), sep = '\t', row.names = F, quote = F)
  
  ##########
  print(file_name)
  print(paste0("共有barcode_id数:", length(unique(test$barcode_id))))
  print(paste0("共有cdr3_id数:", length(unique(test$cdr3_id))))
  print(paste0("去除扩增之后的序列数为：", nrow(test)))
}


# x <- data[, c('sequence_id', 'v_call', 'j_call', 'cdr3', 'UNV', 'barcode', 'CTTGGG', 'UMI', 'barcode_id', 'cdr3_id', 'sequence')] 


setwd('E:/HJY_BCRrep_0628/analysis/HJY_BCRrep_QC/productive/remain/productive/withUMI')

file_list <- list.files(path = './', pattern = '\\.process.tsv$', full.names = F)
for (file_name in file_list){
  file_pre <- sub("\\.tsv$", "", file_name)
  data <- read.table(file_name, sep = '\t', header = T)
  
  barcode_table <- read.table(paste0('./barcode_id/', file_pre, '_clone-pass.tsv'), sep = '\t', header = T)
  barcode_table$barcode_id <- barcode_table$clone_id
  barcode_table <- barcode_table[, c('sequence_id', 'barcode_id')]
  
  cdr3_table <- read.table(paste0('./cdr3_id/', file_pre, '_clone-pass.tsv'), sep = '\t', header = T)
  cdr3_table$cdr3_id <- cdr3_table$clone_id
  cdr3_table <- cdr3_table[, c('sequence_id', 'cdr3_id')]
  
  
  data <- merge(data, barcode_table, by = 'sequence_id', all.x = T)
  data <- merge(data, cdr3_table, by = 'sequence_id', all.x = T)
  
  data <- data %>%
    arrange(barcode_id, cdr3_id)
  
  test <- data %>%
    distinct(barcode_id, cdr3_id, .keep_all = T)
  write.table(test, paste0('no_PCRrepeat/', file_name), sep = '\t', row.names = F, quote = F)
  
  ##########
  print(file_name)
  print(paste0("共有barcode_id数:", length(unique(test$barcode_id))))
  print(paste0("共有cdr3_id数:", length(unique(test$cdr3_id))))
  print(paste0("去除扩增之后的序列数为：", nrow(test)))
}


