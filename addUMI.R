library(stringr)
setwd('E:/HJY_BCRrep_0628/analysis/HJY_BCRrep_QC/productive')


file_list <- list.files(path = './', pattern = '\\.process.tsv$', full.names = F)
for (file_name in file_list){
  data <- read.table(file_name, sep = '\t', header = T)
  
  UNV_true_count <- sum(data$UNV_lv_dis <= 4)
  UNV_CTT_count <- sum(data$UNV_lv_dis <= 4 & data$CTT_lv_dis <= 4)
  ###########################################################
  print(file_name)
  print(paste0("比对成功序列：", nrow(data)))
  print(paste0("有UNV的序列：", UNV_true_count))
  print(paste0("有UNV和CTTGGG的序列：", UNV_CTT_count))
  ###########################################################
  
  data <- data[data$UNV_lv_dis <= 4 & data$CTT_lv_dis <= 4, ]
  data$UMI <- data$barcode
  data$UMI[data$UNV_ham_dis > data$UNV_lv_dis] <- NA
  data$UMI[data$UNV_ham_dis > data$UNV_lv_dis] <- str_match(data$sequence[data$UNV_ham_dis > data$UNV_lv_dis], "GAG(.{14,18}?)CTTGGG")[, 2]
  print(paste0("没有分配到UMI的序列: ", sum(is.na(data$UMI))))
  
  data <- data[!is.na(data$UMI), ]
  write.table(data, paste0('withUMI/', file_name), sep = '\t', row.names = F, quote = F)
  
  barcodeID_table <- data
  barcodeID_table <- barcodeID_table[, c('sequence_id', 'UMI')]
  barcodeID_table$v_call <- 'IGHV3-30'
  barcodeID_table$j_call <- 'IGHJ1'
  barcodeID_table$junction <- barcodeID_table$UMI
  write.table(barcodeID_table, paste0('withUMI/barcode_id/', file_name), sep = '\t', row.names = F, quote = F)
  
  cdr3ID_table <- data
  cdr3ID_table <- cdr3ID_table[c('sequence_id', 'cdr3')]
  cdr3ID_table$v_call <- 'IGHV3-30'
  cdr3ID_table$j_call <- 'IGHJ1'
  cdr3ID_table$junction <- cdr3ID_table$cdr3
  write.table(cdr3ID_table, paste0('withUMI/cdr3_id/', file_name), sep = '\t', row.names = F, quote = F)
}

