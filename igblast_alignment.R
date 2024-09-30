library(stringdist)


setwd('E:/HJY_BCRrep_0628/analysis/HJY_BCRrep_QC')


file_list <- list.files(path = './', pattern = '\\.process.tsv$', full.names = F)
for (file_name in file_list){
  data <- read.table(file_name, sep = '\t', header = T)
  
  UNV <- 'AAGCAGTGGTATCAACGCAGAG'
  CTT <- 'CTTGGG'
  
  # 添加UNV信息
  data$UNV <- substr(data$sequence, 1, 22)
  data$barcode <- substr(data$sequence, 23, 38)
  data$CTTGGG <- substr(data$sequence, 39, 44)
  data$UNV_ham_dis <- stringdist(data$UNV, UNV, method = 'hamming')
  data$UNV_lv_dis <- stringdist(data$UNV, UNV, method = 'lv')
  data$CTT_ham_dis <- stringdist(data$CTTGGG, CTT, method = 'hamming')
  data$CTT_lv_dis <- stringdist(data$CTTGGG, CTT, method = 'lv')
  
  
  # 保留查到VDJ的序列：
  data <- data[data$v_call != '' & data$j_call != '' & data$cdr3 != '' & data$locus == 'IGH', ]
  data <- data[rowSums(is.na(data)) != ncol(data), ]
  
  #######################
  print(file_name)
  print(paste0("比对到的序列数量(v,j,cdr3)为：", nrow(data)))
  #######################
  
  prod <- data[data$productive == 'TRUE', ]
  #######################
  print(paste0("productive的序列有：", nrow(prod)))
  write.table(prod, paste0('productive/', file_name), sep = '\t', row.names = F, quote = F)
  #######################
  
  unprod <- data[data$productive == 'FALSE', ]
  #######################
  print(paste0("unproductive的序列有：", nrow(unprod)))
  write.table(unprod, paste0('unproductive/', file_name), sep = '\t', row.names = F, quote = F)
  #######################
}

rm(data, prod, unprod)


setwd('E:/HJY_BCRrep_0628/analysis/HJY_BCRrep_QC/productive/remain/productive')
file_list <- list.files(path = './', pattern = '\\.process.tsv$', full.names = F)
for (file_name in file_list){
  data <- read.table(file_name, sep = '\t', header = T)
  data$UNV <- substr(data$sequence, 1, 19)
  UNV <- 'AGAGACAGATTGCGCAATG'
  data$UNV_ham_dis <- stringdist(data$UNV, UNV, method = 'hamming')
  data$UNV_lv_dis <- stringdist(data$UNV, UNV, method = 'lv')
  data$barcode <- substr(data$sequence, 20, 29)
  
  data$UMI <- data$barcode
  data <- data[data$UNV_lv_dis <= 5, ]
  print(paste0("分配到UMI的序列有：", nrow(data)))
  # data$UMI[data$UNV_lv_dis > 5] <- NA
  data$UMI[data$UNV_ham_dis > data$UNV_lv_dis] <- NA
  data$UMI[is.na(data$UMI)] <- str_match(data$sequence[is.na(data$UMI)], "GCAATG(.{10})")[,2]
  print(paste0("没有获得UMI的序列有:", sum(is.na(data$UMI))))
  
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
