# to calculate read numbers of different process

# 1. read 1
# 2. merged read
# 3. filter read
# 4. igblast result
# 5. productive
# 6. repeat
# 7. colonal expansion

library(ShortRead)
library(seqinr)

# calculate the number of read 1.
file_path <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/pair-end_read/test_data"
all_file <- list.files(file_path, full.names = TRUE)
R1_file <- all_file[grep("R1", all_file)]

# calculate the number of merged read.
file_path <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/merged_read"
all_file <- list.files(file_path, full.names = TRUE)

# calculate the number of filter read.
file_path <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/filter_read"
all_file <- list.files(file_path, full.names = TRUE)
  #for (i in 1:length(all_file)){
   # print(all_file[i])
   # print(length(readFasta(all_file[i])))
  #}
# calculate the number of productive
file_path <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/IgBlast_result/Parse_IgBlast/Productive/no_pcr"
all_file <- list.files(file_path, full.names = TRUE)
txt_file <- all_file[grep(".txt", all_file)]
for (i in 1:length(txt_file)){
  print(txt_file[i])
  print(length(readLines(txt_file[i]))-1)
}
# calculate the number of after-repear

# calculate the number of after-colonal expansion