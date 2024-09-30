# this program is to write the command line to file for merge use


# locate the folder containing the paired-end read file.
data_path <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/pair-end_read/test_data"
setwd(data_path)
# distinguish read1 & read2
read_1_files <- Sys.glob("*R1*")
read_2_files <- Sys.glob("*R2*")

# specify the output address
output_txt <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/merge_script"
# clear the file before
if (file.exists(output_txt)){
  file.remove(output_txt)
}

# write the command line "cd" to the output address
write(paste("cd", data_path), file = output_txt, append = TRUE)
# generate all the command line once
for (i in 1:length(read_1_files)){
  # generate file name of merged file
  merge_name <- sub("(.*)_R1(.*).fastq", "\\1\\2", read_1_files[i])
  # add the merge filter condition
  merge_condition <- " -v 20 -n 80 -q 20 -j 4"
  # generate the command line of two paired-end read
  # !注意：这里read1和read2应该是一一对应的，否则会报错
  merge_command <- paste("pear", "-f", read_1_files[i], "-r", read_2_files[i], 
                         "-o", paste0("merged_read/", merge_name, ".fastq"),  
                         merge_condition)
  # write these command lines to txt file
  write(merge_command, file = output_txt, append = TRUE)
}

# remove variable and release memory
rm(i, read_1_files, read_2_files, merge_name, merge_condition, data_path, output_txt, merge_command)
