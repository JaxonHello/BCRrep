# this program is to filter all the assembled fastq file among the PEAR result file


# locate folder containing the merge result.
data_path <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/merged_read"
setwd(data_path)
# filter all the assembled read file name
assembled_file <- Sys.glob("*assembled.fastq*")
# loop through all files in folder
file_list <- list.files(data_path)
for (file in file_list){
  # find the assembled file
  if (file %in% assembled_file){
    # keep only the assembled file
    # copy these reads to filter_read folder for filter
    file.copy(file, "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/filter_read")
  }
  else{
    # delete other result from PEAR
    file.remove(file)
  }
}


# remove variable and release memory
rm(assembled_file, data_path, file, file_list)