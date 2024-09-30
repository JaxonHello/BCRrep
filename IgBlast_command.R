# this program is to write the command line for IgBlast


# locate the tool file path
tool_path <- "/Users/jaxonhe/Desktop/BCR_rep/igblast_tool"
# locate the data file path
data_path <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/filter_read"
# specify the output address
output_txt <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/IgBlast_script"
# clear the file before
if (file.exists(output_txt)){
  file.remove(output_txt)
}

# write the command line "cd" to the output address
write(paste("cd", tool_path), file = output_txt, append = TRUE)

# loop all the fasta file in folder
file_list <- list.files(data_path)
for (file in file_list){
  # specify output name
  output_name <- sub("(.*).fasta", "\\1", file)
  # specify the database name
  db_name <- "mouse/mouse_gl"
  # generate the command line for IgBlast
  command <- paste("bin/igblastn -germline_db_V", paste0("database/",  db_name, "_V.fasta"), 
                   "-germline_db_D", paste0("database/",  db_name, "_D.fasta"), "-germline_db_J", 
                   paste0("database/",  db_name, "_J.fasta"), "-organism mouse", "-query", 
                   paste0("myseq/", file), "-auxiliary_data optional_file/mouse_gl.aux",
                   paste0("-show_translation > IgBlast_result/", output_name, ".txt"), "-outfmt 3")
  # write it to txt file
  write(command, file = output_txt, append = TRUE)
}

rm(command, data_path, db_name, file, file_list, output_name, output_txt, tool_path)