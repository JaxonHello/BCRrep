#enter the file location
dir_files <- "C:/Users/Chen-Lab/Desktop/Program/BCR_rep/IgBlast_result/Parse_IgBlast"
setwd(dir_files)
#read the file names
files <- Sys.glob("*.txt")
#create separate folders for each of the group
if(!file.exists("Productive"))
{ dir.create("Productive")}

if(!file.exists("Non_Productive"))
{ dir.create("Non_Productive")}


for(i in 1: length(files)){
  print(files[i])
  Data <- read.delim(files[i],header=T,stringsAsFactors=F)
  #separate and write the productive sequences
  Data_productive <- Data[which(Data$Productive == "Yes" & Data$V.Jframe == "In-frame"),]
  out_productive <- paste("Productive/",files[i], sep ="")
  write.table(Data_productive, out_productive , sep="\t")
  #separate and write the non-productive sequences
  Data_non_productive <- Data[which(Data$Productive == "No"),]
  out_non_productive <- paste("Non_Productive/",files[i], sep ="")
  write.table(Data_non_productive, out_non_productive , sep="\t")
  print("结束")
}

rm(Data, Data_non_productive, Data_productive, dir_files, files, i, out_non_productive, out_productive)