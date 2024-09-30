# this file is to correct read, filter low quality read and convert form to fasta
# BiocManager::install("ShortRead") (downloading may take a lot of time)
library(ShortRead)
library(seqinr)

# part 1: sequence correction: replace the intermittent low quality read (Q<20) to "N"
# part 2: sequence correction: four continous "N" will be trimmed
# part 3: convert the form to ".fasta"
# part 4: remove all the fasta file whose whole "N" > 2%

# define the function read_process(read) to get the filtered read after part1, 2
read_correction <- function(read){
  # split the sequence into a vector
  Nucl_seq <- strsplit(as.character(sread(read)[[1]]), "")[[1]]
  # split the score into a vector
  Score_seq <- strsplit(as.character(quality(read)[[1]]), "")[[1]]
  # loop all the nucleotide with its score
  # origin continous "N" number is 0
  c <- 0
  for (j in 1: length(Score_seq)){
    # convert the ASCII char to Q score number
    n <- as.numeric(charToRaw(Score_seq[j])) - 33
    # trim four continous "N"
    if (c >= 4){
      j <- j-5
      break
    } else if (n<20){
      # find the nucleotide whose Q score < 20, replace it with "N"
      Nucl_seq[j] <- "N"
      c <- c+1
    } else if (n>21){
      c <- 0
    }
  }
  # get the reads after filter
  seq_filtered <- paste(Nucl_seq[1:j], collapse ="")
  return(seq_filtered)
}

# specify the fastq data path
data_path <- getwd()

# define the function to_fasta(file) to write the filter fastq to fasta
to_fasta <- function(file){
  
  # read the fastq file
  fastq_data <- readFastq(file)
  fasta_path <- paste0(sub("(.*).assembled.fastq", '\\1', file), '.fasta')
  print(paste0("经过质量控制前的该文件的总体读数为：", length(fastq_data)))
  
  start_time <- Sys.time()
  for (i in 1:length(fastq_data)){
    # process read with filter
    read_c <- read_correction(fastq_data[i])
    # split the read seq into vector list
    seq_list <- strsplit(read_c, "")[[1]]
    # judge whether the "N" propotion is less that 2%, which is filter part 4
    if (length(which(seq_list == "N"))< 0.02*length(seq_list)){
      write.fasta(sequences = read_c , 
                  names = fastq_data[i]@id[[1]], nbchar = 80, file.out = fasta_path, 
                  open = "a",as.string = TRUE)
    }
  }
  end_time <- Sys.time()
  print(paste0("序列经过filter总共费时为：", (end_time-start_time)))
}


# find all the file with .fastq

file_list <- list.files(path = data_path, pattern = "\\.assembled.fastq$", full.names = TRUE)
for (file in file_list){
  print(file)
  # filter all the task
  to_fasta(file)
}

# remember to remove the variable
rm(data_path, file, file_list)


