# BCR rep analysis pipeline

#### Environment

```shell
source /home/yourname/.bashrc
conda activate BCRrep
cd /your_seq_dir/
```



#### Sequencing Qulity Control Report

```shell
mkdir fastqc_report

fastqc -o ./fastqc_report ./*.fq.gz
```

#### Paired-end Read Merge


```shell
mkdir pear_merge

for file in $(find ./ -type f -name "*_1.fq.gz"); do
filename=$(basename "$file")
filename_pre="${filename%_L1_1.fq.gz}"
echo "-------------------------"
echo "start merge：${filename_pre}"
n_r1=$(zgrep -c '^@' "${filename_pre}_L1_1.fq.gz")
echo "${filename_pre}_L1_1.fq.gz文件序列数为${n_r1}"
n_r2=$(zgrep -c '^@' "${filename_pre}_L1_2.fq.gz")
echo "${filename_pre}_L1_2.fq.gz文件序列数为${n_r2}"
pear -f ${filename_pre}_L1_1.fq.gz -r ${filename_pre}_L1_2.fq.gz -o pear_merge/${filename_pre}_m -v 20 -n 80 -q 20
echo "-------------------------"
done

cd pear_merge
for file in $(find ./ -type f -name "*.assembled.fastq"); do
filename=$(basename "$file")
n_seq=$(zgrep -c '^@' "${filename}")
echo "-------------------------"
echo "start merge：${filename} file has sequence number::${n_seq}"
echo "-------------------------"
done
```

```shell
for file in $(find ./ -type f -name "*_1.fq.gz"); do
filename=$(basename "$file")
filename_pre="${filename%_L1_1.fq.gz}"
echo "-------------------------"
echo "start merge：${filename_pre}"
n_r1=$(zgrep -c '^@' "${filename_pre}_L1_1.fq.gz")
echo "${filename_pre}_L1_1.fq.gz has sequence number: ${n_r1}"
n_r2=$(zgrep -c '^@' "${filename_pre}_L1_2.fq.gz")
echo "${filename_pre}_L1_2.fq.gz has sequence number: ${n_r2}"
pear -f ${filename_pre}_L1_1.fq.gz -r ${filename_pre}_L1_2.fq.gz -o pear_merge_new/${filename_pre}_m -v 20 -n 80
echo "-------------------------"
done
```



#### Quality Control

use `filter.R` scripts to process merge sequence， and do quality control for the low quality base。 the step is as below：

part 1: sequence correction: replace the intermittent low quality read (Q<20) to "N"

part 2: sequence correction: four continous "N" will be trimmed

part 3: remove all the fasta file whose whole "N" > 2%

part 4: convert the form to ".fasta"

```shell
R
```

```R
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
  print(paste0("QC前该文件的总读数为：", length(fastq_data)))
  
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
  print(paste0("序列经过filter总共费时为：", (difftime(end_time, start_time, units = "secs"))))
}


# find all the file with .fastq

file_list <- list.files(path = data_path, pattern = "\\.assembled.fastq$", full.names = TRUE)
for (file in file_list){
  print(paste0('-----------------', file, '-------------------'))
  # filter all the task
  to_fasta(file)
}

# remember to remove the variable
rm(data_path, file, file_list)
```

```shell
mkdir ./QC_fasta
mv *.fasta ./QC_fasta
cd QC_fasta

for file in *.fasta; do
filename=$(basename "$file")
n_seqs=$(grep -o ">" "${filename}" | wc -l)
echo "${filename}经过QC之后的序列有：${n_seqs}"
done
```

#### IgBLAST alignment

```shell
chmod +X ./IgBlast_process.sh
./test.sh -i myseq/HJY_BCRrep_QC/A1_L1_m.fasta -o result/HJY_BCRrep_QC/A1_L1_m.raw.tsv -p result/HJY_BCRrep_QC/A1_L1_m.pro.tsv
```

and for specific directory which has many fasta file, you can use the scripts below：

```shell
dir=Elisa_IgD
for file in $(find ./myseq/${dir} -type f -name "*.fasta"); do
filename=$(basename "$file")
filename_pre="${filename%.fasta}"
echo "--------------------------"
echo "开始查询::"${filename}
./IgBlast_process.sh -d human -i myseq/${dir}/${filename} -o result/${dir}/${filename_pre}.raw.tsv -p result/${dir}/${filename_pre}.process.tsv
echo "--------------------------"
done
```

* after alignment, you can use `has_vjcdr.R`, `addUMI.R`, to get the sepcific barcode and UMI for every read。

####  Get Complete VDJ

the read aligned successfully for all the V, D, cdr3 region is considered as a good BCR read

#### Remove PCR repeat

In the previous step, we have extracted UNV, barcode and short UNV sequence after barcode from the sequenced sequence according to our own library structure. Then, in the real data, due to the influence of sequencing quality and PCR process, the following will appear:

* The beginning is not UNV sequence
* UNV at the beginning and after will be deleted or inserted
* The barcode will be deleted or inserted, resulting in the length not meeting the designed library structure
* Sequencing errors or N in the barcode and UNV structure

For how to remove PCR amplification, the most critical issue is to assign the most appropriate barcode to each sequence.

#### Define Barcode id & cdr3 id

```shell
cd barcode_id

for file in $(find ./ -type f -name "*.tsv"); do
filename=$(basename "$file")
echo "--------------------------"
echo "start requiry::"${filename}
DefineClones.py -d ${filename} --outdir ./ --failed --act set --model ham --norm len --dist 0 --maxmiss 100
echo "--------------------------"
done
```

```shell
cd cdr3_id

for file in $(find ./ -type f -name "*.tsv"); do
filename=$(basename "$file")
echo "--------------------------"
echo "start requiry::"${filename}
DefineClones.py -d ${filename} --outdir ./ --failed --act set --model ham --norm len --dist 0.05 --maxmiss 100
echo "--------------------------"
done
```

#### Define Clones

when remove the PCR repeat, we can now define the clones for BCR
```shell
for file in $(find ./ -type f -name "*.tsv"); do
filename=$(basename "$file")
echo "--------------------------"
echo "start define clones::"${filename}
DefineClones.py -d ${filename} --outdir ./clone_expansion --failed --act set --model ham --norm len --dist 0.15 --maxmiss 100
echo "--------------------------"
done
```

#### DownStream Analysis
<hr>
1. VH gene usagement
2. AA length and hydrophobicity analysis
3. Diversity of VH usagement
4. JSD & KLD
5. VH mutation frequency
#### Draw Treemap

```R
# Treemap
library(treemap)
library(dplyr)

setwd('~/Desktop/clone_expansion_data0628')

file_list <- list.files(path = './', pattern = '\\.process_clone-pass.tsv$', full.names = F)
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
```

