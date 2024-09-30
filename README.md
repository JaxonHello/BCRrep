# BCR rep分析流程

#### 环境配置

```shell
source /home/yzchen_pkuhpc/profile/hjx/.bashrc
conda activate BCRrep
cd /lustre2/yzchen_pkuhpc/profile/hjx/data/HJY_BCRrep_0704/N2413673_80-1550858036_2024-06-30/MS240629-0126/
```



#### 测序质量报告

```shell
mkdir fastqc_report

fastqc -o ./fastqc_report ./*.fq.gz
```

#### 双端merge

由于之前neha处理的时候是通过R语言对每一个文件生成了命令行，但是其实shell里面可以直接做循环的。

```shell
mkdir pear_merge

for file in $(find ./ -type f -name "*_1.fq.gz"); do
filename=$(basename "$file")
filename_pre="${filename%_L1_1.fq.gz}"
echo "-------------------------"
echo "开始merge：${filename_pre}"
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
# 这里的后缀在实际操作的过程中有问题
n_seq=$(zgrep -c '^@' "${filename}")
echo "-------------------------"
echo "开始merge：${filename}的序列有::${n_seq}"
echo "-------------------------"
done
```

```shell
for file in $(find ./ -type f -name "*_1.fq.gz"); do
filename=$(basename "$file")
filename_pre="${filename%_L1_1.fq.gz}"
echo "-------------------------"
echo "开始merge：${filename_pre}"
n_r1=$(zgrep -c '^@' "${filename_pre}_L1_1.fq.gz")
echo "${filename_pre}_L1_1.fq.gz文件序列数为${n_r1}"
n_r2=$(zgrep -c '^@' "${filename_pre}_L1_2.fq.gz")
echo "${filename_pre}_L1_2.fq.gz文件序列数为${n_r2}"
# 由于发现了一些问题，现在取消q 为20 的参数，看看merge的效率。
pear -f ${filename_pre}_L1_1.fq.gz -r ${filename_pre}_L1_2.fq.gz -o pear_merge_new/${filename_pre}_m -v 20 -n 80
echo "-------------------------"
done
```



#### 质量控制

使用`filter.R`脚本处理双端merge序列，对序列进行质量控制(QC)。质量控制的步骤有：

part 1: sequence correction: replace the intermittent low quality read (Q<20) to "N"

part 2: sequence correction: four continous "N" will be trimmed

part 3: remove all the fasta file whose whole "N" > 2%

part 4: convert the form to ".fasta"

由于目前的程序之后结果会在pear_merge文件夹中，所以在R语言代码处理后，还需要将文件移动到QC文件夹中。

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

#### IgBLAST比对

使用的IgBlast是经过改进之后的工具，目前该工具只能查询人类的序列。在igblast工具的子目录下面，它对于单个文件查询的shell命令如下：
```shell
chmod +X ./IgBlast_process.sh
./test.sh -i myseq/HJY_BCRrep_QC/A1_L1_m.fasta -o result/HJY_BCRrep_QC/A1_L1_m.raw.tsv -p result/HJY_BCRrep_QC/A1_L1_m.pro.tsv
```

对于整个文件夹内的所有序列文件，需要通过循环结果查询：

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

* 比对之后的分析顺序是`has_vjcdr.R`, `addUMI.R`,  接着去对应的文件去获取对应的cdr3_id和barcode_id。

#### 比对结果分析

该结果的分析从hpc转到了实验室的server上。

在HJY0628的测序结果中，通过同时有`v`, `j`, `cdr3`的

#### 去除PCR repeat

在上一步我们已经根据我们自己的文库结构从测序的序列中提取出了UNV，barcode和barcode之后的短的UNV序列，然后在真实的数据中，由于测序质量和PCR过程的影响，会出现：

* 开头不是UNV序列
* 开头和之后的UNV会出现删除和插入
* barcode出现删除和插入导致长度不符合设计的文库结构
* barcode和UNV结构中出现测序错误或者N的情况

对于如何去除PCR扩增，最关键的问题就是给每条序列都分配一个最合适的barcode。

#### 定义barcode_id和cdr3_id

```shell
cd barcode_id

for file in $(find ./ -type f -name "*.tsv"); do
filename=$(basename "$file")
echo "--------------------------"
echo "开始查询::"${filename}
DefineClones.py -d ${filename} --outdir ./ --failed --act set --model ham --norm len --dist 0 --maxmiss 100
echo "--------------------------"
done
```

```shell
cd cdr3_id

for file in $(find ./ -type f -name "*.tsv"); do
filename=$(basename "$file")
echo "--------------------------"
echo "开始查询::"${filename}
DefineClones.py -d ${filename} --outdir ./ --failed --act set --model ham --norm len --dist 0.05 --maxmiss 100
echo "--------------------------"
done
```

#### 定义克隆

对于去除PCR repeat之后的序列，我们需要对其进行克隆定义，克隆定义的标准如下：
```shell
for file in $(find ./ -type f -name "*.tsv"); do
filename=$(basename "$file")
echo "--------------------------"
echo "开始定义克隆::"${filename}
DefineClones.py -d ${filename} --outdir ./clone_expansion --failed --act set --model ham --norm len --dist 0.15 --maxmiss 100
echo "--------------------------"
done
```

#### 绘制Treemap

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

#### DownStream Analysis

1. VH gene usagement
2. AA length and hydrophobicity analysis
3. Diversity of VH usagement
4. JSD & KLD
5. VH mutation frequency


