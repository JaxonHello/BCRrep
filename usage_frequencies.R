library(dplyr)
library(ggplot2)

setwd("C:/Users/Chen-Lab/Desktop/Program/BCR_rep/IgBlast_result/Parse_IgBlast/Productive/no_pcr/CE")
file <- "Parse_A2_TTGCTTAGTC_L001_001.txt"
# read the table file, header=T means to read headers.
data <- read.delim(file, header = T, stringsAsFactors = F)
data <- data[,c("VH_top_match", "copies")]

clone <- data %>%
  group_by(VH_top_match) %>%
  # summarize函数用于在每个分组内计算统计信息，它将每个分组的copies列的和放在repeats列
  # calculate the sum of clone counts
  summarize(count_use=sum(copies), type_use=length(copies)) %>%
  mutate(count_freq=count_use/sum(count_use)*100,  type_freq=type_use/sum(type_use)*100) %>%
  arrange(desc(count_freq))

plot_bar <- ggplot(clone, aes(x=reorder(VH_top_match, -count_freq), y=count_freq)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "gene usage frequencies", x = "VH", y = "Percentage(%)") +
  theme(axis.text.x = element_text(angle = 90))

p_line <- ggplot(clone, aes(x = VH_top_match, group = attribute_column)) + 
  geom_line(aes(y = count_freq, color = "Count Frequency")) +
  geom_line(aes(y = type_freq, color = "Type Frequency")) +
  labs(title = "gene usage frequencies", x = "VH", y = "Percentage(%)") +
  scale_color_manual(values = c("Count Frequency" = "blue", "Type Frequency" = "red")) +
  theme(axis.text.x = element_text(angle = 90))

# write.table(clone, "C:/Users/Chen-Lab/Desktop/colon.txt", sep = "\t")
# remove the variables
# file <- "C:/Users/Chen-Lab/Desktop/clone.txt"
rm(data, file)
