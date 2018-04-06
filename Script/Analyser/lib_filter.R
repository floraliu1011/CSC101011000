library(readr)
library(reshape2)
library(ggplot2)
library(data.table)



lib <- fread("~/Desktop/ROP299/phase2/phase2.2_data_analysis/library_sort.csv")
lib[,c(3,4)] <- NULL
lib <- na.omit(lib)
colnames(lib) <-c('cis', 'barcode')
lib <- lib[, .N, by = .(cis, barcode)] 
colnames(lib) <-c('cis', 'barcode', 'count')
separable <- function(bar_seq, lib) {
  df <- lib[barcode == bar_seq][1:2,]
  return ((df$count[1] >= 10 * df$count[2]) && df[2,]$count < 10 )
}
lib_50 <- lib[count > 30]
lib_separable <- lib_50[sapply(lib_50$barcode, separable, lib),]

write.csv()