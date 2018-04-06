# This script takes in a csv file of the format cis-element, barcode, UMI, [spliced isoform]
# And output a csv file containing:
#        (1) the number of inclusion, exclusion, cryptic, total sequence number associated with each 
#            pair of cis-element and barcode
#        (2) the PSI value for each sequence

# import libraries
library(data.table)
library(reshape2)

############# FILL IN THE INFORMATION BELOW #######################
working_directory_path <- "~/Desktop/Github/ROP299"
library_path <- "Output/final_library.csv"
spliced_result_path <- "Output/120_junction_cryptic.csv"
result_path <- "Output/120_junciton_cryptic_result.csv"
lib_count<-F
####################################################################

# import dataset
# choose the csv file that stores the parsed library
setwd(working_directory_path)
lib <- fread(library_path)
# choose the csv file containing spliced parsed splicing result
spliced <- fread(spliced_result_path)


# clean library data table
if(lib_count) {
  lib <- lib[,1:4]
  colnames_final <- c("id", "cis-element", "barcode", "cis_count","exclusion", "inclusion", "cryptic", "unknown")
}else {
  lib_count <- lib[,1:3]
  colnames_final <- c("id", "cis-element", "barcode","exclusion", "inclusion", "cryptic", "unknown")
}

# cast the spliced output table
colnames(spliced) <- c("barcode", "isoform", "UMI")
uniq_spliced <- unique(spliced)
evaluated <- dcast(uniq_spliced[, .N, by=c('barcode', 'isoform')], barcode ~ isoform)
evaluated[is.na(evaluated)] <- 0

raw_spliced_size <- nrow(spliced)
uniq_spliced_size <- nrow(uniq_spliced)
raw_barcode_num <- nrow(evaluated)

# merge two data sets
merged <- merge(lib, evaluated, by='barcode', sort=F)
merged
merged[,c(1,2,3)] <- merged[,c(2,3,1)]
colnames(merged)
colnames(merged) <- colnames_final

merged_barcode_num <- nrow(merged)

# calculate the total and PSI
merged$total <- apply(merged[,4:7], 1, sum)
merged$PSI <- round(merged$inclusion * 100 / (merged$inclusion + merged$exclusion), digits = 2)

basic_stat <-data.frame(raw_spliced_size, uniq_spliced_size, raw_barcode_num, merged_barcode_num)
print(basic_stat)
write.csv(merged, result_path)

