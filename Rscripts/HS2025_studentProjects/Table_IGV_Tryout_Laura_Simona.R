# Create a table with the methylations sites to use in IGV

# Clear R's brain
rm(list=ls())

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)


#import methylation file
data <- readRDS("../resources/methylation_data.rds")

# take only the necessary colomns
jlr_object <- data[c(8,3,4,5,20, 16)]

# create a colomn with the end position
jlr_object$end <- "-"
jlr_object$end <- jlr_object$start

# names for the methylations
jlr_object$name <- paste(data$Name,data$start,data$feature,data$IPDRatio,sep="_")

# use IPDRatio as score
jlr_object$score <- data$IPDRatio

# filter for m6A and IPD > 2.8
jlr_object <- filter(jlr_object, jlr_object$feature=="m6A")
jlr_object <- filter(jlr_object, jlr_object$IPDRatio > 2.8)

#look at IPD ratio
ggplot(data = jlr_object, aes(x = IPDRatio)) + geom_histogram(binwidth = 0.1)

# data file for IGV
jlr_object_2 <- jlr_object[c(1,2,7,8,9,3)]






# save the table as a file
write.table(jlr_object_2,file="methylated_sites_m6A_IPD2.8_tryout.bed",quote = F, row.names = F, col.names = F, sep="\t")


