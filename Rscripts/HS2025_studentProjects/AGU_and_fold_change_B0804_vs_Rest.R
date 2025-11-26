
rm(list=ls())

library(readr)
library(dplyr)
library(stringr)

# Import dataset listing the fold changes in protein expression of B0804 respective to rest (anc,B0403 & B1002)
Samples_B0804_vs_rest <- read_csv("../resources/Samples_B0804_vs_rest.csv", skip = 10)

# select columns with protein_IDs & log2 fold change
sample_selected <- select(Samples_B0804_vs_rest, 'Protein Name', '0804')

# rename column '0804'
sample_selected <- sample_selected %>% rename('log2(fold change)' = '0804')

#Extract AGU number from the 'Protein Name' column
sample_selected$AGU <- str_extract(sample_selected$`Protein Name`, "(?<=\\[protein_id=)[^\\]]+")

# check NAs
sum(is.na(Samples_B0804_vs_rest$`Protein Name`)) # 1 NA
sum(is.na(sample_selected$AGU)) # 135 NAs
which(is.na(sample_selected$AGU))
View(sample_selected[which(is.na(sample_selected$AGU)),])
# -> NAs are Proteins without AGU in Protein Name column

# exclude NAs from sample_selected, cause we only want proteins in df that are identifiable by AGU
sample_selected <- sample_selected[!is.na(sample_selected$AGU),]

sum(is.na(sample_selected)) # 0 -> all NAs removed

sample_selected <- sample_selected[,c('AGU','log2(fold change)')]

# fold changes in brackets indicate that there is at least 1 missing value among the clones for a given protein
# remove rows containing fold changes in brackets, cause they can't be trusted:
sample_selected <- sample_selected %>%
    filter(!grepl("^\\(.*\\)$", `log2(fold change)`))

# there still are rows containing entries "Missing Value" or "Reference Missing" in the column 'log2(fold change)'
# Remove rows containing the word "Missing":
sample_selected <- sample_selected %>%
    filter(!grepl("Missing", `log2(fold change)`))

# convert df to csv file -> upload file in under "Proteins with values/ranks" in STRING for enrichment analysis
write.csv(sample_selected, "AGU_and_fold_change_B0804_vs_Rest.csv", row.names = FALSE)
