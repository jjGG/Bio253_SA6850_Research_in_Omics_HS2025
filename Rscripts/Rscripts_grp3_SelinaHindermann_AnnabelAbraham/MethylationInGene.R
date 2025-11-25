
#look at 6850 only

library(dplyr)
library(ggplot2)
library(tidyr)


methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")
View(methylation_data_EGG)

#look at methylationsites in gene region of RSAU_001268

onlym6A <- methylation_data_EGG %>%
    filter(feature == "m6A")


only6850<- filter(onlym6A, group == 6850)
only6850 <- (filter(only6850, IPDRatio > 2.5))  #filter for only positive IPDRatio
only6850 <- filter(only6850, start > 1368000 & start < 1380000) #focus on a specific region
only6850 <- filter(only6850, strand == "-") #focus on replicate 1

table(only6850$treatment)


only6850 <- select(only6850, start, treatment, motif, description, strand, category, ID, proteinID)

View(only6850)




# look at methylation sites in gene region of RSAU_001770



onlym6A <- methylation_data_EGG %>%
    filter(feature == "m6A")


only6850_2<- filter(onlym6A, group == 6850)
only6850_2 <- (filter(only6850_2, IPDRatio > 2.5))  #filter for only positive IPDRatio
only6850_2 <- filter(only6850_2, start > 1892000 & start < 1898000) #focus on a specific region
only6850_2 <- filter(only6850_2, strand == "+")

table(only6850_2$treatment)


only6850_2 <- select(only6850_2, start, treatment, motif, description, strand, category, ID, proteinID)

View(only6850_2)



