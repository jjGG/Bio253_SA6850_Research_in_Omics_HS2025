

library(dplyr)
library(ggplot2)
library(tidyr)

#look at SB0804

methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")
View(methylation_data_EGG)

#filter feature, group and IPD ratio

onlym6A <- methylation_data_EGG %>%
    filter(feature == "m6A")


onlySB0804 <- filter(onlym6A, group == "SB0804")
onlySB0804 <- (filter(onlySB0804, IPDRatio > 2.5))  #filter for only positive IPDRatio

table(onlySB0804$treatment)


#Create a summary table to see motif presence across treatments

motif_presence <- onlySB0804 %>%
    group_by(start, motif, treatment, description, strand, category, ID, proteinID) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(present = n > 0) %>%
    dplyr::select(-n) %>%
    pivot_wider(
        names_from = treatment,
        values_from = present,
        values_fill = FALSE
    )

View(motif_presence)

table(motif_presence$PASN, motif_presence$TSB)

#Filter for motifs present only in one treatment

PASN_only0804 <- motif_presence %>%
    filter(PASN == TRUE, TSB == FALSE)

View(PASN_only0804)

TSB_only0804 <- motif_presence %>%
    filter(PASN == FALSE, TSB == TRUE)

View(TSB_only0804)



#look at the motifs in the database of SA





#look at 6850

library(dplyr)
library(ggplot2)
library(tidyr)


methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")
View(methylation_data_EGG)

#filter feature, group and IPD ratio

onlym6A <- methylation_data_EGG %>%
    filter(feature == "m6A")


only6850<- filter(onlym6A, group == 6850)
only6850 <- (filter(only6850, IPDRatio > 2.5))  #filter for only positive IPDRatio

table(only6850$treatment)


#Create a summary table to see motif presence across treatments

motif_presence <- only6850 %>%
    group_by(start, motif, treatment, description, strand, category, ID, proteinID) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(present = n > 0) %>%
    dplyr::select(-n) %>%
    pivot_wider(
        names_from = treatment,
        values_from = present,
        values_fill = FALSE
    )

View(motif_presence)

table(motif_presence$PASN, motif_presence$TSB)

#Filter for motifs present only in one treatment

PASN_only6850 <- motif_presence %>%
    filter(PASN == TRUE, TSB == FALSE)

View(PASN_only6850)

TSB_only6850 <- motif_presence %>%
    filter(PASN == FALSE, TSB == TRUE)

View(TSB_only6850)







#look at SB1002

library(dplyr)
library(ggplot2)
library(tidyr)


methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")
View(methylation_data_EGG)

#filter feature, group and IPD ratio

onlym6A <- methylation_data_EGG %>%
    filter(feature == "m6A")


onlySB1002 <- filter(onlym6A, group == "SB1002")
onlySB1002 <- (filter(onlySB1002, IPDRatio > 2.5))  #filter for only positive IPDRatio

table(onlySB1002$treatment)


#Create a summary table to see motif presence across treatments

motif_presence <- onlySB1002 %>%
    group_by(start, motif, treatment, description, strand, category, ID, proteinID) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(present = n > 0) %>%
    dplyr::select(-n) %>%
    pivot_wider(
        names_from = treatment,
        values_from = present,
        values_fill = FALSE
    )

View(motif_presence)

table(motif_presence$PASN, motif_presence$TSB)

#Filter for motifs present only in one treatment

PASN_only1002 <- motif_presence %>%
    filter(PASN == TRUE, TSB == FALSE)

View(PASN_only1002)

TSB_only1002 <- motif_presence %>%
    filter(PASN == FALSE, TSB == TRUE)

View(TSB_only1002)




# try to visualize differences in methylation between treatments for each clone

library(dplyr)
library(tidyr)
library(ggplot2)

methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")

motif_status_for_group <- function(dat, group_val, group_label = NULL) {
    if (is.null(group_label)) group_label <- as.character(group_val)

    dat %>%
        dplyr::filter(
            feature == "m6A",
            group == group_val,
            IPDRatio > 2.5
        ) %>%
        dplyr::group_by(start, motif, treatment, description, strand, category, ID, proteinID) %>%
        dplyr::summarise(n = n(), .groups = "drop") %>%
        dplyr::mutate(present = n > 0) %>%
        dplyr::select(-n) %>%
        tidyr::pivot_wider(
            names_from  = treatment,
            values_from = present,
            values_fill = FALSE
        ) %>%
        dplyr::mutate(
            status = dplyr::case_when(
                PASN & !TSB ~ "PASN",
                !PASN & TSB ~ "TSB",
                TRUE        ~ "None"     # shared or absent discarded
            ),
            group = group_label
        ) %>%
        dplyr::filter(status != "None")   # keep only PASN_only and TSB_only
}

motif_status_all <- dplyr::bind_rows(
    motif_status_for_group(methylation_data_EGG, "SB0804", "SB0804"),
    motif_status_for_group(methylation_data_EGG, 6850,       "6850"),
    motif_status_for_group(methylation_data_EGG, "SB1002",   "SB1002")
)

motif_counts <- motif_status_all %>%
    dplyr::count(group, status)

ggplot(motif_counts, aes(x = group, y = n, fill = status)) +
    geom_col(position = "dodge") +
    labs(
        x = "Group",
        y = "Number of diferences in methylated sites",
        fill = "Treatment",
        title = "Differences in Treatment-specific m6A sites",
    ) +
    theme_minimal(base_size = 14)



#visualiz edifferences in methylatioon only in Genes


library(dplyr)
library(tidyr)
library(ggplot2)

methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")

motif_status_for_group <- function(dat, group_val, group_label = NULL) {
    if (is.null(group_label)) group_label <- as.character(group_val)

    dat %>%
        dplyr::filter(
            feature == "m6A",
            group == group_val,
            category == "gene",
            IPDRatio > 2.5
        ) %>%
        dplyr::group_by(start, motif, treatment, description, strand, category, ID, proteinID) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(present = n > 0) %>%
        dplyr::select(-n) %>%   # <- now explicitly dplyr::select
        tidyr::pivot_wider(
            names_from  = treatment,
            values_from = present,
            values_fill = FALSE
        ) %>%
        dplyr::mutate(
            status = dplyr::case_when(
                PASN & !TSB ~ "PASN",
                !PASN & TSB ~ "TSB",
                TRUE        ~ "None"   # shared/absent discarded
            ),
            group = group_label
        ) %>%
        dplyr::filter(status != "None")   # keep only PASN and TSB
}

# Combine all groups
motif_status_all <- dplyr::bind_rows(
    motif_status_for_group(methylation_data_EGG, "SB0804", "SB0804"),
    motif_status_for_group(methylation_data_EGG, 6850,       "6850"),
    motif_status_for_group(methylation_data_EGG, "SB1002",   "SB1002")
)

# Count per group and status, and force missing combos to 0

motif_counts <- motif_status_all %>%
    dplyr::count(group, status) %>%
    tidyr::complete(
        group  = c("SB0804", "6850", "SB1002"),
        status = c("PASN", "TSB"),
        fill   = list(n = 0)
    )

motif_counts$status <- factor(motif_counts$status, levels = c("PASN", "TSB"))



ggplot(motif_counts, aes(x = group, y = n, fill = status)) +
    geom_col(position = "dodge") +
    scale_fill_manual(
        values = c(
            "PASN" = "forestgreen",   # Grün
            "TSB"  = "steelblue"      # Blau
        )
    ) +
    labs(
        x = "Group",
        y = "Number of differences in methylated sites",
        fill = "Treatment",
        title = "     Differences in m6A modification in Genes Only"
    ) +
    theme_minimal(base_size = 14)
