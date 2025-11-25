
library(dplyr)
library(ggplot2)
library(tidyr)


methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")


#visualize the different features in color

ggplot(df, aes(
    x = IPDRatio,
    y = feature,
    fill = feature
)) +
    geom_density_ridges(
        scale = 0.7,
        rel_min_height = 0,
        alpha = 0.5,
        color = "black"
    ) +
    labs(
        title = "Hill-like Visualization of IPDRatio by Modification Type",
        x = "IPDRatio",
        y = "Modification Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank()
    )




#overlapping plot







ggplot(df, aes(
    x = IPDRatio,
    color = feature,
    fill = feature
)) +
    geom_density(alpha = 0.4, size = 1) +
    labs(
        title = "Overlapping Density Curves of IPDRatio by Modification Type",
        x = "IPDRatio",
        y = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        legend.position = "right",
        panel.grid.minor = element_blank()
    )




#only lok at m6A modifications with IPDRatio > 2.5


df <- methylation_data_EGG %>%
    filter(!is.na(IPDRatio), IPDRatio > 2.5)

ggplot(df, aes(x = IPDRatio)) +
    geom_density(
        aes(color = "m6A", fill = "m6A"),
        alpha = 0.4,
        size  = 1.0
    ) +
    scale_color_manual(values = c("m6A" = "darkgreen")) +
    scale_fill_manual(values  = c("m6A" = "lightgreen")) +
    scale_x_continuous(
        limits = c(2.5, NA),           # start at 2.5
        breaks = seq(2.5, max(df$IPDRatio, na.rm = TRUE), by = 0.5)
    ) +
    labs(
        title = "IPDRatio Density Curve of m6A Modifications",
        x = "IPDRatio",
        y = "Density",
        color = "Feature",
        fill  = "Feature"
    ) +
    theme_minimal(base_size = 14)


#plot median IPDRatio per gene per clone (group) and treatment as boxplots with only m6A sites


library(dplyr)
library(ggplot2)
library(tidyr)


methylation_data_EGG <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")
methylation_data_EGG <- filter(methylation_data_EGG, IPDRatio > 2.5)

library(dplyr)
library(ggplot2)

# 1. Filter to only m6A sites and keep rows with a gene + IPDRatio
df_m6A <- methylation_data_EGG %>%
    filter(feature == "m6A",
           !is.na(IPDRatio),
           !is.na(gene))

# 2. Median IPDRatio per gene per clone (group) and treatment
m6A_gene_medians <- df_m6A %>%
    group_by(group, treatment, gene) %>%
    summarise(
        median_IPDRatio = median(IPDRatio, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    # Extract the strain from clone name (e.g., "6850", "SB0804", "SB1002")
    mutate(strain = sub("_.*", "", group))

# 3. Boxplot for m6A methylation

ggplot(m6A_gene_medians,
       aes(x = group, y = median_IPDRatio, fill = treatment)) +
    geom_boxplot(outlier.size = 0.4) +
    scale_fill_manual(values = c("forestgreen", "steelblue")) +
    facet_wrap(~ strain, scales = "free_x") +
    labs(
        title = "                    Methylation IPDRatios (m6A Only)",
        x = "",
        y = "Average IPDRatio per Gene"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank()
    )




