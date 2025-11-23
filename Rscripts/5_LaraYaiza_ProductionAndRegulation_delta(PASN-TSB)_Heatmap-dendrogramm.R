###BIO253-OMICS - Proteomix-STX-Production in Staphylococus aureus (SA)
###STX- PRODUCTION AND REGULATORY PRTX
###Comparing the difference in Protein Expression within Clones
###by calculating the delat (PASN-TSB) --> and generate a Heatmap + Dendogramm on these changes

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###TODO
### - change titel (better, more intuitive)
### - check axes labelling
### - separate the two strains --> generate 2 separate heatmaps for each strain!
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(readxl)
library(pheatmap)
library(grid)
library(gridExtra)

#### Daten einlesen
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")       # relative folder path!
file <- "resources/LaraYaiza_tabel_STX.xlsx"             # Dateiname Resource
#### richtige Tabelle einlesen
df <- read_excel(file,sheet = "Tabelle2",skip=3)

###Spaltennamen überprüfen:
#names(df)

#Missing values in numeric values umwandlen! --> Missin gmit NA ersetzen!
df[df == "Missing"] <- NA
numeric_cols <- c("crtM", "crtN", "crtP", "crtQ", "crtO",
                  "STX Average", "Growth rate", "H2O2 survival",
                  "sigma-B", "rsbU", "rsbV", "rsbW")

df <- df %>%
    mutate(across(all_of(numeric_cols), ~as.numeric(.)))

##data strucutre check
#str(df)
#head(df)

################################################################################
#### delta = PASN-TSB berechenen für STX-PRODUCTION & STX-REGULATION PrtX ######
################################################################################
delta_df <- df %>%
    pivot_wider(
        names_from = Condition,
        values_from = c(crtM, crtN, crtP, crtQ, crtO,
                        `sigma-B`, rsbU, rsbV, rsbW,
                        `Growth rate`)
    ) %>%
    mutate(
        delta_crtM = crtM_PASN - crtM_TSB,
        delta_crtN = crtN_PASN - crtN_TSB,
        delta_crtP = crtP_PASN - crtP_TSB,
        delta_crtQ = crtQ_PASN - crtQ_TSB,
        delta_crtO = crtO_PASN - crtO_TSB,

        delta_sigmaB = `sigma-B_PASN` - `sigma-B_TSB`,
        delta_rsbU   = rsbU_PASN - rsbU_TSB,
        delta_rsbV   = rsbV_PASN - rsbV_TSB,
        delta_rsbW   = rsbW_PASN - rsbW_TSB,

        delta_growth = `Growth rate_PASN` - `Growth rate_TSB`
    )

### EINDEUTIGE Zeilen-ID erzeugen: Strain_Clone
delta_df <- delta_df %>%
    mutate(Clone_ID = paste(Strain, Clones, sep = "_"))

###Heatmap-Matrix erstellen
heat_data <- delta_df %>%
    select(Clone_ID, starts_with("delta_")) %>%
    column_to_rownames("Clone_ID")

###Spalten mit NA entfernen --> sonst HeatMap verfälscht!????
heat_data_clean <- heat_data %>%
    select(where(~!any(is.na(.))))
##CHECK WHICH GOT "DELETED"
#print("Entfernte Spalten ween NA:")
#print(setdiff(colnames(heat_data),colnames(heat_data_clean)))
###crtQ, crtO, rsbU --> keine delta berechenbar!

###Heatmap erzeugen
pheatmap(
    heat_data_clean,
    scale = "row",
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 11,
    main = "Δ Proteomics + Growth rate (PASN – TSB)"
)


################################################################
###WITHOUT GROWTHRATE - ONLY GENES                           ###
################################################################
heat_data_clean_noGrowth <- heat_data_clean %>%
    select(-delta_growth)

pheatmap(
    heat_data_clean_noGrowth,
    scale = "row",
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 11,
    main = "Δ Proteomics (PASN – TSB, ohne Growth rate)"
)

#############################################################
###MACHT DAS SINN REGULATORY UND PRODUCTION IN EIN MAP???
###ODER SEPARATE PRODUCTION & REGULATORY MAP?

##############################################################
###ONLY STX_PRODUCTION GENES###
##############################################################
heat_stx <- heat_data_clean %>%
    select(delta_crtM, delta_crtN, delta_crtP)

pheatmap(
    heat_stx,
    scale = "row",
    main = "Δ STX-Biosynthesis Proteins"
)

############################################################
###ONLY STX-REGULATORY GENES###
###########################################################
###sigmaB, rsbU, rsbV, rsbW

heat_reg <- heat_data %>%
    select(delta_sigmaB, delta_rsbU, delta_rsbV, delta_rsbW)

# NA-freie Spalten sicherstellen
heat_reg <- heat_reg %>%
    select(where(~ !any(is.na(.))))

pheatmap(
    heat_reg,
    scale = "row",
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 11,
    main = "Δ STX-Regulation (σB & rsb genes)"
)
#############################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#############################################################################
###Regulatory + Production delta(PASN-TSB) + STX-Average
#############################################################################

delta_df <- df %>%
    pivot_wider(
        names_from = Condition,
        values_from = c(crtM, crtN, crtP, crtQ, crtO,
                        `sigma-B`, rsbU, rsbV, rsbW,
                        `Growth rate`)
    ) %>%
    mutate(
        delta_crtM = crtM_PASN - crtM_TSB,
        delta_crtN = crtN_PASN - crtN_TSB,
        delta_crtP = crtP_PASN - crtP_TSB,
        delta_crtQ = crtQ_PASN - crtQ_TSB,
        delta_crtO = crtO_PASN - crtO_TSB,

        delta_sigmaB = `sigma-B_PASN` - `sigma-B_TSB`,
        delta_rsbU   = rsbU_PASN - rsbU_TSB,
        delta_rsbV   = rsbV_PASN - rsbV_TSB,
        delta_rsbW   = rsbW_PASN - rsbW_TSB,

        #delta_growth = `Growth rate_PASN` - `Growth rate_TSB`
    ) %>%
    mutate(Clone_ID = paste(Strain, Clones, sep = "_"))

heat_data <- delta_df %>%
    select(Clone_ID, starts_with("delta_")) %>%
    column_to_rownames("Clone_ID")

heat_data_clean <- heat_data %>%
    select(where(~!any(is.na(.))))

pheatmap(
    heat_data_clean,
    scale = "row",
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 12,
    main = "Protein Abundance Changes PASN vs. TSB"    #titel heatmap
)

### STX baseline heatmap: relative to ancestor per strain

stx_df <- df %>%
    distinct(Strain, Clones, `STX Average`) %>%
    mutate(Clone_ID = paste(Strain, Clones, sep = "_"))

# Calculate strain-wise ancestor reference
ancestor_values <- stx_df %>%
    filter(Clones == "Ancestor") %>%
    select(Strain, ancestor_stx = `STX Average`)

# Join ancestor reference and compute normalized STX
stx_rel <- stx_df %>%
    left_join(ancestor_values, by = "Strain") %>%
    mutate(STX_relative = `STX Average` / ancestor_stx) %>%
    select(Clone_ID, STX_relative) %>%
    column_to_rownames("Clone_ID")

# Heatmap (single column)
pheatmap(
    stx_rel,
    scale = "none",
    cluster_cols = FALSE,      # no column clustering
    cluster_rows = FALSE,      # no row clustering (this removes the dendrogramm!)
    fontsize = 12,
    border_color = NA,
    main = "Phenotypic STX (Relative to Ancestor)"
)
################################################################################
###STX-Average as Heatmap seems a bit an overkill, additional colour scale etc.
###How about as a horizontal number table
###################################################################################
library(tidyverse)
library(ggplot2)

### Desired clone order
clone_order <- c(
    "6850_Ancestor",
    "6850_B0403",
    "6850_B0804",
    "6850_B1002",
    "JE2_Ancestor",
    "JE2_B0604"
)

### Build dataframe
stx_tbl <- df %>%
    distinct(Strain, Clones, `STX Average`) %>%
    mutate(Clone_ID = paste(Strain, Clones, sep = "_")) %>%
    filter(Clone_ID %in% clone_order) %>%
    mutate(Clone_ID = factor(Clone_ID, levels = clone_order))

### PLOT — light grey boxes, no border, white grid spacing
ggplot(stx_tbl, aes(x = Clone_ID, y = 1, label = round(`STX Average`, 2))) +

    # light grey background squares
    geom_tile(fill = "grey95", color = "white", width = 0.95, height = 0.95) +

    # centered black numbers
    geom_text(size = 6, color = "black") +

    # formatting
    scale_y_continuous(expand = c(0,0), limits = c(0.5, 1.5)) +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold")
    ) +

    ggtitle("Baseline STX-Average per Clone")
