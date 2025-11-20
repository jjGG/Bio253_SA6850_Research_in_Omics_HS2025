#### Omics BIO253 - Proteomic STX
#### Matrix with numbers and MULTIPLE ARROW Matrix
###

library(tidyverse)
library(readxl)
library(writexl)

#### Daten einlesen
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")
#path <- "..Bio253_SA6850_Reseach_in_Omics_HS2025/resources"        # relative folder path!
file <- "resources/LaraYaiza_tabel_STX.xlsx"                        # Dateiname Resource
# full_path <- file.path(path, file)
# print(full_path)
# print(file.exists(full_path))

#### richtige Tabelle einlesen
df <- read_excel(file,sheet = "Tabelle2",skip=3)

#### Gen-Spalten definieren
genes <- c("crtM", "crtN", "crtP", "crtQ", "crtO")

#### Sample-Namen generieren (Strain + Clone + Condition)
df <- df %>%
    mutate(Sample = paste(Strain, Clones, Condition, sep = "-"))

#### Missing → NA umwandeln & numerisch machen
df <- df %>%
    mutate(across(all_of(genes), as.character)) %>%   # alles erstmal als text
    mutate(across(all_of(genes), ~ na_if(., "Missing"))) %>%  # Missing -> NA
    mutate(across(all_of(genes), as.numeric))          # zurück zu numeric

#### Matrix mit ZAHLEN direkt erstellen (Samples als Zeilen, Gene als Spalten)
matrix_ready <- df %>%
    select(Sample, all_of(genes)) %>%
    column_to_rownames("Sample")

#### Ausgabe ansehen mit ZAHLEN
matrix_ready

### Matrix mit MEHREREN PFEILEN
arrowify <- function(x, factor = 2) {
    sapply(x, function(v) {
        if (is.na(v)) return("")
        n <- round(abs(v) * factor)
        if (n == 0) return("-")
        if (v > 0) return(paste(rep("↑", n), collapse = ""))
        if (v < 0) return(paste(rep("↓", n), collapse = ""))
    })
}

arrow_matrix <- matrix_ready %>%
    mutate(across(everything(), ~ arrowify(.)))

arrow_matrix

###################################################################
### Erweiterung: Pfeil-Matrix mit sortierter Sample-Reihenfolge ###
### (oben TSB, unten PASN)                                      ###
###################################################################

# 1) Condition-Information extrahieren
sample_info <- df %>%
    select(Sample, Condition)

# 2) Reihenfolge definieren: TSB zuerst, dann PASN / PSAN
sample_order <- sample_info %>%
    mutate(order_group = case_when(
        Condition == "TSB" ~ 1,
        Condition %in% c("PASN", "PSAN") ~ 2,
        TRUE ~ 3
    )) %>%
    arrange(order_group, Sample) %>%
    pull(Sample)

# 3) arrow_matrix in dieser Reihenfolge neu sortieren
arrow_matrix_sorted <- arrow_matrix %>%
    slice(match(sample_order, rownames(arrow_matrix)))

# 4) Ausgabe zeigen
arrow_matrix_sorted

##################################################################################
#################################################################################
##################################################################################
###Saubere Zahlenmatrix delta-Werte - Production & Regulation
### Saubere Zahlenmatrix: Production & Regulatory, inkl. Δ(PASN–TSB)

library(readxl)
library(dplyr)
library(tidyr)

#### Daten einlesen
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")
file <- "resources/LaraYaiza_tabel_STX.xlsx"
df <- read_excel(file, sheet = "Tabelle2", skip = 3)

#### "Missing" → NA und korrekte numerische Konvertierung
df <- df %>%
    mutate(across(
        c(crtQ, crtO, rsbU, rsbV, `sigma-B`, rsbW),
        ~na_if(as.character(.x), "Missing")
    )) %>%
    mutate(across(
        c(crtM, crtN, crtP, crtQ, crtO,
          `sigma-B`, rsbU, rsbV, rsbW),
        ~as.numeric(.x)
    ))


#### ----------------------------------------------
#### 1. Produktionsgene MATRIX
#### ----------------------------------------------

matrix_prod <- df %>%
    unite("Clone_Cond", Strain, Clones, Condition, sep = "-") %>%
    column_to_rownames("Clone_Cond") %>%
    select(crtM, crtN, crtP, crtQ, crtO)

print(matrix_prod)



#### ----------------------------------------------
#### 2. Regulatory MATRIX
#### ----------------------------------------------

matrix_reg <- df %>%
    unite("Clone_Cond", Strain, Clones, Condition, sep = "-") %>%
    column_to_rownames("Clone_Cond") %>%
    select(`sigma-B`, rsbU, rsbV, rsbW)

print(matrix_reg)

#### ----------------------------------------------
#### 3. Combine Production + Regulatory MATRIX
#### ----------------------------------------------

matrix_full <- cbind(matrix_prod, matrix_reg)

print(matrix_full)
#==> THIS WORKS
#WHY DELTA MATRIX NOT WORKING???

#### ----------------------------------------------
#### 4. Δ-Matrix: PASN – TSB für alle Gene
#### ----------------------------------------------

# 1. Rownames in separate Spalten zerlegen
df_full <- matrix_full %>%
    tibble::rownames_to_column("Clone_Cond") %>%
    separate(Clone_Cond, into = c("Strain", "Clone", "Condition"), sep = "-", remove = FALSE)

# 2. PASN/TSB breit machen
df_wide <- df_full %>%
    pivot_wider(
        id_cols = c(Strain, Clone),
        names_from = Condition,
        values_from = c(crtM, crtN, crtP, crtQ, crtO, `sigma-B`, rsbU, rsbV, rsbW)
    )

# 3. Δ = PASN – TSB für jedes Gen berechnen
matrix_delta_full <- df_wide %>%
    mutate(
        dcrtM   = crtM_PASN   - crtM_TSB,
        dcrtN   = crtN_PASN   - crtN_TSB,
        dcrtP   = crtP_PASN   - crtP_TSB,
        dcrtQ   = crtQ_PASN   - crtQ_TSB,
        dcrtO   = crtO_PASN   - crtO_TSB,
        dsigmaB = `sigma-B_PASN` - `sigma-B_TSB`,
        drsbU   = rsbU_PASN   - rsbU_TSB,
        drsbV   = rsbV_PASN   - rsbV_TSB,
        drsbW   = rsbW_PASN   - rsbW_TSB
    ) %>%
    select(Strain, Clone, starts_with("d")) %>%
    unite("Clone_ID", Strain, Clone, sep = "-") %>%
    column_to_rownames("Clone_ID")

print(matrix_delta_full)






