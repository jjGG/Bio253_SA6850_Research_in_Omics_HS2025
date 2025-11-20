### OMICS-2025 — delta_Heatmaps (PASN–TSB)
### HEATMAP FOR THE POSTER!!!!

### Separate HEATMAPS FOR EACH STRAIN JE2 and 6850
### delta(PASN-TA)in each genes/protein expression for PRODUCTION & REGULATORY GENES
### x-Achse= Klones, y-Achse delta-Gennames

library(tidyverse)
library(readxl)
library(pheatmap)
library(grid)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(dplyr)

#### 1) Read data
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")       # relative folder path!
file <- "resources/LaraYaiza_tabel_STX.xlsx"             # Dateiname im Resource folder
#### richtige Tabelle einlesen
df <- read_excel(file,sheet = "Tabelle2",skip=3)

genes <- c("crtM","crtN","crtP","crtQ","crtO",
           "sigma-B","rsbU","rsbV","rsbW")


### 2) Missing → NA & numerisch
df <- df %>%
    mutate(across(all_of(genes), as.character)) %>%
    mutate(across(all_of(genes), ~ na_if(., "Missing"))) %>%
    mutate(across(all_of(genes), as.numeric))

# Beispiel-Sample-Namen (Strain-Clone)
df <- df %>%
    mutate(CloneID = paste(Strain, Clones, sep = "-"))

### 3) PASN & TSB trennen
df_wide <- df %>%
    select(CloneID, Condition, all_of(genes)) %>%
    pivot_wider(
        names_from = Condition,
        values_from = all_of(genes),
        names_sep = "_"
    )

### 4) Delta (PASN – TSB) berechnen
df_delta <- df_wide %>%
    mutate(across(
        .cols = ends_with("PASN"),
        .fns = ~ . - get(sub("_PASN","_TSB", cur_column())),
        .names = "{sub('_PASN','',col)}"
    )) %>%
    select(CloneID, all_of(genes))

# Jetzt ist df_delta: CloneID × Gene mit Δ log2FC

### 5) Matrix für Heatmap
mat <- df_delta %>%
    column_to_rownames("CloneID") %>%
    as.matrix()

### 6) Heatmap OHNE Dendrogramm
pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    main = "Δ Expression (PASN - TSB)",
    fontsize = 12,
    border_color = NA
)
#
# ####WAIT A SECOND
# ####I need for each Gene a BOX with the clones wirten in the -Axes (above or below the box)
# #################################################################################
# ###################################################################################
# #############################################################################
### ### OMICS-2025 — Gene-wise Δ-Expression Boxes (PASN – TSB)
# ### Generates one rectangular heatmap box per gene
# ### With: border, rotated labels, clone labels on top
# ### -------------------------------------------------------------
#
# library(tidyverse)
# library(readxl)
# library(ggplot2)
# library(grid)
#
# ### -------------------------
# ### 1) Load data
# ### -------------------------
# setwd("../Bio253_SA6850_Research_in_Omics_HS2025")       # relative folder path!
# file <- "resources/LaraYaiza_tabel_STX.xlsx"             # Dateiname im Resource folder
# #### richtige Tabelle einlesen
# df <- read_excel(file,sheet = "Tabelle2",skip=3)
#
# genes <- c("crtM","crtN","crtP","crtQ","crtO",
#            "sigma-B","rsbU","rsbV","rsbW")
#
# ### -------------------------
# ### 2) Clean up Missing → NA
# ### -------------------------
#
# df <- df %>%
#     mutate(across(all_of(genes), as.character)) %>%
#     mutate(across(all_of(genes), ~ na_if(., "Missing"))) %>%
#     mutate(across(all_of(genes), as.numeric))
#
# df <- df %>%
#     mutate(CloneID = paste(Strain, Clones, sep = "-"))
#
# ### -------------------------
# ### 3) Create PASN–TSB delta matrix
# ### -------------------------
#
# df_wide <- df %>%
#     select(CloneID, Condition, all_of(genes)) %>%
#     pivot_wider(
#         names_from = Condition,
#         values_from = all_of(genes),
#         names_sep = "_"
#     )
#
# df_delta <- df_wide %>%
#     mutate(across(
#         .cols = ends_with("PASN"),
#         .fns = ~ . - get(sub("_PASN", "_TSB", cur_column())),
#         .names = "{sub('_PASN','',col)}"
#     )) %>%
#     select(CloneID, all_of(genes))
#
# ### -------------------------
# ### 4) Create output folder
# ### -------------------------
#
# dir.create("gene_boxes", showWarnings = FALSE)
#
# ### -------------------------
# ### 5) Generate one box per gene
# ### -------------------------
#
# for (g in genes) {
#
#     df_gene <- df_delta %>%
#         select(CloneID, !!sym(g)) %>%
#         rename(value = !!sym(g))
#
#     outfile <- paste0("gene_boxes/", g, "_box.png")
#
#     ### ---------------------------------------------------------
#     ### CASE 1: all NA → placeholder grey box
#     ### ---------------------------------------------------------
#
#     if (all(is.na(df_gene$value))) {
#
#         png(outfile, width = 700, height = 400)
#         grid.newpage()
#         grid.rect(gp = gpar(fill = "grey90", col = "black", lwd = 2))
#         grid.text(
#             paste0(g, "\n(no data available)"),
#             gp = gpar(fontsize = 18, fontface = "bold")
#         )
#         dev.off()
#
#         message("Generated placeholder box for ", g)
#         next
#     }
#
#     ### ---------------------------------------------------------
#     ### CASE 2: normal heatmap-style tile using ggplot2
#     ### ---------------------------------------------------------
#
#     df_gene$CloneID <- factor(df_gene$CloneID, levels = df_gene$CloneID)
#
#     p <- ggplot(df_gene, aes(x = CloneID, y = g, fill = value)) +
#         geom_tile(color = "black", linewidth = 1.1) +    # border around cells
#         scale_fill_gradient2(
#             low = "#2166ac",
#             mid = "white",
#             high = "#b2182b",
#             midpoint = 0,
#             name = "Δ"
#         ) +
#         scale_x_discrete(position = "top") +             # clone labels above
#         theme_minimal(base_size = 16) +
#         theme(
#             axis.title = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks = element_blank(),
#             axis.text.x = element_text(
#                 angle = 45,
#                 hjust = 0,
#                 vjust = 0,
#                 face = "bold"
#             ),
#             legend.position = "right",
#             panel.grid = element_blank(),
#             plot.margin = margin(10, 10, 10, 10),
#             panel.border = element_rect(
#                 color = "black",
#                 fill = NA,
#                 size = 1.2                        # outer box border
#             )
#         ) +
#         ggtitle(g)
#
#     ggsave(outfile, p, width = 6, height = 3, dpi = 300)
#
#     message("Generated box for ", g)
# }
#
# ### -------------------------
# ### DONE
# ### Boxes saved in: gene_boxes/
# ### browseURL("gene_boxes)
# ### -------------------------

##############################################################################
#########################################################################
###########################################################################
####TODO:
### - make them vector base so when zoomed it not depixeled!
### - scale each row --> increasing contrast
### - force some break/legend
### - we want a GLOBAL SCALE Across all genes (not a per-gene scale)
#######################################################################
#
# ### -------------------------
# ### 6) Generate one box per gene
# ### -------------------------
#
# for (g in genes) {
#
#     df_gene <- df_delta %>%
#         select(CloneID, !!sym(g)) %>%
#         rename(value = !!sym(g))
#
#     outfile <- paste0("gene_boxes/", g, "_box.png")
#
#     # --------------------------------------------------------
#     # CASE 1: whole-box NA (but NOT for crtO)
#     # --------------------------------------------------------
#
#     if (g != "crtO" && all(is.na(df_gene$value))) {
#
#         png(outfile, width = 900, height = 350)
#         grid.newpage()
#
#         grid.rect(gp = gpar(fill = "white", col = "black", lwd = 3))
#         grid.text(g, y = 0.9, gp = gpar(fontsize = 22, fontface="bold"))
#         grid.text("NA", gp = gpar(fontsize = 28, fontface = "bold"))
#
#         dev.off()
#         next
#     }
#
#     # --------------------------------------------------------
#     # CASE 2: normal heatmap, including crtO full NA case
#     # --------------------------------------------------------
#
#     df_gene$CloneID <- factor(df_gene$CloneID, levels = df_gene$CloneID)
#
#     # For crtO: replace all values with NA (ensures tile layout)
#     if (g == "crtO") {
#         df_gene$value <- NA
#     }
#
#     p <- ggplot(df_gene, aes(x = CloneID, y = g, fill = value)) +
#         geom_tile(color = "black", linewidth = 1.0) +
#         scale_fill_gradient2(
#             low = "#2166ac",
#             mid = "white",
#             high = "#b2182b",
#             midpoint = 0,
#             limits = c(-max_abs, max_abs),
#             breaks = legend_breaks,
#             name = "Δ(PASN–TSB)",
#             na.value = "white"
#         ) +
#         geom_text(
#             aes(label = ifelse(is.na(value), "NA", "")),
#             color = "black",
#             fontface = "bold",
#             size = 5
#         ) +
#         scale_x_discrete(position = "top") +
#         theme_minimal(base_size = 18) +
#         theme(
#             axis.title = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks = element_blank(),
#             axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
#             panel.grid = element_blank(),
#             plot.margin = margin(10, 10, 10, 10),
#             panel.border = element_rect(color = "black", fill = NA, size = 1.5)
#         ) +
#         ggtitle(g)
#
#     ggsave(outfile, p, width = 7, height = 3.5, dpi = 300)
# }

###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### TODO - NOT DONE YET!
### add number in middel of each box to one decimal~ gerundet
### LEGENDE make the colour bar horizontal since legend titel is going to be large in order to save space
### Legende " log2(PASN/TSB)" --> whow maximum and minimum and the zero
### labelign does not work in R --> do it with powepoint or in biorender-thingy
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


###############################################################################
### OMICS-2025 — Gene-wise Δ(PASN–TSB) Heatmap Boxes (FINAL FIXED VERSION)
### - Global color scale across all genes
### - Row contrast enhancement
### - NA handling (+ forced NA for crtO)
### - Vector output (PDF)
###############################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(grid)

### ---------------------------------------------------------
### 1. Load data
### ---------------------------------------------------------
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")
file <- "resources/LaraYaiza_tabel_STX.xlsx"
df <- read_excel(file, sheet = "Tabelle2", skip = 3)

genes <- c("crtM","crtN","crtP","crtQ","crtO",
           "sigma-B","rsbU","rsbV","rsbW")

### ---------------------------------------------------------
### 2. Clean + prepare
### ---------------------------------------------------------
df <- df %>%
    mutate(across(all_of(genes), ~na_if(as.character(.), "Missing"))) %>%
    mutate(across(all_of(genes), as.numeric)) %>%
    mutate(CloneID = paste(Strain, Clones, sep = "-"))

### ---------------------------------------------------------
### 3. Compute Δ(PASN – TSB)
### ---------------------------------------------------------
df_wide <- df %>%
    select(CloneID, Condition, all_of(genes)) %>%
    pivot_wider(
        names_from = Condition,
        values_from = all_of(genes),
        names_sep = "_"
    )

df_delta <- df_wide %>%
    mutate(across(
        .cols = ends_with("PASN"),
        .fns  = ~ . - get(sub("_PASN","_TSB",cur_column())),
        .names = "{sub('_PASN','',col)}"
    )) %>%
    select(CloneID, all_of(genes))

### ---------------------------------------------------------
### 4. GLOBAL COLOR SCALE (must come BEFORE plotting!!!)
### ---------------------------------------------------------
all_vals <- unlist(df_delta[, genes])
global_min <- min(all_vals, na.rm = TRUE)
global_max <- max(all_vals, na.rm = TRUE)
global_limit <- max(abs(global_min), abs(global_max))

legend_breaks <- seq(-global_limit, global_limit, by = 1)
message("Global scale: ", -global_limit, " to ", global_limit)

### ---------------------------------------------------------
### 5. Output folder
### ---------------------------------------------------------
dir.create("gene_boxes", showWarnings = FALSE)

### ---------------------------------------------------------
### 6. Generate one box per gene (PDF vector graphics)
### ---------------------------------------------------------
for (g in genes) {

    df_gene <- df_delta %>%
        select(CloneID, value = !!sym(g))

    df_gene$CloneID <- factor(df_gene$CloneID, levels = df_gene$CloneID)

    ### Force crtO to be all NA (requested)
    if (g == "crtO") df_gene$value <- NA

    ### Row contrast enhancement
    row_max <- max(abs(df_gene$value), na.rm = TRUE)
    df_gene$value_scaled <- df_gene$value

    if (!is.infinite(row_max) && row_max > 0) {
        df_gene$value_scaled <- df_gene$value / row_max * global_limit
    }

    ### Plot
    p <- ggplot(df_gene, aes(x = CloneID, y = 1, fill = value_scaled)) +
        geom_tile(color = "black", linewidth = 1.3) +
        geom_text(
            aes(label = ifelse(is.na(value), "NA", "")),
            color = "black", fontface = "bold", size = 5
        ) +
        scale_fill_gradient2(
            low = "#2166ac",
            mid = "white",
            high = "#b2182b",
            midpoint = 0,
            limits = c(-global_limit, global_limit),
            breaks = legend_breaks,
            name = "Δ(PASN–TSB)",
            na.value = "white"
        ) +
        scale_x_discrete(position = "top") +
        scale_y_continuous(expand = c(0,0)) +
        ggtitle(g) +
        theme_minimal(base_size = 18) +
        theme(
            axis.title   = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks   = element_blank(),
            axis.text.x  = element_text(angle = 45, hjust = 0),
            panel.grid   = element_blank(),
            plot.margin  = margin(10, 10, 10, 10),
            panel.border = element_rect(color = "black", fill = NA, size = 1.5)
        )

    ### SAVE VECTOR OUTPUT
    ggsave(paste0("gene_boxes/", g, "_box.pdf"),
           p, width = 7, height = 3.5)

    message("✓ Saved ", g)
}

message("### ALL BOXES GENERATED SUCCESSFULLY ###")

###****************************************************
###############################################################################
### OMICS-2025 — Gene-wise Δ(PASN–TSB) Heatmap Boxes (FINAL)
### - Global color scale
### - Row contrast enhancement
### - Clean legend with fixed breaks
### - Vector output (PDF + SVG)
### - Clean spacing and correct title size
###############################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(grid)

### ---------------------------------------------------------
### 1. Load data
### ---------------------------------------------------------
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")
file <- "resources/LaraYaiza_tabel_STX.xlsx"
df <- read_excel(file, sheet = "Tabelle2", skip = 3)

genes <- c("crtM","crtN","crtP","crtQ","crtO",
           "sigma-B","rsbU","rsbV","rsbW")

### ---------------------------------------------------------
### 2. Clean data
### ---------------------------------------------------------
df <- df %>%
    mutate(across(all_of(genes), ~na_if(as.character(.), "Missing"))) %>%
    mutate(across(all_of(genes), as.numeric)) %>%
    mutate(CloneID = paste(Strain, Clones, sep = "-"))

### ---------------------------------------------------------
### 3. Compute Δ(PASN – TSB)
### ---------------------------------------------------------
df_wide <- df %>%
    select(CloneID, Condition, all_of(genes)) %>%
    pivot_wider(
        names_from = Condition,
        values_from = all_of(genes),
        names_sep = "_"
    )

df_delta <- df_wide %>%
    mutate(across(
        .cols = ends_with("PASN"),
        .fns  = ~ . - get(sub("_PASN","_TSB",cur_column())),
        .names = "{sub('_PASN','',col)}"
    )) %>%
    select(CloneID, all_of(genes))

### ---------------------------------------------------------
### 4. Global symmetric color scale
### ---------------------------------------------------------
all_vals <- unlist(df_delta[ , genes])
global_limit <- max(abs(all_vals), na.rm = TRUE)
legend_breaks <- round(seq(-global_limit, global_limit, by = 1), 2)

### ---------------------------------------------------------
### 5. Output folder
### ---------------------------------------------------------
dir.create("gene_boxes", showWarnings = FALSE)

### ---------------------------------------------------------
### 6. Per-gene plot generation
### ---------------------------------------------------------
for (g in genes) {

    df_gene <- df_delta %>%
        select(CloneID, value = !!sym(g))

    df_gene$CloneID <- factor(df_gene$CloneID, levels = df_gene$CloneID)

    if (g == "crtO") df_gene$value <- NA

    # Row contrast scaling
    row_max <- max(abs(df_gene$value), na.rm = TRUE)
    df_gene$value_scaled <- ifelse(
        is.na(df_gene$value), NA,
        df_gene$value / row_max * global_limit
    )

    p <- ggplot(df_gene, aes(x = CloneID, y = 1, fill = value_scaled)) +
        geom_tile(color = "black", linewidth = 1.2) +

        geom_text(
            aes(label = ifelse(is.na(value), "NA", "")),
            color = "black", fontface = "bold", size = 5
        ) +

        scale_fill_gradient2(
            low = "#2166ac",
            mid = "white",
            high = "#b2182b",
            midpoint = 0,
            limits = c(-global_limit, global_limit),
            breaks = legend_breaks,
            name = "delta(PASN–TSB)",
            na.value = "white"
        ) +

        scale_x_discrete(position = "top") +
        scale_y_continuous(expand = c(0,0)) +

        ggtitle(g) +

        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
            panel.grid = element_blank(),
            plot.margin = margin(5, 10, 5, 10),
            panel.border = element_rect(color = "black", fill = NA, size = 1.5)
        )

    # Save SVG + PDF + PNG
    ggsave(paste0("gene_boxes/", g, "_box.pdf"), p, width = 7, height = 3.5)
    #ggsave(paste0("gene_boxes/", g, "_box.svg"), p, width = 7, height = 3.5)
    #ggsave(paste0("gene_boxes/", g, "_box.png"), p, width = 7, height = 3.5, dpi = 300)
}

###############################################################################
### DONE — Vector files saved in: gene_boxes/
###############################################################################


































###############################################################################
### OMICS-2025 — FINAL GENE BOX SCRIPT
### With:
### - One-box-per-gene (vector format)
### - Global color scale
### - Row contrast (enhanced visibility)
### - Numbers inside each tile (1 decimal)
### - NA labels
### - ONE external horizontal legend ("log2(PASN/TSB)")
###############################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(grid)

### ---------------------------------------------------------
### 1. Load data
### ---------------------------------------------------------
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")
file <- "resources/LaraYaiza_tabel_STX.xlsx"

df <- read_excel(file, sheet = "Tabelle2", skip = 3)

genes <- c("crtM","crtN","crtP","crtQ","crtO",
           "sigma-B","rsbU","rsbV","rsbW")


### ---------------------------------------------------------
### 2. Clean data
### ---------------------------------------------------------
df <- df %>%
    mutate(across(all_of(genes), ~na_if(as.character(.), "Missing"))) %>%
    mutate(across(all_of(genes), as.numeric)) %>%
    mutate(CloneID = paste(Strain, Clones, sep = "-"))


### ---------------------------------------------------------
### 3. Compute Δ(PASN – TSB)
### ---------------------------------------------------------
df_wide <- df %>%
    select(CloneID, Condition, all_of(genes)) %>%
    pivot_wider(
        names_from = Condition,
        values_from = all_of(genes),
        names_sep = "_"
    )

df_delta <- df_wide %>%
    mutate(across(
        .cols = ends_with("PASN"),
        .fns  = ~ . - get(sub("_PASN","_TSB",cur_column())),
        .names = "{sub('_PASN','',col)}"
    )) %>%
    select(CloneID, all_of(genes))


### ---------------------------------------------------------
### 4. Global symmetric color scale for ALL genes
### ---------------------------------------------------------
all_vals <- unlist(df_delta[ , genes])
global_limit <- max(abs(all_vals), na.rm = TRUE)

legend_breaks <- c(-global_limit, 0, global_limit)


### ---------------------------------------------------------
### 5. Output folder
### ---------------------------------------------------------
dir.create("gene_boxes", showWarnings = FALSE)


### ---------------------------------------------------------
### 6. Generate ONE HORIZONTAL LEGEND (SVG + PDF)
### ---------------------------------------------------------
legend_plot <- ggplot(data.frame(x = 1, y = 1, val = 1), aes(x, y, fill = val)) +
    geom_tile() +
    scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        limits = c(-global_limit, global_limit),
        breaks = legend_breaks,
        labels = round(legend_breaks, 2),
        name = "log2(PASN / TSB)"
    ) +
    guides(fill = guide_colourbar(
        title.position = "top",
        direction = "horizontal",
        barwidth = 20,
        barheight = 1.5,
        ticks.colour = "black"
    )) +
    theme_minimal(base_size = 18) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom"
    )

ggsave("gene_boxes/LEGEND_global.pdf", legend_plot, width = 10, height = 2)
ggsave("gene_boxes/LEGEND_global.svg", legend_plot, width = 10, height = 2)



### ---------------------------------------------------------
### 7. Generate the gene boxes (NO legend in each box)
### ---------------------------------------------------------
for (g in genes) {

    df_gene <- df_delta %>% select(CloneID, value = !!sym(g))
    df_gene$CloneID <- factor(df_gene$CloneID, levels = df_gene$CloneID)

    if (g == "crtO") df_gene$value <- NA

    # row contrast scaling
    row_max <- max(abs(df_gene$value), na.rm = TRUE)
    df_gene$value_scaled <- ifelse(
        is.na(df_gene$value), NA,
        df_gene$value / row_max * global_limit
    )

    p <- ggplot(df_gene, aes(x = CloneID, y = 1, fill = value_scaled)) +
        geom_tile(color = "black", linewidth = 1.2) +

        # numbers inside boxes
        geom_text(
            aes(label = ifelse(is.na(value), "NA", sprintf("%.1f", value))),
            fontface = "bold",
            size = 4
        ) +

        scale_fill_gradient2(
            low = "#2166ac",
            mid = "white",
            high = "#b2182b",
            limits = c(-global_limit, global_limit),
            midpoint = 0,
            na.value = "white",
            guide = "none"    # <- REMOVE legend here
        ) +

        scale_x_discrete(position = "top") +
        scale_y_continuous(expand = c(0,0)) +

        ggtitle(g) +

        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1.5),
            plot.margin = margin(5, 10, 5, 10)
        )

    # Save vector + PNG
    ggsave(paste0("gene_boxes/", g, "_box.pdf"), p, width = 7, height = 3.5)
    #ggsave(paste0("gene_boxes/", g, "_box.svg"), p, width = 7, height = 3.5)
    #ggsave(paste0("gene_boxes/", g, "_box.png"), p, width = 7, height = 3.5, dpi = 300)
}


###############################################################################
### DONE — Boxes + ONE GLOBAL LEGEND saved in gene_boxes/
###############################################################################
####max auf 3 und auf -3 min --> statt











