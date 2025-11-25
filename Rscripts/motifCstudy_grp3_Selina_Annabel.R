
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)


methylation_data_EGGnogAnnotated <- readRDS("~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds")

# Filter rows where base_after starts with "C" and strand is "-"

filtered_data <- methylation_data_EGGnogAnnotated[
    !is.na(methylation_data_EGGnogAnnotated$base_after) &
        startsWith(methylation_data_EGGnogAnnotated$base_after, "C") &
        methylation_data_EGGnogAnnotated$strand == "-",
]


# View the new dataset
View(filtered_data)

# Create a dataset containing only the motif column from the filtered rows
motif_dataset <- filtered_data[, "motif", drop = FALSE]

# View the new motif-only dataset
View(motif_dataset)


write.csv(motif_dataset, file = "motif_dataset.csv", row.names = FALSE)


### Compute per-position letter probabilities for 41-mer motifs
### and plot them

## Install/load required packages
packages <- c("dplyr", "tidyr", "ggplot2")
to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(tidyr)
library(ggplot2)

## 1. Extract sequences from motif_dataset (change 'motif' to your column name if needed)
seqs <- as.character(motif_dataset$motif)

## Remove missing entries
seqs <- seqs[!is.na(seqs)]

## Keep only sequences of length 41
seqs <- seqs[nchar(seqs) == 41]

## Basic checks
stopifnot(length(seqs) > 0)
seq_length <- 41
n_seq <- length(seqs)

## 2. Determine alphabet (all distinct letters present)
alphabet <- sort(unique(unlist(strsplit(paste(seqs, collapse = ""), ""))))

## 3. Count letters at each of the 41 positions
counts <- matrix(
    0,
    nrow = seq_length,
    ncol = length(alphabet),
    dimnames = list(
        Position = 1:seq_length,
        Letter   = alphabet
    )
)

for (s in seqs) {
    chars <- strsplit(s, "")[[1]]
    for (pos in seq_len(seq_length)) {
        letter <- chars[pos]
        counts[pos, letter] <- counts[pos, letter] + 1
    }
}

## 4. Convert counts to probabilities
probs <- counts / n_seq  # Position x Letter

## 5. Convert to long data frame for plotting
prob_df <- as.data.frame(probs)
prob_df$Position <- as.numeric(rownames(prob_df))

prob_long <- prob_df %>%
    pivot_longer(
        cols      = all_of(alphabet),
        names_to  = "Letter",
        values_to = "Probability"
    )

## 6. Plot probabilities for each letter at each position
ggplot(prob_long, aes(x = Position, y = Probability, color = Letter)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = 1:seq_length) +
    labs(
        title = "Per-position letter probabilities across 41-letter motifs",
        x = "Position in motif",
        y = "Probability"
    ) +
    theme_minimal()




#Make the Plot with big Letters

install.packages("ggseqlogo", dependencies = TRUE)

library(tidyverse)
library(ggseqlogo)

# Load your data
methylation_data_EGGnogAnnotated <- readRDS(
    "~/Bio253_Research_in_Omics_HS2025/resources/methylation_data_EGGnogAnnotated.rds"
)

# Filter rows of interest
filtered_data <- methylation_data_EGGnogAnnotated %>%
    filter(
        !is.na(base_after),
        startsWith(base_after, "C"),
        strand == "-"
    )

# Extract motifs
seqs <- filtered_data$motif %>%
    as.character() %>%
    discard(is.na) %>%        # remove NA
    keep(~ nchar(.) == 41)    # keep only 41-mers

# Safety check
stopifnot(length(seqs) > 0)

# Plot a sequence logo with large letters
ggseqlogo(seqs, method = "prob") +
    theme_classic(base_size = 18) +
    labs(
        title = "Sequence Logo of Motifs",
        x = "Position",
        y = "Probability")






