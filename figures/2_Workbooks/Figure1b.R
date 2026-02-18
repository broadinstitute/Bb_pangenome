########## Figure 1b ##########
# This notebook calculates the probabilities of plasmid presence with Clopper-Pearson CIs
# based on the presence / absence data within our sampled isolates.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(binom)
library(purrr)
library(RColorBrewer)
library(scales)

# input files
plasmid_file = "../../output/genotyping/replicons/best_hits_1000bp_v11.csv"
metadata_file = "../0_Data/2_Processed/MetadataFull.csv"
plasmidbinary_file = "../0_Data/2_Processed/plasmids_binary_long.csv" # generated in `Figure1a.ipynb`


# import and clean plasmid data
plasmids <- read.csv(plasmid_file)
colnames(plasmids)[colnames(plasmids) == "assembly_id"] <- "Isolate"
plasmids$Isolate <- gsub("GCF_\\d+\\.\\d+_", "", plasmids$Isolate)
plasmids$Isolate <- gsub("_genomic", "", plasmids$Isolate)
plasmids <- unique(plasmids[, c("Isolate", "plasmid_name")])
plasmids$value <- 1 
plasmids_wide <- reshape(plasmids, 
                         idvar = "Isolate", 
                         timevar = "plasmid_name", 
                         direction = "wide")
row.names(plasmids_wide) <- NULL
plasmids_wide[is.na(plasmids_wide)] <- 0
colnames(plasmids_wide) <- gsub("value\\.", "", colnames(plasmids_wide))
plasmids_wide[, -1] <- lapply(plasmids_wide[, -1], as.integer)
head(plasmids_wide)

# import and clean metadata
metadata <- read.csv(metadata_file,
                     header=TRUE,
                     sep=",")
metadata$Strain <- metadata$Strain_ID
metadata <- metadata %>%
  pivot_longer(
    cols = c(Strain_ID, Alias, Assembly_ID),
    names_to = "Name_Type",
    values_to = "Isolate"
  ) %>%
  as.data.frame()
head(metadata)

# combine plasmids with metadata
df <- plasmids_wide %>%
  left_join(metadata[,c("Isolate", "RST", "OspC")], by = "Isolate")
head(df)

# pull the total number of isolates
n_isolates <- dim(df)[1]

# count the number of occurrences of each plasmid across the sampled isolates
df_count <- df %>%
  summarise(across(starts_with("lp") | starts_with("cp"), sum)) %>%
  as.data.frame()
df_count <- df_count %>%
  pivot_longer(cols = everything(), 
               names_to = "Plasmid", 
               values_to = "Count")
df_count <- df_count[order(df_count$Count, decreasing = TRUE), ]
df_count$n <- n_isolates

# compute Clopper-Pearson CIs 
df_count_cp <- df_count %>%
  mutate(
    Proportion = Count / n, 
    CI = map2(Count, n, ~ binom.test(.x, .y, conf.level = 0.95)$conf.int),
    LB = map_dbl(CI, ~ .x[1]),
    UB = map_dbl(CI, ~ .x[2])
  ) %>%
  select(Plasmid, Proportion, LB, UB) %>%
  as.data.frame()
df_count_cp$Plasmid <- factor(df_count_cp$Plasmid, levels = df_count_cp$Plasmid[order(df_count_cp$Proportion)])

# pull plasmid color assignments from a previously generated file
plasmids <- read.csv(plasmidbinary_file,
                     header=TRUE,
                     sep=",")
plasmid_colors <- plasmids %>%
  select(c("plasmid", "color")) %>%
  distinct()
color_vector <- setNames(plasmid_colors$color, plasmid_colors$plasmid)

# define function for plotting 
errorbar_plot <- function(df, title) {
  ggplot(df, aes(x = Proportion, y = Plasmid, color = Plasmid)) +
    geom_errorbar(
      aes(xmin = LB, xmax = UB),
      position = position_dodge(0.75),
      width = 0.5
    ) +
    scale_color_manual(values = color_vector) +
    geom_point(position = position_dodge(0.75)) +
    ggtitle(title) +
    ylab("Plasmid") +
    xlab("Probability of Presence") +
    theme_bw() + 
    theme(legend.position="none") 
}

# generate Figure 1b
errorbar_plot(df_count_cp, "Probabilities of Plasmid Presence (Clopper-Pearson CI)")