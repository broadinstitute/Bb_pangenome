########## Figure S2 ##########
# This notebook calculates confidence intervals for the probability of plasmid presence by genotype.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(binom)
library(purrr)
library(RColorBrewer)
library(scales)
library(patchwork)

# input files
plasmid_file = "/Users/rl275/Projects/bb_longread/Manuscript/0_Data/0_Raw/best_hits_1000bp_v11.csv"
metadata_file = "/Users/rl275/Projects/bb_longread/Manuscript/0_Data/2_Processed/MetadataFull.csv"

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

# count the total number of isolates for each RST
rst_counts <- df %>%
  count(RST) %>% 
  mutate(RST=as.character(RST)) %>%
  bind_rows(summarise(., RST = "Total", n = sum(n)))
rst_counts

# count the number of occurrences of each plasmid across each RST
df_rst_count <- df %>%
  group_by(RST) %>%
  summarise(across(starts_with("lp") | starts_with("cp"), sum)) %>%
  as.data.frame()
df_rst_count <- melt(df_rst_count, id.vars = "RST", variable.name = "Plasmid", value.name = "Count")
df_rst_count <- dcast(df_rst_count, Plasmid ~ RST, value.var = "Count")
df_rst_count$Total <- rowSums(df_rst_count[, -1]) 
df_rst_count <- df_rst_count[order(df_rst_count$Total, decreasing = TRUE), ]
df_rst_count$Plasmid <- factor(df_rst_count$Plasmid, levels = df_rst_count$Plasmid)
row.names(df_rst_count) <-  NULL
head(df_rst_count)

# convert counts to proportions
df_rst_prop <- df_rst_count
for (rst in rst_counts$RST) {
  df_rst_prop[[as.character(rst)]] <- df_rst_count[[as.character(rst)]] / rst_counts$n[rst == rst_counts$RST]
}
head(df_rst_prop)

# compute Clopper-Pearson CIs 
df_rst_cp <- df_rst_count %>%
  pivot_longer(cols = -Plasmid, names_to = "RST", values_to = "Count") %>%
  left_join(rst_counts, by = "RST") %>% 
  mutate(
    Proportion = Count / n,
    CI = map2(Count, n, ~ binom.test(.x, .y, conf.level = 0.95)$conf.int),
    LB = map_dbl(CI, ~ .x[1]),
    UB = map_dbl(CI, ~ .x[2])
  ) %>%
  select(Plasmid, RST, Proportion, LB, UB) %>%
  as.data.frame()

# add plasmid type
df_rst_cp <- df_rst_cp %>%
  mutate(
    Type = if_else(grepl("^cp", Plasmid), "Circular",
                   if_else(grepl("^lp", Plasmid), "Linear", NA_character_))
  )

# plotting function for RST
errorbar_plot_rst <- function(df, title) {
  ggplot(df, aes(x = Proportion, y = Plasmid, colour = RST)) +
    geom_errorbar(
      aes(xmin = LB, xmax = UB),
      position = position_dodge(0.75),
      width = 0.5
    ) +
    geom_point(position = position_dodge(0.75)) +
    ggtitle(title) +
    scale_color_manual(name="RST", values = c("#FF0000", "#00FF00", "#0000FF", "grey35"), rst_counts$RST) +
    ylab("Plasmid") +
    xlab("Probability of Presence") +
    theme_bw()
}

# generate plots for RST
p1 <- errorbar_plot_rst(df_rst_cp[df_rst_cp$Type == "Circular", ], 
                  "Circular Plasmids, RST")
p2 <- errorbar_plot_rst(df_rst_cp[df_rst_cp$Type == "Linear", ], 
                  "Linear Plasmids, RST")


# count the total number of isolates for each OspC type
ospc_counts <- df %>%
  count(OspC) %>%
  bind_rows(summarise(., OspC = "Total", n = sum(n)))
ospc_counts

# count the number of occurrences of each plasmid across each OspC type
df_ospc_count <- df %>%
  group_by(OspC) %>%
  summarise(across(starts_with("lp") | starts_with("cp"), sum)) %>%
  as.data.frame()
df_ospc_count <- melt(df_ospc_count, id.vars = "OspC", variable.name = "Plasmid", value.name = "Count")
df_ospc_count <- dcast(df_ospc_count, Plasmid ~ OspC, value.var = "Count")
df_ospc_count$Total <- rowSums(df_ospc_count[, -1]) 
df_ospc_count <- df_ospc_count[order(df_ospc_count$Total, decreasing = TRUE), ]
df_ospc_count$Plasmid <- factor(df_ospc_count$Plasmid, levels = df_ospc_count$Plasmid)
row.names(df_ospc_count) <-  NULL
head(df_ospc_count)

# convert counts to proportions
df_ospc_prop <- df_ospc_count
for (ospc in ospc_counts$OspC) {
  df_ospc_prop[[as.character(ospc)]] <- df_ospc_count[[as.character(ospc)]] / ospc_counts$n[ospc == ospc_counts$OspC]
}
head(df_ospc_prop)

# compute Clopper-Pearson CIs 
df_ospc_cp <- df_ospc_count %>%
  pivot_longer(cols = -Plasmid, names_to = "OspC", values_to = "Count") %>%
  left_join(ospc_counts, by = "OspC") %>%
  mutate(
    Proportion = Count / n,
    CI = map2(Count, n, ~ binom.test(.x, .y, conf.level = 0.95)$conf.int),
    LB = map_dbl(CI, ~ .x[1]),
    UB = map_dbl(CI, ~ .x[2])
  ) %>%
  select(Plasmid, OspC, Proportion, LB, UB) %>%
  as.data.frame()

# add plasmid type
df_ospc_cp <- df_ospc_cp %>%
  mutate(
    Type = if_else(grepl("^cp", Plasmid), "Circular",
                   if_else(grepl("^lp", Plasmid), "Linear", NA_character_))
  )

# assign colors by OspC type
ospc_types <- unique(df_ospc_cp$OspC)
ospc_types_letters <- ospc_types[ospc_types != "Total"]
ospc_colors = hue_pal()(length(ospc_types_letters))
df_ospc_colors <- data.frame(
  OspC = ospc_types_letters,
  Color = ospc_colors
)

df_ospc_colors <- rbind(df_ospc_colors, data.frame(OspC = "Total", Color = "grey35"))
df_ospc_colors

# plotting function for OspC
errorbar_plot_ospc <- function(df, title, color_df) {
  ggplot(df, aes(x = Proportion, y = Plasmid, colour = OspC)) +
    geom_errorbar(
      aes(xmin = LB, xmax = UB),
      position = position_dodge(0.75),
      width = 0.5
    ) +
    geom_point(position = position_dodge(0.75)) +
    ggtitle(title) +
    scale_color_manual(values = setNames(
      df_ospc_colors$Color, df_ospc_colors$OspC 
    )) +
    ylab("Plasmid") +
    xlab("Probability of Presence") +
    theme_bw()
}

df_selected <- df_ospc_cp %>%
  filter(OspC %in% c("A", "K", "Total"))

# generate plots for OspC
p3 <- errorbar_plot_ospc(df_selected[df_selected$Type == "Circular", ], 
                 "Circular Plasmids, OspC Type")
p4 <- errorbar_plot_ospc(df_selected[df_selected$Type == "Linear", ], 
                 "Linear Plasmids, OspC Type")

# output panelled plot
(p1 + p2 + p3 + p4) +
  plot_layout(ncol = 2)
