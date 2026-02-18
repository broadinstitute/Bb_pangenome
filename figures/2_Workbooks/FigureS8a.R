########## Figure S8a ##########
# This notebook generates boxplots comparing genome length by RST and OspC type.

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(gridExtra)

# input files
metadata_file = "../0_Data/2_Processed/MetadataFull.csv"
plasmid_file = "../0_Data/0_Raw/best_hits_1000bp_v11.csv"

# process and clean metadata
metada <- read.csv(metadata_file,
                   header=TRUE,
                   sep=",")
metada$Strain <- metada$Strain_ID
metada <- metada %>%
  pivot_longer(
    cols = c(Strain_ID, Alias, Assembly_ID),
    names_to = "Name_Type", 
    values_to = "Isolate" 
  )
metada$Isolate <- gsub("GCF_\\d+\\.\\d+_", "", metada$Isolate)
metada$Isolate <- gsub("_genomic", "", metada$Isolate)

# import and clean plasmid annotations
annotations <- read.csv(plasmid_file,
                        header=TRUE,
                        sep=",")
annotations$assembly_id <- gsub("GCF_\\d+\\.\\d+_", "", annotations$assembly_id)
annotations$assembly_id <- gsub("_genomic", "", annotations$assembly_id)
annotations_condensed <- annotations %>%
  group_by(assembly_id) %>%
  summarise(total_contig_len = sum(contig_len), .groups = "drop") %>%
  rename(Isolate=assembly_id) %>%
  left_join(metada, by="Isolate") %>%
  mutate(RST = as.factor(RST)) %>%
  as.data.frame()
annotations_condensed

# set legend / formatting parameters
plottitlesize <- 24
axistitlesize <- 16
axistextsize <- 14
legendtitlesize <- 16
legendtextsize <- 14

# genome length by RST
box_rst <- ggplot(annotations_condensed, aes(x = RST, y = total_contig_len)) +
  geom_boxplot(aes(fill=RST),outlier.shape = NA) +
  geom_dotplot(
    aes(x=RST),
    fill="black", binaxis = 'y', stackdir = 'center', dotsize = 0.35, stackratio = 1.5
  ) +
  geom_text_repel(
    data = annotations_condensed %>%
      group_by(RST) %>%
      filter(total_contig_len > quantile(total_contig_len, 0.75) + 1.5 * IQR(total_contig_len) |
               total_contig_len < quantile(total_contig_len, 0.25) - 1.5 * IQR(total_contig_len)),
    aes(label = Strain, color=RST), 
    vjust = -1, size = 3.5, show.legend=FALSE
  ) +
  labs(
    title = "Genome Length by RST",
    x = "RST",
    y = "Genome Length"
  ) +
  scale_fill_manual(
    values = c("1" = "#FF0000", "2" = "#00FF00", "3" = "#0000FF")  
  ) +
  scale_color_manual(
    values = c("1" = "#FF0000", "2" = "#00FF00", "3" = "#0000FF") 
  ) +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c(1,2),c(2,3),c(1,3))) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = plottitlesize, face = "bold"),      
    axis.title = element_text(size = axistitlesize),                  
    axis.text = element_text(size = axistextsize),               
    legend.title = element_text(size = legendtitlesize),   
    legend.text = element_text(size = legendtextsize)  
  )

# genome length by OspC type
ospc_types <-sort(unique(lipo_counts$OspC))
ospc_colors = hue_pal()(length(ospc_types))

df_ospc_colors <- data.frame(
  OspC = ospc_types,
  Color = ospc_colors
)
ospc_colors <- setNames(df_ospc_colors$Color, df_ospc_colors$OspC)

box_ospc <- ggplot(annotations_condensed, aes(x = OspC, y = total_contig_len)) +
  geom_boxplot(aes(fill=OspC),outlier.shape = NA) + 
  geom_dotplot(
    aes(OspC),
    fill="black", binaxis = 'y', stackdir = 'center', dotsize = 0.35, stackratio = 1.5
  ) +
  geom_text_repel(
    data = annotations_condensed %>%
      group_by(OspC) %>%
      filter(total_contig_len > quantile(total_contig_len, 0.75) + 1.5 * IQR(total_contig_len) |
               total_contig_len < quantile(total_contig_len, 0.25) - 1.5 * IQR(total_contig_len)),
    aes(label = Strain, color=OspC), 
    vjust = -1, size = 3.5, show.legend=FALSE
  ) +
  labs(
    title = "Genome Length by OspC Type",
    x = "OspC",
    y = "Genome Length"
  ) +
  scale_fill_manual(
    values = ospc_colors
  ) +
  scale_color_manual(
    values = ospc_colors
  ) +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = list(c("A","K"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = plottitlesize, face = "bold"),      
    axis.title = element_text(size = axistitlesize),                  
    axis.text = element_text(size = axistextsize),               
    legend.title = element_text(size = legendtitlesize),   
    legend.text = element_text(size = legendtextsize)  
  )

# display boxplots side by side
grid.arrange(box_rst, box_ospc, ncol = 2, widths = c(0.3, 0.7))




