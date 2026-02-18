########## Figure 1a ##########
# This notebook constructs the typing scheme and plasmid presence phylogenetic tree shown in Figure 1.
# It also contains code for Figure S8c.

library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(gridExtra)
library(ggstar)
library(ggnewscale)
library(ggpubr)
library(tidytree)
library(treeio)
library(phytools)
library(phyloseq)
library(tidyr)
library(dplyr)
library(reshape2)
library(TDbook)
library(stringr)
library(gtools)
library(viridis)
library(ggrepel)
library(scales)
library(RColorBrewer)

##### Figure 1a #####
# input files
tree_file = "../../output/results/v9/no_merge_paralogs/fasttree/fasttree.newick"
metadata_file = "../0_Data/2_Processed/MetadataFull.csv"
plasmidbinary_file = "../0_Data/2_Processed/plasmids_binary_long.csv" # generated in `Figure1a.ipynb`

# set legend / formatting parameters
ncol=6
ts_width = 0.00075
ts_offset = 0.115
kh = 0.25
kw = 0.6

# midpoint root the tree file and clean node labels
tr <- midpoint.root(read.tree(tree_file))
tr$tip.label <- gsub("GCF_\\d+\\.\\d+_", "", tr$tip.label)
tr$tip.label <- gsub("_genomic", "", tr$tip.label)

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
names(metada)[names(metada) == "Isolate"] <- "ID"
old_isolates <- c(metada$id)
metadata <- metada %>%
  select(c("ID", "Location", "Strain", 
           "MLST",
           "RST", "Imputed",
           "OspC", "WGS", "BAPS",
           "Accessory", "Core", "DBScan"))
metadata$Location <- factor(metadata$Location, levels=c("Connecticut", "Rhode Island", "New York", "Wisconsin", "Slovenia"))
metadata$MLST <- factor(metadata$MLST)
metadata$RST <- factor(metadata$RST)
metadata$Imputed <- factor(metadata$Imputed)
metadata$OspC <- factor(metadata$OspC)
metadata$WGS <- factor(metadata$WGS, levels=c("A", "B.1", "B.2", "C", ""))
metadata$BAPS <- factor(metadata$BAPS)
metadata$Accessory <- factor(metadata$Accessory)
metadata$Core <- factor(metadata$Core)
metadata$DBScan <- factor(metadata$DBScan)


# construct the base tree
p <- ggtree(tr,
            size=0.5)
p <- p %<+% metadata

# color nodes by RST
p1 <-p +
  geom_tippoint(mapping=aes(colour=RST, shape=Imputed),
                size=3,
                stroke=0,
                show.legend=FALSE
  ) +
  scale_color_manual(
    name="RST",
    values=c("#FF0000", "#00FF00", "#0000FF"),
    guide=guide_legend(override.aes = list(label = "\u2022", size = 3), keywidth=0.3, keyheight=0.3, ncol=3, order=1)
  ) +
  scale_size_manual(
    name = "Imputed",
    values = c(16,17,18),
    guide = guide_legend(
      override.aes = list(size = 3), 
      keywidth = 0.3, 
      keyheight = 0.3, 
      ncol = , 
      order = 1
    )
  ) +
  new_scale_color()

# annotate nodes with labels colored by OspC type
p2 <- p1 +
  geom_tiplab(mapping=aes(label=Strain, colour=OspC),
              align=TRUE,
              linetype=3,
              size=3,
              linesize=0.2,
              offset=0.001,
              show.legend=FALSE) +
  scale_color_discrete(
    name="OspC",
  ) +
  new_scale_color()

# add layer for RST
p3 <-p2 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=RST),
    width=ts_width,
    offset=0.75
  ) +
  scale_fill_manual(
    name="RST",
    values=c("#FF0000", "#00FF00", "#0000FF"),
    guide=guide_legend(keywidth=0.5, keyheight=0.3, ncol=3, order=1)
  ) +
  theme(
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")
  )

# add layer for OspC
p4 <-p3 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=OspC),
    width=ts_width,
    offset=ts_offset
  ) +
  scale_fill_discrete(
    name="OspC",
    guide=guide_legend(keywidth=0.5, keyheight=0.3, ncol=ncol, order=3)
  ) +
  theme(
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")
  )

# add layer for MLST type
p5 <-p4 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=MLST),
    width=ts_width,
    offset=ts_offset
  ) +
  scale_fill_discrete(name="MLST",
                    guide=guide_legend(keywidth=0.5, keyheight=0.3, ncol=ncol, order=4)
  ) +
  theme(
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")
  )

# add layer for BAPS
p6 <-p5 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=BAPS),
    width=ts_width,
    offset=ts_offset*1.75
  ) +
  scale_fill_viridis(name="BAPS",
                     option="mako",
                     discrete=TRUE,
                     direction=-1,
                    guide=guide_legend(keywidth=0.5, keyheight=0.3, ncol=ncol, order=5)
  ) +
  theme(
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")
  )

# add layer for PopPUNK Core
p7 <-p6 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Core),
    width=ts_width,
    offset=ts_offset
  ) +
  scale_fill_viridis(name="Core",
                     option="mako",
                     discrete=TRUE,
                     direction=-1,
                     guide=guide_legend(keywidth=0.5, keyheight=0.3, ncol=ncol, order=6)
  ) +
  theme(
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")
  )

# add layer for PopPUNK Accessory
p8 <-p7 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Accessory),
    width=ts_width,
    offset=ts_offset
  ) +
  scale_fill_viridis(name="Accessory",
                     option="mako",
                     discrete=TRUE,
                     direction=-1,
                     guide=guide_legend(keywidth=0.5, keyheight=0.3, ncol=ncol, order=7)
  ) +
  theme(
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")
  )

# add layer for DBScan
p9 <-p8 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=DBScan),
    width=ts_width,
    offset=ts_offset
  ) +
  scale_fill_viridis(name="DBScan",
                     option="mako",
                     discrete=TRUE,
                     direction=-1,
                     guide=guide_legend(keywidth=0.5, keyheight=0.3, ncol=ncol, order=8)
  ) +
  theme(
    legend.title=element_text(size=10), 
    legend.text=element_text(size=8),
    legend.spacing.y = unit(0.02, "cm")
  )

# import and process plasmid presence file
plasmids <- read.csv(plasmidbinary_file,
                     header=TRUE,
                     sep=",")
names(plasmids)[names(plasmids) == "Isolate"] <- "ID"
plasmids$plasmid <- factor(plasmids$plasmid, levels=unique(plasmids$plasmid))
plasmids$presence <- factor(plasmids$presence, levels=unique(plasmids$presence))

plasmid_colors <- plasmids %>%
  select(c("plasmid", "color")) %>%
  distinct()
plasmid_colors <- plasmid_colors[order(plasmid_colors$plasmid, decreasing=FALSE),"color"]
plasmid_colors

# add heatmap depicting plasmid presence
p10 <- p9 +
  new_scale_fill() +
  geom_fruit(
    data=plasmids,
    geom=geom_tile,
    mapping=aes(x=plasmid, y=ID, fill=presence),
    color="white",
    pwidth=2,
    offset=ts_offset*0.75
  ) +
  scale_fill_manual(
    name="Plasmid Composition",
    values=plasmid_colors,
    na.translate=FALSE,
    guide=guide_legend(keywidth=kw,
                       keyheight=kh,
                       ncol=3,
                       order=8
    )
  ) +
  theme(plot.title=element_text(size=24, face="bold"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9),) +
  ggtitle("Typing Schemes and Plasmid Presence for Long-Read Isolates")
p10

# NB: save the resulting image as a pdf that is 20in x 12in

##### Figure S8C #####
# set formatting parameters
plottitlesize <- 24
axistitlesize <- 16
axistextsize <- 14
legendtitlesize <- 16
legendtextsize <- 14

# calculate plasmid counts for each isolate
plasmid_counts <- plasmids %>%
  group_by(ID) %>%
  summarise(Num_Plasmids = n(), .groups = "drop") %>%
  left_join(metadata, by="ID")

# plasmid counts by RST
box_rst <- ggplot(plasmid_counts, aes(x = RST, y = Num_Plasmids)) +
  geom_boxplot(aes(fill=RST),outlier.shape = NA) +
  geom_dotplot(
    aes(x=RST),
    fill="black", binaxis = 'y', stackdir = 'center', dotsize = 0.35, stackratio = 1.5
  ) +
  geom_text_repel(
    data = plasmid_counts %>%
      group_by(RST) %>%
      filter(Num_Plasmids > quantile(Num_Plasmids, 0.75) + 1.5 * IQR(Num_Plasmids) |
               Num_Plasmids < quantile(Num_Plasmids, 0.25) - 1.5 * IQR(Num_Plasmids)),
    aes(label = Strain, color=RST), 
    vjust = -1, size = 3.5, show.legend=FALSE
  ) +
  labs(
    title = "Plasmid Count by RST",
    x = "RST",
    y = "Plasmid Count"
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

# plasmid counts by OspC type
ospc_types <-sort(unique(plasmid_counts$OspC))
ospc_colors = hue_pal()(length(ospc_types))
df_ospc_colors <- data.frame(
  OspC = ospc_types,
  Color = ospc_colors
)
ospc_colors <- setNames(df_ospc_colors$Color, df_ospc_colors$OspC)

box_ospc <- ggplot(plasmid_counts, aes(x = OspC, y = Num_Plasmids)) +
  geom_boxplot(aes(fill=OspC),outlier.shape = NA) + 
  geom_dotplot(
    aes(OspC),
    fill="black", binaxis = 'y', stackdir = 'center', dotsize = 0.35, stackratio = 1.5
  ) +
  geom_text_repel(
    data = plasmid_counts %>%
      group_by(OspC) %>%
      filter(Num_Plasmids > quantile(Num_Plasmids, 0.75) + 1.5 * IQR(Num_Plasmids) |
               Num_Plasmids < quantile(Num_Plasmids, 0.25) - 1.5 * IQR(Num_Plasmids)), 
    aes(label = Strain, color=OspC), 
    vjust = -1, size = 3.5, show.legend=FALSE
  ) +
  labs(
    title = "Plasmid Count by OspC Type",
    x = "OspC",
    y = "Plasmid Count"
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
