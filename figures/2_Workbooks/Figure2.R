########## Figure 2 ##########
# This notebook generates various phylogenetic trees with heatmaps displaying
# gene group presence according to encoding replicon.
# In addition to Figure 2, Figures S5, S7, and S8b are produced.

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
library(dplyr)
library(tidyr)
library(reshape2)
library(TDbook)
library(stringr)
library(gtools)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(shiny)
library(ggiraph)

# input files
tree_file = "../../output/results/v9/no_merge_paralogs/intermediates/fasttree/fasttree_panaroo_v1.newick"
metadata_file = "../0_Data/2_Processed/MetadataFull.csv"
fullgenes_file = "../0_Data/2_Processed/GeneTree_Ordered.csv" # generated in `Figure2.ipynb`
thresholdedgenes_file = "../0_Data/2_Processed/GeneTree_Ordered_Thresh.csv" # generated in `Figure2.ipynb`

# set legend / formatting parameters
size_legendtitle = 16
size_legendtext=12
kw=1
kh=1

# midpoint root the tree file and clean node labels
tr <- midpoint.root(read.tree(tree_file))
tr$tip.label <- gsub("GCF_\\d+\\.\\d+_", "", tr$tip.label)

# process and clean metadata
tr$tip.label <- gsub("_genomic", "", tr$tip.label)
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
           "Accessory", "Core"))
metadata$MLST <- factor(metadata$MLST)
metadata$RST <- factor(metadata$RST)
metadata$Imputed <- factor(metadata$Imputed)
metadata$OspC <- factor(metadata$OspC)
metadata$WGS <- factor(metadata$WGS, levels=c("A", "B.1", "B.2", "C", ""))
metadata$BAPS <- factor(metadata$BAPS)
metadata$Accessory <- factor(metadata$Accessory)
metadata$Core <- factor(metadata$Core)

# construct the base tree
p <- ggtree(tr,size=0.5) 
p <- p %<+% metadata

# color nodes by RST
p <-p +
  geom_tippoint(mapping=aes(colour=RST),
                size=3,
                stroke=0
  ) +
  scale_color_manual(
    name="RST",
    values=c("#FF0000", "#00FF00", "#0000FF"),
    guide=guide_legend(override.aes = list(label = "\u2022", size = 5), keywidth=kw, keyheight=kh, ncol=1, order=1)
  ) +
  theme(
    legend.title=element_text(size=size_legendtitle),
    legend.text=element_text(size=size_legendtext),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  new_scale_color()

# annotate nodes with labels colored by OspC type
p <- p +
  geom_tiplab(mapping=aes(label=Strain, colour=OspC),
              align=TRUE,
              linetype=3,
              size=3,
              linesize=0.5,
              offset=0.002,
              show.legend=TRUE) +
  scale_color_discrete(
    name="OspC",
    guide=guide_legend(keywidth=kw, keyheight=kh, ncol=2, order=2)
  ) +
  theme(
    legend.title=element_text(size=size_legendtitle),
    legend.text=element_text(size=size_legendtext),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  new_scale_color()
 

##### Full Gene Tree ##### 
genes <- read.csv(fullgenes_file,
                   header=TRUE,
                   sep=",")
 
names(genes)[names(genes) == "assembly"] <- "ID"

all_replicons = unique(genes$replicon_name)
best_replicons = unique(genes$best_replicon)
non_best_replicons <- setdiff(all_replicons, best_replicons)
combined_replicons <- c(best_replicons, non_best_replicons)
 
genes$replicon_name <- factor(genes$replicon_name, levels=combined_replicons)
genes$gene <- factor(genes$gene, levels=unique(genes$gene))
 
colors <- genes %>%
   select(c("replicon_name", "color")) %>%
   distinct()
colors <- colors %>%
  arrange(ifelse(replicon_name == "Unclassified", 1, 0))
plasmid_colors <- setNames(colors$color, colors$replicon_name)
 
# add heatmap with gene presence by replicon
p_genes <- p +
   geom_fruit(
     data=genes,
     geom=geom_tile,
     mapping=aes(x=gene, y=ID, fill=replicon_name),
     pwidth=5,
     offset=0.4
   ) +
   scale_fill_manual(
     name="Plasmid Composition",
     values=plasmid_colors,
     na.translate=FALSE,
     guide=guide_legend(keywidth=0.5,
                        keyheight=0.5,
                        ncol=2,
                        order=3
     )
   ) +
   theme(plot.title=element_text(size=24, face="bold"),
         legend.title=element_text(size=12),
         legend.text=element_text(size=10)) +
   ggtitle("Borrelia Phylogeny and Gene Presence by Plasmid") +
   new_scale_fill()
 p_genes
 

 
##### Figure S8B #####
# set formatting parameters
plottitlesize <- 24
axistitlesize <- 16
axistextsize <- 14
legendtitlesize <- 16
legendtextsize <- 14
 
# subset the full gene file for lipoproteins and aggregate at the isolate level
genes_liposubset <- genes[genes$lipo == "Lipoprotein", ]
lipo_counts <- genes_liposubset %>%
   group_by(ID) %>%
   summarise(Num_Lipos = n(), .groups = "drop") %>%
   left_join(metadata, by="ID")
lipo_counts
 
# lipoprotein count by RST
box_rst <- ggplot(lipo_counts, aes(x = RST, y = Num_Lipos)) +
   geom_boxplot(aes(fill=RST),outlier.shape = NA) +
   geom_dotplot(
     aes(x=RST),
     fill="black", binaxis = 'y', stackdir = 'center', dotsize = 0.35, stackratio = 1.5
   ) +
   geom_text_repel(
     data = lipo_counts %>%
       group_by(RST) %>%
       filter(Num_Lipos > quantile(Num_Lipos, 0.75) + 1.5 * IQR(Num_Lipos) |
                Num_Lipos < quantile(Num_Lipos, 0.25) - 1.5 * IQR(Num_Lipos)),
     aes(label = Strain, color=RST), 
     vjust = -1, size = 3.5, show.legend=FALSE
   ) +
   labs(
     title = "Lipoprotein Count by RST",
     x = "RST",
     y = "Lipoprotein Count"
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
 
# lipoprotein count by OspC type
ospc_types <-sort(unique(lipo_counts$OspC))
ospc_colors = hue_pal()(length(ospc_types))
df_ospc_colors <- data.frame(
   OspC = ospc_types,
   Color = ospc_colors
)
ospc_colors <- setNames(df_ospc_colors$Color, df_ospc_colors$OspC)
 
box_ospc <- ggplot(lipo_counts, aes(x = OspC, y = Num_Lipos)) +
   geom_boxplot(aes(fill=OspC),outlier.shape = NA) +
   geom_dotplot(
     aes(OspC),
     fill="black", binaxis = 'y', stackdir = 'center', dotsize = 0.35, stackratio = 1.5
   ) +
   geom_text_repel(
     data = lipo_counts %>%
       group_by(OspC) %>%
       filter(Num_Lipos > quantile(Num_Lipos, 0.75) + 1.5 * IQR(Num_Lipos) |
                Num_Lipos < quantile(Num_Lipos, 0.25) - 1.5 * IQR(Num_Lipos)),
     aes(label = Strain, color=OspC), 
     vjust = -1, size = 4, show.legend=FALSE
   ) +
   labs(
     title = "Lipoprotein Count by OspC Type",
     x = "OspC",
     y = "Lipoprotein Count"
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
 
##### Plasmid-specific gene trees (e.g. Figure S5) #####
add_heatmap <- function(base_tree, replicon_name) {
  # filter down to just the genes that are present on the specified replicon
  subset <- genes[genes$replicon == replicon_name, ]
  
  # if a replicon has only one gene, skip the tree generation
  if (dplyr::n_distinct(subset$gene) < 2) {
    message(sprintf("Skipping %s: fewer than 2 unique genes.", replicon_name))
    return(base_tree)
  }

  # reshape data to a wide format
  mat <- subset |>
    dplyr::select(ID, gene) |>
    dplyr::distinct() |>        
    dplyr::mutate(x = 1) |>
    tidyr::pivot_wider(names_from = gene,
                       values_from = x,
                       values_fill = 0)
    m <- as.matrix(mat[, -1])
  
  # perform hierarchical clustering on columns
  gene_dist <- dist(t(m), method = "euclidean")
  gene_hc <- hclust(gene_dist, method = "average")
  gene_order <- gene_hc$labels[gene_hc$order]
  subset$gene <- factor(subset$gene, levels = gene_order)
  
  # set plasmid color assignments
  plasmid_colors <- subset %>%
    dplyr::select(c("replicon_name", "color")) %>%
    distinct()
  plasmid_colors <- plasmid_colors[order(plasmid_colors$replicon_name), "color"]
  
  # append plasmid-specific heatmap data to tree
  base_tree +
    geom_fruit(
      data=subset,
      geom=geom_tile,
      mapping=aes(x=gene, y=ID, fill=replicon_name),
      pwidth=3,
      offset=0.6
    ) +
    scale_fill_manual(
      name="Plasmid Composition",
      values=plasmid_colors,
      na.translate=FALSE,
      guide=guide_legend(keywidth=0.5,
                         keyheight=0.5,
                         order=3)
    ) +
    theme(plot.title=element_text(size=16, face="bold"),
          legend.title=element_text(size=12),
          legend.text=element_text(size=10)) +
    ggtitle(sprintf("Gene Presence in %s for Long-Read Isolates", replicon_name)) +
    new_scale_fill()
  
  # save figure
  file_name = paste("../3_Figures/Plasmid Gene Trees/gene_", replicon_name, ".png", sep="")
  ggsave(file_name, width=3000, height=2250, units="px", dpi=300)
}

# iterate through plasmids to output heatmaps
for (r in unique(genes$replicon_name)){
  print(r)
  add_heatmap(p, r)
}
 
##### Figure 2 #####
# import thresholded tree input data and filter to non-chromosomal genes
genesthresh <- read.csv(thresholdedgenes_file,
                        header=TRUE,
                        sep=",")
names(genesthresh)[names(genesthresh) == "assembly"] <- "ID"
genesthresh <- genesthresh %>% filter(best_replicon != "chromosome")
genesthresh$replicon_name <- factor(genesthresh$replicon_name, levels=combined_replicons)
genesthresh$gene <- factor(genesthresh$gene, levels=unique(genesthresh$gene))

# build tree with all replicons other than the chromosome
p_genesthresh <- p +
  geom_fruit(
    data=genesthresh,
    geom=geom_tile,
    mapping=aes(x=gene, y=ID, fill=replicon_name),
    pwidth=25,
    offset=1.75,
  ) +
  scale_fill_manual(
    name="Plasmid Composition",
    values=plasmid_colors,
    na.translate=FALSE,
    guide=guide_legend(keywidth=kw,
                       keyheight=kw,
                       ncol=2,
                       order=3
    )
  ) +
  scale_x_ggtree() + 
  theme(plot.title=element_text(size=24, face="bold"),
        legend.title=element_text(size=size_legendtitle),
        legend.text=element_text(size=size_legendtext)) +
  ggtitle("Borrelia Phylogeny and Gene Presence by Plasmid") +
  new_scale_fill()
p_genesthresh

##### File S1 (larger version with x-axis labels) #####
p_genesthresh <- p +
  geom_fruit(
    data=genesthresh,
    geom=geom_tile,
    mapping=aes(x=gene, y=ID, fill=replicon_name),
    pwidth=25,
    offset=1.75,
    axis.params = list(axis="x",
                       text.size=2,
                       text.angle=90,
                       hjust=0,
                       vjust=0.5)
  ) +
  scale_fill_manual(
    name="Plasmid Composition",
    values=plasmid_colors,
    na.translate=FALSE,
    guide=guide_legend(keywidth=kw/2,
                       keyheight=kw/2,
                       ncol=2,
                       order=3
    )
  ) +
  scale_x_ggtree() + 
  theme(plot.title=element_text(size=75, face="bold"),
        legend.title=element_text(size=40),
        legend.text=element_text(size=30)) +
  ggtitle("Borrelia Phylogeny and Gene Presence by Plasmid") +
  new_scale_fill()
p_genesthresh

##### Figure S7 #####
# import data and filter to non-chromosomal lipoproteins
lipos <- read.csv(thresholdedgenes_file,
                  header=TRUE,
                  sep=",")
names(lipos)[names(lipos) == "assembly"] <- "ID"
lipos <- lipos %>% filter(best_replicon != "chromosome")
lipos <- lipos %>% filter(lipo == "Lipoprotein")
lipos$replicon_name <- factor(lipos$replicon_name, levels=combined_replicons)
lipos$gene <- factor(lipos$gene, levels=unique(lipos$gene))

# build tree with lipoproteins that are not typically located on the chromosome
p_lipos <- p +
  geom_fruit(
    data=lipos,
    geom=geom_tile,
    mapping=aes(x=gene, y=ID, fill=replicon_name,
    ),
    pwidth=5,
    offset=0.4,
    axis.params = list(axis="x",
                       text.size=2.5,
                       text.angle=90,
                       hjust=0,
                       vjust=0.5)
  ) +
  scale_fill_manual(
    name="Plasmid Composition",
    values=plasmid_colors,
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       ncol=2,
                       order=3
    )
  ) +
  theme(plot.title=element_text(size=24, face="bold"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) +
  ggtitle("Borrelia Phylogeny and Lipoprotein Presence by Plasmid") +
  new_scale_fill()
p_lipos

