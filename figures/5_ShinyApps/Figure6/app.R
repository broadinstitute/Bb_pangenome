########## Figure 6 - Interactive ##########
# This code is for an RShiny app which produces an interactive version of Figure 6.
# The user can filter by replicon and / or gene prevalence threshold.

options(renv.config.ignore.Rprofile = TRUE)

library(shiny)
library(shinyWidgets)
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ggiraph)
library(dplyr)
library(tidyr)
library(phytools)
library(showtext)


# input files (all in the same folder as the app script)
tree_file = "core_gene_alignment.nwk"
metadata_file = "Lemieux2023_TableS2.tsv"
scaffolds_file = "/Users/rl275/Projects/bb_longread/Manuscript/0_Data/2_Processed/ScaffoldTree.csv"

# set legend / formatting parameters
size_legendtitle = 16
size_legendtext = 12
kw=1
kh=1

# midpoint root the tree file and clean node labels
tr <- midpoint.root(read.tree(tree_file))
tr$tip.label <- gsub("GCF_\\d+\\.\\d+_", "", tr$tip.label)
tr$tip.label <- gsub("_genomic", "", tr$tip.label)

# process and clean metadata
metada <- read.csv(metadata_file,
               header=TRUE,
               sep="\t")
names(metada)[names(metada) == "sample_name"] <- "ID"
metada$Disseminated[is.na(metada$Disseminated)] <- ""
metadata <- metada %>%
  select(c("ID", "OspC_Type", 
           "RST_Type",
           "MLST", "Disseminated"))
metadata$Strain <- metadata$ID
metadata$MLST <- factor(metadata$MLST)
metadata$RST_Type <- factor(metadata$RST_Type)
metadata$OspC_Type <- factor(metadata$OspC_Type)
metadata$Disseminated <- factor(metadata$Disseminated, levels=c("L", "D", ""), labels=c("Localized", "Disseminated", ""))

# import gene tree data
p <- ggtree(tr,size=0.25) 
p <- p %<+% metadata

# color nodes by RST
p <- p +
  geom_tippoint(mapping=aes(colour=RST_Type),
                size=0.25,
                stroke=0
  ) +
  scale_color_manual(
    name="RST",
    values=c("#FF0000", "#00FF00", "#0000FF"),
    guide=guide_legend(override.aes = list(label = "\u2022", size = 3), keywidth=kw/2, keyheight=kh/2, ncol=3, order=1)
  ) +
  theme(
    legend.title=element_text(size=size_legendtitle/2),
    legend.text=element_text(size=size_legendtext/2),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  new_scale_color()

# annotate nodes with labels colored by OspC type
p <- p +
  geom_tiplab(mapping=aes(label=Strain, colour=OspC_Type),
              align=TRUE,
              linetype=1,
              size=0.5,
              linesize=0.1,
              offset=0.005,
              show.legend=FALSE) +
  scale_color_discrete(
    name="OspC"
  ) +
  theme(
    legend.title=element_text(size=size_legendtitle/2),
    legend.text=element_text(size=size_legendtext/2),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  new_scale_color()

# add layer for OspC
p <- p +  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=OspC_Type),
    width=0.0015,
    offset=1.5
  ) +
  scale_fill_discrete(name="OspC",
                      guide=guide_legend(keywidth=kw/2, keyheight=kh/2, ncol=3, order=2)
  ) +
  theme(
    legend.title=element_text(size=size_legendtitle/2), 
    legend.text=element_text(size=size_legendtext/2),
    legend.spacing.y = unit(0.02, "cm")
  )

# add layer for dissemination status
p_diss <- p +  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Disseminated),
    width=0.0015,
    offset=0.25
  ) +
  scale_fill_manual(name="Disseminated",
                    values=c("lightgreen", "darkgreen", "white"),
                    guide=guide_legend(keywidth=kw/2, keyheight=kh/2, ncol=3, order=3)
  ) +
  theme(
    legend.title=element_text(size=size_legendtitle/2), 
    legend.text=element_text(size=size_legendtext/2),
    legend.spacing.y = unit(0.02, "cm")
  )


# import scaffold tree data
genes <- read.csv(scaffolds_file,
                  header=TRUE,
                  sep=",")

names(genes)[names(genes) == "assembly"] <- "ID"

gene_pa <- genes[,c("gene", "ID")]
gene_pa <- gene_pa %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = gene,
    values_from = presence,
    values_fill = list(presence = 0)
  ) %>% 
  as.data.frame()
rownames(gene_pa) <- gene_pa$ID
gene_pa$ID <- NULL

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

# set up UI for shiny app
ui <- fluidPage(
  tags$head(
    tags$style(HTML(paste0("
    #reps .checkbox { margin-bottom: 1px !important; }
    #reps .checkbox label { font-size: ", size_legendtext, "px !important; }
    #reps .checkbox input[type='checkbox'] {
      transform: scale(0.8) !important;
      margin-right: 1px !important;
    }
  ")))
  ),
  titlePanel("Borrelia Phylogeny and Gene Presence by Plasmid ~ Scaffolds"),
  # checkbox menu for selecting replicons
  sidebarLayout(
    sidebarPanel(
      width = 2,
      checkboxGroupInput(
        "reps", "Replicons",
        choices = sort(unique(genes$replicon_name)),
        selected = unique(genes$replicon_name)
      )
    ),
    mainPanel(
      width = 10,
      girafeOutput("gtree", height = "900px")
    )
  ),
  # slider for setting gene prevalence threshold
  fluidRow(
    column(
      width = 12,
      br(),
      sliderTextInput(
        inputId = "prev_thresh",
        label = "Gene Prevalence Threshold",
        choices = sprintf("%.2f", seq(0, 1, 0.1)),
        selected = "0.00",
        grid = TRUE,
        width = "100%"
      )
    )
  )
)

# set up the interactive plot
server <- function(input, output, session) {
  # filter data based on the selected replicons and threshold value
  output$gtree <- renderGirafe({
    genes_filtered <- genes |>
      filter(
        replicon_name %in% input$reps,
        prevalence >= input$prev_thresh
      ) |>
      mutate(gene = factor(gene, levels = unique(gene))) |>
      mutate(
        genesafe = gsub("['\"]", "_", gene)
      )
  # add heatmap to the tree
  p_genes <- p_diss +  new_scale_fill() +
      geom_fruit(
        data=genes_filtered,
        geom=geom_tile_interactive,
        mapping=aes(x=gene,
                    y=ID,
                    fill=replicon_name,
                    tooltip=gene,
                    data_id = genesafe
        ),
        pwidth=10,
        offset=4
      ) +
      scale_fill_manual(
        name="Plasmid",
        values = plasmid_colors,
        na.translate = FALSE,
        guide = guide_legend(
          keywidth = kw/2,
          keyheight = kw/2,
          ncol = 2,
          order = 3
        )
      ) +
      scale_x_ggtree() +
      theme(plot.title=element_text(size=24, face="bold"),
            legend.title=element_text(size=size_legendtitle/2),
            legend.text=element_text(size=size_legendtext/2)) +
      new_scale_fill()
  # interactive widget area 
  girafe(
      ggobj = p_genes,
      options = list(
        opts_selection(type = "single", only_shiny = FALSE),
        opts_zoom(max = 20),
        opts_toolbar(saveaspng = TRUE)
      )
    )
  })
}

# run shiny app
shinyApp(ui, server)

