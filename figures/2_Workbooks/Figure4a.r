########## Figure 4a ##########
# This notebook constructs the ternary plot shown in Figure 4a.
# It also produces Figure S11.

library(dplyr)
library(ggplot2)
library(ggtern)

# input file
ternary_file = "../0_Data/2_Processed/ternary_data.csv"

# import data
tern <- read.csv(ternary_file, header=TRUE, sep=",")
tern <- tern[order(tern$gene), ]

# set up a way to plot lines from edge midpoints to the triangle's center
edge_to_center <- data.frame(
  x = c(0.5, 0.5, 0), 
  y = c(0.5, 0, 0.5), 
  z = c(0, 0.5, 0.5), 
  xend = c(1/3, 1/3, 1/3),
  yend = c(1/3, 1/3, 1/3),
  zend = c(1/3, 1/3, 1/3)
)

##### Figure 4a #####
ggtern(tern, aes(x = RST2, y = RST1, z = RST3)) +
  geom_point(aes(size=gene, fill=gene, stroke=surface_lipo),
             alpha=0.75,
             shape=21,
             position= position_jitter_tern(x=0.025, y=0.025, z=0.025),
             show.legend=FALSE) +
  scale_color_manual(values="black")+
  scale_discrete_manual(aesthetics="stroke", values=c(0,0.75))+
  scale_fill_manual(values=tern$color)+
  scale_size_manual(values=tern$size/5+1) +
  guides(color = "none",
         size="none")  +
  theme_showarrows() +
  theme(tern.panel.background = element_rect(fill = "grey75", color = NA),
        plot.title=element_text(size=28, face="bold", hjust=0.5),
        tern.axis.title.T = element_text(size = 24), 
        tern.axis.title.L = element_text(size = 24),
        tern.axis.title.R = element_text(size = 24),
        tern.axis.arrow.text.T = element_text(size = 22),
        tern.axis.arrow.text.L = element_text(size = 22),
        tern.axis.arrow.text.R = element_text(size = 22),
        tern.axis.text.T  = element_text(size = 20),
        tern.axis.text.L  = element_text(size = 20),
        tern.axis.text.R  = element_text(size = 20)
        )+
  geom_segment(data = edge_to_center,
               aes(x = x, y = y, z = z, xend = xend, yend = yend, zend = zend),
               color = "#808080", linetype = "dashed", size = 0.5)


##### Figure S11 #####
ggtern(tern, aes(x = PC2, y = PC1, z = PC3)) +
  geom_point(aes(size=gene, fill=gene, stroke=surface_lipo),
             alpha=0.75,
             shape=21,
             position= position_jitter_tern(x=0.025, y=0.025, z=0.025),
             show.legend=FALSE) +
  scale_color_manual(values="black")+
  scale_discrete_manual(aesthetics="stroke", values=c(0,0.75))+
  scale_fill_manual(values=tern$color)+
  scale_size_manual(values=tern$size/5+1) +
  labs(title = "Gene Group Distributions by 3D PCA") +
  guides(color = "none",
         size="none")  +
  theme_showarrows() +
  theme(tern.panel.background = element_rect(fill = "grey75", color = NA),
        plot.title=element_text(size=28, face="bold", hjust=0.5),
        tern.axis.title.T = element_text(size = 24),
        tern.axis.title.L = element_text(size = 24),
        tern.axis.title.R = element_text(size = 24),
        tern.axis.arrow.text.T = element_text(size = 22),
        tern.axis.arrow.text.L = element_text(size = 22),
        tern.axis.arrow.text.R = element_text(size = 22),
        tern.axis.text.T  = element_text(size = 20),
        tern.axis.text.L  = element_text(size = 20),
        tern.axis.text.R  = element_text(size = 20)
  )+
  geom_segment(data = edge_to_center,
               aes(x = x, y = y, z = z, xend = xend, yend = yend, zend = zend),
               color = "#808080", linetype = "dashed", size = 0.5)