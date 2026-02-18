#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com
##to reproduce the ordination plot based on Sorensen distances indicated in fig. 3c
#required packages
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMRplus","magick",
"ComplexHeatmap","DECIPHER","philentropy","ggtern","venn")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

set.seed(123)

shape.list1<-c(16,15,17)

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
#	legend.position="none"
	)

	mat=read.delim( "gc_table.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(t(mat))
	OTU=otu_table(mat, taxa_are_rows=T) 
#	tab2=read.delim( "gc_summary.txt", sep="\t", row.names=1, header=T)
#	DT.tab2=data.table(tab2[,-28], keep.rownames=T, key="gene_cluster_id")
#	DT2=unique(DT.tab2, by="gene_cluster_id")
#	DT3=DT2[DT2$num_genomes_gene_cluster_has_hits > "84"]
	tax=read.delim("GC_annotation.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
	rownames(sd)<-sd$ID
	sd<-(sd[-1])
	SD=sample_data(sd)

	physeq=phyloseq(OTU, TAXA, SD) 
#	break()
	dist_sor=distance(t(mat), method = "sorensen")
	colnames(dist_sor)=rownames(dist_sor)=rownames(t(mat))
	dist_Sor=as.dist(dist_sor)

	Ord=phyloseq::ordinate( physeq, "PCoA" , dist_Sor )
	
	pJC<-plot_ordination( physeq, Ord, color="OspC", shape="RST")
	pJC$layers<-pJC$layers[-1]
	
	p1=pJC+geom_point(size=4, alpha=0.7)+
	scale_shape_manual(values=shape.list1)+ 
	theme_bw()+theme_new

#	pdf("Figure4A.pdf", useDingbats=FALSE)
	print(p1)
#	dev.off()
