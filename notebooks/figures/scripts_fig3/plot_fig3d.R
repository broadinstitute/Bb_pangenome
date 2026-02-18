#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com
##to plot core vs accessory distances fig. 3d
#required packages
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMRplus","magick",
"ComplexHeatmap","DECIPHER","philentropy","ggtern","venn")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(123)

color.1<-c("#2bb065","#de8c00","#c77cff","#47b3da","#ff64b0","#ef5656")
color.2<-c("#b02b76","#2bb065","#47b3da","#ef5656")
shape.1<-c(15,17,16,19,13)
shape.2<-c(16,15,17)

theme_new1 <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1, size=10, color="black")
	)
theme_new2 <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none"
	)
theme_change <- theme(
#	legend.position="none",
	plot.background = element_blank( ),
	panel.grid.minor = element_blank( ),
	panel.grid.major = element_blank( ),
	)

# upload and prepare phyloseq objects***
	mat=read.delim( "gc_table.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(t(mat))
	OTU=otu_table(mat, taxa_are_rows=T) 
#	tab2=read.delim( "gc_summary.txt", sep="\t", row.names=1, header=T)
#	DT.tab2=data.table(tab2[,-28], keep.rownames=T, key="gene_cluster_id")
#	DT2=unique(DT.tab2, by="gene_cluster_id")
#	DT3=DT2[DT2$num_genomes_gene_cluster_has_hits > "84"]
#	tax=read.delim("GC_annotation.txt", sep="\t", row.names=1, header=1)
	tax=read.delim("LD83_gc_annotation.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
	rownames(sd)<-sd$ID
	sd<-(sd[-1])
	SD=sample_data(sd)


#	TREE=read.tree("SGCs_LD83_rooted.tre")
	dend <- ReadDendrogram("SGCs_LD83_rooted.tre")  
	dist_core <- cophenetic(dend)

	mat.core=as.matrix(dist_core)
	DT.core=data.table(reshape2::melt(mat.core), keep.rownames=T, key="rn")
	DT.core$Var3=paste(DT.core$Var1, DT.core$Var2,sep=".")
	setnames(DT.core, "value", "coph.dist")
	setkey(DT.core, "Var3")


	physeq=phyloseq(OTU, TAXA, SD) 
#	dist_JC=phyloseq::distance( physeq, "sor" , Binary=TRUE)

	physeq.filter = filter_taxa(physeq, function(x) sum(x) < (0.95*length(x)) , TRUE)

#	dist.acc=phyloseq::distance( physeq.filter, "jaccard" , Binary=TRUE)
	
	mat.filter=as(otu_table(physeq.filter), "matrix")
	dist.acc=distance(t(mat.filter), method = "jaccard")
	colnames(dist.acc)=rownames(dist.acc)=rownames(t(mat.filter))
	dist.Acc=as.dist(dist.acc)

	mat.acc=as.matrix(dist.Acc)
	DT.acc=data.table(reshape2::melt(mat.acc), keep.rownames=T, key="rn")
	DT.acc$Var3=paste(DT.acc$Var1, DT.acc$Var2,sep=".")
	setkey(DT.acc, "Var3")
	setnames(DT.acc, "value", "jc.dist")

	DT=merge(DT.core, DT.acc, by.x="Var3", by.y="Var3")
	DT<-DT[DT$coph.dist!="0" & DT$jc.dist!="0"]

	sd2=sd[,c(1,5,7,9,10)]
	DT.sd=data.table(sd2, keep.rownames=T, key="rn")
	
	DT=merge(DT, DT.sd, by.x="Var1.x", by.y="rn")
	setnames(DT, c("BAPS_Core_level1","OspC","RST","core","accessory"), 
	c("BAPSv1","OspCv1","RSTv1","corev1","accv1"))

	DT=merge(DT, DT.sd, by.x="Var2.y", by.y="rn")
	setnames(DT, c("BAPS_Core_level1","OspC","RST","core","accessory"),
	c("BAPSv2","OspCv2","RSTv2","corev2","accv2"))

	DT$rst=ifelse(DT$RSTv1=="RST1" & DT$RSTv2=="RST1", "RST1",
	ifelse(DT$RSTv1=="RST2" & DT$RSTv2=="RST2","RST2", 
	ifelse(DT$RSTv1=="RST3" & DT$RSTv2=="RST3", "RST3","interRST" )))

	DT$RST=ifelse(DT$RSTv1=="RST1" & DT$RSTv2=="RST1", "RST1",
	ifelse(DT$RSTv1=="RST2" & DT$RSTv2=="RST2","RST2", 
	ifelse(DT$RSTv1=="RST3" & DT$RSTv2=="RST3", "RST3",
	ifelse(DT$RSTv1 %in% c("RST1","RST2") & DT$RSTv2 %in% c("RST1","RST2"), "RST12", 
	ifelse(DT$RSTv1 %in% c("RST1","RST3") & DT$RSTv2 %in% c("RST1","RST3"), "RST13","RST23")))))

	pp.dist=ggplot(DT, aes(coph.dist, jc.dist ))
	
	
	p2=pp.dist+geom_jitter(aes(color=RST),alpha=.5,width = 2e-4)+
	theme_bw()+xlab("core genome divergence (cophenetic distances)")+ 
	ylab("accessory genome divergence (Jaccard distances)")+
	geom_smooth(method="glm", color="darkblue")+
	scale_color_manual(values=color.1)+
	theme(legend.position="top")+
	theme_new2

#	pdf("Figure4B.pdf",useDingbats=FALSE)
	print(p2)
#	dev.off()
