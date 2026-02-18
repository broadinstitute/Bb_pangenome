#R scripts written by M. Amine Hassani - hassani.medamine@gmail.com
##to plot full average nucleotide identity heatmap - fig. 3a
#required packages 
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMRplus","magick",
"ComplexHeatmap","DECIPHER","philentropy", "corrplot","agricolae")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c("#b02b76","#ffa500","#00b2b2","#155832","#4f7b95")
shape.1<-c(15,17,16,19,13)

list.1<-c("ESI425H","ESI36H","ASM1913465v1","ASM4079078v1","ASM4079076v1","ASM336729v1",
"URI34H","UNY193P","UNY203P","UCT35H","UNY1128P","UNY1083P","UWI248P","UWI247P",
"UCT113H","URI88H","UNY208P","UNY1090P","ASM4079073v1","UWI263P","ASM4079071v1",
"UCT92H","UCT32H","URI47H","UNY149P","UCT110H","URI46H","URI44H","ASM4079079v1",
"URI117H","URI103H","UCT30H","URI86H","URI112H","URI118H","UNY1038P","ASM2466217v1",
"URI33H","URI48H","UNY990P","UWI283P","XYZ459H","UNY172P","UNY1085P","ASM4079075v1",
"B418PP","UCT50H","UNY169P","ASM4079074v1","URI36H","URI56H","B500PP","B331PP","ASM2466215v1",
"PFhe_I_PB_Ill_cons","ESI403H","ESI26H","UCT96H","UCT124H","UNY1032P","URI101H","UCT29H",
"URI87H","URI107H","URI93H","URI102H","UCT109H","URI39H","URI89H","URI41H","URI42H","URI40H",
"URI111H","UCT31H","URI120H","URI91H","ASM4079077v1","ASM2466219v1","ASM4079080v1","ASM215150v1",
"ASM215148v1","ASM215146v1")

set.seed(123)
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1, size=10, color="black")
	)
# upload and prepare phyloseq objects***
#	mat1=read.delim("ANI_pi_LD83.txt", sep="\t", row.names=1, header=T)
	mat2=read.delim("ANI_fpi_LD83.txt", sep="\t", row.names=1, header=T)

#	mat=as.matrix(mat1)
	mat=as.matrix(mat2)

	mat[mat < 0.95] <- NA
	mat=mat-0.97
	
#	pdf("LD83_SGCs_rooted.pdf")
	TREE=read.tree("SGCs_LD83_rooted.tre")
#	dev.off()

	dend <- ReadDendrogram("SGCs_LD83_rooted.tre")
	
	hm1=Heatmap(mat, name = "ANI", use_raster=F,show_column_dend=T,col=c("#ffbaba","#a70000"), row_dend_width= unit(2,"cm"), row_names_gp=gpar(fontsize=5), cluster_rows=dend, show_column_names=F,cluster_columns=dend, width = unit(13, "cm"), height=unit(13,"cm"),rect_gp = gpar(col = "white", lwd = 2),na_col = "white")
	
	hm2=Heatmap(mat, name = "ANI", use_raster=F,show_column_dend=T,col=c("#6565f7","#de8a0d","#d42643"), row_dend_width= unit(2,"cm"), row_names_gp=gpar(fontsize=5), cluster_rows=dend, show_column_names=F,cluster_columns=dend, width = unit(13, "cm"), height=unit(13,"cm"),rect_gp = gpar(col = "white", lwd = 2),na_col = "white")
	
#	pdf("ANI_heatmap_plot.pdf", useDingbats=FALSE)
	hm2
#	dev.off()
	
	sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
	rownames(sd)<-sd$ID
	sd<-(sd[-1])
	sd2=sd[,c(1,5,7,9,10)]


	hm3=Heatmap(as.matrix(sd2[,1]), cluster_rows=dend, width = unit(3, "cm"), use_raster=F, row_names_gp=gpar(fontsize=5)  )
	
#	pdf("ANI_heatmap_meta_plot.pdf", useDingbats=FALSE)
	hm3
#	dev.off()

	hm1+hm2+hm3
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
