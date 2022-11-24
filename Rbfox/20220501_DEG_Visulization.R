library("Rsubread")
library("Rsamtools")
library("pheatmap")
library("ggplot2")
library("edgeR")
library("limma")
library("biomaRt")
library("GOstats")
library("DESeq2")
library("xlsx")

rm(list=ls())

######load and merge data from rawcounts, rpkm, and DESeq2 results#####
AllgeneRPKM<-read.csv("20220408_RPKM_Allgene_Annotation.csv")
rownames(AllgeneRPKM)<-AllgeneRPKM$entrezgene_id

DESeq2All<-read.csv("20220328_DEseq2_MUTvsCtrl_All_Annotation.csv")
rownames(DESeq2All)<-DESeq2All$entrezgene_id
  
Alllist<-merge(x = AllgeneRPKM,y = DESeq2All,by='row.names', all=TRUE)
colnames(Alllist)

Alllist<-Alllist[,-c(1,2,21:29)]


####pheatmap and PCA DEGs###
Alllist.noNA<-Alllist[!is.na(Alllist$log2FoldChange) & !is.na(Alllist$padj),]
DEGs<-Alllist.noNA[abs(Alllist.noNA$log2FoldChange) >1 & Alllist.noNA$padj<0.05,]
dim(DEGs)
colnames(DEGs)
DEGs$ctrlmeanRpkm<-apply(DEGs[,13:15],1,mean)

DEGs[order((DEGs[,"ctrlmeanRpkm"]),decreasing=FALSE),]

ColAnno<-data.frame("Group"=c(rep("Ctrl",3),rep("Mut",3)))
rownames(ColAnno)<-colnames(DEGs[,13:18])
ColannoColor<-list(Group=c(Ctrl="#d5d4d4", Mut="#ff84fd"))  
plotheatmap<-pheatmap(DEGs[,13:18], treeheight_col = 20,
         cluster_rows = T,border_color = "gray",show_rownames = F,
         cluster_cols = T,scale = "row",cellwidth = 20,treeheight_row = 0,
         cellheight = 0.15,display_numbers = F,annotation_col = ColAnno,
         annotation_colors = ColannoColor,
         color = colorRampPalette(c("#276acc", "white", "#f45864"))(50))

ggsave(filename = "20220501_DEGs_heatmap.pdf",plot = plotheatmap,
       device = "pdf",width = 5,height = 12)

# transpose the data to have variables (genes) as columns
data_for_PCA <- t(DEGs[,13:18])
dim(data_for_PCA)
## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned
# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")
## calculate MDS
mds <- cmdscale(dist(data_for_PCA))
#Samples representation
plot(mds[,1], -mds[,2], type="p", xlab=paste("Dimension1_",eig_pc[1],sep=""), 
     ylab=paste("Dimension2_",eig_pc[2],sep=""), main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 

#####Volcano plot with DEGs result####
# add a column of NAs
Alllist.noNA$diffexpressed <- "No"
# set "UP" 
Alllist.noNA$diffexpressed[Alllist.noNA$log2FoldChange > 1 & Alllist.noNA$padj < 0.05] <- "Up"
# set "DOWN"
Alllist.noNA$diffexpressed[Alllist.noNA$log2FoldChange < -1  & Alllist.noNA$padj < 0.05] <- "Down"
sum(Alllist.noNA$diffexpressed=="No")
sum(Alllist.noNA$diffexpressed=="Up")
sum(Alllist.noNA$diffexpressed=="Down")

pVol<-ggplot(Alllist.noNA,aes(x=log2FoldChange,y=-log10(padj),col=diffexpressed))+
  geom_point(alpha=12,size=2,pch=16) +
  #geom_point(x=-Allgenes[Allgenes$X=="22433",]$logFC,
  #            y=-log10(Allgenes[Allgenes$X=="22433",]$adj.P.Val),color="red",size=3)+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  scale_x_continuous(limits = c(-10,10))+
 scale_y_continuous(limits = c(0,20))
# geom_text_repel(
#  data = Allgenes[Allgenes$X == c("22433"),],
# aes(label = "Xbp1"),
# size = 5,color=c("black")
# )

ggsave(filename = "20220501_VolcanoPlot.pdf",plot = pVol,
       device = "pdf",width = 6,height = 8)


#####Chaperon gene expression######
ChaperonGenes<-read.xlsx("ChaperonGenes.xlsx",sheetIndex = 1,header = T)
rownames(ChaperonGenes)<-ChaperonGenes$gene
ChaperonGenesExpr<-AllgeneRPKM[AllgeneRPKM$mgi_symbol %in% ChaperonGenes$gene,]
rownames(ChaperonGenesExpr)<-ChaperonGenesExpr$mgi_symbol
colnames(ChaperonGenesExpr)
ChaperonGenesExpr<-merge(ChaperonGenesExpr,ChaperonGenes,by=0,all=T)
rownames(ChaperonGenesExpr)<-ChaperonGenesExpr$mgi_symbol

RowAnno<-data.frame("Family"=ChaperonGenesExpr$family)
rownames(RowAnno)<-rownames(ChaperonGenesExpr)
ColAnno<-data.frame("Group"=c(rep("Ctrl",3),rep("Mut",3)))
rownames(ColAnno)<-colnames(ChaperonGenesExpr[,15:20])
ColannoColor<-list(Group=c(Ctrl="#d5d4d4", Mut="#ff84fd"))  

ChaperonGeneplot<-pheatmap(ChaperonGenesExpr[,15:20], treeheight_col = 20,
         cluster_rows = F,border_color = NA,show_rownames = T,
         cluster_cols = F,scale = "row",cellwidth = 20,treeheight_row = 20,
         cellheight = 8,display_numbers = F,annotation_row = RowAnno,
         annotation_col=ColAnno,annotation_colors = ColannoColor,
         color = colorRampPalette(c("#276acc", "white", "#f45864"))(50))

ggsave(filename = "20220501_ChaperonGeneplot.pdf",plot = ChaperonGeneplot,
       device = "pdf",width = 6,height = 8)

