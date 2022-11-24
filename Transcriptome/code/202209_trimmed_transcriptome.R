library("Rsubread")
library("Rsamtools")
library("pheatmap")
library("ggplot2")
library("edgeR")
library("limma")
library("biomaRt")
library("GOstats")
library("matrixStats")
####load featureCounte result####
counts<-read.csv("./Transcriptome/data/featureCounts/20220925_featureCounts_onTrimmed.csv",sep=',')

dim(counts)
counts[1:5,]
colnames(counts)<-sub(pattern = ".bam",replacement = "",x=sub(pattern = "..result.subread_alignment.",replacement = "",x=colnames(counts)))
counts[1:5,]
colnames(counts)<-c("Geneid","Chr","Start","End","Strand","Length",
                    "MPP4_3","MPP3_1","MPP1_3","MPP2_2","MPP3_3",
                    "HSC_1","HSC_2","MPP1_4","MPP4_1","MPP1_1",
                    "HSC_4","MPP2_3","MPP3_2","MPP2_1","HSC_3","MPP4_2")
####Add gene name to count table####
DEGgene_id <- as.character(counts$Geneid)
length(DEGgene_id)

mart<- useDataset("mmusculus_gene_ensembl",
                  useMart("ENSEMBL_MART_ENSEMBL"))
detags.IDs <- getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name","entrezgene_id","strand","transcript_length"),
  values=DEGgene_id,
  mart= mart)

dim(detags.IDs)
head(detags.IDs)
# remove duplicates
detags.IDs.matrix<-detags.IDs[-which(duplicated(detags.IDs$entrezgene_id)),]
dim(detags.IDs.matrix)
head(detags.IDs.matrix)
# select genes of interest only
rownames(detags.IDs.matrix)<-detags.IDs.matrix$entrezgene_id
rownames(counts)<-counts$Geneid
dim(counts)
# join the two tables
counts.kept<-counts[rownames(counts) %in% detags.IDs.matrix$entrezgene_id,]
dim(counts.kept)
counts.annot <- merge(detags.IDs.matrix,counts.kept,by='row.names', all=TRUE)
dim(counts.annot)
head(counts.annot)
# check the annotated table
head(counts.annot )
dim(counts.annot )
nrow(counts.annot )
write.csv(counts.annot ,file="./Transcriptome/result/20220926_counts_Annotation.csv")
colnames(counts.annot)
rownames(counts.annot)<-counts.annot$entrezgene_id
####calculate RPKM####
yRpkm<-rpkm(y =counts.annot[,14:29],
            gene.length = counts.annot$Length,
            log = F )
colnames(yRpkm)<-paste(colnames(yRpkm),"_RPKM",sep="")
counts.annot.RPKM <- merge(counts.annot,yRpkm,by=0, all=TRUE)
colnames(counts.annot.RPKM)
counts.annot.RPKM<-counts.annot.RPKM[,-c(1,2)]
rownames(counts.annot.RPKM)<-counts.annot.RPKM$entrezgene_id
write.csv(counts.annot.RPKM,file="./Transcriptome/result/20220926_counts_Annotation_RPKM.csv")

counts.annot.RPKM[counts.annot.RPKM$mgi_symbol %in% c("Rbfox1","Rbfox2","Rbfox3"),]


####load experiment design####
experiment_design <- read.table(file ="./Transcriptome/data/ExperimentDesign.txt",header =T )
# set the rownames to the sampleID to allow for ordering
rownames(experiment_design) <- experiment_design$SampleID
# order the design following the counts sample order
experiment_design.ord <- experiment_design[colnames(counts.annot.RPKM)[13:28],]
# look at the design
experiment_design.ord
# list the ordered samples for future use
samples <- as.character(experiment_design.ord$SampleID)
# create factors for future plotting
group<-factor(experiment_design.ord$condition)
group

####Basic QC plots#####
# density plot of raw read counts (log10)
pdf(file="./Transcriptome/result/20220926_Raw_read_counts_per_gene.density.pdf")
logcounts <- log(counts.annot.RPKM[,13],10) 
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.45),
     xlab="Raw read counts per gene (log10)",
     ylab="Density")
for (s in 14:28){
  logcounts <- log(counts.annot.RPKM[,s],10) 
  d <- density(logcounts)
  lines(d)
}
dev.off()

## box plots of raw read counts (log10)
pdf(file="./Transcriptome/result/20220926_Raw_read_counts_per_gene.boxplot.pdf")
logcounts <- log((counts.annot.RPKM[,29:44]+1),10)
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
axis(2)
axis(1,at=c(1:length(samples)),labels=colnames(logcounts),las=2,cex.axis=0.8)
dev.off()

# select data for the 100 most highly expressed genes
select <- order(rowMeans(counts.annot.RPKM[,29:44]), decreasing=TRUE)[1:1000]
highexprgenes_counts <-counts.annot.RPKM[select,29:44]
# heatmap with sample name on X-axis
pdf(file="./Transcriptome/result/20220926_High1000_expr_genes_heatmap.pdf")
heatmap(as.matrix(highexprgenes_counts), col=topo.colors(50), margin=c(10,6))
dev.off()

# select data for the 500 most highly variable expressed genes
logcounts <- cpm(counts.annot.RPKM[,13:28],log=TRUE)
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
#colnames(highly_variable_lcpm)<- group
# plot
pdf(file="./Transcriptome/result/20220926_High_var2000_genes.heatmap.group.pdf")
heatmap(highly_variable_lcpm, col = topo.colors(50), margin=c(10,6))
dev.off()
 
Vari500<-pheatmap(highly_variable_lcpm[order(highly_variable_lcpm[,1]),],
                  cale="row",show_rownames = F,cutree_rows = 3,
                  cellheight = 0.3,cellwidth = 20,
                  cluster_rows = T,cluster_cols = T,
                  border_color=NA,
                  color =colorRampPalette(c("blue","white","red"))(1000))
ggsave(plot = Vari500,filename = "./Transcriptome/result/20220926_Vari1000_Expr.pdf",
       device = "pdf",width = 8,height=40)

# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highly_variable_lcpm)
dim(data_for_PCA)
## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned
# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
pdf(file="./Transcriptome/result/20220926_High_var_genes.PCA_PropExplainedV.pdf")
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")
dev.off()
## calculate MDS
mds <- cmdscale(dist(data_for_PCA))
#Samples representation
pdf(file="./Transcriptome/result/20220926_High_var_genes.PCA.pdf")
plot(mds[,1], -mds[,2], type="p", xlab=paste("Dimension1_",eig_pc[1],sep=""), 
     ylab=paste("Dimension2_",eig_pc[2],sep=""), main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
dev.off()

####PCA############
# select data for the 100 most highly expressed genes
select <- order(rowMeans(counts.annot.RPKM[,13:28]), decreasing=TRUE)[1:100]
highexprgenes_counts <- counts.annot.RPKM[,13:28][select,]
# annotate the data with condition group as labels
#colnames(highexprgenes_counts)<- group
# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)
dim(data_for_PCA)
## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned
# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
pdf(file="./Transcriptome/result/20220926_High_expr_genesPCA_PropExplainedVariance.pdf")
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")
dev.off()
## calculate MDS
mds <- cmdscale(dist(data_for_PCA))
#Samples representation
pdf(file="./Transcriptome/result/20220926_High_expr_genesPCA.pdf")
plot(mds[,1], -mds[,2], type="p", xlab=paste("Dimension1_",eig_pc[1],sep=""), 
     ylab=paste("Dimension2_",eig_pc[2],sep=""), main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
dev.off()


####DESeq2####
cts<-counts.annot.RPKM[,13:28]
coldata<-experiment_design
i<-1
for (i in 1:nrow(coldata)){
  if (coldata$condition[i]!="HSC"){
    coldata$condition[i]<-"MPPs"
  }
}
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds$condition <- factor(dds$condition, levels = c("HSC","MPPs"))


###DEGs
dds <- DESeq(dds)
resultsNames(dds) 
res <- results(dds,alpha =0.05,cooksCutoff = FALSE,independentFiltering = F)
summary(res)
sum(res$padj <= 0.05& abs(res$log2FoldChange)>1, na.rm=TRUE)
write.csv(x = res,file = "./Transcriptome/result/20220926_DESeq2_HSCvsMPPs_All.csv")

########Gene Annotation#######
DEseq2All<-read.csv("./Transcriptome/result/20220926_DESeq2_HSCvsMPPs_All.csv",header = T)
rownames(DEseq2All)<-DEseq2All$X
DEseq2All<-res
DEseq2All_id<- as.character(rownames(DEseq2All))
length(DEseq2All_id)
DEseq2All["12462",]

detags.IDs <- getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name","entrezgene_id","strand","transcript_length"),
  values=DEseq2All_id,
  mart= mart)

dim(detags.IDs)
head(detags.IDs)
# remove duplicates
detags.IDs.matrix<-detags.IDs[-which(duplicated(detags.IDs$entrezgene_id)),]
dim(detags.IDs.matrix)
head(detags.IDs.matrix)
# select genes of interest only
rownames(detags.IDs.matrix)<-detags.IDs.matrix$entrezgene_id
#entrez_genes.annot <- detags.IDs.matrix[as.character(DEGentrez_genes),]
# join the two tables
DEseq2All.kept<-DEseq2All[rownames(DEseq2All) %in% detags.IDs.matrix$entrezgene_id,]
DEseq2All.annot<-merge(as.data.frame(detags.IDs.matrix), as.data.frame(DEseq2All.kept), by='row.names', all=TRUE)

# check the annotated table
head(DEseq2All.annot )
dim(DEseq2All.annot )
nrow(DEseq2All.annot )
write.csv(DEseq2All.annot ,file="./Transcriptome/result/20220926_DEseq2_MUTvsCtrl_All_Annotation.csv")

DEseq2All.annot.up<-DEseq2All.annot[DEseq2All.annot$padj < 0.05 &DEseq2All.annot$log2FoldChange>0,]
write.csv(DEseq2All.annot.up ,file="./Transcriptome/result/20220926_DEseq2_MUTvsCtrl_Up_Annotation.csv")
DEseq2All.annot.down<-DEseq2All.annot[DEseq2All.annot$padj < 0.05 &DEseq2All.annot$log2FoldChange<0,]
write.csv(DEseq2All.annot.down ,file="./Transcriptome/result/20220926_DEseq2_MUTvsCtrl_Down_Annotation.csv")

DEseq2All.annot.noNA<-DEseq2All.annot[rowSums(is.na(DEseq2All.annot)) == 0 ,]
DEseq2All.annot.noNA.sig<-DEseq2All.annot.noNA[DEseq2All.annot.noNA$padj < 0.05 & abs (DEseq2All.annot.noNA$log2FoldChange)>1,]
rownames(DEseq2All.annot.noNA.sig)<-DEseq2All.annot.noNA.sig$entrezgene_id

heatmapgenes<-counts.annot.RPKM[rownames(DEseq2All.annot.noNA.sig),29:44]
dim(heatmapgenes)

colnames(heatmapgenes)<-stringr::str_sub(string = colnames(heatmapgenes),end=-6)

p11<-pheatmap(heatmapgenes,
             cluster_cols = T,treeheight_col = 40,
             cluster_rows = T,treeheight_row = 20,
             border_color = NA,show_rownames = F,
             scale = "row",display_numbers = F,
             color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))

ColAnno<-data.frame("Cell"=stringr::str_sub(colnames(heatmapgenes),end=-3))
rownames(ColAnno)<-colnames(heatmapgenes)
ann_colors <- list(
Cell = c(HSC="#f8766d",MPP1="#7fcd00",MPP2="#ffa33b",
        MPP3="#00d5fa",MPP4="#3dc1ff"))

p11col_dend<-p11[[2]]

p11col_dend <-dendextend::rotate(p11col_dend, order =c("HSC_4" ,"HSC_1","HSC_2","HSC_3",
                                            "MPP1_1", "MPP1_4", "MPP1_3",
                                            "MPP2_3" ,"MPP2_2","MPP2_1",
                                            "MPP3_1","MPP3_2","MPP3_3",
                                            "MPP4_3","MPP4_1","MPP4_2"))


p22<-pheatmap(heatmapgenes, 
             cluster_cols=as.hclust(p11col_dend),cutree_rows = 2,
             treeheight_col = 40,cutree_cols = 2,
             cluster_rows = T,treeheight_row = 10,border_color = NA,show_rownames = F,
             scale = "row",display_numbers = F,
            annotation_col=ColAnno,
            annotation_colors = ann_colors,
             color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))

ggsave(filename = "./Transcriptome/result/20221001_heatmap_byDEGSs.pdf",
       plot=p22,device = "pdf",width = 10,height = 10)


###RBP volcano plot####
###import RBPDB###
RBPDB<-read.csv("RBPDB_v1.3.1_proteins_mouse_2012-11-21.csv",header = F)
head(RBPDB)
ENSid<-RBPDB$V2
###Convert EnsemblID to Gene names###
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
geneName<-getBM(attributes = c("ensembl_gene_id","mgi_symbol", "description"),filters = "ensembl_gene_id",
                values =ENSid,mart = mart)
head(geneName)
###Add gene names to original file###
RBPDB_with_MgiSymbol<-merge(x = RBPDB,y = geneName,by.x ="V2" ,by.y ="ensembl_gene_id")
head(RBPDB_with_MgiSymbol)
sort(RBPDB_with_MgiSymbol$mgi_symbol)
write.table(x = RBPDB_with_MgiSymbol,file = "RBPDB_with_MgiSymbol.csv",sep="\t")
#check how merge works
RBPDB_with_MgiSymbol[RBPDB_with_MgiSymbol$mgi_symbol %in% c("Rbfox1","Rbfox2","Rbfox3"),]

###load RPKM data generated from Trummph's 2014 CSC paper###
RBPDB_with_MgiSymbol<-read.csv("./Transcriptome/data/RBPDB_with_MgiSymbol.csv",sep='\t')
head(RBPDB_with_MgiSymbol)
colnames(RBPDB_with_MgiSymbol)
RPKM<-counts.annot.RPKM
head(RPKM)
dim(RPKM)

RBPDB_RPKM<-RPKM[RPKM$mgi_symbol %in% RBPDB_with_MgiSymbol$mgi_symbol,]
RBPDB_RPKM<-RBPDB_RPKM[rownames(RBPDB_RPKM)!="19652",]
dim(RBPDB_RPKM)
colnames(RBPDB_RPKM)
head(RBPDB_RPKM[,-c(1:28)])
rownames(RBPDB_RPKM)<-RBPDB_RPKM$mgi_symbol
pheatmap(RBPDB_RPKM[,c(29:44)],scale = "row",cluster_row = T,
         treeheight_row = 0,cluster_cols = F,show_colnames = T,
         show_rownames = F)

###load rowcounts data###
rownames(DEseq2All.annot)<-DEseq2All.annot$entrezgene_id
RBPDB_RPKM_DEGs<-DEseq2All.annot[DEseq2All.annot$mgi_symbol %in% RBPDB_with_MgiSymbol$mgi_symbol,]
dim(RBPDB_RPKM_DEGs)
RBPDB_RPKM_DEGs<-RBPDB_RPKM_DEGs[!duplicated(RBPDB_RPKM_DEGs$mgi_symbol),]
rownames(RBPDB_RPKM_DEGs)<-RBPDB_RPKM_DEGs$mgi_symbol
rownames(RBPDB_RPKM_DEGs)
RBPDB_RPKM_DEGs<-as.data.frame(RBPDB_RPKM_DEGs)
RBPDB_RPKM_DEGs["Rbfox2",]
# RBPDB_RPKM_DEGs$ensembl_gene_id<-rownames(RBPDB_RPKM_DEGs)
# RBPDB_RPKM_DEGs<-merge(RBPDB_RPKM_DEGs,RBPDB_RPKM,
#                        by.x ="ensembl_gene_id",by.y="ensembl_gene_id" )

# add a column of NAs
RBPDB_RPKM_DEGs$diffexpressed <- "No"
# set "UP" 
RBPDB_RPKM_DEGs$diffexpressed[RBPDB_RPKM_DEGs$log2FoldChange < (-1) & RBPDB_RPKM_DEGs$padj < 0.05] <- "Up"
# set "DOWN"
RBPDB_RPKM_DEGs$diffexpressed[RBPDB_RPKM_DEGs$log2FoldChange >1 & RBPDB_RPKM_DEGs$padj < 0.05] <- "Down"
sum(RBPDB_RPKM_DEGs$diffexpressed =="No")
sum(RBPDB_RPKM_DEGs$diffexpressed =="Up")
sum(RBPDB_RPKM_DEGs$diffexpressed =="Down")

volplot<-ggplot(data = as.data.frame(RBPDB_RPKM_DEGs),
                mapping = aes(x=-(log2FoldChange),y=-log10(padj),col=diffexpressed))+
  geom_point(alpha=12,size=2,pch=16)+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",panel.grid=element_blank())+
  scale_color_manual(values=c( "#9a9a9a", "#f45864","#276acc"),
                     labels=c(paste("No_",sum(RBPDB_RPKM_DEGs$diffexpressed =="No"),sep=""), 
                              paste("Up_",sum(RBPDB_RPKM_DEGs$diffexpressed =="Up"),sep=""),
                              paste("Down_",sum(RBPDB_RPKM_DEGs$diffexpressed =="Down"),sep="")))+
  ggrepel::geom_label_repel(aes(label=ifelse(abs(log2FoldChange)>1 & padj<0.05,
                                    as.character(mgi_symbol),'')))
ggsave(filename = "./Transcriptome/result/20221001_RBPDB_DEGsVolcano.pdf",plot = volplot,
       device = "pdf",width = 10,height = 10)
write.csv(x = RBPDB_RPKM_DEGs,file = "./Transcriptome/result/20221001_RBPDB_DEGs.csv")


####Extract Rbfox RPKMs####
RbfoxFamily<-counts.annot.RPKM[counts.annot.RPKM$mgi_symbol %in% c("Rbfox1", "Rbfox2","Rbfox3"),]
write.csv(x = RbfoxFamily,file = "./Transcriptome/result/20221001_RbfoxFamily_RPKM.csv")





####DESeq2 between HSC and MPP1####
cts1<-counts.annot.RPKM[,c(15,18,19,20,22,23,27)]
coldata1<-experiment_design[experiment_design$condition %in% c("MPP1","HSC"),]
j<-1
for (j in 1:nrow(coldata1)){
  if (coldata1$condition[j]!="HSC"){
    coldata1$condition[j]<-"MPP1"
  }
}
library("DESeq2")
dds1 <- DESeqDataSetFromMatrix(countData = cts1,
                              colData = coldata1,
                              design = ~ condition)
dds1$condition <- factor(dds1$condition, levels = c("HSC","MPP1"))


###DEGs
dds1 <- DESeq(dds1)
resultsNames(dds1) 
res1 <- results(dds1,alpha =0.05,cooksCutoff = FALSE,independentFiltering = F)
summary(res1)
sum(res1$padj <= 0.05& abs(res1$log2FoldChange)>1, na.rm=TRUE)
write.csv(x = res1,file = "./Transcriptome/result/20220926_DESeq2_HSCvsMPP1_All.csv")

########Gene Annotation#######
#DEseq2All<-read.csv("./Transcriptome/result/20220926_DESeq2_HSCvsMPP1_All.csv",header = T)
# DEseq2HSC_vs_MPP1<-res1
# rownames(DEseq2HSC_vs_MPP1)<-DEseq2HSC_vs_MPP1$X
DEseq2HSC_vs_MPP1<-res1
DEseq2HSC_vs_MPP1_id<- as.character(rownames(DEseq2HSC_vs_MPP1))
length(DEseq2HSC_vs_MPP1_id)
DEseq2HSC_vs_MPP1["12462",]

detags.IDs1 <- getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name","entrezgene_id","strand","transcript_length"),
  values=DEseq2HSC_vs_MPP1_id,
  mart= mart)

dim(detags.IDs1)
head(detags.IDs1)
# remove duplicates
detags.IDs.matrix1<-detags.IDs1[-which(duplicated(detags.IDs1$entrezgene_id)),]
dim(detags.IDs.matrix1)
head(detags.IDs.matrix1)
# select genes of interest only
rownames(detags.IDs.matrix1)<-detags.IDs.matrix1$entrezgene_id
#entrez_genes.annot <- detags.IDs.matrix[as.character(DEGentrez_genes),]
# join the two tables
DEseq2HSC_vs_MPP1.kept<-DEseq2HSC_vs_MPP1[rownames(DEseq2HSC_vs_MPP1) %in% detags.IDs.matrix1$entrezgene_id,]
DEseq2HSC_vs_MPP1.annot<-merge(as.data.frame(detags.IDs.matrix1), as.data.frame(DEseq2HSC_vs_MPP1.kept), by='row.names', all=TRUE)

# check the annotated table
head(DEseq2HSC_vs_MPP1.annot )
dim(DEseq2HSC_vs_MPP1.annot )
nrow(DEseq2HSC_vs_MPP1.annot )


DEseq2HSC_vs_MPP1.annot.noNA<-DEseq2HSC_vs_MPP1.annot[rowSums(is.na(DEseq2HSC_vs_MPP1.annot)) == 0 ,]
DEseq2HSC_vs_MPP1.annot.noNA.sig<-DEseq2HSC_vs_MPP1.annot.noNA[DEseq2HSC_vs_MPP1.annot.noNA$padj < 0.05 & abs (DEseq2HSC_vs_MPP1.annot.noNA$log2FoldChange)>1,]
rownames(DEseq2HSC_vs_MPP1.annot.noNA.sig)<-DEseq2HSC_vs_MPP1.annot.noNA.sig$entrezgene_id

heatmapgenes1<-counts.annot.RPKM[rownames(DEseq2HSC_vs_MPP1.annot.noNA.sig),c(29:44)]
dim(heatmapgenes1)

colnames(heatmapgenes1)<-stringr::str_sub(string = colnames(heatmapgenes1),end=-6)

p111<-pheatmap(heatmapgenes1,
              cluster_cols = T,treeheight_col = 40,
              cluster_rows = T,treeheight_row = 20,
              border_color = NA,show_rownames = F,
              scale = "row",display_numbers = F,
              color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))

ColAnno1<-data.frame("Cell"=stringr::str_sub(colnames(heatmapgenes1),end=-3))
rownames(ColAnno1)<-colnames(heatmapgenes1)
ann_colors1 <- list(
  Cell = c(HSC="#f8766d",MPP1="#7fcd00",MPP2="#ffa33b",
           MPP3="#00d5fa",MPP4="#3dc1ff"))

p111col_dend<-p111[[2]]

p111col_dend <-dendextend::rotate(p111col_dend, order =c("HSC_4" ,"HSC_1","HSC_2","HSC_3",
                                                       "MPP1_1", "MPP1_4", "MPP1_3",
                                                       "MPP2_3" ,"MPP2_2","MPP2_1",
                                                       "MPP3_1","MPP3_2","MPP3_3",
                                                       "MPP4_3","MPP4_1","MPP4_2"))


p122<-pheatmap(heatmapgenes1, 
              cluster_cols=as.hclust(p111col_dend),cutree_rows = 2,
              treeheight_col = 40,cutree_cols = 2,
              cluster_rows = T,treeheight_row = 10,border_color = NA,show_rownames = F,
              scale = "row",display_numbers = F,
              annotation_col=ColAnno1,
              annotation_colors = ann_colors1,
              color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))

ggsave(filename = "./Transcriptome/result/20221018_heatmap_byDEGSs_HSCvsMPP1.pdf",
       plot=p122,device = "pdf",width = 10,height = 10)


####20221019  top 300 varible genes between HSC and MPP1####
colnames(counts.annot.RPKM)
rownames(counts.annot.RPKM)
counts.annot.RPKM$HSCmean<-rowMeans(as.matrix(counts.annot.RPKM[,c(34,35,39,43)]))
counts.annot.RPKM$MPP1mean<-rowMeans(as.matrix(counts.annot.RPKM[,c(31,36,38)]))
counts.annot.RPKM$Allmean<-rowMeans(as.matrix(counts.annot.RPKM[,c(29:44)]))
#counts.annot.RPKM.tophalf<-counts.annot.RPKM[rowMeans(as.matrix(counts.annot.RPKM[,c(29:44)]))>median(counts.annot.RPKM$Allmean),]
HSCvsMPP1_Sds<-rowSds(as.matrix(counts.annot.RPKM[,c("HSCmean","MPP1mean")]))

#HSCvsMPP1_Sds<-rowSds(as.matrix(counts.annot.RPKM[,c(34,35,39,43,31,36,38)]))
select <- order(HSCvsMPP1_Sds, 
                decreasing=TRUE)[1:200]
heatmapgenes3 <-counts.annot.RPKM[select,c(29:44)]
dim(heatmapgenes3)

p311<-pheatmap(heatmapgenes3,
               cluster_cols = T,treeheight_col = 40,
               cluster_rows = T,treeheight_row = 20,
               border_color = NA,show_rownames = F,
               scale = "row",display_numbers = F,
               color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))

ColAnno3<-data.frame("Cell"=stringr::str_sub(colnames(heatmapgenes3),end=-8))
rownames(ColAnno3)<-colnames(heatmapgenes3)
ann_colors3 <- list(
  Cell = c(HSC="#f8766d",MPP1="#7fcd00",MPP2="#ffa33b",
           MPP3="#00d5fa",MPP4="#3dc1ff"))

p311col_dend<-p311[[2]]

p311col_dend <-dendextend::rotate(p311col_dend, order =paste(c("HSC_4" ,"HSC_1","HSC_2","HSC_3",
                                                         "MPP1_1", "MPP1_4", "MPP1_3",
                                                         "MPP2_3" ,"MPP2_2","MPP2_1",
                                                         "MPP3_1","MPP3_2","MPP3_3",
                                                         "MPP4_3","MPP4_1","MPP4_2"),"_RPKM",sep=""))


p322<-pheatmap(heatmapgenes3, 
               cluster_cols=as.hclust(p311col_dend),cutree_rows = 2,
               treeheight_col = 40,cutree_cols = 2,
               cluster_rows = T,treeheight_row = 10,border_color = NA,show_rownames = F,
               scale = "row",display_numbers = F,
               annotation_col=ColAnno3,
               annotation_colors = ann_colors3,
               color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))

ggsave(filename = "./Transcriptome/result/20221019_heatmap_byTop300Varible_HSCvsMPP1.pdf",
       plot=p322,device = "pdf",width = 10,height = 10)

