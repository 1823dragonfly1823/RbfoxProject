library("Rsubread")
library("Rsamtools")
library("pheatmap")
library("ggplot2")
library("edgeR")
library("limma")
library("biomaRt")
library("GOstats")
library("RColorBrewer")
library("Mus.musculus") 
library("DESeq2")

######
load("20220321_to_featureCounts.RData")
experiment_design <- read.table(file ="ExperimentDesign.txt",header =T )
# set the rownames to the sampleID to allow for ordering
rownames(experiment_design) <- experiment_design$SampleID
# order the design following the counts sample order
experiment_design.ord <- experiment_design[colnames(Rawcounts$counts),]
# look at the design
experiment_design.ord
# list the ordered samples for future use
samples <- as.character(experiment_design.ord$SampleID)
# create factors for future plotting
group<-factor(experiment_design.ord$condition)
group

dim(Rawcounts$counts)

geneid <- rownames(Rawcounts$counts) 

mart<- useDataset("mmusculus_gene_ensembl",
                  useMart("ENSEMBL_MART_ENSEMBL"))
genes <- getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name","entrezgene_id","strand","transcript_length"),
  values=geneid,
  mart= mart)
dim(genes)
head(genes)
# remove duplicates
genes.matrix<-genes[-which(duplicated(genes$entrezgene_id)),]
dim(genes.matrix)
head(genes.matrix)
# select genes of interest only
rownames(genes.matrix)<-genes.matrix$entrezgene_id
head(rownames(genes.matrix))
# join the two tables
Rawcounts_counts.kept<-Rawcounts$counts[rownames(Rawcounts$counts) %in% genes.matrix$entrezgene_id,]
Rawcounts_counts.kept.annot <- merge(genes.matrix,Rawcounts_counts.kept,by='row.names', all=TRUE)
# check the annotated table
head(Rawcounts_counts.kept.annot )
dim(Rawcounts_counts.kept.annot )
nrow(Rawcounts_counts.kept.annot )
write.csv(Rawcounts_counts.kept.annot ,file="20220408_Rawcount_Allgene_Annotation.csv")

####CalculateRPKM####
yRpkm<-rpkm(y =Rawcounts_counts.kept.annot[,c(8:13)], 
            normalized.lib.sizes=T,
            gene.length = Rawcounts_counts.kept.annot$transcript_length,log = F )
colnames(yRpkm)<-c("Ctrl1_RPKM","Ctrl2_RPKM","Ctrl3_RPKM","Mut1_RPKM","Mut2_RPKM","Mut3_RPKM")

apply(yRpkm,2,median)
yRpkm.annot <- cbind(Rawcounts_counts.kept.annot,yRpkm)
rownames(yRpkm.annot)<-yRpkm.annot$entrezgene_id
write.csv(yRpkm.annot[,-1],file="20220408_RPKM_Allgene_Annotation.csv")

dim(yRpkm.annot)

ctrlmeanRpkm<-apply(yRpkm.annot[,14:16],1,mean)
mutmeanRpkm<-apply(yRpkm.annot[,17:19],1,mean)

Rawcount.filter<-Rawcounts_counts.kept.annot[ctrlmeanRpkm>median(ctrlmeanRpkm) | mutmeanRpkm>median(mutmeanRpkm),]
rownames(Rawcount.filter)<-Rawcount.filter$entrezgene_id
write.csv(Rawcount.filter[,-1],file="20220408_Rawcount_Filteredgene_Annotation.csv")
dim(Rawcount.filter)

####edgeR-limma-voom####
myfilteredcounts <- Rawcount.filter[,c(8:13)]
y <- DGEList(myfilteredcounts,lib.size = colSums(myfilteredcounts),
             norm.factors = calcNormFactors(myfilteredcounts),
             samples = experiment_design$SampleID,
             group = experiment_design$condition)
dim(y$counts)
# check your samples grouping
experiment_design.ord[colnames(y$counts),]$condition == group
# create design matrix for limma
design <- model.matrix(~0+group)
# substitute "group" from the design column names
#colnames(design)<- gsub("group","",colnames(design))
colnames(design)<-levels(group)
# check your design matrix
design

# calculate normalization factors between libraries
y <- calcNormFactors(y,method = "TMM")

#Unsupervised clustering of samples
lcpm <- cpm(y, log=TRUE)
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="Sample groups")

# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix <- makeContrasts(Ctrl-Mut,levels=design)
cont.matrix 

# normalize the read counts with 'voom' function
v <- voom(y,design,plot=T)
# extract the normalized read counts
counts.voom <- v$E

# fit linear model for each gene given a series of libraries
fit <- lmFit(v,design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit <- contrasts.fit(fit, cont.matrix)
# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit <- eBayes(fit)
plotSA(fit)
options(digits=3)
# check the output fit
dim(fit)
summary(decideTests(fit))
# set adjusted pvalue threshold and log fold change threshold
mypval=0.05
myfc=1
# get the coefficient name for the comparison  of interest
colnames(fit$coefficients)
mycoef<-colnames(fit$coefficients)
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit,coef=mycoef)
# get the full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])
limma.res.annotation<-merge(yRpkm.annot[rownames(limma.res),-c(1)],limma.res,by='row.names', all=TRUE)
rownames(limma.res.annotation)<-limma.res.annotation$entrezgene_id
write.csv(limma.res.annotation,file="20220408_DEG_Allgenes_edgeR.csv")

# get significant DE genes only (adjusted p-value < mypval)
limma.res.annotation.p0.05 <- limma.res.annotation[limma.res.annotation$adj.P.Val<0.05,]
dim(limma.res.annotation.p0.05)
# write limma output table for significant genes into a tab delimited file
write.csv(limma.res.annotation.p0.05,file="20220408_DEG_pvalue_0.05_edgeR.csv")

####GenerateGSEAtablefromedgeR####
x<-DGEList(myfilteredcounts,lib.size = colSums(myfilteredcounts),
           norm.factors = calcNormFactors(myfilteredcounts),
           samples = experiment_design$SampleID,
           group = experiment_design$condition)
x <- calcNormFactors(x)
GSEA_table <- cpm(x, log=TRUE)
colnames(GSEA_table)<-c("WT1","WT2","WT3","MUT1","MUT2","MUT3")
dim(genes)
geneid <- rownames(GSEA_table) 
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL"), 
                keytype="ENTREZID")
rownames(genes)<-genes$ENTREZID
GSEA_table<-merge(GSEA_table,genes,by='row.names',all=T)
GSEA_table<-GSEA_table[!is.na(GSEA_table$SYMBOL),]
rownames(GSEA_table)<-GSEA_table$ENTREZID
write.csv(GSEA_table,"20220408_GSEA_table_limma.csv")

####Gene Set Enrichment from edgeR##########
# Define list of genes of interest (DE genes - EntrezGene IDs)
DEGEntrezIDs <- as.character(limma.res.annotation.p0.05$entrezgene_id)
# Define the universe
AllIDs <- getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name","entrezgene_id",'strand'),
  values= rownames(Rawcount.filter),
  mart= mart)
universeids <- as.character(unique(AllIDs$entrezgene_id))
length(universeids)
# define the p-value cut off for the hypergeometric test
hgCutoff <- 0.05
params <- new("GOHyperGParams",annotation="org.Mm.eg.db",
              geneIds=DEGEntrezIDs,universeGeneIds=universeids,
              ontology="BP",
              pvalueCutoff=hgCutoff,testDirection="over")
#  Run the test
mm <- hyperGTest(params)
# Check results
mm
## Get the p-values of the test
mm.pv <- pvalues(mm)
## Adjust p-values for multiple test (FDR)
mm.pv.fdr <- p.adjust(mm.pv,'fdr')
## select the GO terms with adjusted p-value less than the cut off
sigGO.ID <- names(mm.pv.fdr[mm.pv.fdr < hgCutoff])
length(sigGO.ID)
# get table from HyperG test result
df <- summary(mm)
# keep only significant GO terms in the table
GOannot.table <- df[df[,1] %in% sigGO.ID,]
head(GOannot.table)

# Create text report of the significantly over-represented GO terms
write.csv(GOannot.table,file="20220408_GOterms_OverRep_BP_limmavoom.csv")
# Create html report of all over-represented GO terms
htmlReport(mm, file="20220408_GOterms_OverRep_BP_limmavoom.html")

##################################################################
####DESeq2######
cts<-Rawcount.filter[,c(8:13)]
coldata<-experiment_design
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds$condition <- factor(dds$condition, levels = c("Ctrl","Mut"))
###DEGs
dds <- DESeq(dds)
resultsNames(dds) 
res <- results(dds,alpha =0.05,cooksCutoff = FALSE,independentFiltering = F)
summary(res)
sum(res$padj >= 0.05& res$log2FoldChange<0, na.rm=TRUE)
res.annonatation<-yRpkm.annot[rownames(res),]
DEG.annotation<-merge(as.data.frame(res),res.annonatation,by="row.names",all=T)
write.csv(DEG.annotation,file="20220408_DEG_Allgenes_DEseq2.csv")

DEG.annotation.pvalue_0.05<-DEG.annotation[DEG.annotation$padj < 0.05 ,]
write.csv(DEG.annotation.pvalue_0.05 ,file="220220408_DEG_pvalue_0.05_DEseq2.csv")
DEG.annotation.up<-DEG.annotation[DEG.annotation$padj < 0.05 & DEG.annotation$log2FoldChange>0,]
write.csv(DEG.annotation.up ,file="20220408_DEG_Upgenes_DEseq2.csv")
DEG.annotation.down<-DEG.annotation[DEG.annotation$padj < 0.05 & DEG.annotation$log2FoldChange<0,]
write.csv(DEG.annotation.down ,file="20220408_DEG_Downgenes_DEseq2.csv")

####Gene Set Enrichment##########
# Define list of genes of interest (DE genes - EntrezGene IDs)
resSigUpEntrezIDs <- as.character(DEG.annotation.pvalue_0.05$entrezgene_id)
# Define the universe
AllIDs <- getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name","entrezgene_id",'strand'),
  values= rownames(Rawcount.filter),
  mart= mart)
universeids <- as.character(unique(AllIDs$entrezgene_id))
length(universeids)

# define the p-value cut off for the hypergeometric test
hgCutoff <- 0.05
params <- new("GOHyperGParams",annotation="org.Mm.eg.db",
              geneIds=resSigUpEntrezIDs,universeGeneIds=universeids,
              ontology="BP",
              pvalueCutoff=hgCutoff,testDirection="over")
#  Run the test
mm <- hyperGTest(params)
# Check results
mm
## Get the p-values of the test
mm.pv <- pvalues(mm)
## Adjust p-values for multiple test (FDR)
mm.pv.fdr <- p.adjust(mm.pv,'fdr')
## select the GO terms with adjusted p-value less than the cut off
sigGO.ID <- names(mm.pv.fdr[mm.pv.fdr < hgCutoff])
length(sigGO.ID)

# get table from HyperG test result
df <- summary(mm)
# keep only significant GO terms in the table
GOannot.table <- df[df[,1] %in% sigGO.ID,]
head(GOannot.table)

# Create text report of the significantly over-represented GO terms
write.csv(GOannot.table,file="20220408_GOterms_OverRep_BP_DEseq2.csv")
# Create html report of all over-represented GO terms
htmlReport(mm, file="20220408_GOterms_OverRep_BP_DEseq2.html")

####generateGSEA table from DESEQ2####
DEseq2GCT<-DEG.annotation
head(DEseq2GCT)
DEseq2GCT$fcsign <- sign(DEseq2GCT$log2FoldChange)
DEseq2GCT$logP<--log10(DEseq2GCT$pvalue)
DEseq2GCT$metric<- DEseq2GCT$logP/DEseq2GCT$fcsign
DEseq2GCT<-DEseq2GCT[,c("mgi_symbol", "metric")]
DEseq2GCT.filtered <- na.omit(DEseq2GCT)
head(DEseq2GCT.filtered)
write.csv(DEseq2GCT.filtered,file="20220408_GSEA_table_DEseq2.csv",quote=F)
