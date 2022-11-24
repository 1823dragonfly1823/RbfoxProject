library("Seurat")
library("biomaRt")
library("ggplot2")
library("dplyr")

#setwd("./Downloads/GSE81682_HTSeq_counts/")

####Import data####
GSE81682_HTSeq_counts<-read.delim("./data/GSE81682_HTSeq_counts.txt")
dim(GSE81682_HTSeq_counts)
GSE81682_HTSeq_counts[1:4,1:4]

table(sub(x = colnames(GSE81682_HTSeq_counts),pattern = "_.*",replacement = ""))
# ####Convert Enzembl ID in the original data into Gene symbol####
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
geneName<-getBM(attributes = c("ensembl_gene_id","mgi_symbol", "description"),
                filters = "ensembl_gene_id",
                values =GSE81682_HTSeq_counts$ID,
                mart = mart)
dim(geneName)

#####Add gene names to original file####
GSE81682_HTSeq_counts_with_MgiSymbol<-merge(x = geneName,y= GSE81682_HTSeq_counts,by.x ="ensembl_gene_id",by.y="ID")
head(GSE81682_HTSeq_counts_with_MgiSymbol)
GSE81682_HTSeq_counts_with_MgiSymbol[1:4,1:5]
GSE81682_HTSeq_counts_with_MgiSymbol<-GSE81682_HTSeq_counts_with_MgiSymbol[!duplicated(GSE81682_HTSeq_counts_with_MgiSymbol$mgi_symbol),]
rownames(GSE81682_HTSeq_counts_with_MgiSymbol)<-GSE81682_HTSeq_counts_with_MgiSymbol$mgi_symbol
GSE81682<-GSE81682_HTSeq_counts_with_MgiSymbol[,-c(1,3)]
GSE81682[1:5,1:10]
colnames(GSE81682)<-sub(pattern = "\\_","_",colnames(GSE81682))
rownames(GSE81682)<-sub(pattern = "\\_","_",rownames(GSE81682))
rownames(GSE81682)<-sub(pattern = "\\.","_",rownames(GSE81682))
rownames(GSE81682)<-sub(pattern = "\\%","_",rownames(GSE81682))
rownames(GSE81682)<-sub(pattern = "\\(","_",rownames(GSE81682))
rownames(GSE81682)<-sub(pattern = "\\)","_",rownames(GSE81682))
GSE81682<-GSE81682[rownames(GSE81682)!="",]
write.csv(x =GSE81682,file = "./result/20220924_GSE81682_for_seurat.csv" )

####create Seurat object####
GSE81682_seurat<-CreateSeuratObject(counts = GSE81682[,-c(1)])
GSE81682_seurat@meta.data$orig.ident<-stringr::str_sub(string = rownames(GSE81682_seurat@meta.data),
                                                       end = -5)
View(GSE81682_seurat@meta.data)
# Add number of UMIs per gene for each cell to metadata
GSE81682_seurat$log10GenesPerUMI <- log10(GSE81682_seurat$nFeature_RNA)/log10(GSE81682_seurat$nCount_RNA)
max(GSE81682_seurat$log10GenesPerUMI[GSE81682_seurat$log10GenesPerUMI!="NaN"])

##Extract mito genes from features list
mitogenes<-grep('^mt-',GSE81682_seurat@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)

# Compute percent mito ratio
GSE81682_seurat$mitoRatio <- PercentageFeatureSet(object = GSE81682_seurat, pattern = "^mt-")
GSE81682_seurat$mitoRatio <- GSE81682_seurat@meta.data$mitoRatio/100
max(GSE81682_seurat$mitoRatio[GSE81682_seurat$mitoRatio!="NaN"])

#compute ribosomal gene ratio
Ribogenes<-grep('^Rp[sl][[:digit:]]',GSE81682_seurat@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)
GSE81682_seurat$RiboRatio <- PercentageFeatureSet(object = GSE81682_seurat, pattern = "^Rp[sl][[:digit:]]")
GSE81682_seurat$RiboRatio <- GSE81682_seurat@meta.data$RiboRatio/100

# Rename columns
GSE81682_seurat@meta.data<-dplyr::rename (GSE81682_seurat@meta.data,seq_folder=orig.ident,
                          nUMI = nCount_RNA,
                          nGene = nFeature_RNA)
GSE81682_seurat@meta.data$sample<-GSE81682_seurat@meta.data$seq_folder
View(GSE81682_seurat@meta.data)

nUMIplot<-ggplot(GSE81682_seurat@meta.data[GSE81682_seurat@meta.data$nUMI>100,],
                 aes(sample,nUMI))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=500,color="red",size=1)+
  annotate(geom="text", x=0.5, y=500, label="500",
           size=6,color="blue")

nGeneplot<-ggplot(GSE81682_seurat@meta.data[GSE81682_seurat@meta.data$nUMI>100,],
                  aes(sample,nGene))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  scale_y_continuous()+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=250,color="red",size=1)+
  annotate(geom="text", x=0.5, y=250, label="250",
           size=6,color="blue")

nlog10GenesPerUMIplot<-ggplot(GSE81682_seurat@meta.data[GSE81682_seurat@meta.data$nUMI>100,],
                              aes(sample,log10GenesPerUMI))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.8,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.8, label="0.8",
           size=6,color="blue")

mitoRatioplot<-ggplot(GSE81682_seurat@meta.data[GSE81682_seurat@meta.data$nUMI>100,],
                      aes(sample,mitoRatio))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.2,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.2, label="0.2",
           size=6,color="blue")

RiboRatioplot<-ggplot(GSE81682_seurat@meta.data[GSE81682_seurat@meta.data$nUMI>100,],
                      aes(sample,RiboRatio))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.1,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.1, label="0.1",
           size=6,color="blue")

QCplot<-(nUMIplot+nGeneplot+nlog10GenesPerUMIplot+mitoRatioplot+RiboRatioplot)
ggsave("./result/20220924_QCplot.pdf",QCplot,"pdf",width=8,height=6)

# Filter#######
##Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat<-subset(x = GSE81682_seurat, 
                           subset= ((nUMI >= 200000) & 
                                      (nGene >= 4000) &
                                      (log10GenesPerUMI > 0.6) &
                                      (mitoRatio < 0.1) &
                                      (RiboRatio < 0.1)))

table(filtered_seurat@meta.data$sample)
filtered_combinednumber<-data.frame(filtered_seurat@meta.data$sample)

#SCTransform
options(future.globals.maxSize = 4000 * 1024^2)
##perform the cell cycle scoring and sctransform on all samples
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
filtered_seurat <- CellCycleScoring(filtered_seurat, g2m.features=stringr::str_to_title((tolower(cc.genes$g2m.genes))), 
                                        s.features=stringr::str_to_title((tolower(cc.genes$s.genes))),set.ident=T)
filtered_seurat <- SCTransform(filtered_seurat, 
                                   vars.to.regress = c("mitoRatio",
                                                       "S.Score", 
                                                       "G2M.Score",
                                                       "RiboRatio"),
                                   variable.features.n=3000)

seurat_integrated<-filtered_seurat

###UMAP visualization
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated,
                            ndims.print = 1:50,nfeatures.print = 5)

# Plot PCA
PCAPlot(seurat_integrated) 
ElbowPlot(seurat_integrated,ndims = 50)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:15,
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,reduction = "umap",group.by = "sample")


# Run TSNE
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:15,
                             reduction = "pca")
seurat_integrated <- RunTSNE(seurat_integrated,
                             dims = 1:15,tsne.method = "FIt-SNE",
                             reduction = "pca",perplexity=50)
# Plot TSNE                             
DimPlot(seurat_integrated,reduction = "tsne",
        group.by = "sample",pt.size = 1)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:15)
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.8,1))
View(seurat_integrated@meta.data)  

# Assign identity of clusters
Idents(object = seurat_integrated) <- "SCT_snn_res.1"

# Plot the UMAP
DimPlot(seurat_integrated,
                 reduction = "umap",
                 label = T,repel = T,
                 label.size = 5,
                 pt.size = 0.5)

DimPlot(seurat_integrated,
                 reduction = "tsne",
                 label = T,repel = T,
                 label.size = 5,
                 pt.size = 1)

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "seq_folder")) %>%
  dplyr::count(ident, seq_folder) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)
table(seurat_integrated@meta.data$seq_folder)

# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio","RiboRatio")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

DefaultAssay(seurat_integrated) <- "RNA"
Vlnplot1<-VlnPlot(object = seurat_integrated, 
                 features = c("Kitl","Ly6a","Procr","Trpc6", "Cd34","Flt3", "Slamf1","Cd48",
                              "Hoxb5","Fgd5",
                              "Mpl","Hba-a1","Smim1",
                              "Cd3e","Itgam","Cd19",
                     "Ly6g","Ly76",
                    "Fcgr2b","Fcgr3","Cd34"
                    ),
                 log=T,ncol = 5)
ggsave(filename = "./result/20221001_VlnPlot_1.pdf",plot = Vlnplot1,
       device = "pdf",
       width=20,height = 20,limitsize = F)

VlnPlot(object = seurat_integrated, 
        features = c("Cd3","Cd11b","Cd19",
                     "Ly6g","Ly76","Kitl",
                     "Sca1","Flt3",
                     "Slamf1","Cd34", "Cd48"
                     ),
        log=T,ncol = 2)
DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "SCT_snn_res.1"
HSCmarker<-FeaturePlot(object = seurat_integrated,pt.size = 1.5,
            cols = c("grey","#f8766d"),min.cutoff = 'q0',max.cutoff = 'q90',
            features =c("Slamf1","Cd48","Procr","Mpl", "Hoxb5","Fgd5","Mecom",
                        "Cd34","Flt3"),ncol = 3,label = T,label.size = 6)
ggsave(filename ="./result/20221001_featureplot_HSCmarker.pdf",plot = HSCmarker,device = "pdf", 
        width =25,height = 20)
HSCmarkertsne<-FeaturePlot(object = seurat_integrated,pt.size = 1.5,reduction = "tsne",
                       cols = c("grey","#f8766d"),min.cutoff = 'q0',max.cutoff = 'q90',
                       features =c("Slamf1","Cd48","Procr","Mpl", "Hoxb5","Fgd5","Mecom",
                                   "Cd34","Flt3"),ncol = 3,label = T,label.size = 6)
ggsave(filename ="./result/20221001_featureplot_HSCmarker_tSNE.pdf",plot = HSCmarkertsne,device = "pdf", 
       width =25,height = 20)
HSCmarkerVln<-VlnPlot(object = seurat_integrated, 
        features = c("Slamf1","Cd48","Procr","Mpl", "Hoxb5","Fgd5","Mecom",
                     "Cd34","Flt3"),pt.size = 1,
        log=T,ncol = 3
        # ,cols=c("#f8766d","#cd9600", "#00bfc4",
        #                       "#7cae00","#00a9ff","#ff61cc")
        )
ggsave(filename ="./result/20221001_Vlnplot_HSCmarker.pdf",plot = HSCmarkerVln,device = "pdf", 
       width =25,height = 20)

DefaultAssay(seurat_integrated) <- "RNA"
j<-0
for (j in (0:5)){
  conserved_markers <- FindMarkers(seurat_integrated, 
                                            ident.1 = j,
                                            verbose = T,only.pos=T) 
  write.csv(x=conserved_markers,file = paste("./result/cluster---",j,".csv",sep=''))
}


####extract cell ID####
metaData<-seurat_integrated@meta.data
All_CellID<-rownames(metaData)
HSC_CellID<-rownames(metaData[metaData$SCT_snn_res.0.3=='0',])
LMPP_CellID<-rownames(metaData[metaData$SCT_snn_res.0.3=='1',])
CMP_CellID<-rownames(metaData[metaData$SCT_snn_res.0.3=='2',])
GMP_CellID<-rownames(metaData[metaData$SCT_snn_res.0.3=='3',])
MEP_CellID<-rownames(metaData[metaData$SCT_snn_res.0.3=='4',])
MPP_CellID<-rownames(metaData[metaData$SCT_snn_res.0.3=='5',])

length(HSC_CellID)
table(metaData$SCT_snn_res.1)

write.csv(x = All_CellID,file = "./result/20221001_CellID_postFilter.csv")

####Load all samples from GEO###
AllSample_GEO<-read.csv("./data/sample.csv")
glimpse(AllSample_GEO)
# in sample.csv file
# rows 1-852, HSPC
# rows 853-1068, LT-HSC
# rows 1069-1920, Pro
# rows 1921-2772 HSPC_Rep
# rows 2773-2988, LT-HSC_Rep
# rows 2989-3840, Pro_Rep
CellID<-"NA"
CellID[1:852]<-sub(pattern =".*HSPC",replacement = "HSPC",x=AllSample_GEO[1:852,]$Title)
CellID[853:1068]<-sub(pattern =".*LT-HSC",replacement = "LT.HSC",x=AllSample_GEO[853:1068,]$Title)
CellID[1069:1920]<-sub(pattern =".*Prog",replacement = "Prog",x=AllSample_GEO[1069:1920,]$Title)
CellID[1921:2772]<-sub(pattern =".*HSPC_Rep",replacement = "HSPC",x=AllSample_GEO[1921:2772,]$Title)
CellID[2773:2988]<-sub(pattern =".*LT-HSC_Rep",replacement = "LT.HSC",x=AllSample_GEO[2773:2988,]$Title)
CellID[2989:3840]<-sub(pattern =".*Prog_Rep",replacement = "Prog",x=AllSample_GEO[2989:3840,]$Title)

tail(CellID)
length(CellID)

AllSample_GEO$CellID<-CellID
##Filter cells####
Samples_postFilter<-AllSample_GEO[AllSample_GEO$CellID %in% All_CellID,]
dim(Samples_postFilter)
colnames(Samples_postFilter)
HSPC_Samples_postFilter<-Samples_postFilter[c(1:696,1617:2312),]
LT.HSC_Samples_postFilter<-Samples_postFilter[c(697:864,2313:2480),]
Prog_Samples_postFilter<-Samples_postFilter[c(865:1616,2481:3232),]

Prog_Samples_postFilter$CellID[c(611:625,1363:1377)][14:30]
write.table(x =(Prog_Samples_postFilter[c(611:625,1363:1377),-(13)][14:30,]),quote = F,row.names = F,file = "20220802Prog_Samples_postFilter.txt",sep=',' )

write.table(x =(Samples_postFilter[,-(13)]),quote = F,row.names = F,file = "Samples_postFilter.txt",sep=',' )
write.table(x =(HSPC_Samples_postFilter[,-(13)]),quote = F,row.names = F,file = "HSPC_Samples_postFilter.txt",sep=',' )
write.table(x =(LT.HSC_Samples_postFilter[,-(13)]),quote = F,row.names = F,file = "LT.HSC_Samples_postFilter.txt",sep=',' )
write.table(x =(Prog_Samples_postFilter[,-(13)]),quote = F,row.names = F,file = "Prog_Samples_postFilter.txt",sep=',' )
write.table(x =(LT.HSC_Samples_postFilter[,-(13)][c(107,275),]),quote = F,row.names = F,file = "LT.HSC_Samples_postFilter_137.txt",sep=',' )

##create matrix for HSPC cells####
HSPC_countit_out_folder<-"./data/countit_out/HSPC/"
HSPC_filelist<-list.files(HSPC_countit_out_folder)
i<-1
cell.alt3<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/alt3.count.txt",sep=''),sep='\t')
cell.alt5<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/alt5.count.txt",sep=''),sep='\t')
cell.cass<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/cass.count.txt",sep=''),sep='\t')
cell.iret<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/iret.count.txt",sep=''),sep='\t')
cell.mutx<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/mutx.count.txt",sep=''),sep='\t')
cell.taca<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/taca.count.txt",sep=''),sep='\t')
exonTags<-data.frame("ID"=c(cell.alt3$name,cell.alt5$name,cell.cass$name,
                            cell.iret$name,cell.mutx$name,cell.taca$name))
exonTags[length(cell.alt3$name)+1,]
i<-1
for (i in (1:696)){
  cell.alt3<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/alt3.count.txt",sep=''),sep='\t')
  cell.alt5<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/alt5.count.txt",sep=''),sep='\t')
  cell.cass<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/cass.count.txt",sep=''),sep='\t')
  cell.iret<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/iret.count.txt",sep=''),sep='\t')
  cell.mutx<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/mutx.count.txt",sep=''),sep='\t')
  cell.taca<-read.csv(paste(HSPC_countit_out_folder,HSPC_filelist[i],"/taca.count.txt",sep=''),sep='\t')
  
  exonTags$exon<-c(cell.alt3$exonTags,cell.alt5$exonTags,cell.cass$exonTags,
                   cell.iret$junctionTags,cell.mutx$X5.ExonTags,
                   cell.taca$exonTags)
  colnames(exonTags)[ncol(exonTags)]<-noquote(sub(pattern="-",replacement = ".",
                                                  x=sub(pattern = "_countit_out",replacement = "",x=HSPC_filelist[i])))
}

##create matrix for LT.HSC cells####
LT.HSC_countit_out_folder<-"./data/countit_out/LT.HSC/"
LT.HSC_filelist<-list.files(LT.HSC_countit_out_folder)
type<-c("alt3","alt5","cass","iret","mutx","taca")
i<-1
for (i in (1:168)){
  cell.alt3<-read.csv(paste(LT.HSC_countit_out_folder,LT.HSC_filelist[i],"/alt3.count.txt",sep=''),sep='\t')
  cell.alt5<-read.csv(paste(LT.HSC_countit_out_folder,LT.HSC_filelist[i],"/alt5.count.txt",sep=''),sep='\t')
  cell.cass<-read.csv(paste(LT.HSC_countit_out_folder,LT.HSC_filelist[i],"/cass.count.txt",sep=''),sep='\t')
  cell.iret<-read.csv(paste(LT.HSC_countit_out_folder,LT.HSC_filelist[i],"/iret.count.txt",sep=''),sep='\t')
  cell.mutx<-read.csv(paste(LT.HSC_countit_out_folder,LT.HSC_filelist[i],"/mutx.count.txt",sep=''),sep='\t')
  cell.taca<-read.csv(paste(LT.HSC_countit_out_folder,LT.HSC_filelist[i],"/taca.count.txt",sep=''),sep='\t')
  exonTags$exon<-c(cell.alt3$exonTags,cell.alt5$exonTags,cell.cass$exonTags,
                   cell.iret$junctionTags,cell.mutx$X5.ExonTags,
                   cell.taca$exonTags)
  colnames(exonTags)[ncol(exonTags)]<-noquote(sub(pattern="-",replacement = ".",
                                                  x=sub(pattern = "_countit_out",replacement = "",x=LT.HSC_filelist[i])))
}

##create matrix for Prog cells####
Prog_countit_out_folder<-"./data/countit_out/Prog/"
Prog_filelist<-list.files(Prog_countit_out_folder)
i<-1
for (i in (1:752)){
  cell.alt3<-read.csv(paste(Prog_countit_out_folder,Prog_filelist[i],"/alt3.count.txt",sep=''),sep='\t')
  cell.alt5<-read.csv(paste(Prog_countit_out_folder,Prog_filelist[i],"/alt5.count.txt",sep=''),sep='\t')
  cell.cass<-read.csv(paste(Prog_countit_out_folder,Prog_filelist[i],"/cass.count.txt",sep=''),sep='\t')
  cell.iret<-read.csv(paste(Prog_countit_out_folder,Prog_filelist[i],"/iret.count.txt",sep=''),sep='\t')
  cell.mutx<-read.csv(paste(Prog_countit_out_folder,Prog_filelist[i],"/mutx.count.txt",sep=''),sep='\t')
  cell.taca<-read.csv(paste(Prog_countit_out_folder,Prog_filelist[i],"/taca.count.txt",sep=''),sep='\t')
  
  exonTags$exon<-c(cell.alt3$exonTags,cell.alt5$exonTags,cell.cass$exonTags,
                   cell.iret$junctionTags,cell.mutx$X5.ExonTags,
                   cell.taca$exonTags)
  colnames(exonTags)[ncol(exonTags)]<-noquote(sub(pattern="-",replacement = ".",
                                                  x=sub(pattern = "_countit_out",replacement = "",x=Prog_filelist[i])))
}

write.table(x = colnames(exonTags),file = "./result/20220924_SplicingMatrixColnames.txt")

###create splicing seurat####
splicing_seurat<-CreateSeuratObject(counts =exonTags[,-1],project = "splicing")
counts <- GetAssayData(object = splicing_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- matrixStats::rowMaxs(as.matrix(counts)) >= 20 & 
  (Matrix::rowSums(nonzero[,1:696]) >= (696*0.5) |
  Matrix::rowSums(nonzero[,697:864]) >= (168*0.5) |
  Matrix::rowSums(nonzero[,865:1616]) >= (752*0.5))

filtered_counts <- counts[keep_genes, ]
dim(filtered_counts)
filtered_splicing_seurat <- CreateSeuratObject(filtered_counts, 
                                        meta.data = splicing_seurat@meta.data)

filtered_splicing_seurat <- SCTransform(filtered_splicing_seurat,
                               variable.features.n=3000)

###UMAP visualization
# Run PCA
filtered_splicing_seurat <- RunPCA(object = filtered_splicing_seurat,
                            ndims.print = 1:50,nfeatures.print = 5)
# Plot PCA
PCAPlot(filtered_splicing_seurat) 
ElbowPlot(filtered_splicing_seurat,ndims = 50)
# Run UMAP
filtered_splicing_seurat <- RunUMAP(filtered_splicing_seurat, 
                             dims = 1:50,
                             reduction = "pca")
# Plot UMAP                             
DimPlot(filtered_splicing_seurat,reduction = "umap",group.by = "orig.ident")

# Run TSNE
filtered_splicing_seurat <- RunTSNE(filtered_splicing_seurat, 
                             dims = 1:50,
                             reduction = "pca")
filtered_splicing_seurat <- RunTSNE(filtered_splicing_seurat,
                             dims = 1:50,tsne.method = "FIt-SNE",
                             reduction = "pca",perplexity=10)
# Plot TSNE                             
DimPlot(filtered_splicing_seurat,reduction = "tsne",
        group.by = "orig.ident",pt.size = 1)

# Determine the K-nearest neighbor graph
filtered_splicing_seurat <- FindNeighbors(object = filtered_splicing_seurat, 
                                   dims = 1:50)
# Determine the clusters for various resolutions                                
filtered_splicing_seurat <- FindClusters(object = filtered_splicing_seurat,
                                  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.2,1.35,1.5,1.75,2,2.5,3))


##plot with splicing pattern cluster####
Idents(object = filtered_splicing_seurat) <- "SCT_snn_res.1"
pSplicing<-DimPlot(filtered_splicing_seurat,
        reduction = "tsne",
        label = T,repel = T,
        label.size = 6,
        pt.size = 1)+
  ggtitle("Exon number cluster")

View(filtered_splicing_seurat@meta.data)  
dim(seurat_integrated@meta.data)
rownames(filtered_splicing_seurat@meta.data)
##plot with transcriptome cluster####

filtered_splicing_seurat@meta.data$TranscriptomeCluster<-
  seurat_integrated@meta.data[sub(pattern = "_",replacement ="-" ,x=rownames(filtered_splicing_seurat@meta.data)),]$SCT_snn_res.1
  
View(filtered_splicing_seurat@meta.data)
Idents(object = filtered_splicing_seurat) <-"TranscriptomeCluster"
pTranscriptome<-DimPlot(filtered_splicing_seurat,
                   reduction = "tsne",
                   label = T,repel = T,
                   label.size = 6,
                   pt.size = 1)+
  ggtitle("Transcriptome cluster")
##plot with sorting index####
Idents(object = filtered_splicing_seurat) <- "orig.ident"
pIndex<-DimPlot(filtered_splicing_seurat,
        reduction = "tsne",
        label = T,repel = T,
        label.size = 6,
        pt.size = 1,
        cols=c("#7cae00","#f8766d","#00a9ff"))+
  ggtitle("Cell sorting index cluster")

pSplicingUMAP<-pIndex+pTranscriptome+pSplicing
# ggsave(filename = "20220803_pSplicingUMAP.pdf",plot = pSplicingUMAP,
#         device = "pdf",width = 20,height = 6)

###add splicing pattern into seurat_integrated@meta.data matadata####
transcriptomeMetadata<-seurat_integrated@meta.data
seurat_integrated@meta.data$SplicingCluster<-
  filtered_splicing_seurat@meta.data[sub(pattern = "-",replacement = "_",x=rownames(transcriptomeMetadata)),]$SCT_snn_res.0.5
View(seurat_integrated@meta.data)
##plot with transcriptome cluster####
Idents(object = seurat_integrated) <- "SCT_snn_res.1"
pTranscriptome0<-DimPlot(seurat_integrated,
                        reduction = "tsne",
                        label = T,repel = T,
                        label.size = 6,
                        pt.size = 1
                        # ,cols=c("#f8766d","#cd9600", "#00bfc4",
                        #        "#7cae00","#00a9ff","#ff61cc")
                        )+
  ggtitle("Transcriptome cluster")

##plot with splicing pattern cluster####
Idents(object = seurat_integrated) <- "SplicingCluster"
pSplicing0<-DimPlot(seurat_integrated,
                   reduction = "tsne",
                   label = T,repel = T,
                   label.size = 6,
                   pt.size = 1,
                   cols=c("#cd9600","#7cae00","#f8766d","#00a9ff",
                          "#00bfc4","#ff61cc"))+
  ggtitle("Exon number cluster")

View(seurat_integrated@meta.data)  

##plot with sorting index####
Idents(object = seurat_integrated) <- "seq_folder"
pIndex0<-DimPlot(seurat_integrated,
                reduction = "tsne",
                label = T,repel = T,
                label.size = 6,
                pt.size = 1,
                cols=c("#7cae00","#f8766d","#00a9ff"))+
  ggtitle("Cell sorting index cluster")

pTranscriptomeTSNE<-pIndex0+pTranscriptome0+pSplicing0
# ggsave(filename = "./result/20220924_pTranscriptomeTSNE.pdf",plot = pTranscriptomeTSNE,
#        device = "pdf",width = 20,height = 6)

####Mpl_Exon9&10_PSI####
seurat_integrated@assays$SplicingRNA<-filtered_splicing_seurat@assays$RNA
DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "SCT_snn_res.0.3"
DefaultAssay(seurat_integrated) <- "SplicingRNA"
Mpl_Exon_PSI<-FeaturePlot(object = seurat_integrated,features = "35115",pt.size = 1.5,
            cols = c("grey","#f8766d"),min.cutoff = 'q0',max.cutoff = 'q80',
            label = T,label.size = 6)+
  ggtitle("Mpl_Exon9&10_JunctionCount")
ggsave(filename = "Mpl_Exon9&10_JunctionCount.pdf",plot =Mpl_Exon_PSI,
       device ="pdf")

####using countit_diff PSI results####
##LT.HSC
LT.HSC_countit_diff_folder<-"./data/countit_out/LT.HSC_countit_diff/"
LT.HSC_diff_filelist<-list.files(LT.HSC_countit_diff_folder)
i<-1
LT.HSC_diff.alt3<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/alt3.csv",sep=''),sep='\t')
LT.HSC_diff.alt5<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/alt5.csv",sep=''),sep='\t')
LT.HSC_diff.cass<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/cass.csv",sep=''),sep='\t')
LT.HSC_diff.iret<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/iret.csv",sep=''),sep='\t')
LT.HSC_diff.mutx<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/mutx.csv",sep=''),sep='\t')
LT.HSC_diff.taca<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/taca.csv",sep=''),sep='\t')

colnames(LT.HSC_diff.taca)
diff.dataframe<-data.frame("ID"=c(cell.alt3$name,cell.alt5$name,cell.cass$name,
                            cell.iret$name,cell.mutx$name,cell.taca$name))
diff.dataframe[length(LT.HSC_diff.alt3$name)+1,]
i<-1
for (i in (1:(168/2))){
  LT.HSC_diff.alt3<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/alt3.csv",sep=''),sep='\t')
  LT.HSC_diff.alt5<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/alt5.csv",sep=''),sep='\t')
  LT.HSC_diff.cass<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/cass.csv",sep=''),sep='\t')
  LT.HSC_diff.iret<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/iret.csv",sep=''),sep='\t')
  LT.HSC_diff.mutx<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/mutx.csv",sep=''),sep='\t')
  LT.HSC_diff.taca<-read.csv(paste(LT.HSC_countit_diff_folder,LT.HSC_diff_filelist[i],"/taca.csv",sep=''),sep='\t')
  
  diff.dataframe$exon_I1<-c(LT.HSC_diff.alt3[,12],LT.HSC_diff.alt5[,12],LT.HSC_diff.cass[,11],
                            LT.HSC_diff.iret[,11],LT.HSC_diff.mutx[,11],LT.HSC_diff.taca[,11])
  diff.dataframe$exon_I2<-c(LT.HSC_diff.alt3[,13],LT.HSC_diff.alt5[,13],LT.HSC_diff.cass[,12],
                            LT.HSC_diff.iret[,12],LT.HSC_diff.mutx[,12],LT.HSC_diff.taca[,12])
  
  colnames(diff.dataframe)[ncol(diff.dataframe)-1]<-noquote(substr(x=(sub(pattern = "I_g1.",replacement = "",
                                                                        x =colnames(LT.HSC_diff.alt3)[12])),
                                                                 start = 1,stop = 10))
  colnames(diff.dataframe)[ncol(diff.dataframe)]<-noquote(substr(x=(sub(pattern = "I_g2.",replacement = "",
                                                                          x =colnames(LT.HSC_diff.alt3)[13])),
                                                                   start = 1,stop = 10))
}

##Prog
Prog_countit_diff_folder<-"./data/countit_out/Prog_countit_diff/"
Prog_diff_filelist<-list.files(Prog_countit_diff_folder)

i<-1
for (i in (1:(752/2))){
  Prog_diff.alt3<-read.csv(paste(Prog_countit_diff_folder,Prog_diff_filelist[i],"/alt3.csv",sep=''),sep='\t')
  Prog_diff.alt5<-read.csv(paste(Prog_countit_diff_folder,Prog_diff_filelist[i],"/alt5.csv",sep=''),sep='\t')
  Prog_diff.cass<-read.csv(paste(Prog_countit_diff_folder,Prog_diff_filelist[i],"/cass.csv",sep=''),sep='\t')
  Prog_diff.iret<-read.csv(paste(Prog_countit_diff_folder,Prog_diff_filelist[i],"/iret.csv",sep=''),sep='\t')
  Prog_diff.mutx<-read.csv(paste(Prog_countit_diff_folder,Prog_diff_filelist[i],"/mutx.csv",sep=''),sep='\t')
  Prog_diff.taca<-read.csv(paste(Prog_countit_diff_folder,Prog_diff_filelist[i],"/taca.csv",sep=''),sep='\t')
  
  diff.dataframe$exon_I1<-c(Prog_diff.alt3[,12],Prog_diff.alt5[,12],Prog_diff.cass[,11],
                            Prog_diff.iret[,11],Prog_diff.mutx[,11],Prog_diff.taca[,11])
  diff.dataframe$exon_I2<-c(Prog_diff.alt3[,13],Prog_diff.alt5[,13],Prog_diff.cass[,12],
                            Prog_diff.iret[,12],Prog_diff.mutx[,12],Prog_diff.taca[,12])
  
  colnames(diff.dataframe)[ncol(diff.dataframe)-1]<-noquote(substr(x=(sub(pattern = "I_g1.",replacement = "",
                                                                          x =colnames(Prog_diff.alt3)[12])),
                                                                   start = 1,stop = 8))
  colnames(diff.dataframe)[ncol(diff.dataframe)]<-noquote(substr(x=(sub(pattern = "I_g2.",replacement = "",
                                                                        x =colnames(Prog_diff.alt3)[13])),
                                                                 start = 1,stop = 8))
}
dim(diff.dataframe)
levels(diff.dataframe[,-1])

##HSPC
HSPC_countit_diff_folder<-"./data/countit_out/HSPC_countit_diff/"
HSPC_diff_filelist<-list.files(HSPC_countit_diff_folder)

i<-1
for (i in (1:(696/2))){
  HSPC_diff.alt3<-read.csv(paste(HSPC_countit_diff_folder,HSPC_diff_filelist[i],"/alt3.csv",sep=''),sep='\t')
  HSPC_diff.alt5<-read.csv(paste(HSPC_countit_diff_folder,HSPC_diff_filelist[i],"/alt5.csv",sep=''),sep='\t')
  HSPC_diff.cass<-read.csv(paste(HSPC_countit_diff_folder,HSPC_diff_filelist[i],"/cass.csv",sep=''),sep='\t')
  HSPC_diff.iret<-read.csv(paste(HSPC_countit_diff_folder,HSPC_diff_filelist[i],"/iret.csv",sep=''),sep='\t')
  HSPC_diff.mutx<-read.csv(paste(HSPC_countit_diff_folder,HSPC_diff_filelist[i],"/mutx.csv",sep=''),sep='\t')
  HSPC_diff.taca<-read.csv(paste(HSPC_countit_diff_folder,HSPC_diff_filelist[i],"/taca.csv",sep=''),sep='\t')
  
  diff.dataframe$exon_I1<-c(HSPC_diff.alt3[,12],HSPC_diff.alt5[,12],HSPC_diff.cass[,11],
                            HSPC_diff.iret[,11],HSPC_diff.mutx[,11],HSPC_diff.taca[,11])
  diff.dataframe$exon_I2<-c(HSPC_diff.alt3[,13],HSPC_diff.alt5[,13],HSPC_diff.cass[,12],
                            HSPC_diff.iret[,12],HSPC_diff.mutx[,12],HSPC_diff.taca[,12])
  
  colnames(diff.dataframe)[ncol(diff.dataframe)-1]<-noquote(substr(x=(sub(pattern = "I_g1.",replacement = "",
                                                                          x =colnames(HSPC_diff.alt3)[12])),
                                                                   start = 1,stop = 8))
  colnames(diff.dataframe)[ncol(diff.dataframe)]<-noquote(substr(x=(sub(pattern = "I_g2.",replacement = "",
                                                                        x =colnames(HSPC_diff.alt3)[13])),
                                                                 start = 1,stop = 8))
}
tail(diff.dataframe)
diff.dataframe[1:6,1:6]
##Replace NA with 0
diff.dataframe1<-diff.dataframe
diff.dataframe1[is.na(diff.dataframe1)]<-0
##Create Seurat object
diff.splicing_seurat1 <- CreateSeuratObject(diff.dataframe1[,-1])
diff.splicing_seurat1<-diff.splicing_seurat1
View(diff.splicing_seurat1@assays$RNA@data)
diff.splicing_seurat1 <- SCTransform(diff.splicing_seurat1,
                                     variable.features.n=5000)
# Run PCA
diff.splicing_seurat1 <- RunPCA(object = diff.splicing_seurat1,
                                ndims.print = 1:50,nfeatures.print = 5)
# Plot PCA
PCAPlot(diff.splicing_seurat1) 
ElbowPlot(diff.splicing_seurat1,ndims = 50)

# Run UMAP
diff.splicing_seurat1 <- RunUMAP(diff.splicing_seurat1, 
                                 dims = 1:10,
                                 reduction = "pca")
# Plot UMAP                             
DimPlot(diff.splicing_seurat1,
        reduction = "umap",group.by = "orig.ident")
# Run TSNE
diff.splicing_seurat1 <- RunTSNE(diff.splicing_seurat1, 
                                 dims = 1:50,
                                 reduction = "pca")
#diff.splicing_seurat1 <- RunTSNE(diff.splicing_seurat1,
#                                         dims = 1:50,tsne.method = "FIt-SNE",
#                                         reduction = "pca",perplexity=10)
# Plot TSNE                             
DimPlot(diff.splicing_seurat1,reduction = "tsne",
        group.by = "orig.ident",pt.size = 1)

# Determine the K-nearest neighbor graph
diff.splicing_seurat1 <- FindNeighbors(object = diff.splicing_seurat1, 
                                       dims = 1:50)
# Determine the clusters for various resolutions                                
diff.splicing_seurat1 <- FindClusters(object = diff.splicing_seurat1,
                                      resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5,2,3,5))

View(diff.splicing_seurat1@meta.data)
unique(diff.splicing_seurat1@meta.data$SCT_snn_res.0.4)
Idents(object = diff.splicing_seurat1) <- "SCT_snn_res.1"
DimPlot(diff.splicing_seurat1,
        reduction = "umap",
        label = T,repel = T,
        label.size = 6,
        pt.size = 1
        # cols=c("#cd9600","#7cae00","#f8766d","#00a9ff","#00bfc4","#ff61cc")
)+
  ggtitle("PSI cluster")

rownames(diff.splicing_seurat1@meta.data)


###add splicing pattern into seurat_integrated@meta.data matadata####
transcriptomeMetadata<-seurat_integrated@meta.data
seurat_integrated@meta.data$SplicingCluster<-
  filtered_splicing_seurat@meta.data[sub(pattern = "-",replacement = "_",x=rownames(transcriptomeMetadata)),]$SCT_snn_res.0.5
seurat_integrated@meta.data$PSICluster<-
  diff.splicing_seurat1@meta.data[sub(pattern = "-",replacement = "_",x=rownames(transcriptomeMetadata)),]$SCT_snn_res.0.6

View(seurat_integrated@meta.data)
##plot with transcriptome cluster####
Idents(object = seurat_integrated) <- "SCT_snn_res.0.3"
pTranscriptome0<-DimPlot(seurat_integrated,
                         reduction = "tsne",
                         label = T,repel = T,
                         label.size = 6,
                         pt.size = 1,
                         cols=c("#f8766d","#cd9600", "#00bfc4",
                                "#7cae00","#00a9ff","#ff61cc"))+
  ggtitle("Transcriptome cluster")

##plot with splicing pattern cluster####
Idents(object = seurat_integrated) <- "SplicingCluster"
pSplicing0<-DimPlot(seurat_integrated,
                    reduction = "tsne",
                    label = T,repel = T,
                    label.size = 6,
                    pt.size = 1,
                    cols=c("#cd9600","#7cae00","#f8766d","#00a9ff",
                           "#00bfc4","#ff61cc"))+
  ggtitle("Exon number cluster")

##plot with sorting index####
Idents(object = seurat_integrated) <- "seq_folder"
pIndex0<-DimPlot(seurat_integrated,
                 reduction = "tsne",
                 label = T,repel = T,
                 label.size = 6,
                 pt.size = 1,
                 cols=c("#7cae00","#f8766d","#00a9ff"))+
  ggtitle("Cell sorting index cluster")

##plot with PSI####
Idents(object = seurat_integrated) <- "PSICluster"
pPSI0<-DimPlot(seurat_integrated,
               reduction = "tsne",
               label = T,repel = T,
               label.size = 6,
               pt.size = 1)+
  ggtitle("PSI cluster")

pTSNE<-pIndex0+pTranscriptome0+pSplicing0+pPSI0

ggsave(filename = "./result/20220924_tsne_withPSI.pdf",plot = pTSNE,
       device ="pdf",width = 10,height=10 )

####Assign cell type identity###
#cell type annotation downloaded from http://blood.stemcells.cam.ac.uk/single_cell_atlas.html##
all_cell_type<-read.csv(file = "./data/all_cell_types.txt",sep='\t')
rownames(all_cell_type)<-sub(pattern = "LT-HSC",replacement = "LT.HSC",x = rownames(all_cell_type))
#rownames(all_cell_type)<-sub(pattern = "_",replacement = "-",x = rownames(all_cell_type))

TypesofCells_in_SeuratIntegrated<-all_cell_type[rownames(seurat_integrated@meta.data),
                                                c("LTHSC","LMPP","MPP","CMP","MEP","GMP")]
dim(TypesofCells_in_SeuratIntegrated)
#colnames(TypesofCells_in_SeuratIntegrated)<-sub(pattern = "_broad",replacement = "",x = colnames(TypesofCells_in_SeuratIntegrated))
sum(rowSums(TypesofCells_in_SeuratIntegrated)==0)
TypesofCells_in_SeuratIntegrated$celltype<-NA

k<-1;j<-1
for (k in 1:nrow(TypesofCells_in_SeuratIntegrated)) {
  for (j in (1:6)){
    if (TypesofCells_in_SeuratIntegrated[k,j]==1) {
      TypesofCells_in_SeuratIntegrated[k,"celltype"]<-colnames(TypesofCells_in_SeuratIntegrated)[j]
    }
  }
}

seurat_integrated@meta.data$celltype<-
  TypesofCells_in_SeuratIntegrated[rownames(seurat_integrated@meta.data),"celltype"]

View(seurat_integrated@meta.data)
TypesofCells_in_SeuratIntegrated["LT.HSC-018","celltype"]
pcelltype<-DimPlot(object =seurat_integrated,group.by = "celltype",
        reduction = "tsne",cols = c())+
  ggtitle("Cell type")


ggsave(filename = "./result/20221001_celltype.pdf",
       device = "pdf",plot = pcelltype)



###add cell type to diff.splicing_seurat1@meta.data###
diff.splicing_seurat1@meta.data$celltype<-
  TypesofCells_in_SeuratIntegrated[sub(pattern ="_" ,replacement ="-" ,x = rownames(diff.splicing_seurat1@meta.data)),"celltype"]
View(diff.splicing_seurat1@meta.data)

DimPlot(diff.splicing_seurat1,
        reduction = "tsne",
        label = T,repel = T,
        label.size = 6,
        pt.size = 1
        # cols=c("#cd9600","#7cae00","#f8766d","#00a9ff","#00bfc4","#ff61cc")
)+  ggtitle("PSI cluster")

DimPlot(diff.splicing_seurat1,
        reduction = "tsne",group.by = "celltype")

####20221001 Reanalysis transcriptome cluster per celltype plot#####
seurat_integrated1<-seurat_integrated
Idents(object = seurat_integrated1) <- "SCT_snn_res.1"
#rename transcriptome cluster name and plot transcriptome
seurat_integrated1<-RenameIdents(object = seurat_integrated1,
                                 '0'="LMPP",
                                 '1'="HSC",
                                 '2'="MEP",
                                 '3'="MEP",
                                 '4'="MPP",
                                 '5'="CMP/GMP",
                                 '6'="CMP/GMP",
                                 '7'="CMP/GMP",
                                 '8'="MEP",
                                 '9'="NA")
##color value
# "#3dc1ff", LMPP
# "#f8766d", HSC
# "#ffa33b", MEP
# "#7fcd00", MPP
# "#00d5fa", CMP/GMP
# "gray" NA
pTranscriptome<-DimPlot(seurat_integrated1,
        reduction = "tsne",
        label = T,repel = T,
        label.size = 6,
        pt.size = 1,
        cols=c("#3dc1ff","#f8766d","#ffa33b","#7fcd00","#00d5fa","gray"))+
ggtitle("transcriptome cluster")

DefaultAssay(seurat_integrated1) <- "RNA"

Vlnplot1<-VlnPlot(object = seurat_integrated1,
                  features = c("Kit","Ly6a","Procr","Trpc6", "Cd34","Flt3", "Slamf1","Cd48",
                               "Hoxb5","Fgd5",
                               "Mpl","Hba-a1","Smim1",
                               "Cd3e","Itgam","Cd19",
                               "Ly6g","Ly76",
                               "Fcgr2b","Fcgr3","Cd34"
                  ),cols=c("#3dc1ff","#f8766d","#ffa33b","#7fcd00","#00d5fa","gray"),
                  log=T,ncol = 5)
ggsave(filename  = "./result/20221001_VlnPlot_1.pdf",plot = Vlnplot1,
       device = "pdf",
       width=20,height = 20,limitsize = F)

VlnPlot(object = seurat_integrated1, 
        features = c("Cd3","Cd11b","Cd19",
                     "Ly6g","Ly76","Kitl",
                     "Sca1","Flt3",
                     "Slamf1","Cd34", "Cd48"
        ),
        log=T,ncol = 2)

HSCmarker<-FeaturePlot(object = seurat_integrated1,pt.size = 1.5,
                       cols = c("grey","#f8766d"),min.cutoff = 'q0',max.cutoff = 'q90',
                       features =c("Slamf1","Cd48","Procr","Mpl", "Hoxb5","Fgd5","Mecom",
                                   "Cd34","Flt3"),ncol = 3,label = T,label.size = 6)
ggsave(filename ="./result/20221001_featureplot_HSCmarkerUMAP.pdf",plot = HSCmarker,device = "pdf", 
       width =25,height = 20)
HSCmarkertsne<-FeaturePlot(object = seurat_integrated1,pt.size = 1.5,reduction = "tsne",
                           cols = c("grey","#f8766d"),min.cutoff = 'q0',max.cutoff = 'q90',
                           features =c("Slamf1","Cd48","Procr","Mpl", "Hoxb5","Fgd5","Mecom",
                                       "Cd34","Flt3"),ncol = 3,label = T,label.size = 6)
ggsave(filename ="./result/20221001_featureplot_HSCmarker_tSNE.pdf",plot = HSCmarkertsne,device = "pdf", 
       width =25,height = 20)
HSCmarkerVln<-VlnPlot(object = seurat_integrated1, 
                      features = c("Slamf1","Cd48","Procr","Mpl", "Hoxb5","Fgd5","Mecom",
                                   "Cd34","Flt3"),pt.size = 1,
                      cols=c("#3dc1ff","#f8766d","#ffa33b","#7fcd00","#00d5fa","gray"),
                      log=T,ncol = 3
                     )
ggsave(filename ="./result/20221001_Vlnplot_HSCmarker.pdf",plot = HSCmarkerVln,device = "pdf", 
       width =25,height = 20)

#plot cell type
Idents(object = seurat_integrated1) <- "celltype"

pCelltype<-DimPlot(RenameIdents(seurat_integrated1,
                     "LTHSC"="HSC",
                     "LMPP"="LMPP",
                     "MPP"="MPP",
                     "CMP"="CMP/GMP",
                     "GMP"="CMP/GMP",
                     "MEP"="MEP",
                     "NA"="NA"),
        reduction = "tsne",
        label = T,repel = T,
        label.size = 6,
        pt.size = 1,
        cols=c("#f8766d","#3dc1ff","#7fcd00","#00d5fa","#ffa33b","gray"))+
  ggtitle("cell type cluster")
#plot exon junction number 
seurat_integrated1@meta.data$SplicingCluster<-
  filtered_splicing_seurat@meta.data[rownames(seurat_integrated1@meta.data),]$SCT_snn_res.1.5
Idents(object = seurat_integrated1) <- "SplicingCluster"

pExonJunctionNumber<-DimPlot(RenameIdents(object = seurat_integrated1,
                                         '0'="LMPP",
                                         '1'="HSC",
                                         '2'="MEP",
                                         '3'="CMP/GMP",
                                         '4'="MEP",
                                         '5'="MPP",
                                         '6'="MPP",
                                         '7'="CMP/GMP",
                                         '8'="NA"),
        reduction = "tsne",
        label = T,repel = T,
        label.size = 6,
        pt.size = 1,
        cols=c("#3dc1ff","#f8766d","#ffa33b","#00d5fa","#7fcd00","gray"))+
  ggtitle("Exon junction number cluster")

##plot sorting index##
Idents(seurat_integrated1)<-"seq_folder"
pSortingIndex<-DimPlot(object = seurat_integrated1,reduction = "tsne",
                             label = T,repel = T,
                             label.size = 6,
                             pt.size = 1,
                       cols=c("#7fcd00","#f8766d","#00d5fa")
                            )+
  ggtitle("Seq folder cluster")

pTSNE_20221001<-pSortingIndex + pCelltype + pTranscriptome + pExonJunctionNumber
ggsave(filename = "./result/20221001_pTSNE.pdf",plot = pTSNE_20221001,device = "pdf",
      height = 20,width = 20)

Fig1a<-FeaturePlot(object = seurat_integrated1,pt.size = 1.5,reduction = "tsne",
                                  cols = c("grey","#f8766d"),min.cutoff = 'q0',max.cutoff = 'q90',
                                  features =c("Procr","Mpl",
                                              "Cd34","Flt3"),ncol = 2,label = T,label.size = 6)
FIg1b<- pTranscriptome + pExonJunctionNumber

ggsave(filename ="./result/20221001_Fig1a.pdf",plot = Fig1a,device = "pdf", 
       width =20,height = 20)
ggsave(filename ="./result/20221001_Fig1b.pdf",plot = FIg1b,device = "pdf", 
       width =20,height =10)
