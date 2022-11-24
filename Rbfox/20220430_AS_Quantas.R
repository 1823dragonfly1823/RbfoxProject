library("dplyr")


##Load expression data
Expr<-read.table("/Users/mac data/Data/Seq/NextSeq/Rbfox/202110/AnalysisResults/20220122_huijuan/CtrlvsMut.expr.diff.nogroup4.txt",header = T)
rownames(Expr)<-Expr$gene.ids
dim(Expr)

##Filter out genes whose RPKM is lower than median of all genes in both group
RPKM.Ctrl.median<-median(Expr$RPKM.Ctrl.)
RPKM.Mut.median<-median(Expr$RPKM.Mut.)
Expr.filtered<-Expr[Expr$RPKM.Ctrl.>RPKM.Ctrl.median | Expr$RPKM.Mut.>RPKM.Mut.median,]
##Filter out top 25% genes 
RPKM.Ctrl.top25<-top_frac(Expr,0.25,RPKM.Ctrl.)
RPKM.Mut.top25<-top_frac(Expr,0.25,RPKM.Mut.)
Expr.top25<-Expr[Expr$gene.ids %in% RPKM.Ctrl.top25$gene.ids | Expr$gene.ids %in% RPKM.Mut.top25$gene.ids,]

##Load Splicing data_alt3
alt3<-read.table("alt3_dataset.diff.txt",header = T)
head(alt3)
dim(alt3)
alt3Geneid<-sub("//.*","",alt3$gene)
alt3$geneid<-alt3Geneid

alt3.filtered<-alt3[alt3$geneid %in% Expr.filtered$gene.ids,]
alt3.filtered.sig<-alt3.filtered[abs(alt3.filtered$dI_g1_vs_g2)>0.1 & alt3.filtered$FDR <0.05,]
write.csv(alt3.filtered.sig,"CtrlvsMut.alt3.diff_medianfiltered_dI0.1_FDR0.05.csv")

alt3.top25<-alt3[alt3$geneid %in% Expr.top25$gene.ids,]
alt3.top25.sig<-alt3.top25[abs(alt3.top25$dI_g1_vs_g2)>0.1 & alt3.top25$FDR <0.05,]
write.csv(alt3.top25.sig,"CtrlvsMut.alt3.diff_top25_dI0.1_FDR0.05.csv")

##Load Splicing data_alt5
alt5<-read.table("alt5_dataset.diff.txt",header = T)
head(alt5)
dim(alt5)
alt5Geneid<-sub("//.*","",alt5$gene)
alt5$geneid<-alt5Geneid

alt5.filtered<-alt5[alt5$geneid %in% Expr.filtered$gene.ids,]
alt5.filtered.sig<-alt5.filtered[abs(alt5.filtered$dI_g1_vs_g2)>0.1 & alt5.filtered$FDR <0.05,]
write.csv(alt5.filtered.sig,"CtrlvsMut.alt5.diff_medianfiltered_dI0.1_FDR0.05.csv")

alt5.top25<-alt5[alt5$geneid %in% Expr.top25$gene.ids,]
alt5.top25.sig<-alt5.top25[abs(alt5.top25$dI_g1_vs_g2)>0.1 & alt5.top25$FDR <0.05,]
write.csv(alt5.top25.sig,"CtrlvsMut.alt5.diff_top25_dI0.1_FDR0.05.csv")

##Load Splicing data_iret
iret<-read.table("iret_dataset.diff.txt",header = T)
head(iret)
dim(iret)
iretGeneid<-sub("//.*","",iret$gene)
iret$geneid<-iretGeneid

iret.filtered<-iret[iret$geneid %in% Expr.filtered$gene.ids,]
iret.filtered.sig<-iret.filtered[abs(iret.filtered$dI_g1_vs_g2)>0.1 & iret.filtered$FDR <0.05,]
write.csv(iret.filtered.sig,"CtrlvsMut.iret.diff_medianfiltered_dI0.1_FDR0.05.csv")

iret.top25<-iret[iret$geneid %in% Expr.top25$gene.ids,]
iret.top25.sig<-iret.top25[abs(iret.top25$dI_g1_vs_g2)>0.1 & iret.top25$FDR <0.05,]
write.csv(iret.top25.sig,"CtrlvsMut.iret.diff_top25_dI0.1_FDR0.05.csv")

##Load Splicing data_cass
cass<-read.table("cass_dataset.diff.txt",header = T)
head(cass)
dim(cass)
cassGeneid<-sub("//.*","",cass$gene)
cass$geneid<-cassGeneid

cass.filtered<-cass[cass$geneid %in% Expr.filtered$gene.ids,]
cass.filtered.sig<-cass.filtered[abs(cass.filtered$dI_g1_vs_g2)>0.1 & cass.filtered$FDR <0.05,]
write.csv(cass.filtered.sig,"CtrlvsMut.cass.diff_medianfiltered_dI0.1_FDR0.05.csv")

cass.top25<-cass[cass$geneid %in% Expr.top25$gene.ids,]
cass.top25.sig<-cass.top25[abs(cass.top25$dI_g1_vs_g2)>0.1 & cass.top25$FDR <0.05,]
write.csv(cass.top25.sig,"CtrlvsMut.cass.diff_top25_dI0.1_FDR0.05.csv")

##Load Splicing data_taca
taca<-read.table("taca_dataset.diff.txt",header = T)
head(taca)
dim(taca)
tacaGeneid<-sub("//.*","",taca$gene)
taca$geneid<-tacaGeneid

taca.filtered<-taca[taca$geneid %in% Expr.filtered$gene.ids,]
taca.filtered.sig<-taca.filtered[abs(taca.filtered$dI_g1_vs_g2)>0.1 & taca.filtered$FDR <0.05,]
write.csv(taca.filtered.sig,"CtrlvsMut.taca.diff_medianfiltered_dI0.1_FDR0.05.csv")

taca.top25<-taca[taca$geneid %in% Expr.top25$gene.ids,]
taca.top25.sig<-taca.top25[abs(taca.top25$dI_g1_vs_g2)>0.1 & taca.top25$FDR <0.05,]
write.csv(taca.top25.sig,"CtrlvsMut.taca.diff_top25_dI0.1_FDR0.05.csv")


##Load Splicing data_mutx
mutx<-read.table("mutx_dataset.diff.txt",header = T)
head(mutx)
dim(mutx)
mutxGeneid<-sub("//.*","",mutx$gene)
mutx$geneid<-mutxGeneid

mutx.filtered<-mutx[mutx$geneid %in% Expr.filtered$gene.ids,]
mutx.filtered.sig<-mutx.filtered[abs(mutx.filtered$dI_g1_vs_g2)>0.1 & mutx.filtered$FDR <0.05,]
write.csv(mutx.filtered.sig,"CtrlvsMut.mutx.diff_medianfiltered_dI0.1_FDR0.05.csv")

mutx.top25<-mutx[mutx$geneid %in% Expr.top25$gene.ids,]
mutx.top25.sig<-mutx.top25[abs(mutx.top25$dI_g1_vs_g2)>0.1 & mutx.top25$FDR <0.05,]
write.csv(mutx.top25.sig,"CtrlvsMut.mutx.diff_top25_dI0.1_FDR0.05.csv")

