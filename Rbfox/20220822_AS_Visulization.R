library("dplyr")
library("pheatmap")
library("ggplot2")
library("xlsx")
library("tidyr")

rm(list=ls())
setwd("/Users/mac data/Data/Seq/NextSeq/Rbfox/20220321_trimmed/Olego-Quantas/20220822_Visulization_Coverage20_Diff0.1_FDR0.05")
#####alt3 results####
alt3<-read.table("../alt3_dataset.diff.txt",header = T)
head(alt3)
dim(alt3)
alt3Geneid<-sub(".*//","",alt3$gene)
alt3$geneid<-alt3Geneid
alt3.noNA<-alt3[!is.na(alt3$dI_g1_vs_g2) & 
                  !is.na(alt3$FDR) &
                  alt3$coverage>20,]
dim(alt3.noNA)

#####Volcano plot with noNA result####
# add a column of NAs
alt3.noNA$diffexpressed <- "No"
# set "UP" 
alt3.noNA$diffexpressed[alt3.noNA$dI_g1_vs_g2 > 0.1 & alt3.noNA$FDR < 0.05] <- "Down"
# set "DOWN"
alt3.noNA$diffexpressed[alt3.noNA$dI_g1_vs_g2 < -0.1 & alt3.noNA$FDR < 0.05] <- "Up"
sum(alt3.noNA$diffexpressed=="No")
sum(alt3.noNA$diffexpressed=="Up")
sum(alt3.noNA$diffexpressed=="Down")

pAlt3<-ggplot(alt3.noNA,aes(x=-dI_g1_vs_g2,y=-log10(FDR),col=diffexpressed))+
  geom_point(alpha=12,size=1,pch=16) +
  #geom_point(x=-Allgenes[Allgenes$X=="22433",]$logFC,
  #            y=-log10(Allgenes[Allgenes$X=="22433",]$adj.P.Val),color="red",size=3)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"),
                     labels=c(paste("Down_",sum(alt3.noNA$diffexpressed=="Down"),sep=""),
                              paste("No_",sum(alt3.noNA$diffexpressed=="No"),sep=""), 
                              paste("Up_",sum(alt3.noNA$diffexpressed=="Up"),sep="")))+ 
  theme(legend.position = "bottom",panel.grid=element_blank())+
  ggtitle("alt3")

ggsave(filename = "20220822_alt3.pdf",plot = pAlt3,device = "pdf",width = 3,height = 4)

#####alt5 results####
alt5<-read.table("../alt5_dataset.diff.txt",header = T)
head(alt5)
dim(alt5)
alt5Geneid<-sub(".*//","",alt5$gene)
alt5$geneid<-alt5Geneid
alt5.noNA<-alt5[!is.na(alt5$dI_g1_vs_g2) &
                  !is.na(alt5$FDR) &
                  alt5$coverage>20,]
dim(alt5.noNA)
#####Volcano plot with noNA result####
# add a column of NAs
alt5.noNA$diffexpressed <- "No"
# set "UP" 
alt5.noNA$diffexpressed[alt5.noNA$dI_g1_vs_g2 > 0.1 & alt5.noNA$FDR < 0.05] <- "Down"
# set "DOWN"
alt5.noNA$diffexpressed[alt5.noNA$dI_g1_vs_g2 < -0.1 & alt5.noNA$FDR < 0.05] <- "Up"
sum(alt5.noNA$diffexpressed=="No")
sum(alt5.noNA$diffexpressed=="Up")
sum(alt5.noNA$diffexpressed=="Down")

palt5<-ggplot(alt5.noNA,aes(x=-dI_g1_vs_g2,y=-log10(FDR),col=diffexpressed))+
  geom_point(alpha=12,size=1,pch=16) +
  #geom_point(x=-Allgenes[Allgenes$X=="22433",]$logFC,
  #            y=-log10(Allgenes[Allgenes$X=="22433",]$adj.P.Val),color="red",size=3)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"),
                     labels=c(paste("Down_",sum(alt5.noNA$diffexpressed=="Down"),sep=""),
                              paste("No_",sum(alt5.noNA$diffexpressed=="No"),sep=""), 
                              paste("Up_",sum(alt5.noNA$diffexpressed=="Up"),sep="")))+ 
  theme(legend.position = "bottom",panel.grid=element_blank())+
  ggtitle("alt5")

ggsave(filename = "20220822_alt5.pdf",plot = palt5,device = "pdf",width = 3,height = 4)

#####cass results####
cass<-read.table("../cass_dataset.diff.txt",header = T)
head(cass)
dim(cass)
cassGeneid<-sub(".*//","",cass$gene)
cass$geneid<-cassGeneid
cass.noNA<-cass[!is.na(cass$dI_g1_vs_g2) 
                & !is.na(cass$FDR)
                & cass$coverage>20,]
dim(cass.noNA)
#####Volcano plot with noNA result####
# add a column of NAs
cass.noNA$diffexpressed <- "No"
# set "UP" 
cass.noNA$diffexpressed[cass.noNA$dI_g1_vs_g2 > 0.1 & cass.noNA$FDR < 0.05] <- "Down"
# set "DOWN"
cass.noNA$diffexpressed[cass.noNA$dI_g1_vs_g2 < -0.1 & cass.noNA$FDR < 0.05] <- "Up"
sum(cass.noNA$diffexpressed=="No")
sum(cass.noNA$diffexpressed=="Up")
sum(cass.noNA$diffexpressed=="Down")

pcass<-ggplot(cass.noNA,aes(x=-dI_g1_vs_g2,y=-log10(FDR),col=diffexpressed))+
  geom_point(alpha=12,size=1,pch=16) +
  #geom_point(x=-Allgenes[Allgenes$X=="22433",]$logFC,
  #            y=-log10(Allgenes[Allgenes$X=="22433",]$adj.P.Val),color="red",size=3)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"),
                     labels=c(paste("Down_",sum(cass.noNA$diffexpressed=="Down"),sep=""),
                              paste("No_",sum(cass.noNA$diffexpressed=="No"),sep=""), 
                              paste("Up_",sum(cass.noNA$diffexpressed=="Up"),sep="")))+ 
  theme(legend.position = "bottom",panel.grid=element_blank())+
  ggtitle("cass")

ggsave(filename = "20220822_cass.pdf",plot = pcass,device = "pdf",width = 3,height = 4)

#####iret results####
iret<-read.table("../iret_dataset.diff.txt",header = T)
head(iret)
dim(iret)
iretGeneid<-sub(".*//","",iret$gene)
iret$geneid<-iretGeneid
iret.noNA<-iret[!is.na(iret$dI_g1_vs_g2) 
                & !is.na(iret$FDR)
                &iret$coverage>20,]
dim(iret.noNA)
#####Volcano plot with noNA result####
# add a column of NAs
iret.noNA$diffexpressed <- "No"
# set "UP" 
iret.noNA$diffexpressed[iret.noNA$dI_g1_vs_g2 > 0.1 & iret.noNA$FDR < 0.05] <- "Down"
# set "DOWN"
iret.noNA$diffexpressed[iret.noNA$dI_g1_vs_g2 < -0.1 & iret.noNA$FDR < 0.05] <- "Up"
sum(iret.noNA$diffexpressed=="No")
sum(iret.noNA$diffexpressed=="Up")
sum(iret.noNA$diffexpressed=="Down")

piret<-ggplot(iret.noNA,aes(x=-dI_g1_vs_g2,y=-log10(FDR),col=diffexpressed))+
  geom_point(alpha=12,size=1,pch=16) +
  #geom_point(x=-Allgenes[Allgenes$X=="22433",]$logFC,
  #            y=-log10(Allgenes[Allgenes$X=="22433",]$adj.P.Val),color="red",size=3)+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"),
                     labels=c(paste("Down_",sum(iret.noNA$diffexpressed=="Down"),sep=""),
                              paste("No_",sum(iret.noNA$diffexpressed=="No"),sep=""), 
                              paste("Up_",sum(iret.noNA$diffexpressed=="Up"),sep="")))+ 
  theme(legend.position = "bottom",panel.grid=element_blank())+
  ggtitle("iret")

ggsave(filename = "20220822_iret.pdf",plot = piret,device = "pdf",width = 3,height = 4)


#####mutx results####
mutx<-read.table("../mutx_dataset.diff.txt",header = T)
head(mutx)
dim(mutx)
mutxGeneid<-sub(".*//","",mutx$gene)
mutx$geneid<-mutxGeneid
mutx.noNA<-mutx[!is.na(mutx$dI_g1_vs_g2) 
                & !is.na(mutx$FDR)
                &mutx$coverage>20,]
dim(mutx.noNA)
#####Volcano plot with noNA result####
# add a column of NAs
mutx.noNA$diffexpressed <- "No"
# set "UP" 
mutx.noNA$diffexpressed[mutx.noNA$dI_g1_vs_g2 > 0.1 & mutx.noNA$FDR < 0.05] <- "Down"
# set "DOWN"
mutx.noNA$diffexpressed[mutx.noNA$dI_g1_vs_g2 < -0.1 & mutx.noNA$FDR < 0.05] <- "Up"
sum(mutx.noNA$diffexpressed=="No")
sum(mutx.noNA$diffexpressed=="Up")
sum(mutx.noNA$diffexpressed=="Down")

pmutx<-ggplot(mutx.noNA,aes(x=-dI_g1_vs_g2,y=-log10(FDR),col=diffexpressed))+
  geom_point(alpha=12,size=1,pch=16) +
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"),
                     labels=c(paste("Down_",sum(mutx.noNA$diffexpressed=="Down"),sep=""),
                              paste("No_",sum(mutx.noNA$diffexpressed=="No"),sep=""), 
                              paste("Up_",sum(mutx.noNA$diffexpressed=="Up"),sep="")))+ 
  theme(legend.position = "bottom",panel.grid=element_blank())+
  ggtitle("mutx")

ggsave(filename = "20220822_mutx.pdf",plot = pmutx,device = "pdf",width = 3,height = 4)

#####taca results####
taca<-read.table("../taca_dataset.diff.txt",header = T)
head(taca)
dim(taca)
tacaGeneid<-sub(".*//","",taca$gene)
taca$geneid<-tacaGeneid
taca.noNA<-taca[!is.na(taca$dI_g1_vs_g2) 
                & !is.na(taca$FDR)
                & taca$coverage>20,]
dim(taca.noNA)
label.x<--taca.noNA[taca.noNA$gene %in% "17480//Mpl",]$dI_g1_vs_g2
label.y<--log10(taca.noNA[taca.noNA$gene %in% "17480//Mpl",]$FDR)

#####Volcano plot with noNA result####
# add a column of NAs
taca.noNA$diffexpressed <- "No"
# set "UP" 
taca.noNA$diffexpressed[taca.noNA$dI_g1_vs_g2 > 0.1 & taca.noNA$FDR < 0.05] <- "Down"
# set "DOWN"
taca.noNA$diffexpressed[taca.noNA$dI_g1_vs_g2 < -0.1 & taca.noNA$FDR < 0.05] <- "Up"
sum(taca.noNA$diffexpressed=="No")
sum(taca.noNA$diffexpressed=="Up")
sum(taca.noNA$diffexpressed=="Down")

ptaca<-ggplot(taca.noNA,aes(x=-dI_g1_vs_g2,y=-log10(FDR),col=diffexpressed))+
  geom_point(alpha=12,size=1,pch=16) +
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"),
                     labels=c(paste("Down_",sum(taca.noNA$diffexpressed=="Down"),sep=""),
                              paste("No_",sum(taca.noNA$diffexpressed=="No"),sep=""), 
                              paste("Up_",sum(taca.noNA$diffexpressed=="Up"),sep="")))+ 
  theme(legend.position = "bottom",panel.grid=element_blank())+
  ggtitle("taca")+
 annotate("text",x=label.x,y=label.y,label='Mpl',
             color="#ff84fd",
             size=7 , angle=0, fontface="bold")+
  geom_point(x=label.x,y=label.y,color="#ff84fd",size=1)


ggsave(filename = "20220822_taca_MPL.pdf",plot = ptaca,device = "pdf",width = 3,height = 4)

Allevent<-rbind(alt3[,c("gene","type","name","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")],
                alt5[,c("gene","type","name","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")],
                cass[,c("gene","type","name","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")],
                iret[,c("gene","type","name","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")],
                mutx[,c("gene","type","name","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")],
                taca[,c("gene","type","name","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")])
write.csv(x =Allevent, file = "20220521bind.csv")


# 
# Allevent.noNA<-Allevent[!is.na(Allevent$I_g1.Control.) & !is.na(Allevent$I_g2.Mutant.),]
# 
# alt3_sig_events<-alt3[alt3$coverage>20 & abs(alt3$dI_g1_vs_g2)>0.1 & alt3$FDR<0.05,][,c("gene","type","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")]
# alt5_sig_events<-alt5[alt5$coverage>20 & abs(alt5$dI_g1_vs_g2)>0.1 & alt5$FDR<0.05,][,c("gene","type","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")]
# cass_sig_events<-cass[cass$coverage>20 & abs(cass$dI_g1_vs_g2)>0.1 & cass$FDR<0.05,][,c("gene","type","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")]
# iret_sig_events<-iret[iret$coverage>20 & abs(iret$dI_g1_vs_g2)>0.1 & iret$FDR<0.05,][,c("gene","type","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")]
# mutx_sig_events<-mutx[mutx$coverage>20 & abs(mutx$dI_g1_vs_g2)>0.1 & mutx$FDR<0.05,][,c("gene","type","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")]
# taca_sig_events<-taca[taca$coverage>20 & abs(taca$dI_g1_vs_g2)>0.1 & taca$FDR<0.05,][,c("gene","type","I_g1.Control.","I_g2.Mutant.","coverage","dI_g1_vs_g2","FDR")]
# All_sig_events<-rbind(alt3_sig_events,alt5_sig_events,cass_sig_events,
#                       iret_sig_events,mutx_sig_events,taca_sig_events)
# dim(All_sig_events)
# All_sig_events$Geneid<-sub(".*//","",All_sig_events$gene)
# SigGens<-All_sig_events$Geneid[!duplicated(All_sig_events$Geneid)]
# write.csv(file = "202207719_significant_spliced_genes_Coverage20_Diff0.1_FDR0.5.csv",x=SigGens)

Allevent$GeneID<-sub(".*//",'',Allevent$gene)
head(Allevent)
Allevent_sig<-Allevent[abs(Allevent$dI_g1_vs_g2)>0.1 & Allevent$FDR<0.05,]
head(Allevent_sig)
sum(Allevent_sig$coverage>20)
p101<-pheatmap(mat = Allevent_sig[,3:4],show_rownames = F,cellwidth = 50,cellheight = 0.1,
         color = colorRampPalette(c("#00118d", "white","#be9a44"))(1000))
ColAnno101<-data.frame("group"=c("Ctrl","Mut"))
rownames(ColAnno101)<-colnames(Allevent_sig[,3:4])

clusDesignation <- cutree(as.hclust(p110[[1]]), 2)
clusDesignation[clusDesignation==1]
clusDesignation[clusDesignation==3]
####density plot####
Allevent..long<-pivot_longer(data = Allevent,
  cols = c("I_g1.Control.","I_g2.Mutant." ),
  names_to = "Group",
  names_prefix = " ",
  values_to = "PSI",
  values_drop_na = TRUE)

##define scale_y_log2 function
library(scales)
scale_y_log2 <- function (...) 
{
  scale_y_continuous(..., trans = log2_trans())
}


AllDensity<-ggplot(Allevent.noNA.long,aes(x=PSI,color=Group))+
  geom_density(alpha=0.3)+ 
  scale_y_log10()+
  theme(axis.line = element_line(colour = "black", size = 0.3))+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(panel.grid = element_blank())

ggsave(filename = "20220521_AllEventDesnityPLot.pdf",plot = AllDensity,
       device = "pdf",width = 4,height = 4)


###20221018 Venn plot between DifferEvents of HSC_vs_MPPS in Trummp 2012 CSC and DifferEvents of wt_vs_Rbfox2del HSCs
Allevent$uniqueIndex<-paste("Num",c(1:nrow(Allevent)),sep="")
rownames(Allevent)<-Allevent$uniqueIndex
DiffeEvents<-Allevent[Allevent$coverage>20 & abs(Allevent$dI_g1_vs_g2)>0.1 & Allevent$FDR<0.05,]

write.csv(x =DiffeEvents ,file = "./20221018_DiffEvents_Coverage20_FDR0.05_dI_0.1.csv")
