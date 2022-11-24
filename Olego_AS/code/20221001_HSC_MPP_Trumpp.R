library("tidyr")
library("ggplot2")
library("pheatmap")
library("matrixStats")
library("seriation")
library("dendextend")
library("edgeR")
library("topGO")
library("dplyr")

####Load results from quantas pepline####
dir<-"./Olego_AS/data/countit_result/countit_diff/"
filelist<-list.files(dir)
groups<-filelist[stringr::str_detect(filelist,"alt3")]
groups<-sub(pattern ="_alt3.csv",replacement = "",x = groups)
###Load HSC vs MPPs###
HSC_MPPs.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPPs_cass.csv',header = T)
HSC_MPPs.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPPs_taca.csv',header = T)
HSC_MPPs.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPPs_alt3.csv',header = T)
HSC_MPPs.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPPs_alt5.csv',header = T)
HSC_MPPs.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPPs_mutx.csv',header = T)
HSC_MPPs.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPPs_iret.csv',header = T)

HSC_MPPs<-rbind(HSC_MPPs.cass[,c(1:15)],
                HSC_MPPs.taca[,c(1:15)],
                HSC_MPPs.alt5[,c(1:9,11:16)],
                HSC_MPPs.alt3[,c(1:9,11:16)],
                HSC_MPPs.mutx[,c(1:15)],
                HSC_MPPs.iret[,c(1:15)])
dim(HSC_MPPs)

HSC_MPPs$uniqueIndex<-paste("Num",c(1:nrow(HSC_MPPs)),sep="")
rownames(HSC_MPPs)<-HSC_MPPs$uniqueIndex
HSC_MPPs.sig<-HSC_MPPs[HSC_MPPs$coverage>20 &
                         abs(HSC_MPPs$dI_g1_vs_g2)>0.1&
                         HSC_MPPs$FDR<0.05,]

write.csv(x = HSC_MPPs.sig,file = "./Olego_AS/result/20221018_HSC_vs_MPPs_DiffSplicingEvents.csv")

###Load HSC vs MPP1###
HSC_MPP1.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP1_cass.csv',header = T)
HSC_MPP1.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP1_taca.csv',header = T)
HSC_MPP1.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP1_alt3.csv',header = T)
HSC_MPP1.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP1_alt5.csv',header = T)
HSC_MPP1.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP1_mutx.csv',header = T)
HSC_MPP1.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP1_iret.csv',header = T)

HSC_MPP1<-rbind(HSC_MPP1.cass[,c(1:15)],
                HSC_MPP1.taca[,c(1:15)],
                HSC_MPP1.alt5[,c(1:9,11:16)],
                HSC_MPP1.alt3[,c(1:9,11:16)],
                HSC_MPP1.mutx[,c(1:15)],
                HSC_MPP1.iret[,c(1:15)])
dim(HSC_MPP1)

HSC_MPP1$uniqueIndex<-paste("Num",c(1:nrow(HSC_MPP1)),sep="")
rownames(HSC_MPP1)<-HSC_MPP1$uniqueIndex
HSC_MPP1.sig<-HSC_MPP1[HSC_MPP1$coverage>20 &
                         abs(HSC_MPP1$dI_g1_vs_g2)>0.1&
                         HSC_MPP1$FDR<0.05,]

write.csv(x = HSC_MPP1.sig,file = "./Olego_AS/result/20221018_HSC_vs_MPP1_DiffSplicingEvents.csv")

Wt_Rb2delHSC<-read.csv("/Users/mac\ data/Data/Seq/NextSeq/Rbfox/20220321_trimmed/Olego-Quantas/20220822_Visulization_Coverage20_Diff0.1_FDR0.05/20221018_DiffEvents_Coverage20_FDR0.05_dI_0.1.csv")

colnames(HSC_MPP1.sig)
Wt_Rb2delHSC.HSChigh<-Wt_Rb2delHSC[Wt_Rb2delHSC$dI_g1_vs_g2>0,]
Wt_Rb2delHSC.HSClow<-Wt_Rb2delHSC[Wt_Rb2delHSC$dI_g1_vs_g2<0,]
HSC_MPP1.sig.HSChigh<-HSC_MPP1.sig[HSC_MPP1.sig$dI_g1_vs_g2<0,]
HSC_MPP1.sig.HSClow<-HSC_MPP1.sig[HSC_MPP1.sig$dI_g1_vs_g2>0,]

x<-list(HSC_MPP1.sig.HSChigh$uniqueIndex,
        HSC_MPP1.sig.HSClow$uniqueIndex,
        Wt_Rb2delHSC.HSChigh$uniqueIndex,
        Wt_Rb2delHSC.HSClow$uniqueIndex)
library("ggVennDiagram")
ggVennDiagram(x, label_alpha = 0,
              category.names = c("HSC_MPP1.HSChigh" , "HSC_MPP1.HSClow " ,
                                 "Wt_Rb2delHSC.HSChigh", "Wt_Rb2delHSC.HSClow"))

###Load HSC vs MPP2###
HSC_MPP2.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP2_cass.csv',header = T)
HSC_MPP2.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP2_taca.csv',header = T)
HSC_MPP2.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP2_alt3.csv',header = T)
HSC_MPP2.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP2_alt5.csv',header = T)
HSC_MPP2.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP2_mutx.csv',header = T)
HSC_MPP2.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP2_iret.csv',header = T)

HSC_MPP2<-rbind(HSC_MPP2.cass[,c(1:15)],
                HSC_MPP2.taca[,c(1:15)],
                HSC_MPP2.alt5[,c(1:9,11:16)],
                HSC_MPP2.alt3[,c(1:9,11:16)],
                HSC_MPP2.mutx[,c(1:15)],
                HSC_MPP2.iret[,c(1:15)])
dim(HSC_MPP2)

###Load HSC vs MPP3###
HSC_MPP3.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP3_cass.csv',header = T)
HSC_MPP3.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP3_taca.csv',header = T)
HSC_MPP3.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP3_alt3.csv',header = T)
HSC_MPP3.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP3_alt5.csv',header = T)
HSC_MPP3.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP3_mutx.csv',header = T)
HSC_MPP3.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP3_iret.csv',header = T)

HSC_MPP3<-rbind(HSC_MPP3.cass[,c(1:15)],
                HSC_MPP3.taca[,c(1:15)],
                HSC_MPP3.alt5[,c(1:9,11:16)],
                HSC_MPP3.alt3[,c(1:9,11:16)],
                HSC_MPP3.mutx[,c(1:15)],
                HSC_MPP3.iret[,c(1:15)])
dim(HSC_MPP3)

###Load HSC vs MPP4###
HSC_MPP4.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP4_cass.csv',header = T)
HSC_MPP4.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP4_taca.csv',header = T)
HSC_MPP4.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP4_alt3.csv',header = T)
HSC_MPP4.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP4_alt5.csv',header = T)
HSC_MPP4.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP4_mutx.csv',header = T)
HSC_MPP4.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_MPP4_iret.csv',header = T)

HSC_MPP4<-rbind(HSC_MPP4.cass[,c(1:15)],
                HSC_MPP4.taca[,c(1:15)],
                HSC_MPP4.alt5[,c(1:9,11:16)],
                HSC_MPP4.alt3[,c(1:9,11:16)],
                HSC_MPP4.mutx[,c(1:15)],
                HSC_MPP4.iret[,c(1:15)])
dim(HSC_MPP4)

######HSC_2_MPP1_4######
HSC_2_MPP1_4.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_2_MPP1_4_cass.csv',header = T)
HSC_2_MPP1_4.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_2_MPP1_4_taca.csv',header = T)
HSC_2_MPP1_4.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_2_MPP1_4_alt3.csv',header = T)
HSC_2_MPP1_4.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_2_MPP1_4_alt5.csv',header = T)
HSC_2_MPP1_4.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_2_MPP1_4_mutx.csv',header = T)
HSC_2_MPP1_4.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_2_MPP1_4_iret.csv',header = T)

HSC_2_MPP1_4<-rbind(HSC_2_MPP1_4.cass[,c(1:15)],
                    HSC_2_MPP1_4.taca[,c(1:15)],
                    HSC_2_MPP1_4.alt5[,c(1:9,11:16)],
                    HSC_2_MPP1_4.alt3[,c(1:9,11:16)],
                    HSC_2_MPP1_4.mutx[,c(1:15)],
                    HSC_2_MPP1_4.iret[,c(1:15)])
# HSCvsMPP1.noNA<-HSCvsMPP1[!is.na(HSCvsMPP1$I_g1.HSC.) & !is.na(HSCvsMPP1$I_g2.MPP1.),]
# HSCvsMPP1.noNA.long<-pivot_longer(data = HSCvsMPP1.noNA,
#                              cols = c("I_g1.HSC.","I_g2.MPP1." ),
#                              names_to = "Group",
#                              names_prefix = " ",
#                              values_to = "PSI",
#                              values_drop_na = TRUE)
dim(HSC_2_MPP1_4)
HSC_2_MPP1_4[c(35,64,7645,34567),1]
colnames(HSC_2_MPP1_4)

######HSC_3_MPP4_2######
HSC_3_MPP4_2.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_3_MPP4_2_cass.csv',header = T)
HSC_3_MPP4_2.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_3_MPP4_2_taca.csv',header = T)
HSC_3_MPP4_2.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_3_MPP4_2_alt3.csv',header = T)
HSC_3_MPP4_2.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_3_MPP4_2_alt5.csv',header = T)
HSC_3_MPP4_2.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_3_MPP4_2_mutx.csv',header = T)
HSC_3_MPP4_2.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_3_MPP4_2_iret.csv',header = T)

HSC_3_MPP4_2<-rbind(HSC_3_MPP4_2.cass[,c(1:15)],
                    HSC_3_MPP4_2.taca[,c(1:15)],
                    HSC_3_MPP4_2.alt5[,c(1:9,11:16)],
                    HSC_3_MPP4_2.alt3[,c(1:9,11:16)],
                    HSC_3_MPP4_2.mutx[,c(1:15)],
                    HSC_3_MPP4_2.iret[,c(1:15)])

######HSC_4_MPP2_3######
HSC_4_MPP2_3.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_4_MPP2_3_cass.csv',header = T)
HSC_4_MPP2_3.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_4_MPP2_3_taca.csv',header = T)
HSC_4_MPP2_3.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_4_MPP2_3_alt3.csv',header = T)
HSC_4_MPP2_3.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_4_MPP2_3_alt5.csv',header = T)
HSC_4_MPP2_3.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_4_MPP2_3_mutx.csv',header = T)
HSC_4_MPP2_3.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/HSC_4_MPP2_3_iret.csv',header = T)

HSC_4_MPP2_3<-rbind(HSC_4_MPP2_3.cass[,c(1:15)],
                    HSC_4_MPP2_3.taca[,c(1:15)],
                    HSC_4_MPP2_3.alt5[,c(1:9,11:16)],
                    HSC_4_MPP2_3.alt3[,c(1:9,11:16)],
                    HSC_4_MPP2_3.mutx[,c(1:15)],
                    HSC_4_MPP2_3.iret[,c(1:15)])

######MPP1_3_MPP2_2######
MPP1_3_MPP2_2.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP1_3_MPP2_2_cass.csv',header = T)
MPP1_3_MPP2_2.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP1_3_MPP2_2_taca.csv',header = T)
MPP1_3_MPP2_2.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP1_3_MPP2_2_alt3.csv',header = T)
MPP1_3_MPP2_2.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP1_3_MPP2_2_alt5.csv',header = T)
MPP1_3_MPP2_2.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP1_3_MPP2_2_mutx.csv',header = T)
MPP1_3_MPP2_2.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP1_3_MPP2_2_iret.csv',header = T)

MPP1_3_MPP2_2<-rbind(MPP1_3_MPP2_2.cass[,c(1:15)],
                     MPP1_3_MPP2_2.taca[,c(1:15)],
                     MPP1_3_MPP2_2.alt5[,c(1:9,11:16)],
                     MPP1_3_MPP2_2.alt3[,c(1:9,11:16)],
                     MPP1_3_MPP2_2.mutx[,c(1:15)],
                     MPP1_3_MPP2_2.iret[,c(1:15)])

######MPP3_2_MPP2_1######
MPP3_2_MPP2_1.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_2_MPP2_1_cass.csv',header = T)
MPP3_2_MPP2_1.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_2_MPP2_1_taca.csv',header = T)
MPP3_2_MPP2_1.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_2_MPP2_1_alt3.csv',header = T)
MPP3_2_MPP2_1.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_2_MPP2_1_alt5.csv',header = T)
MPP3_2_MPP2_1.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_2_MPP2_1_mutx.csv',header = T)
MPP3_2_MPP2_1.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_2_MPP2_1_iret.csv',header = T)

MPP3_2_MPP2_1<-rbind(MPP3_2_MPP2_1.cass[,c(1:15)],
                     MPP3_2_MPP2_1.taca[,c(1:15)],
                     MPP3_2_MPP2_1.alt5[,c(1:9,11:16)],
                     MPP3_2_MPP2_1.alt3[,c(1:9,11:16)],
                     MPP3_2_MPP2_1.mutx[,c(1:15)],
                     MPP3_2_MPP2_1.iret[,c(1:15)])

######MPP3_3_HSC_1######
MPP3_3_HSC_1.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_3_HSC_1_cass.csv',header = T)
MPP3_3_HSC_1.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_3_HSC_1_taca.csv',header = T)
MPP3_3_HSC_1.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_3_HSC_1_alt3.csv',header = T)
MPP3_3_HSC_1.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_3_HSC_1_alt5.csv',header = T)
MPP3_3_HSC_1.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_3_HSC_1_mutx.csv',header = T)
MPP3_3_HSC_1.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP3_3_HSC_1_iret.csv',header = T)

MPP3_3_HSC_1<-rbind(MPP3_3_HSC_1.cass[,c(1:15)],
                    MPP3_3_HSC_1.taca[,c(1:15)],
                    MPP3_3_HSC_1.alt5[,c(1:9,11:16)],
                    MPP3_3_HSC_1.alt3[,c(1:9,11:16)],
                    MPP3_3_HSC_1.mutx[,c(1:15)],
                    MPP3_3_HSC_1.iret[,c(1:15)])

######MPP4_1_MPP1_1######
MPP4_1_MPP1_1.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_1_MPP1_1_cass.csv',header = T)
MPP4_1_MPP1_1.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_1_MPP1_1_taca.csv',header = T)
MPP4_1_MPP1_1.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_1_MPP1_1_alt3.csv',header = T)
MPP4_1_MPP1_1.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_1_MPP1_1_alt5.csv',header = T)
MPP4_1_MPP1_1.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_1_MPP1_1_mutx.csv',header = T)
MPP4_1_MPP1_1.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_1_MPP1_1_iret.csv',header = T)

MPP4_1_MPP1_1<-rbind(MPP4_1_MPP1_1.cass[,c(1:15)],
                     MPP4_1_MPP1_1.taca[,c(1:15)],
                     MPP4_1_MPP1_1.alt5[,c(1:9,11:16)],
                     MPP4_1_MPP1_1.alt3[,c(1:9,11:16)],
                     MPP4_1_MPP1_1.mutx[,c(1:15)],
                     MPP4_1_MPP1_1.iret[,c(1:15)])

######MPP4_3_MPP3_1######
MPP4_3_MPP3_1.cass<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_3_MPP3_1_cass.csv',header = T)
MPP4_3_MPP3_1.taca<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_3_MPP3_1_taca.csv',header = T)
MPP4_3_MPP3_1.alt3<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_3_MPP3_1_alt3.csv',header = T)
MPP4_3_MPP3_1.alt5<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_3_MPP3_1_alt5.csv',header = T)
MPP4_3_MPP3_1.mutx<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_3_MPP3_1_mutx.csv',header = T)
MPP4_3_MPP3_1.iret<-read.table('./Olego_AS/data/countit_result/countit_diff/MPP4_3_MPP3_1_iret.csv',header = T)

MPP4_3_MPP3_1<-rbind(MPP4_3_MPP3_1.cass[,c(1:15)],
                     MPP4_3_MPP3_1.taca[,c(1:15)],
                     MPP4_3_MPP3_1.alt5[,c(1:9,11:16)],
                     MPP4_3_MPP3_1.alt3[,c(1:9,11:16)],
                     MPP4_3_MPP3_1.mutx[,c(1:15)],
                     MPP4_3_MPP3_1.iret[,c(1:15)])

####Combine all individual samples###
HSC_MPP_AllIndividual<-cbind(HSC_MPPs,HSC_MPP1[,c(10:15)],HSC_MPP2[,c(10:15)],
                             HSC_MPP3[,c(10:15)],HSC_MPP4[,c(10:15)],
                             HSC_2_MPP1_4[,c(11,12)],HSC_3_MPP4_2[,c(11,12)],
                             HSC_4_MPP2_3[,c(11,12)],MPP1_3_MPP2_2[,c(11,12)],
                             MPP3_2_MPP2_1[,c(11,12)],MPP3_3_HSC_1[,c(11,12)],
                             MPP4_1_MPP1_1[,c(11,12)],MPP4_3_MPP3_1[,c(11,12)])
dim(HSC_MPP_AllIndividual)

HSC_MPP_AllIndividual.nona<-na.omit(HSC_MPP_AllIndividual)
colnames(HSC_MPP_AllIndividual.nona)

HSC_MPPs_Sds<-rowSds(as.matrix(HSC_MPP_AllIndividual.nona[,c(12,11)]))
HSC_MPP1_Sds<-rowSds(as.matrix(HSC_MPP_AllIndividual.nona[,c(12,17)]))
HSC_MPP2_Sds<-rowSds(as.matrix(HSC_MPP_AllIndividual.nona[,c(12,23)]))
HSC_MPP3_Sds<-rowSds(as.matrix(HSC_MPP_AllIndividual.nona[,c(12,29)]))
HSC_MPP4_Sds<-rowSds(as.matrix(HSC_MPP_AllIndividual.nona[,c(12,35)]))
Sds_table<-cbind(HSC_MPP1_Sds,HSC_MPP2_Sds,HSC_MPP3_Sds,HSC_MPP4_Sds)
#select <- order(apply(Sds_table,1,max), 
 #               decreasing=TRUE)[1:1000]
select <- order(HSC_MPP1_Sds, 
                decreasing=TRUE)[1:300]
highvarigenes_counts <-HSC_MPP_AllIndividual.nona[select,]
colnames(highvarigenes_counts)[40:55]<-sub(pattern = ".*\\.",replacement = "",
                                           x = stringr::str_sub(colnames(highvarigenes_counts[,40:55]),end=-2))

p1<-pheatmap(highvarigenes_counts[,40:55],
         cluster_cols = T,treeheight_col = 40,
         cluster_rows = T,treeheight_row = 20,
         border_color = NA,show_rownames = F,
         scale = "row",display_numbers = F,
         color = colorRampPalette(c("blue", "white", "red"))(50))

RowAnno<-data.frame("type"=highvarigenes_counts$type)
rownames(RowAnno)<-rownames(highvarigenes_counts)

ColAnno<-data.frame("Cell"=stringr::str_sub(colnames(highvarigenes_counts[,40:55]),end=-3))
rownames(ColAnno)<-colnames(highvarigenes_counts[,40:55])
ann_colors <- list(
  Cell = c(HSC="#f8766d",MPP1="#7fcd00",MPP2="#ffa33b",
           MPP3="#00d5fa",MPP4="#3dc1ff"),
  type=c(alt3="#ff9289",alt5="#cda1ff",cass="#cabd00",
         iret="#ff81d1",mutx="#ff86ff",taca="#00dcc2"))

p1col_dend<-p1[[2]]
p1col_dend <- rotate(p1col_dend, order =c("HSC_1","HSC_4","HSC_2","HSC_3",
                                            "MPP1_3" ,"MPP1_1" ,"MPP1_4",
                                            "MPP2_3","MPP2_1","MPP2_2",
                                            "MPP3_1","MPP3_2","MPP3_3",
                                            "MPP4_3","MPP4_1","MPP4_2"))
p1row_dend<-p1[[1]]
p1row_dend <- dendextend::rotate(p1row_dend, order =rev(p1row_dend[["labels"]]))

p2<-pheatmap(highvarigenes_counts[,40:55],
             cluster_cols=as.hclust(p1col_dend),cutree_rows = 2,
             treeheight_col = 40,cutree_cols = 2,
             cluster_rows = as.hclust(p1row_dend),treeheight_row = 10,border_color = NA,show_rownames = F,
             scale = "row",display_numbers = F,annotation_row = RowAnno,
             annotation_col=ColAnno,annotation_colors = ann_colors,
             color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))
ggsave(filename = "./Olego_AS/result/20221001_HSC_MPP_MaxIndividualSd_top300_heatmap.pdf",
       plot=p2,device = "pdf",width = 10,height = 10)

clusDesignation <- cutree(as.hclust(p1[[1]]), 2)
clusDesignation[clusDesignation==1]
clusDesignation[clusDesignation==2]

highvarigenes_counts$clusDesignation<-clusDesignation
highvarigenes_counts$GeneGene<-sub(".*//","",highvarigenes_counts$gene)

write.csv(x =highvarigenes_counts,file = "./Olego_AS/result/20221001_HSC_MPP_MaxIndividualSd_top300.csv" )
write.csv(x =HSC_MPP_AllIndividual.nona,file = "./Olego_AS/result/20221001_HSC_MPPall.noNA.csv" )



####20221019_HSC_MPP1_sig events###
colnames(HSC_MPP_AllIndividual.nona)
HSC_MPP1_sig<-HSC_MPP_AllIndividual.nona[abs(HSC_MPP_AllIndividual.nona[,19])>0.1 & HSC_MPP_AllIndividual.nona[,21]<0.05 & HSC_MPP_AllIndividual.nona[,16]>20,] 

select1<- rownames(HSC_MPP1_sig)
highvarigenes_counts1 <-HSC_MPP_AllIndividual.nona[select1,]
colnames(highvarigenes_counts1)[40:55]<-sub(pattern = ".*\\.",replacement = "",
                                            x = stringr::str_sub(colnames(highvarigenes_counts1[,40:55]),end=-2))

p11<-pheatmap(highvarigenes_counts1[,40:55],
              cluster_cols = T,treeheight_col = 40,
              cluster_rows = T,treeheight_row = 20,
              border_color = NA,show_rownames = F,
              scale = "row",display_numbers = F,
              color = colorRampPalette(c("blue", "white", "red"))(50))

RowAnno1<-data.frame("type"=highvarigenes_counts1$type)
rownames(RowAnno1)<-rownames(highvarigenes_counts1)

ColAnno1<-data.frame("Cell"=stringr::str_sub(colnames(highvarigenes_counts1[,40:55]),end=-3))
rownames(ColAnno1)<-colnames(highvarigenes_counts1[,40:55])
ann_colors1 <- list(
  Cell = c(HSC="#f8766d",MPP1="#7fcd00",MPP2="#ffa33b",
           MPP3="#00d5fa",MPP4="#3dc1ff"),
  type=c(alt3="#ff9289",alt5="#cda1ff",cass="#cabd00",
         iret="#ff81d1",mutx="#ff86ff",taca="#00dcc2"))

p11col_dend<-p11[[2]]
p11col_dend <- rotate(p11col_dend, order =c("HSC_1","HSC_4","HSC_2","HSC_3",
                                            "MPP1_3" ,"MPP1_1" ,"MPP1_4",
                                            "MPP2_3","MPP2_1","MPP2_2",
                                            "MPP3_1","MPP3_2","MPP3_3",
                                            "MPP4_3","MPP4_1","MPP4_2"))
p11row_dend<-p11[[1]]
p11row_dend <- dendextend::rotate(p11row_dend, order =rev(p11row_dend[["labels"]]))

p12<-pheatmap(highvarigenes_counts1[,40:55],
              cluster_cols=as.hclust(p11col_dend),cutree_rows = 2,
              treeheight_col = 40,cutree_cols = 2,
              cluster_rows = as.hclust(p11row_dend),treeheight_row = 10,border_color = NA,show_rownames = F,
              scale = "row",display_numbers = F,annotation_row = RowAnno1,
              annotation_col=ColAnno1,annotation_colors = ann_colors1,
              color = colorRampPalette(c("#00118d", "white","#be9a44"))(50))
ggsave(filename = "./Olego_AS/result/20221019_HSC_MPP1_SigEvent_heatmap.pdf",
       plot=p12,device = "pdf",width = 10,height = 10)
