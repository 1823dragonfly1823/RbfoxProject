 README file



 GSE81682_HTSeq_counts and AllSample_GEO table 
 1920 cells in total
 3840 sequencing files in total
 number i and (1920+i) are from the sample sample, but sequenced twice

 HSPC    LT.HSC   Prog 
   852        216    852 

post filter (subset= ((nUMI >= 200000) & 
                                      (nGene >= 4000)&
                                      (log10GenesPerUMI>0.6)&
                                      (mitoRatio < 0.1)&
                                      (RiboRatio<0.1)))
  HSPC LT.HSC   Prog 
   696    168    752 


   LT-HSC: Procr, Mpl, Hobx5, Fdg5, Ctnaal1
MEP: Gata1, Gypa


0: HSC
1: LMPP
2: CMP
3: GMP
4: MEP
5: MPP

in sample.csv file
rows 1-852, HSPC
rows 853-1068, LT-HSC
rows 1069-1920, Prog
rows 1921-2772 HSPC_Rep
rows 2773-2988, LT-HSC_Rep
rows 2989-3840, Prog_Rep

In /ocean/projects/tra110037p/lgao3/result/GSE81682 folder  (the sra cache folder limit is 20GB, emtpy this folder frequently otherwize there will be error)
../../software/sradownloader/sradownloader  --outdir XXXXXXX 'LT.HSC_Samples_postFilter.txt'  --noena

Olego and quantas:
#!/bin/bash

cd /ocean/projects/tra110037p/lgao3/result/GSE81682

files=(`ls ./FastqFiles/LT.HSC | grep fastq.gz `)
for i in {0..5}
do
        echo ${files[i]}
        echo ${files[i+168]}
        b=`basename "${files[i]}"  _Mus_musculus_RNA-Seq.fastq.gz`
        echo $b
        #searchstring="Long-term_Haematopoietic_stem_cell_"
        #catname=${b%*$searchstring}
        catname=${b##*Long-term_Haematopoietic_stem_cell_}
        echo $catname
        cat "./FastqFiles/LT.HSC/${files[i]}" "./FastqFiles/LT.HSC/${files[i+168]}" > "./FastqFiles/LT.HSC/${catname}.fastq.gz"

/ocean/projects/tra110037p/lgao3/olego/olego -v -t 32 -r /ocean/projects/tra110037p/lgao3/olego/models/mm.cfg -j /ocean/projects/tra110037p/lgao3/IndexAnnotationfiles/mm10.intron.hmr.bed -o "./FastqFiles/LT.HSC/${catname}.sam"  /ocean/projects/tra110037p/lgao3/IndexAnnotationfiles/GRCm38.primary_assembly.genome.fa.gz "./FastqFiles/LT.HSC/${catname}.fastq.gz"

perl /ocean/projects/tra110037p/lgao3/olego/sam2bed.pl --uniq -v "./FastqFiles/LT.HSC/${catname}.sam" "./FastqFiles/LT.HSC/${catname}.unique.bed"

perl /ocean/projects/tra110037p/lgao3/Quantas/quantas-master/countit/summarize_splicing_wrapper.pl -c ./"${catname}.cache" -v -big  -conf /ocean/projects/tra110037p/lgao3/IndexAnnotationfiles/mm10.conf -dbkey mm10 -cass -taca -alt5 -alt3 -mutx -iret "./FastqFiles/LT.HSC/${catname}.unique.bed" "./FastqFiles/LT.HSC/${catname}_countit_out"

done


evi1 (HSC)
Flt3cre
cd41


Pro fastq file canceled because of limited space
to continue, run 
../../software/sradownloader/sradownloader --outdir ./FastqFiles/Prog/ Prog_Samples_postFilter_Continued.txt --noena
(marked 20220719_10:39am)

count fule number
ls | wc -l

R extract string: 
sub(" xxx.*", "", x)                 # Extract characters before pattern
sub(".*xxx ", "", x)                 # Extract characters after pattern


LT-HSC_214.fastq.gz
LT-HSC_215.fastq.gz
LT-HSC_216.fastq.gz

Prog_quantas: 0..124
Prog_quantas_1: 125..249
Prog_quantas_2: 250..374
Prog_quantas_3: 375..499
Prog_quantas_4: 500..624
Prog_quantas_5: 625..751

HSPC_quantas_3: 0..44
HSPC_quantas_4: 45..88
HSPC_quantas_5: 89..132
HSPC_quantas_6: 133..177
HSPC_quantas_7: 178..222
HSPC_quantas_8: 223..268
slurm-10091172.out-----slurm-10091177.out

LT_HSC:
slurm-10051795.out: LT-HSC_001-LT-HSC_214
slurm-10082099.out: LT-HSC_215-LT-HSC_216 
HSPC: 
slurm-10065941.out: HSPC_488-HSPC_797 correct
slurm-10065773.out: HSPC_001-HSPC_212 correct

Pro:
slurm-10083216.out: Pro_002-Pro_140
slurm-10083217.out: Pro_144-Prog_277
slurm-10083218.out: Pro_292-Prog_417
slurm-10083219.out: Pro_433-Prog_573
slurm-10083220.out: Pro_574-Prog_693
slurm-10083221.out: Pro_718-Prog_852



Prog_140.fastq.gz
Prog_278.fastq.gz
Prog_419.fastq.gz
Prog_695.fastq.gzcd

HSPC_Samples_postFilter[-c(1:177,401:647,697:873,1097:1343),]

HSPC_quantas_101  0..44
HSPC_quantas_102  45..89
HSPC_quantas_103  90..135
HSPC_quantas_104  136..179
HSPC_quantas_105  180..225
HSPC_quantas_106 226..271



 slurm-10121957.out: HSPC_425--430
  slurm-10120403.out: Pro_702--709  717


  HSPC: correct 424 + HSPC_not_correct_counotit_REDO20220801 272



  20220803: 
  create quantas countit_diff group.conf files
  
  #!/bin/bash
  
files=(`ls ./countit_out | grep countit_out`)
n=(`ls ./countit_out | grep countit_out | wc -l`)
echo $n
for i in {0..694..2}
do
        echo $i
        b1=`basename "${files[i]}"  _countit_out`
        b2=`basename "${files[i+1]}"  _countit_out`
        echo $b1
        echo $b2
        printf '%s\t%s\n' ./countit_out/${files[i]} $b1 > "./groupconf/${b1}_${b2}.group.conf"
        printf '%s\t%s\n' ./countit_out/${files[i+1]}  $b2 >>"./groupconf/${b1}_${b2}.group.conf"

done


Run quantas to analyze diff splicing event
#!/bin/bash
grouplist=(`ls ./groupconf | grep .group.conf`)

for j in {0..347}
do
        b=`basename "${grouplist[j]}"  .group.conf`
        echo $b
        mkdir "./countit_diff/${b}"
        type=("cass" "taca" "alt3" "alt5" "mutx" "iret")
        for i in "${type[@]}"
        do
        perl /ocean/projects/tra110037p/lgao3/Quantas/quantas-master/countit/test_splicing_diff.pl -type $i -v --min-cov 20 --id2gene2symbol /ocean/projects/tra110037p/lgao3/IndexAnnotationfiles/mm10/Mm.seq.all.AS.chrom.can.id2gene2symbol "./groupconf/${b}.group.conf" "./countit_diff/${b}/${i}.csv"
        done
done
