setwd("/scratch/Arun/Projects/lectin_indel_ase/meta5")
library(tidyr)
library(stringr)

filenames <- read.table(file = "ref_105_208_228_892_bams.list")
filenames$V1 <- gsub("_ref.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
filenames$V1 <- gsub("^161","136" ,filenames$V1)
filenames$V1 <- gsub("^208","161" ,filenames$V1)
filenames$V1 <- gsub("^228","208" ,filenames$V1)
filenames$V1 <- gsub("^280","217" ,filenames$V1)
filenames$V1 <- gsub("^358","228" ,filenames$V1)
filenames$V1 <- gsub("^386","280" ,filenames$V1)
filenames$V1 <- gsub("^406","350" ,filenames$V1)
filenames$V1 <- gsub("^409","386" ,filenames$V1)
filenames$V1 <- gsub("^437","406" ,filenames$V1)
filenames$V1 <- gsub("^486","409" ,filenames$V1)
filenames$V1 <- gsub("^509","486" ,filenames$V1)
filenames$V1 <- gsub("^517","509" ,filenames$V1)
filenames$V1 <- gsub("^563","517" ,filenames$V1)
filenames$V1 <- gsub("^774","787" ,filenames$V1)

allele_counts <- read.table(file = "HP_ref_105_208_228_892_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 34,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])

new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('ref',filenames[i,], sep='_'),paste('alt',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full <- as.data.frame(t(cbind(new2)))
colnames(full) <- "SNP"
full$SNP <- as.numeric(as.character(full$SNP))
full$Sample <- rownames(full)
full <- full %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))

full2 <- full
full2$Sample <- paste(full2$DGRP,full2$Source,full2$Allele, sep = "_")
full2$Sample2 <- paste(full2$DGRP,full2$Allele, sep = "_")

refplot <- ggplot(full2[full2$Allele == "ref",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  ggtitle("Reference allele")
altplot <- ggplot(full2[full2$Allele == "alt",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  ggtitle("437 allele")

grid.arrange(refplot,altplot, ncol = 1)


ase_indels <- read.csv(file = "ase_indels.csv", stringsAsFactors = F)
ase_indels$X7bp_ins_1 <- gsub("^A$","D",ase_indels$X7bp_ins_1)
ase_indels$X7bp_ins_1 <- gsub("^P$","A",ase_indels$X7bp_ins_1)
ase_indels$X21bp_ins_3 <- gsub("^A$","D",ase_indels$X21bp_ins_3)
ase_indels$X21bp_ins_3 <- gsub("^P$","A",ase_indels$X21bp_ins_3)
ase_indels$X8bp_ins_2 <- gsub("^P$","D",ase_indels$X8bp_ins_2)
ase_indels2 <- ase_indels
ase_indels$indel_hap <- paste(ase_indels$X7bp_ins_1,ase_indels$X8bp_ins_2,ase_indels$X21bp_ins_3)
ase_indels <- as.data.frame(cbind(ase_indels$line,ase_indels$indel_hap))
colnames(ase_indels)<- c("DGRP","indel_hap")

for (r in 1:nrow(ase_indels2)){
  if (as.character(ase_indels2$snp.at.pos.3.in.21bp.indel[r]) %in% "G"){
    ase_indels2$X21bp_ins_3[r] <- "P(G)"
  }
  if (as.character(ase_indels2$snp.at.pos.3.in.21bp.indel[r]) %in% "A"){
    ase_indels2$X21bp_ins_3[r] <- "P(A)"
  }
}

ase_indels2$indel_hap <- paste(ase_indels2$X7bp_ins_1,ase_indels2$X8bp_ins_2,ase_indels2$X21bp_ins_3)
ase_indels2 <- as.data.frame(cbind(ase_indels2$line,ase_indels2$indel_hap))
colnames(ase_indels2)<- c("DGRP","indel_hap")


full_ref <- full[full$Allele == "ref",]
full_alt <- full[full$Allele == "alt",]
full_af <- full_alt
full_af$AF <- full_ref$SNP/(full_alt$SNP+full_ref$SNP)
full_af <- left_join(full_af,ase_indels, by = "DGRP")
afplot <- ggplot(full_af,aes(x=indel_hap, y=AF, color = Source) )+
  geom_boxplot()+
  ggtitle("Reference AF")
afplot

full_ref <- full[full$Allele == "ref",]
full_alt <- full[full$Allele == "alt",]
full_af <- full_alt
full_af$AF <- full_ref$SNP/(full_alt$SNP+full_ref$SNP)
full_af2 <- left_join(full_af,ase_indels2, by = "DGRP")
afplot2 <- ggplot(full_af2,aes(x=indel_hap, y=AF, color = Source) )+
  geom_boxplot()+
  ggtitle("Reference AF")
afplot2


##Technical replicate

full_af_T1 <- full_af[full_af$TechRep == "T1",]
full_af_T2 <- full_af[full_af$TechRep == "T2",]
full_af_T1$sample <- paste(full_af_T1$DGRP,full_af_T1$Source,full_af_T1$BioRep, sep = "_")
full_af_T2$sample <- paste(full_af_T2$DGRP,full_af_T2$Source,full_af_T2$BioRep, sep = "_")

full_af_T <- inner_join(full_af_T1,full_af_T2,by="sample")
ggplot(full_af_T,aes(x=AF.x,y=AF.y,color=Source.x))+
  geom_point()+
  #facet_grid(cols = vars(Allele.x))+
  xlab("Technical replicate 1")+
  ylab("Technical replicate 2")+
  geom_abline(intercept=0,slope=1)

##Biological replicate

full_B1_T1 <- full[full$BioRep == "B1" & full$TechRep == "T1"& full$Source == "cDNA",]
full_B1_T2 <- full[full$BioRep == "B1" & full$TechRep == "T2"& full$Source == "cDNA",]
full_B2_T1 <- full[full$BioRep == "B2" & full$TechRep == "T1"& full$Source == "cDNA",]
full_B2_T2 <- full[full$BioRep == "B2" & full$TechRep == "T2"& full$Source == "cDNA",]

full_B1 <- full_B1_T1
full_B1$SNP <- full_B1_T1$SNP + full_B1_T2$SNP
full_B2 <- full_B2_T1
full_B2$SNP <- full_B2_T1$SNP + full_B2_T2$SNP

full_B1_ref <- full_B1[full_B1$Allele == "ref",]
full_B1_alt <- full_B1[full_B1$Allele == "alt",]
full_B1_af <- full_B1_alt
full_B1_af$AF <- full_B1_ref$SNP/(full_B1_alt$SNP+full_B1_ref$SNP)
full_B1_af$sample <- paste(full_B1_af$DGRP,full_B1_af$Source, sep = "_")

full_B2_ref <- full_B2[full_B2$Allele == "ref",]
full_B2_alt <- full_B2[full_B2$Allele == "alt",]
full_B2_af <- full_B2_alt
full_B2_af$AF <- full_B2_ref$SNP/(full_B2_alt$SNP+full_B2_ref$SNP)
full_B2_af$sample <- paste(full_B2_af$DGRP,full_B2_af$Source, sep = "_")

full_af_B <- inner_join(full_B1_af,full_B2_af,by="sample")
ggplot(full_af_B,aes(x=AF.x,y=AF.y,color=Source.x))+
  geom_point()+
  #facet_grid(cols = vars(Allele.x))+
  xlab("Biological replicate 1")+
  ylab("Biological replicate 2")+
  geom_abline(intercept=0,slope=1)

##expression variation
full_B1_T1 <- full[full$BioRep == "B1" & full$TechRep == "T1" & full$Source == "cDNA",]
full_B1_T2 <- full[full$BioRep == "B1" & full$TechRep == "T2"& full$Source == "cDNA",]
full_B2_T1 <- full[full$BioRep == "B2" & full$TechRep == "T1"& full$Source == "cDNA",]
full_B2_T2 <- full[full$BioRep == "B2" & full$TechRep == "T2"& full$Source == "cDNA",]

full_T1_gdna <- full[full$BioRep == "B1" & full$TechRep == "T1" & full$Source == "gDNA",]
full_T2_gdna <- full[full$BioRep == "B1" & full$TechRep == "T2"& full$Source == "gDNA",]

full_allsum <- full_B1_T1
full_allsum$SNP <- full_B1_T1$SNP + full_B1_T2$SNP + full_B2_T1$SNP + full_B2_T2$SNP
full_allsum_gdna <- full_T1_gdna
full_allsum_gdna$SNP <- full_T1_gdna$SNP + full_T2_gdna$SNP
full_allsum <- rbind(full_allsum,full_allsum_gdna)

full_allsum_ref <- full_allsum[full_allsum$Allele == "ref",]
full_allsum_alt <- full_allsum[full_allsum$Allele == "alt",]
full_allsum_af <- full_allsum_alt
full_allsum_af$AF <- full_allsum_ref$SNP/(full_allsum_alt$SNP+full_allsum_ref$SNP)
ggplot(full_allsum_af,aes(x=Source, y=AF) )+
  geom_point()+
  ggtitle("Total expression variation")

#View(full_allsum_af)

full_allsum_af2 <- left_join(full_allsum_af,ase_indels2, by = "DGRP")
full_allsum_af2 <- full_allsum_af2[full_allsum_af2$Source == "cDNA",]
full_allsum_af2$Expressed <- ""
for (r in 1:nrow(full_allsum_af2)){
  if (full_allsum_af2$AF[r] > 0.3){
    full_allsum_af2$Expressed[r] <- "E"
  }
  if (full_allsum_af2$AF[r] < 0.3){
    full_allsum_af2$Expressed[r] <- "NE"
  }
}

ggplot(full_allsum_af2, aes(x=indel_hap, fill = Expressed))+
  geom_histogram(stat="count")

qpcr <- read.csv(file= "qpcr_dgrp_lectin_exp.csv", header = T)
colnames(qpcr) <-c("DGRP", "Fold.induction")
qpcr$DGRP <- as.character(qpcr$DGRP)
full_allsum_af2_qpcr <- left_join(full_allsum_af2,qpcr,by = "DGRP")
ggplot(full_allsum_af2_qpcr, aes(y=AF,x=Fold.induction,color=indel_hap))+
  geom_point()

full_allsum_af2 <- left_join(full_allsum_af,ase_indels, by = "DGRP")
full_allsum_af2 <- full_allsum_af2[full_allsum_af2$Source == "cDNA",]
full_allsum_af2$Expressed <- ""
for (r in 1:nrow(full_allsum_af2)){
  if (full_allsum_af2$AF[r] > 0.3){
    full_allsum_af2$Expressed[r] <- "E"
  }
  if (full_allsum_af2$AF[r] < 0.3){
    full_allsum_af2$Expressed[r] <- "NE"
  }
}

ggplot(full_allsum_af2, aes(x=indel_hap, fill = Expressed))+
  geom_histogram(stat="count")+
  xlab("Indel haplotype")

dgrp_snps <- read.csv(file = "lectin_promoter_and_coding.csv")
dgrp_snps <- dgrp_snps[dgrp_snps$pos > 3717772,]
dgrp_snps <- dgrp_snps[dgrp_snps$pos < 3718772,]
dgrp_snps <- dgrp_snps[!dgrp_snps$pos == 3718040,]
rownames(dgrp_snps) <- dgrp_snps$id
dgrp_snps <- dgrp_snps[-c(1:9)]
dgrp_snps <- dgrp_snps[!dgrp_snps$line_437 == dgrp_snps$line_892,]
dgrp_snps <- as.data.frame(t(dgrp_snps))
rownames(dgrp_snps) <- gsub("line_","",rownames(dgrp_snps))
dgrp_snps$DGRP <- rownames(dgrp_snps)

full_allsum_af2 <- inner_join(full_allsum_af2,dgrp_snps,by = "DGRP")
for (i in 11:ncol(full_allsum_af2)){
  newt <- table(full_allsum_af2$Expressed,full_allsum_af2[,i])
  if (ncol(newt) == 3 ){
    newt <- newt[,-c(1)]
  }
  #newt <- prop.table(newt)
  #print(fisher.test(newt)$p.value)
  if (fisher.test(newt)$p.value < 0.05){
    #print(newt)
  }
}

full_allsum_af2$indel_hap <- as.character(full_allsum_af2$indel_hap)
full_allsum_af3 <- separate(full_allsum_af2, col = "indel_hap", into = c("7 BP Indel","8 BP Indel","21 BP Indel"), sep = " ", remove = F)
full_allsum_af3$`7 BP Indel` <- gsub("^A$","P",full_allsum_af3$`7 BP Indel`)
full_allsum_af3$`7 BP Indel` <- gsub("^D$","A",full_allsum_af3$`7 BP Indel`)
full_allsum_af3$`21 BP Indel` <- gsub("^A$","P",full_allsum_af3$`21 BP Indel`)
full_allsum_af3$`21 BP Indel` <- gsub("^D$","A",full_allsum_af3$`21 BP Indel`)
full_allsum_af3$`8 BP Indel` <- gsub("^D$","P",full_allsum_af3$`8 BP Indel`)

full_allsum_af4 <- cbind(full_allsum_af3[3],full_allsum_af3[12],full_allsum_af3[13:ncol(full_allsum_af3)],full_allsum_af3[11],full_allsum_af3[10])
full_allsum_af4 <- melt(full_allsum_af4, id.vars = c("DGRP","Expressed"), variable_name = "Position")
full_allsum_af4$Position <- gsub("2L_","",full_allsum_af4$Position)
full_allsum_af4 <- full_allsum_af4[!full_allsum_af4$value == "-",]
full_allsum_af4$value <- gsub("0","Reference allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("2","Alternate allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("^A$","Reference allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("^P$","Alternate allele",full_allsum_af4$value)

ggplot(data = full_allsum_af4) + 
  geom_bar(mapping = aes(x = factor(Position, levels = unique(Position)), fill = Expressed), position = "stack")+
  facet_grid(rows = vars(value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Number of lines")+
  xlab("Variant")
  

##exp association
for (i in 14:ncol(full_allsum_af3)){
  newt <- table(full_allsum_af3[,13],full_allsum_af3[,i])
  if (ncol(newt) == 3 ){
    newt <- newt[,-c(1)]
  }
  newt <- prop.table(newt,1)
  #print(fisher.test(newt)$p.value)
  if ((newt[1,1]==1 & newt[2,2] ==1)|(newt[1,2]==1 & newt[2,1] ==1)){
    print(colnames(full_allsum_af3[i]))
  }
}
##7bp association
for (i in 14:ncol(full_allsum_af3)){
  newt <- table(full_allsum_af3[,10],full_allsum_af3[,i])
  if (ncol(newt) == 3 ){
    newt <- newt[,-c(1)]
  }
  newt <- prop.table(newt,1)
  #print(fisher.test(newt)$p.value)
  if ((newt[1,1]==1 & newt[2,2] ==1)|(newt[1,2]==1 & newt[2,1] ==1)){
    print(colnames(full_allsum_af3[i]))
  }
}
##8bp association
for (i in 14:ncol(full_allsum_af3)){
  newt <- table(full_allsum_af3[,11],full_allsum_af3[,i])
  if (ncol(newt) == 3 ){
    newt <- newt[,-c(1)]
  }
  newt <- prop.table(newt,1)
  #print(fisher.test(newt)$p.value)
  if ((newt[1,1]==1 & newt[2,2] ==1)|(newt[1,2]==1 & newt[2,1] ==1)){
    print(colnames(full_allsum_af3[i]))
  }
}
##21bp association
for (i in 14:ncol(full_allsum_af3)){
  newt <- table(full_allsum_af3[,12],full_allsum_af3[,i])
  if (ncol(newt) == 3 ){
    newt <- newt[,-c(1)]
  }
  newt <- prop.table(newt,1)
  #print(fisher.test(newt)$p.value)
  if ((newt[1,1]==1 & newt[2,2] ==1)|(newt[1,2]==1 & newt[2,1] ==1)){
    print(colnames(full_allsum_af3[i]))
  }
}

for (i in 10:12){
  newt <- table(full_allsum_af3[,13],full_allsum_af3[,i])
  if (ncol(newt) == 3 ){
    newt <- newt[,-c(1)]
  }
  newt <- prop.table(newt,1)
  #print(fisher.test(newt)$p.value)
  if ((newt[1,1]==1 & newt[2,2] ==1)|(newt[1,2]==1 & newt[2,1] ==1)){
    print(colnames(full_allsum_af3[i]))
  }
}


##stop codons dgrp

individuals <- read.csv(file = "TableS1_individuals.csv")
individuals <- individuals[,c(1,4,10,11)]
datades <- read.csv(file = "TableS2_populations.csv")
individuals_des <- left_join(individuals,datades, by =c("Population"="Population.ID") )

anno <- read.csv(file = "lectin_anno_all.csv")
anno$all <- paste(anno$Pos,anno$Alt,sep=",")
anno_pos <- cbind(anno[1],anno[3])
anno <- anno[-c(1:2)]

dpgplectinaltvariants <- read.csv(file = "dpgplectinaltvariants.csv", check.names = FALSE)
dpgplectinaltvariants <- dpgplectinaltvariants[!names(dpgplectinaltvariants) %in% "ZW184"]
dpgplectinaltvariants$all <- paste(dpgplectinaltvariants$Pos,dpgplectinaltvariants$Alt,sep=",")
dpgplectinaltvariants <- right_join(anno,dpgplectinaltvariants, by =c("all") )
dpgplectinaltvariants <- right_join(anno_pos,dpgplectinaltvariants, by =c("Pos") )

for (r in 1:nrow(dpgplectinaltvariants)){
  if (is.na(dpgplectinaltvariants$ann.y[r])){
    dpgplectinaltvariants$ann.y[r] <- dpgplectinaltvariants$ann.x[r]
  }
}
dpgplectinaltvariants <- dpgplectinaltvariants[-c(3)]
colnames(dpgplectinaltvariants)[2] <- "ann"
rownames(dpgplectinaltvariants) <- paste(dpgplectinaltvariants$ann,make.names(dpgplectinaltvariants$Pos, unique=TRUE),sep="_")
dpgplectinaltvariants <- dpgplectinaltvariants[!duplicated(dpgplectinaltvariants$Pos),]
dpgplectinaltvariants <- dpgplectinaltvariants[-c(1:5)]
dpgplectinaltvariants <- as.data.frame(t(dpgplectinaltvariants))
samplenames <- as.data.frame(rownames(dpgplectinaltvariants))
colnames(samplenames) <- "Sample"
sranames_meta_des <- inner_join(samplenames,individuals_des, by = c("Sample"="Stock.ID"))
dpgplectinaltvariants$Sample <- rownames(dpgplectinaltvariants)
dpgplectinaltvariants <- right_join(sranames_meta_des,dpgplectinaltvariants, by =c("Sample") )
dpgplectinaltvariants[,16:ncol(dpgplectinaltvariants)] <- apply(dpgplectinaltvariants[,17:ncol(dpgplectinaltvariants)],2,as.numeric)

##check stopcodons in dgrp
fulldgrp <- dpgplectinaltvariants[dpgplectinaltvariants$Data.Group %in% c("DGRP"),]
fulldgrp <- fulldgrp[!is.na(fulldgrp$Country),]
locus2 <- fulldgrp[,16:ncol(fulldgrp)]
locus3 <- locus2[,grep("stop_gained",colnames(locus2))]
locus3 <- cbind(fulldgrp[1:15],locus3)
dgrp_stopcodon <- cbind(locus3[1],locus3[21])
dgrp_stopcodon$Sample <- gsub("RAL-","",dgrp_stopcodon$Sample)
colnames(dgrp_stopcodon)[1] <- "DGRP"

full_allsum_af5 <- left_join(full_allsum_af2,dgrp_stopcodon,by = "DGRP")

ggplot(full_allsum_af5, aes(x=as.factor(stop_gained_X3717487), fill = Expressed))+
  geom_histogram(stat="count")+
  xlab("Premature stop codon")

ggplot(full_allsum_af5[full_allsum_af5$stop_gained_X3717487 == 0,], aes(x=indel_hap, fill = Expressed))+
  geom_histogram(stat="count")+
  xlab("Premature stop codon")

##alt_136_161_217_280_371_409_427_584_787_820_859
filenames <- read.table(file = "alt_136_161_217_280_371_409_427_584_787_820_859_bams.list")
filenames$V1 <- gsub("_alt.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
filenames$V1 <- gsub("^161","136" ,filenames$V1)
filenames$V1 <- gsub("^208","161" ,filenames$V1)
filenames$V1 <- gsub("^228","208" ,filenames$V1)
filenames$V1 <- gsub("^280","217" ,filenames$V1)
filenames$V1 <- gsub("^358","228" ,filenames$V1)
filenames$V1 <- gsub("^386","280" ,filenames$V1)
filenames$V1 <- gsub("^406","350" ,filenames$V1)
filenames$V1 <- gsub("^409","386" ,filenames$V1)
filenames$V1 <- gsub("^437","406" ,filenames$V1)
filenames$V1 <- gsub("^486","409" ,filenames$V1)
filenames$V1 <- gsub("^509","486" ,filenames$V1)
filenames$V1 <- gsub("^517","509" ,filenames$V1)
filenames$V1 <- gsub("^563","517" ,filenames$V1)
filenames$V1 <- gsub("^774","787" ,filenames$V1)
allele_counts <- read.table(file = "HP_alt_136_161_217_280_371_409_427_584_787_820_859_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 34,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])
new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('ref',filenames[i,], sep='_'),paste('alt',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full_136_161_217_280_371_409_427_584_787_820_859 <- as.data.frame(t(cbind(new2)))
colnames(full_136_161_217_280_371_409_427_584_787_820_859) <- "SNP"
full_136_161_217_280_371_409_427_584_787_820_859$SNP <- as.numeric(as.character(full_136_161_217_280_371_409_427_584_787_820_859$SNP))
full_136_161_217_280_371_409_427_584_787_820_859$Sample <- rownames(full_136_161_217_280_371_409_427_584_787_820_859)
full_136_161_217_280_371_409_427_584_787_820_859 <- full_136_161_217_280_371_409_427_584_787_820_859 %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))

##alt_350_386_406_627_486
filenames <- read.table(file = "alt_350_386_406_627_486_bams.list")
filenames$V1 <- gsub("_alt.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
allele_counts <- read.table(file = "HP_alt_350_386_406_627_486_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 34,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])
new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('ref',filenames[i,], sep='_'),paste('alt',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full_350_386_406_627_486 <- as.data.frame(t(cbind(new2)))
colnames(full_350_386_406_627_486) <- "SNP"
full_350_386_406_627_486$SNP <- as.numeric(as.character(full_350_386_406_627_486$SNP))
full_350_386_406_627_486$Sample <- rownames(full_350_386_406_627_486)
full_350_386_406_627_486 <- full_350_386_406_627_486 %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))

##alt_509
filenames <- read.table(file = "alt_509_bams.list")
filenames$V1 <- gsub("_alt.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
filenames$V1 <- gsub("^161","136" ,filenames$V1)
filenames$V1 <- gsub("^208","161" ,filenames$V1)
filenames$V1 <- gsub("^228","208" ,filenames$V1)
filenames$V1 <- gsub("^280","217" ,filenames$V1)
filenames$V1 <- gsub("^358","228" ,filenames$V1)
filenames$V1 <- gsub("^386","280" ,filenames$V1)
filenames$V1 <- gsub("^406","350" ,filenames$V1)
filenames$V1 <- gsub("^409","386" ,filenames$V1)
filenames$V1 <- gsub("^437","406" ,filenames$V1)
filenames$V1 <- gsub("^486","409" ,filenames$V1)
filenames$V1 <- gsub("^509","486" ,filenames$V1)
filenames$V1 <- gsub("^517","509" ,filenames$V1)
filenames$V1 <- gsub("^563","517" ,filenames$V1)
filenames$V1 <- gsub("^774","787" ,filenames$V1)
allele_counts <- read.table(file = "HP_alt_509_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 34,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])
new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('ref',filenames[i,], sep='_'),paste('alt',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full_509 <- as.data.frame(t(cbind(new2)))
colnames(full_509) <- "SNP"
full_509$SNP <- as.numeric(as.character(full_509$SNP))
full_509$Sample <- rownames(full_509)
full_509 <- full_509 %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))

##alt_517
filenames <- read.table(file = "alt_517_bams.list")
filenames$V1 <- gsub("_alt.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
allele_counts <- read.table(file = "HP_alt_517_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 34,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])
new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('ref',filenames[i,], sep='_'),paste('alt',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full_517 <- as.data.frame(t(cbind(new2)))
colnames(full_517) <- "SNP"
full_517$SNP <- as.numeric(as.character(full_517$SNP))
full_517$Sample <- rownames(full_517)
full_517 <- full_517 %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))

##alt_822
filenames <- read.table(file = "alt_822_bams.list")
filenames$V1 <- gsub("_alt.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
filenames$V1 <- gsub("^161","136" ,filenames$V1)
filenames$V1 <- gsub("^208","161" ,filenames$V1)
filenames$V1 <- gsub("^228","208" ,filenames$V1)
filenames$V1 <- gsub("^280","217" ,filenames$V1)
filenames$V1 <- gsub("^358","228" ,filenames$V1)
filenames$V1 <- gsub("^386","280" ,filenames$V1)
filenames$V1 <- gsub("^406","350" ,filenames$V1)
filenames$V1 <- gsub("^409","386" ,filenames$V1)
filenames$V1 <- gsub("^437","406" ,filenames$V1)
filenames$V1 <- gsub("^486","409" ,filenames$V1)
filenames$V1 <- gsub("^509","486" ,filenames$V1)
filenames$V1 <- gsub("^517","509" ,filenames$V1)
filenames$V1 <- gsub("^563","517" ,filenames$V1)
filenames$V1 <- gsub("^774","787" ,filenames$V1)
allele_counts <- read.table(file = "HP_alt_822_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 34,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])
new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('ref',filenames[i,], sep='_'),paste('alt',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full_822 <- as.data.frame(t(cbind(new2)))
colnames(full_822) <- "SNP"
full_822$SNP <- as.numeric(as.character(full_822$SNP))
full_822$Sample <- rownames(full_822)
full_822 <- full_822 %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))


##alt_437, mapping bias
filenames <- read.table(file = "alt_437_bams.list")
filenames$V1 <- gsub("_alt.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
filenames$V1 <- gsub("^161","136" ,filenames$V1)
filenames$V1 <- gsub("^208","161" ,filenames$V1)
filenames$V1 <- gsub("^228","208" ,filenames$V1)
filenames$V1 <- gsub("^280","217" ,filenames$V1)
filenames$V1 <- gsub("^358","228" ,filenames$V1)
filenames$V1 <- gsub("^386","280" ,filenames$V1)
filenames$V1 <- gsub("^406","350" ,filenames$V1)
filenames$V1 <- gsub("^409","386" ,filenames$V1)
filenames$V1 <- gsub("^437","406" ,filenames$V1)
filenames$V1 <- gsub("^486","409" ,filenames$V1)
filenames$V1 <- gsub("^509","486" ,filenames$V1)
filenames$V1 <- gsub("^517","509" ,filenames$V1)
filenames$V1 <- gsub("^563","517" ,filenames$V1)
filenames$V1 <- gsub("^774","787" ,filenames$V1)
allele_counts <- read.table(file = "HP_alt_437_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 34,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])
new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('alt',filenames[i,], sep='_'),paste('ref',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full_437 <- as.data.frame(t(cbind(new2)))
colnames(full_437) <- "SNP"
full_437$SNP <- as.numeric(as.character(full_437$SNP))
full_437$Sample <- rownames(full_437)
full_437 <- full_437 %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))


plot(full[full$Allele == "ref",]$SNP,full_437[full_437$Allele == "ref",]$SNP, xlab = "D.mel reference", ylab = "DGRP 437 ref", main = "Reference allele")
abline(1,1)
plot(full[full$Allele == "alt",]$SNP,full_437[full_437$Allele == "alt",]$SNP, xlab = "D.mel reference", ylab = "DGRP 437 ref", main = "437 allele")
abline(1,1)

full$sample <- rownames(full)
full_437$sample <- rownames(full_437)
full_mapbias <- inner_join(full,full_437,by="sample")
ggplot(full_mapbias,aes(x=SNP.x,y=SNP.y,color=Source.x))+
  geom_point()+
  facet_grid(cols = vars(Allele.x))+
  xlab("D. melanogaster reference")+
  ylab("437 reference")+
  geom_abline(intercept=0,slope=1)

plot(full$SNP,full_136_161_217_280_371_409_427_584_787_820_859$SNP, xlab = "D.mel reference", ylab = "DGRP 161 ref")
abline(1,1)
plot(full$SNP,full_136_161_217_280_371_409_427_584_787_820_859$SNP, xlab = "D.mel reference", ylab = "DGRP 161 ref")
abline(1,1)
plot(full$SNP,full_350_386_406_627_486$SNP, xlab = "D.mel reference", ylab = "DGRP 386 ref")
abline(1,1)
plot(full$SNP,full_509$SNP, xlab = "D.mel reference", ylab = "DGRP 509 ref")
abline(1,1)
plot(full$SNP,full_517$SNP, xlab = "D.mel reference", ylab = "DGRP 517 ref")
abline(1,1)
plot(full$SNP,full_822$SNP, xlab = "D.mel reference", ylab = "DGRP 822 ref")
abline(1,1)


## stop codon allele - 3717487

filenames <- read.table(file = "ref_105_208_228_892_bams.list")
filenames$V1 <- gsub("_ref.*","",filenames$V1)
filenames$V1 <- gsub("_meta_5.*B","_meta5_B",filenames$V1)
filenames$V1 <- gsub("^161","136" ,filenames$V1)
filenames$V1 <- gsub("^208","161" ,filenames$V1)
filenames$V1 <- gsub("^228","208" ,filenames$V1)
filenames$V1 <- gsub("^280","217" ,filenames$V1)
filenames$V1 <- gsub("^358","228" ,filenames$V1)
filenames$V1 <- gsub("^386","280" ,filenames$V1)
filenames$V1 <- gsub("^406","350" ,filenames$V1)
filenames$V1 <- gsub("^409","386" ,filenames$V1)
filenames$V1 <- gsub("^437","406" ,filenames$V1)
filenames$V1 <- gsub("^486","409" ,filenames$V1)
filenames$V1 <- gsub("^509","486" ,filenames$V1)
filenames$V1 <- gsub("^517","509" ,filenames$V1)
filenames$V1 <- gsub("^563","517" ,filenames$V1)
filenames$V1 <- gsub("^774","787" ,filenames$V1)

allele_counts <- read.table(file = "HP_ref_105_208_228_892_lectin_ase.vcf")
allele_counts <- allele_counts[allele_counts$V2 == 419,]
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)])
fix <- as.data.frame(allele_counts[,1:8])

new2 <- ""
for (i in 1:ncol(gt)){
  new <- gt  %>% separate(colnames(gt)[i], c("GT","AD","DP","GQ","PL"), sep = ":")
  new <- new  %>% separate(AD, c(paste('ref',filenames[i,], sep='_'),paste('alt',filenames[i,], sep='_')), sep = ",")
  refname <- paste('ref',filenames[i,],sep='_')
  altname <- paste('alt',filenames[i,],sep='_')
  new2 <- cbind(new2,new[refname],new[altname])
}
new2 <- new2[,-1]
full <- as.data.frame(t(cbind(new2)))
colnames(full) <- "SNP"
full$SNP <- as.numeric(as.character(full$SNP))
full$Sample <- rownames(full)
full <- full %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep"))

full2 <- full
full2$Sample <- paste(full2$DGRP,full2$Source,full2$Allele, sep = "_")
full2$Sample2 <- paste(full2$DGRP,full2$Allele, sep = "_")

refplot <- ggplot(full2[full2$Allele == "ref",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  ggtitle("Reference allele")
altplot <- ggplot(full2[full2$Allele == "alt",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  ggtitle("Stop gained allele")

grid.arrange(refplot,altplot, ncol = 1)
