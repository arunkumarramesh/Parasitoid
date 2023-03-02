library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(dplyr)
library(plyr)
library(reshape2)
library(multcomp)

###### Process vcf generated from D. mel reference, same as DGRP-892 #####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(paste(getwd(),"/vcfs",sep=""))
filenames <- read.table(file = "ref_105_208_228_892_bams.list") ## read in DGRP sample names
### file names were incorrectly copied over from lab notebook, correcting them here
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
allele_counts <- read.table(file = "HP_ref_105_208_228_892_lectin_ase.vcf") ## read in vcf file containing allele depths
allele_counts <- allele_counts[allele_counts$V2 == 34,] ### this is the variant being studied for ASE
gt <- as.data.frame(allele_counts[,10:ncol(allele_counts)]) ## genotype data
fix <- as.data.frame(allele_counts[,1:8]) ## meta data
### process genotype matrix to get ref and alt allele depths for variant of interest
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
full <- full %>% separate(Sample, c("Allele","DGRP","Source","Primer","BioRep","TechRep")) ### separate allele depths into different categories
full2 <- full
full2$Sample <- paste(full2$DGRP,full2$Source,full2$Allele, sep = "_") ### sample name with cdna or gdna source
full2$Sample2 <- paste(full2$DGRP,full2$Allele, sep = "_") ### sample name with just line and allele

############## ref (892) and alt (437) allele counts plot #############
## first only plot for 892 and 437
isolate <- full2[full2$DGRP %in% "892",]
isolate$Allele <- gsub("ref","892",isolate$Allele)
isolate$Allele <- gsub("alt","437",isolate$Allele)
isolate2 <- isolate %>% 
  dplyr::group_by(Allele,Source)%>% 
  dplyr::summarise(Depth=mean(SNP),ci=(sd(SNP)/sqrt(3))*qt(p=0.05/2, df=2,lower.tail=F))

t.test(isolate[isolate$Source %in% "gDNA",]$SNP~isolate[isolate$Source %in% "gDNA",]$Allele)
t.test(isolate[isolate$Source %in% "cDNA",]$SNP~isolate[isolate$Source %in% "cDNA",]$Allele)

ann_text <- data.frame(x = c(1.5,1.51),y = 540,lab = c("****",""), Source = c("cDNA","gDNA"),Allele = "437")
ann_text2 <- data.frame(x = c(1.5,NA),y = 520,lab = c("****",""), Source = c("cDNA","gDNA"),Allele = "437")

pdf(file="../plots/437_892_ase.pdf",height=3,width = 2.5)
ggplot(isolate2,aes(x=Allele, y=Depth,fill=Allele) )+
  geom_col()+
  scale_fill_manual(values = c("red","blue"))+
  geom_errorbar(aes(ymin=Depth-ci, ymax=Depth+ci), width=0.3)+ 
  geom_jitter(data=isolate,aes(x=Allele,y=SNP),width=0.4,height=0,size=1, alpha=0.5)+
  facet_wrap(~Source)+ 
  geom_text(data=ann_text,aes(x=x,y=y,label=lab))+
  geom_segment(data=ann_text2,aes(x=x-0.4,xend=x+0.4,y=y,yend=y))+
  theme(legend.position = "none",panel.border = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 13)) +
  theme(strip.text = element_text(size = 13))
dev.off()

## now plotting for all DGRPs

refplot <- ggplot(full2[full2$Allele == "ref",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  ggtitle("Reference allele")+
  ylab("Allele depth")
altplot <- ggplot(full2[full2$Allele == "alt",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  ggtitle("437 (alternate) allele")+
  ylab("Allele depth")
#grid.arrange(refplot,altplot, ncol = 1)
snpdepthplot <- ggplot(full2[full2$Allele == "ref",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  ylab("Allele depth")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

################ Group DGRPs by insertion haplotype ###################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ase_indels <- read.csv(file = "ase_indels.csv", stringsAsFactors = F) ### document containing Sanger seq calls upstream indels in DGRPs
colnames(ase_indels)[1] <- "line"
ase_indels$X7bp_ins_1 <- gsub("^A$","D",ase_indels$X7bp_ins_1)
ase_indels$X7bp_ins_1 <- gsub("^P$","I",ase_indels$X7bp_ins_1)
ase_indels$X21bp_ins_3 <- gsub("^A$","D",ase_indels$X21bp_ins_3)
ase_indels$X21bp_ins_3 <- gsub("^P$","I",ase_indels$X21bp_ins_3)
ase_indels$X8bp_ins_2 <- gsub("^P$","I",ase_indels$X8bp_ins_2)
ase_indels$X8bp_ins_2 <- gsub("^A$","D",ase_indels$X8bp_ins_2)
ase_indels2 <- ase_indels
ase_indels$indel_hap <- paste(ase_indels$X7bp_ins_1,ase_indels$X8bp_ins_2,ase_indels$X21bp_ins_3) ### indel haplotype for DGRP lines
ase_indels <- as.data.frame(cbind(ase_indels$line,ase_indels$indel_hap))
colnames(ase_indels)<- c("DGRP","indel_hap")


full2_forplot <- full2
full2_forplot <- left_join(full2_forplot,ase_indels,by="DGRP")
pdf("plots/ase_ins_hap.pdf", height = 3, width = 7)
ggplot(full2_forplot[full2_forplot$Allele == "ref",],aes(x=DGRP, y=SNP, color = Source) )+
  geom_boxplot()+
  facet_grid(~indel_hap, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x") +  # Let the width of facets vary and force all bars to have the same width.

  ylab("Allele depth")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

## stats to check if DDI is different from DDD

forstats <- full2_forplot[full2_forplot$Source %in% "cDNA" & full2_forplot$Allele %in% "ref",]
colnames(forstats)[1:2] <- c("ref","alt")
forstats$alt <- full2_forplot[full2_forplot$Source %in% "cDNA" & full2_forplot$Allele %in% "alt",1]
model_ase <- glm(data=forstats,cbind(ref,alt)~indel_hap,family=quasibinomial())

## check for times differneces between allele depths for haplotypes
forstats2 <- forstats
forstats2$Prop <- forstats2$ref/(forstats2$ref+forstats2$alt)
  
forstats2 %>%
  dplyr::group_by(indel_hap) %>%
  dplyr::summarise(depth=mean(Prop))

## to plot proportions instead

propplot <- full2_forplot[full2_forplot$Allele %in% "ref",]
propplot$Prop <- propplot$SNP/(propplot$SNP + full2_forplot[full2_forplot$Allele %in% "alt",]$SNP)

pdf("plots/ase_ins_hap_prop.pdf", height = 3.5, width = 5)
ggplot(propplot,aes(x=DGRP, y=Prop, color = Source) )+
  geom_boxplot()+
  facet_grid(~indel_hap, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x") +  # Let the width of facets vary and force all bars to have the same width.
  
  ylab("Proportion of non DGRP-437 allele")+
  xlab("Inbred Line (DGRP)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = c(0.85, 0.4),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=13)) +
  theme(axis.title = element_text(size = 14)) +
  theme(legend.title=element_text(size=13)) +
  theme(legend.text=element_text(size=13)) +
  theme(strip.text.x = element_text(size = 14))
dev.off()

################ Group DGRPs by indel haplotype ###################
##convert indel abs/pres into derived and ancestral states
##ase_indels.csv - Sanger seq calls of indels for DGRPs
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ase_indels <- read.csv(file = "ase_indels.csv", stringsAsFactors = F) ### document containing Sanger seq calls upstream indels in DGRPs
colnames(ase_indels)[1] <- "line"
### ancestral states were identified through alignments to multiple closely related drosophila species
ase_indels$X7bp_ins_1 <- gsub("^A$","D",ase_indels$X7bp_ins_1)
ase_indels$X7bp_ins_1 <- gsub("^P$","A",ase_indels$X7bp_ins_1)
ase_indels$X21bp_ins_3 <- gsub("^A$","D",ase_indels$X21bp_ins_3)
ase_indels$X21bp_ins_3 <- gsub("^P$","A",ase_indels$X21bp_ins_3)
ase_indels$X8bp_ins_2 <- gsub("^P$","D",ase_indels$X8bp_ins_2)
ase_indels2 <- ase_indels
ase_indels$indel_hap <- paste(ase_indels$X7bp_ins_1,ase_indels$X8bp_ins_2,ase_indels$X21bp_ins_3) ### indel haplotype for DGRP lines
ase_indels <- as.data.frame(cbind(ase_indels$line,ase_indels$indel_hap))
colnames(ase_indels)<- c("DGRP","indel_hap")

##genotype accounting for SNP in 21 BP indel
for (r in 1:nrow(ase_indels2)){
  if (as.character(ase_indels2$snp.at.pos.3.in.21bp.indel[r]) %in% "G"){
    ase_indels2$X21bp_ins_3[r] <- "(G)"
  }
  if (as.character(ase_indels2$snp.at.pos.3.in.21bp.indel[r]) %in% "A"){
    ase_indels2$X21bp_ins_3[r] <- "(A)"
  }
}
ase_indels2$indel_hap <- paste(ase_indels2$X7bp_ins_1,ase_indels2$X8bp_ins_2,ase_indels2$X21bp_ins_3) ### indel haplotype for DGRP lines with snp 21 BP indel
ase_indels2 <- as.data.frame(cbind(ase_indels2$line,ase_indels2$indel_hap))
colnames(ase_indels2)<- c("DGRP","indel_hap")

##frequency of reference allele by indel haplotype
full_ref <- full[full$Allele == "ref",]
full_alt <- full[full$Allele == "alt",]
full_af <- full_alt
full_af$AF <- full_ref$SNP/(full_alt$SNP+full_ref$SNP) ### frequency of ref (non 437) allele 
full_af <- left_join(full_af,ase_indels, by = "DGRP")
afplot <- ggplot(full_af,aes(x=indel_hap, y=AF, color = Source) )+
  geom_boxplot()+
  ylab("Allele Frequency")+
  xlab("Indel haplotype")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


##frequency of reference allele by indel haplotype, accounting for 21 BP SNP
full_ref <- full[full$Allele == "ref",]
full_alt <- full[full$Allele == "alt",]
full_af <- full_alt
full_af$AF <- full_ref$SNP/(full_alt$SNP+full_ref$SNP) ### frequency of ref (non 437) allele 
full_af2 <- left_join(full_af,ase_indels2, by = "DGRP")
afplot2 <- ggplot(full_af2,aes(x=indel_hap, y=AF, color = Source) )+
  geom_boxplot()+
  ylab("Allele Frequency")+
  xlab("Indel haplotype")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

############# Comparison of Technical replicates ############# 
## there were two technical replicates
## split data from two technical replicates
full_af_T1 <- full_af[full_af$TechRep == "T1",]
full_af_T2 <- full_af[full_af$TechRep == "T2",]
full_af_T1$sample <- paste(full_af_T1$DGRP,full_af_T1$Source,full_af_T1$BioRep, sep = "_")
full_af_T2$sample <- paste(full_af_T2$DGRP,full_af_T2$Source,full_af_T2$BioRep, sep = "_")
## combine technical replicates by sample
full_af_T <- inner_join(full_af_T1,full_af_T2,by="sample")
plottechrep <- ggplot(full_af_T,aes(x=AF.x,y=AF.y,color=Source.x))+
  geom_point()+
  #facet_grid(cols = vars(Allele.x))+
  xlab("AF for Technical replicate 1")+
  ylab("AF for Technical replicate 2")+
  geom_abline(intercept=0,slope=1)+ 
  labs(color = "Source")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title = element_blank()) +
  theme(
    legend.position = c(.4, .9),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

############# Comparison of Biological replicates #################
## there were two biological replicates for cDNA
## split data from two technical and biological replicates
full_B1_T1 <- full[full$BioRep == "B1" & full$TechRep == "T1"& full$Source == "cDNA",]
full_B1_T2 <- full[full$BioRep == "B1" & full$TechRep == "T2"& full$Source == "cDNA",]
full_B2_T1 <- full[full$BioRep == "B2" & full$TechRep == "T1"& full$Source == "cDNA",]
full_B2_T2 <- full[full$BioRep == "B2" & full$TechRep == "T2"& full$Source == "cDNA",]
### combine data for two technical replicates for a given biological replicate, separatly
full_B1 <- full_B1_T1
full_B1$SNP <- full_B1_T1$SNP + full_B1_T2$SNP
full_B2 <- full_B2_T1
full_B2$SNP <- full_B2_T1$SNP + full_B2_T2$SNP
full_B1_ref <- full_B1[full_B1$Allele == "ref",]
full_B1_alt <- full_B1[full_B1$Allele == "alt",]
full_B1_af <- full_B1_alt
full_B1_af$AF <- full_B1_ref$SNP/(full_B1_alt$SNP+full_B1_ref$SNP) ### frequency of ref (non 437) allele for biological replicate 1
full_B1_af$sample <- paste(full_B1_af$DGRP,full_B1_af$Source, sep = "_")
full_B2_ref <- full_B2[full_B2$Allele == "ref",]
full_B2_alt <- full_B2[full_B2$Allele == "alt",]
full_B2_af <- full_B2_alt
full_B2_af$AF <- full_B2_ref$SNP/(full_B2_alt$SNP+full_B2_ref$SNP) ### frequency of ref (non 437) allelefor biological replicate 2
full_B2_af$sample <- paste(full_B2_af$DGRP,full_B2_af$Source, sep = "_")
full_af_B <- inner_join(full_B1_af,full_B2_af,by="sample")
plotbiorep <- ggplot(full_af_B,aes(x=AF.x,y=AF.y,color=Source.x))+
  geom_point()+
  #facet_grid(cols = vars(Allele.x))+
  xlab("AF for Biological replicate 1")+
  ylab("AF for Biological replicate 2")+
  geom_abline(intercept=0,slope=1)+ 
  labs(color = "Source")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(
    legend.position = c(.4, .9),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

plottechbiorep <- plot_grid(plottechrep,plotbiorep,labels = "AUTO")

############# Total expression variation #################
## this section compares frequency non-437  allele in cDNA and gDNA
### first split biological and technical replicates and gDNA
full_B1_T1 <- full[full$BioRep == "B1" & full$TechRep == "T1" & full$Source == "cDNA",]
full_B1_T2 <- full[full$BioRep == "B1" & full$TechRep == "T2"& full$Source == "cDNA",]
full_B2_T1 <- full[full$BioRep == "B2" & full$TechRep == "T1"& full$Source == "cDNA",]
full_B2_T2 <- full[full$BioRep == "B2" & full$TechRep == "T2"& full$Source == "cDNA",]
full_T1_gdna <- full[full$BioRep == "B1" & full$TechRep == "T1" & full$Source == "gDNA",]
full_T2_gdna <- full[full$BioRep == "B1" & full$TechRep == "T2"& full$Source == "gDNA",]
## now combine data for biological replicates
full_allsum <- full_B1_T1
full_allsum$SNP <- full_B1_T1$SNP + full_B1_T2$SNP + full_B2_T1$SNP + full_B2_T2$SNP
## now combine gDNA data from technical replicates
full_allsum_gdna <- full_T1_gdna
full_allsum_gdna$SNP <- full_T1_gdna$SNP + full_T2_gdna$SNP
full_allsum <- rbind(full_allsum,full_allsum_gdna)
full_allsum_ref <- full_allsum[full_allsum$Allele == "ref",]
full_allsum_alt <- full_allsum[full_allsum$Allele == "alt",]
full_allsum_af <- full_allsum_alt
full_allsum_af$AF <- full_allsum_ref$SNP/(full_allsum_alt$SNP+full_allsum_ref$SNP) ### calculate frequency of non-437 allele
## plot of non-437 af in cDNA and gDNA. gDNA at 0.5 as expected, and clear bimodality for cDNA, indicating either an allele is expressed or it isn't
ggplot(full_allsum_af,aes(x=Source, y=AF) )+
  geom_point()+
  ggtitle("Total expression variation")

########### Indel haplotype by expression #############
##haplotypes accounting for SNP in 21 BP indel
full_allsum_af2 <- left_join(full_allsum_af,ase_indels2, by = "DGRP")
full_allsum_af2 <- full_allsum_af2[full_allsum_af2$Source == "cDNA",]
## given bimodality of allele expression, samples were classified as either having induced lectin following infection and those that did not induce the gene
full_allsum_af2$Induced <- ""
for (r in 1:nrow(full_allsum_af2)){
  if (full_allsum_af2$AF[r] > 0.3){
    full_allsum_af2$Induced[r] <- "Induced"
  }
  if (full_allsum_af2$AF[r] < 0.3){
    full_allsum_af2$Induced[r] <- "Not Induced"
  }
}
### plot of indel haplotype and induction of lectin-24a accounting for snp in 21bp indel
ggplot(full_allsum_af2, aes(x=indel_hap, fill = Induced))+
  geom_histogram(stat="count")+
  xlab("Indel haplotype")+
  ylab("Number of lines")+
  theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


##haplotypes for three indels only
full_allsum_af2 <- left_join(full_allsum_af,ase_indels, by = "DGRP")
full_allsum_af2 <- full_allsum_af2[full_allsum_af2$Source == "cDNA",]
## given bimodality of allele expression, samples were classified as either having induced lectin following infection and those that did not induce the gene
full_allsum_af2$Induced <- ""
for (r in 1:nrow(full_allsum_af2)){
  if (full_allsum_af2$AF[r] > 0.3){
    full_allsum_af2$Induced[r] <- "Induced"
  }
  if (full_allsum_af2$AF[r] < 0.3){
    full_allsum_af2$Induced[r] <- "Not Induced"
  }
}
### plot of indel haplotype and induction of lectin-24a
indelhaplotypeplot <- ggplot(full_allsum_af2, aes(x=indel_hap, fill = Induced))+
  geom_histogram(stat="count")+
  xlab("Indel haplotype")+
  ylab("Number of lines")+
  theme(legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("plots/indel haplotype expression.pdf",height = 1.8, width = 4)
indelhaplotypeplot
dev.off()

##chi square test for association on indel haplotype with capacity to induce lectin-24A expression
chisq.test(table(full_allsum_af2$indel_hap,full_allsum_af2$Induced))

### Comparison of ASE estimates (DGRP-Test x DGRP-437) and qPCR estimates for DGRP-Test lines #####
### read in upstream lectin-24A indel calls for DGRPs and convert to ancestral and derived calls
qpcr_dgrp_indels <- read.csv(file= "qpcr_dgrp_indels.csv", header = T)
qpcr_dgrp_indels$X7bp_ins_1 <- gsub("A","D",qpcr_dgrp_indels$X7bp_ins_1)
qpcr_dgrp_indels$X7bp_ins_1 <- gsub("P","A",qpcr_dgrp_indels$X7bp_ins_1)
qpcr_dgrp_indels$X8bp_ins_2 <- gsub("P","D",qpcr_dgrp_indels$X8bp_ins_2)
qpcr_dgrp_indels$X21bp_ins_3 <- gsub("A","D",qpcr_dgrp_indels$X21bp_ins_3)
qpcr_dgrp_indels$X21bp_ins_3 <- gsub("P","A",qpcr_dgrp_indels$X21bp_ins_3)
## accounting for snp in 21 bp indel
qpcr_dgrp_indels$snp.at.pos.3.in.21bp.indel <- gsub("A","(A)",qpcr_dgrp_indels$snp.at.pos.3.in.21bp.indel)
qpcr_dgrp_indels$snp.at.pos.3.in.21bp.indel <- gsub("G","(G)",qpcr_dgrp_indels$snp.at.pos.3.in.21bp.indel)
qpcr_dgrp_indels$X21bp_ins_snp <- paste(qpcr_dgrp_indels$X21bp_ins_3,qpcr_dgrp_indels$snp.at.pos.3.in.21bp.indel,sep = "")
qpcr_dgrp_indels$X21bp_ins_snp <- gsub("DNA","D",qpcr_dgrp_indels$X21bp_ins_snp)
qpcr_dgrp_indels$haplotype <- paste(qpcr_dgrp_indels$X7bp_ins_1,qpcr_dgrp_indels$X8bp_ins_2,qpcr_dgrp_indels$X21bp_ins_3) ### indel haplotype only
qpcr_dgrp_indels$haplotype2 <- paste(qpcr_dgrp_indels$X7bp_ins_1,qpcr_dgrp_indels$X8bp_ins_2,qpcr_dgrp_indels$X21bp_ins_snp)  ### indel haplotype with snp in 21 bp indel
qpcr_dgrp_indels <- cbind(qpcr_dgrp_indels[1],qpcr_dgrp_indels[7:8])
qpcr <- read.csv(file= "qpcr_dgrp.csv", header = T) ## read in qPCR data
qpcr_mod <- as.data.frame(cbind(qpcr[qpcr$Sample %in% "Treatment",1],qpcr[qpcr$Sample %in% "Treatment",]$delta_ct-qpcr[qpcr$Sample %in% "Control",]$delta_ct,qpcr[qpcr$Sample %in% "Control",]$delta_ct)) ## get relevant columns, delta ct for control and infected samples and log2fc following infection
colnames(qpcr_mod) <- c("line","log2foldchange","deltaCTlectin")
qpcr_indels <- left_join(qpcr_dgrp_indels, qpcr_mod, by = "line") ## combine ASE and qPCR data
qpcr_indels$ASE <- mapvalues(qpcr_indels$haplotype, from = c("A A A","A D A","D D A","D A D","D A A"),to= c("Induced","Induced","Induced","Not Induced","Not Induced")) ## convert haplotypes into induced or not induced using results from previous section
#model <- lm(data = qpcr_indels, deltaCTlectin  ~ factor(haplotype))
#summary(glht(model, linfct = mcp("factor(haplotype)" = "Tukey")), test = adjusted("holm"))
### compare logFC for induced and non-induced haplotypes as designated by ASE
qpcrfc <- ggplot(data = qpcr_indels, aes(x=haplotype,y=log2foldchange))+
  geom_violin()+
  geom_point(aes(colour = ASE))+
  xlab("Indel haplotype")+
  ylab(expression(paste("Log(2)FC ",italic('Lectin-24A'))))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
### compare delta CT in control samples, looking at basal expression in DGRP lines, for induced and non-induced haplotypes as designated by ASE
deltactlectin <- ggplot(data = qpcr_indels, aes(x=haplotype,y=-deltaCTlectin))+
  geom_violin()+
  geom_point(aes(colour = ASE))+
  xlab("Indel haplotype")+
  ylab(expression(paste("No infection ",Delta,'CT ', italic('Lectin-24A'))))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(breaks=seq(2, 13, 2),limits = c(2,13))
### compare delta CT in infected samples, looking at induced expression in DGRP lines, for induced and non-induced haplotypes as designated by ASE
qpcr_mod2 <- as.data.frame(cbind(qpcr[qpcr$Sample %in% "Treatment",1],qpcr[qpcr$Sample %in% "Treatment",]$delta_ct-qpcr[qpcr$Sample %in% "Control",]$delta_ct,qpcr[qpcr$Sample %in% "Treatment",]$delta_ct))
colnames(qpcr_mod2) <- c("line","log2foldchange","deltaCTlectin")
qpcr_indels2 <- left_join(qpcr_dgrp_indels, qpcr_mod2, by = "line")
qpcr_indels2$ASE <- mapvalues(qpcr_indels2$haplotype, from = c("A A A","A D A","D D A","D A D","D A A"),to= c("Induced","Induced","Induced","Not Induced","Not Induced"))
deltactlectin2 <- ggplot(data = qpcr_indels2, aes(x=haplotype,y=-deltaCTlectin))+
  geom_violin()+
  geom_point(aes(colour = ASE))+
  xlab("Indel haplotype")+
  ylab(expression(paste(Delta,'CT Lectin-24A')))+
  ggtitle("Infected samples")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  scale_y_continuous(breaks=seq(2, 13, 2),limits = c(2,13))

## a different value to plot delta CT in infected and no infection samples for lectin-24A and no reverse transciptase controls
qpcr_melt <- qpcr
qpcr_melt$delta_ct_lectin <-  qpcr_melt$ct_rpl - qpcr_melt$ct_lectin ## calculate delta ct by subtracting ct of lectin from rpl
qpcr_melt$delta_ct_lectin_nort <- qpcr_melt$ct_lectin_no_rt - qpcr_melt$ct_rpl ## calculate delta ct by subtracting ct of no rt control from rpl
summary(aov(data = qpcr_melt, delta_ct_lectin_nort~line+Sample)) ## test to see any difference in no RT controls by line or infection, none detected

qpcr_melt <- qpcr_melt[-c(3:6)]
qpcr_melt <- melt(qpcr_melt,id.vars = c("line","Sample"))
qpcr_melt <- qpcr_melt[order(-qpcr_melt$line),]
colnames(qpcr_melt)[3:4] <- c("Gene","CT")
qpcr_melt$Gene <- mapvalues(qpcr_melt$Gene,from=c("delta_ct_lectin","delta_ct_lectin_nort"),to=c("Lectin-24A","NRT"))
qpcr_melt$Sample <-  mapvalues(qpcr_melt$Sample,from=c("Control","Treatment"),to=c("No infection","Infection"))
qpcr_deltactdgrp <- ggplot(qpcr_melt,aes(x=-CT,y=as.factor(line),color=Gene))+
  geom_point()+
  facet_grid(.~Sample)+
  xlab(expression(paste(Delta,"CT")))+
  ylab("DGRP line")+
  theme(legend.title = element_blank())+
  scale_colour_manual(values = c("orange","blue"))

write.csv(qpcr_melt[(qpcr_melt$Sample %in% "No infection") & (qpcr_melt$Gene %in% "Lectin-24A"), ],file="no_infection_dgrp_lectin.csv")  ## write out data for no infection lectin-24a to see baseline expression

### plot of DGRP-437 and DGRP-892 fold change following infection, highlighting the different for key lines that were previously summerised by haplotype
ressus_qpcr <- ggplot(data = qpcr_indels2[qpcr_indels2$line %in% c(437,892),],aes(x=as.factor(line),y=log2foldchange))+
  geom_bar(stat = "identity") +
  xlab("DGRP line")+
  ylab(expression(paste("Log(2)FC ",italic('Lectin-24A'))))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black"))

dgrppqcr <- plot_grid(qpcr_deltactdgrp,ressus_qpcr,rel_widths = c(3,1),labels = "AUTO",ncol=2)
qpcr_ase <- plot_grid(deltactlectin,qpcrfc,ncol=2,labels = c("C","D"),rel_widths = c(1.4,1))
qpcr_ase2 <- plot_grid(dgrppqcr,qpcr_ase,ncol=1,rel_heights = c(1.1,1))

qpcr_model <- qpcr_indels
qpcr_model$haplotype <- as.factor(qpcr_model$haplotype)
### analysis of variance to test effect of haplotype on delta ct  and log2FC for lectin-24A following infection
model <- aov(data = qpcr_model, deltaCTlectin~haplotype)
model2 <- aov(data = qpcr_model, log2foldchange~haplotype)

##new plot to compare qPCR induction data with ASE data
qpcr_melt_new <- qpcr
qpcr_melt_new$CT <- (qpcr_melt_new$ct_rpl - qpcr_melt_new$ct_lectin) - ( qpcr_melt_new$ct_rpl - qpcr_melt_new$ct_lectin_no_rt) ## this metric subtracts the delta CT lectin from delta CT for NRT 
qpcr_dgrp_indels <- read.csv(file= "qpcr_dgrp_indels.csv", header = T)
qpcr_dgrp_indels$X7bp_ins_1 <- gsub("A","D",qpcr_dgrp_indels$X7bp_ins_1)
qpcr_dgrp_indels$X7bp_ins_1 <- gsub("P","I",qpcr_dgrp_indels$X7bp_ins_1)
qpcr_dgrp_indels$X8bp_ins_2 <- gsub("P","I",qpcr_dgrp_indels$X8bp_ins_2)
qpcr_dgrp_indels$X8bp_ins_2 <- gsub("A","D",qpcr_dgrp_indels$X8bp_ins_2)
qpcr_dgrp_indels$X21bp_ins_3 <- gsub("A","D",qpcr_dgrp_indels$X21bp_ins_3)
qpcr_dgrp_indels$X21bp_ins_3 <- gsub("P","I",qpcr_dgrp_indels$X21bp_ins_3)
## accounting for snp in 21 bp indel
qpcr_dgrp_indels$haplotype <- paste(qpcr_dgrp_indels$X7bp_ins_1,qpcr_dgrp_indels$X8bp_ins_2,qpcr_dgrp_indels$X21bp_ins_3) ### indel haplotype only
qpcr_melt_new <- left_join(qpcr_melt_new,qpcr_dgrp_indels,by="line")

## test significance of CT values
model <- lm(data=qpcr_melt_new,CT ~ as.factor(haplotype))
tukey.test <- TukeyHSD(aov(model),  conf.level=0.95)
plot(tukey.test , las=1 , col="brown") # Tuckey test representation

labeldat <- data.frame(haplotype=names(table(qpcr_melt_new$haplotype)),x=1,xend=c(4,3,2,4,5),y=c(5,8,14,14,14),yend=c(5,8,14,14,14) )
labeldat2 <- data.frame(haplotype=names(table(qpcr_melt_new$haplotype)),x=c(2.5,2,1.5,2.5,3),y=c(6,9,15,15,15),lab=c("a","bc","c","c","c") )

qpcr_melt_new %>%
  dplyr::group_by(haplotype) %>%
  dplyr::summarise(CT = mean(CT))


pdf("plots/qpcr_ins_hap.pdf", height = 3, width = 7)
ggplot(qpcr_melt_new,aes(x=as.factor(line), y=CT) )+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(~haplotype, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x") +  # Let the width of facets vary and force all bars to have the same width.
  
  ylab(expression('Log'[2]*' relative '*italic('Lectin-24A')*' expression'))+
  xlab("DGRP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_segment(data=labeldat,aes(x=x,xend=xend,y=y,yend=yend))+
  geom_text(data=labeldat2,aes(x=x,y=y,label=lab))
dev.off()


pdf("plots/qpcr_ins_hap_ase_lines.pdf", height = 3, width = 7)
ggplot(qpcr_melt_new[qpcr_melt_new$line %in% full2_forplot$DGRP,],aes(x=as.factor(line), y=CT) )+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(~haplotype, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x") +  # Let the width of facets vary and force all bars to have the same width.
  
  ylab(expression('Log'[2]*' relative '*italic('Lectin-24A')*' expression'))+
  xlab("DGRP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

### check which contrasts are significantly different
summary(model)
summary(model2)
summary(glht(model, linfct = mcp(haplotype = "Tukey")), test = adjusted("holm"))


### analysis of variance to test effect of induced and non-induced haplotypes on delta ct  and log2FC for lectin-24A following infection
t.test(data = qpcr_model, deltaCTlectin~ASE)
t.test(data = qpcr_model, log2foldchange~ASE)


################# Lectin upstream SNPs and ASE status ###################
## this part checks if any snp in 1KB upstream region of lectin-24a is fully associated with induction as noted by ASE
dgrp_snps <- read.csv(file = "lectin_promoter_and_coding.csv") ## file with DGRP Freeze 2.0 calls for upstream SNPs for lectin-24a
## get upstream 1 kb region
dgrp_snps <- dgrp_snps[dgrp_snps$pos > 3717772,]
dgrp_snps <- dgrp_snps[dgrp_snps$pos < 3718772,]
dgrp_snps <- dgrp_snps[!dgrp_snps$pos == 3718040,] ### remove 8bp indel called in DGRP freeze 2.0
rownames(dgrp_snps) <- dgrp_snps$id
dgrp_snps <- dgrp_snps[-c(1:9)]
dgrp_snps <- as.data.frame(t(dgrp_snps))
rownames(dgrp_snps) <- gsub("line_","",rownames(dgrp_snps))
dgrp_snps$DGRP <- rownames(dgrp_snps)
full_allsum_af2 <- inner_join(full_allsum_af2,dgrp_snps,by = "DGRP") ### combine snp data with ase induction data
full_allsum_af2$indel_hap <- as.character(full_allsum_af2$indel_hap)
### convert haplotype to individual indel calls
full_allsum_af3 <- separate(full_allsum_af2, col = "indel_hap", into = c("7 BP Indel","8 BP Indel","21 BP Indel"), sep = " ", remove = F)
full_allsum_af3$`7 BP Indel` <- gsub("^A$","P",full_allsum_af3$`7 BP Indel`)
full_allsum_af3$`7 BP Indel` <- gsub("^D$","A",full_allsum_af3$`7 BP Indel`)
full_allsum_af3$`21 BP Indel` <- gsub("^A$","P",full_allsum_af3$`21 BP Indel`)
full_allsum_af3$`21 BP Indel` <- gsub("^D$","A",full_allsum_af3$`21 BP Indel`)
full_allsum_af3$`8 BP Indel` <- gsub("^D$","P",full_allsum_af3$`8 BP Indel`)
full_allsum_af3 <- full_allsum_af3[full_allsum_af3$`21 BP Indel` %in% "P",] ## only keep lines where 21 BP is present, since we already know it controls expression
full_allsum_af4 <- cbind(full_allsum_af3[3],full_allsum_af3[12],full_allsum_af3[13:ncol(full_allsum_af3)],full_allsum_af3[11],full_allsum_af3[10])
full_allsum_af4 <- melt(full_allsum_af4, id.vars = c("DGRP","Induced"))
colnames(full_allsum_af4)[3] <- "Position"
full_allsum_af4$Position <- gsub("2L_","",full_allsum_af4$Position)
full_allsum_af4 <- full_allsum_af4[!full_allsum_af4$value == "-",]
### convert calls to reference and alternate for upstream snps and indels
full_allsum_af4$value <- gsub("0","Reference allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("2","Alternate allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("^A$","Reference allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("^P$","Alternate allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("A//)$","Reference allele",full_allsum_af4$value)
full_allsum_af4$value <- gsub("G//)$","Alternate allele",full_allsum_af4$value)
## plot of variant association with lectin-24a induction
variant_ase <- ggplot(data = full_allsum_af4) + 
  geom_bar(mapping = aes(x = factor(Position, levels = unique(Position)), fill = Induced), position = "stack")+
  facet_grid(rows = vars(value))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Number of lines")+
  xlab("Variant")+
  ylab("Number of lines")+
  theme(legend.title = element_blank())+theme(
    legend.position = c(.65, 1),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))


pdf("plots/snps_exp.pdf",height=4,width=7)
variant_ase
dev.off()
################# mapping to alternate references ###################
## point of this section is to test for mapping biases, using different references to assess lectin-24A induction. Repeating many of the steps above
## references were generated by substituting snps for alternate alleles prior to mapping and each section, starting with "alt" refers to a different reference containing alleles indicate in DGRP lines specified

##alt_136_161_217_280_371_409_427_584_787_820_859
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(paste(getwd(),"/vcfs",sep=""))
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
full$sample <- rownames(full)
full_437$sample <- rownames(full_437)
full_mapbias <- inner_join(full,full_437,by="sample")
full_mapbias$Allele.x <- gsub("ref","Unmodified reference allele",full_mapbias$Allele.x)
full_mapbias$Allele.x <- gsub("alt","Unmodified alternate allele",full_mapbias$Allele.x)
plotmappingbias <- ggplot(full_mapbias,aes(x=SNP.x,y=SNP.y,color=Source.x))+
  geom_point()+
  facet_grid(cols = vars(Allele.x))+
  xlab("Allele depth from unmodified reference")+
  ylab("Allele depth from modified reference")+
  geom_abline(intercept=0,slope=1)+ 
  labs(color = "Source")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


########## Combine Supplementary Plots ############
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
pdf("plots/ase_supp.pdf", height = 7, width = 7)
plot_grid(plottechbiorep,plotmappingbias,ncol=1,labels = c("","C"),rel_heights = c(1,0.9))
dev.off()
