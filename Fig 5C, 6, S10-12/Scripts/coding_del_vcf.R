## the only purpuse of the script is to create a vcf file for the 165 BP lectin-24a coding deletion. I used the dgn_lectin.R script and remove steps involving sample filtering. Many of the comments refer to the other script and are not related to this one. 
##coding deletion calls were made by checking for consistent NA values for SNPs that occur in the insertion wild type sequence. If a sample has NA for all 9 snps that occur in this region but is called for other surrounding SNPs, then it is designated as having the coding deletion

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyr)
library(dplyr)
library(stringr)
library(adegenet)
library(plyr)
library(hierfstat)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(car)
library(lme4)
library(data.table)
library(sjmisc)
library(phylotools)
library(ggpubr)
library(ape)
library(Biostrings)
library(ggtree)
library(phangorn)
library(cowplot)
library(colorspace)
library(RColorBrewer)
library(seqinr)
library(grid)
library(LDheatmap)
library(genetics)
library(varhandle)
library(VariantAnnotation)
library(snpStats)
library(rehh)
library(gggenes)
library(pBrackets)
library(pheatmap)
library(ggplotify)

### Aim is to study molecular evolution of lectin-24A in DPGP samples

#### 1. Grouping DPGP populations by region  ####

B_POP <- as.character(read.table("../Input/B_POP")[[1]])
CO_POP <- as.character(read.table("../Input/CO_POP")[[1]])	 
EA_POP <- as.character(read.table("../Input/EA_POP")[[1]])	 
EB_POP <- as.character(read.table("../Input/EB_POP")[[1]])	 
ED_POP <- as.character(read.table("../Input/ED_POP")[[1]])	 
EF_POP <- as.character(read.table("../Input/EF_POP")[[1]])	 
EG_POP <- as.character(read.table("../Input/EG_POP")[[1]])	 
ER_POP <- as.character(read.table("../Input/ER_POP")[[1]])	 
FR_POP <- as.character(read.table("../Input/FR_POP")[[1]])	 
GA_POP <- as.character(read.table("../Input/GA_POP")[[1]])	 
GU_POP <- as.character(read.table("../Input/GU_POP")[[1]])	 
I_POP <- as.character(read.table("../Input/I_POP")[[1]])	 
N_POP <- as.character(read.table("../Input/N_POP")[[1]])	 
NG_POP <- as.character(read.table("../Input/NG_POP")[[1]])  
RAL_POP <- as.character(read.table("../Input/RAL_POP")[[1]])		 
RG_POP <- as.character(read.table("../Input/RG_POP")[[1]])	 
SB_POP <- as.character(read.table("../Input/SB_POP")[[1]])	 
SD_POP <- as.character(read.table("../Input/SD_POP")[[1]])	 
SF_POP <- as.character(read.table("../Input/SF_POP")[[1]])	 
SP_POP <- as.character(read.table("../Input/SP_POP")[[1]])	 
T_POP <- as.character(read.table("../Input/T_POP")[[1]])	 
UG_POP <- as.character(read.table("../Input/UG_POP")[[1]])	 
W_POP <- as.character(read.table("../Input/W_POP")[[1]])	 
ZI_POP <- as.character(read.table("../Input/ZI_POP")[[1]])	 
ZS_POP <- as.character(read.table("../Input/ZS_POP")[[1]])	 
ZW_POP <- as.character(read.table("../Input/ZW_POP")[[1]])


### grouping populations by broader regions
ocenia <- T_POP
asia <- B_POP
america <- c(I_POP,RAL_POP,W_POP)
europe_north_africa <- c(N_POP,FR_POP,EG_POP)
southern_africa <- c(SB_POP,SD_POP,SF_POP,SP_POP,ZI_POP,ZS_POP,ZW_POP)
central_africa <- c(CO_POP,RG_POP,GA_POP)
west_africa <- c(GU_POP,NG_POP)
east_africa <- c(UG_POP,EA_POP,EB_POP,ED_POP,EF_POP,ER_POP)

### creating a list of regions
cont <- list(ocenia,asia,america,europe_north_africa,southern_africa,central_africa,west_africa,east_africa)
names(cont) <- c("Ocenia","Asia","North America","Europe & North Africa","Southern Africa","Central Africa","West Africa","East Africa")

#### 3. Extract genotypes for lectin-24A SNPs in DPGP samples from vcf  ####

dataChunk <- read.table("../Input/CRISP_dpgp_lectin24a_annot.vcf")  ### variants in DPGP were called jointly for all samples genome-wide using CRISP, vcf with lectin region was previously subsetted
dataChunk <- dataChunk[dataChunk$V2 < 3718142,] ### restrict to upstream lectin region
dataChunk <- dataChunk[!duplicated(dataChunk$V2),] ### remove duplicated positons
dataChunk <- dataChunk[!dataChunk$V9 == "GT:GQ:DP:ADf:ADr",] ### remove SNPs without any genotypes in samples
gt <- as.data.frame(dataChunk[,10:ncol(dataChunk)])  ### genotype data
fix <- as.data.frame(dataChunk[,1:8]) ### meta data

### shorten sample names
names_files <- read.table("../Input/file.list",stringsAsFactors = FALSE)
for (n in 1:nrow(names_files)) {
  names_files[n,] <- gsub("_paired_wasp_mapped_reorder_readgroup_mq20_sort_marked_realign.bam","",names_files[n,])
}
names_files <- as.character(names_files$V1)
colnames(gt) <- names_files

### process vcf files to extract the reference and alternate calls for each sample
new2 <- ""
for (i in 1:ncol(gt)){
  gt_sub <- gt[i] ### per sample unprocessed genotype
  new <- gt_sub %>% separate(colnames(gt_sub), into =c("MLAC","GQ","DP","ADf","ADr","ADb"), sep = ":") ### MLAC contains genotypes called by CRISP
  new <- new %>% separate(MLAC, into = c("ref","gt"),sep="/")  ### after / is genotype of sample, expecting single call in homozygous genomes
  if(nrow(new[new$ref != 0,]) > 0){
    new[new$ref != 0,]$gt <- NA
  }
  genosample <- as.data.frame(new$gt)
  colnames(genosample) <- names_files[i] ### add SRA names to samples
  genosample[as.numeric(new$GQ ) < 10 | as.numeric(new$GQ ) %in% NA,] <- NA  ### if genotype quality less than 10 or NA, make genotype call NA
  new2 <- cbind(new2,genosample)
}
new2 <- new2[,-1] ### remove empty row
full <- cbind (fix,new2) ### add genotype and meta data
full <- full[full$V8 %like% "VT=SNV;",]  ### only keep SNPs

### convert gt calls to numeric
full[,9:ncol(full)] <- apply(full[,9:ncol(full)],2,as.character)
full[,9:ncol(full)] <- apply(full[,9:ncol(full)],2,as.numeric)


#### 4. Get variants annotations and variant names  ####

### Get variants SNPEff annotations from info column in vcf 
### also get protein names for coding region SNPs and DNA based names for the remainder, names obtained from SnpEff

ann=""
ann_allsnps=""
ann_dna_snpname=""
ann_prot_snpname=""
for (i in 1:nrow(full)){
  altalleles <- str_split(full[i,5],",")[[1]]
  ### extract annotation for biallelic SNPs
  if(length(altalleles) < 2){
    ann[i] <- strsplit(full$V8[i], split ="\\|")[[1]][2]
    ann_dna_snpname[i] <- strsplit(full$V8[i], split ="\\|")[[1]][10]
    ann_prot_snpname[i] <- strsplit(full$V8[i], split ="\\|")[[1]][11]
    ann_allsnps <- rbind(ann_allsnps,c(full$V2[i],altalleles,ann[i]))
  }
  ### extract annotation for multiallelic SNPs
  if(length(altalleles) > 1){
    annwrk2 <- strsplit(full$V8[i], split ="=|\\|,")[[1]]
    annwrk2 <- annwrk2[14:length(annwrk2)]
    annwrk2 <- as.data.frame(str_split_fixed(annwrk2,"\\|",12))
    annwrk2 <- annwrk2[1:11]
    annwrk2 <- annwrk2[!duplicated(annwrk2$V1),]
    annwrk2 <- annwrk2[annwrk2$V1 %in% na.omit(altalleles),]
    ### if all alleles have same annotation, use that one
    if(length(table(annwrk2$V2))<1){
      
      ann[i] <- annwrk2$V2[1]
      ann_dna_snpname[i] <- strsplit(full$V8[i], split ="\\|")[[1]][10]
      ann_prot_snpname[i] <- strsplit(full$V8[i], split ="\\|")[[1]][11]
      ann_allsnps <- rbind(ann_allsnps,c(full$V2[i],annwrk2$V2[1:2]))
    } else{
      ### if one of the alternate alleles is at much higher frequency, extract annotation for that SNP
      gencount <- table(as.numeric(full[i,9:ncol(full)]))
      gencount <- gencount[-c(which(names(gencount) %in% 0))]
      if (max(prop.table(gencount))>0.79){
        ann[i] <- annwrk2[as.numeric(names(gencount)[gencount == max(gencount)]),]$V2
        ann_dna_snpname[i] <- strsplit(full$V8[i], split ="\\|")[[1]][10]
        ann_prot_snpname[i] <- strsplit(full$V8[i], split ="\\|")[[1]][11]
      } else {
        ### if there are multiple alleles at similar freq with different annotations, code annotation as NA for that site
        ann[i] <- NA
      }
      ann_allsnps <- rbind(ann_allsnps,as.matrix(unname(cbind(full$V2[i],annwrk2[1:2]))))
    }
    
  }
}

### use protein based names if it exists, else use dna based names for variants
ann_snpname <- ann_prot_snpname
for (n in 1:length(ann_snpname)){
  if(sjmisc::is_empty(ann_snpname[n])){
    ann_snpname[n] <- ann_dna_snpname[n]
  } 
}
names(ann_snpname) <- paste(full$V1,full$V2,ann,sep="_")


#### 5. Add variants annotations and sample and population metadata to SNP matrix  ####

full <- full[!duplicated(full$V2),] ### remove duplicated sites if exists
full_snps <- full[9:ncol(full)] ### get genotype calls only
rownames(full_snps) <- paste(full$V1,full$V2,ann,sep="_") ## rownames combination of chr and pos
full_snps$ann <- ann ## annotation column for genotype matrix
full_snps$Pos <- full$V2 ## positions column for genotype matrix
full_snps <- full_snps[!is.na(full_snps$ann),]  ### remove SNPs without any annotation, mainly multiallelic SNPs with different effects

### process genotype matrix, add run names
full_snps <- full_snps[-c((ncol(full_snps)-1):ncol(full_snps))]  ### get genotype matrix
full_snps <- as.data.frame(t(full_snps)) #transpose genotype matrix
full_snps$Run <- rownames(full_snps) ### add run column using rownames
sranames <- full_snps$Run ### separate vector for run names
sranames <- as.data.frame(sranames)
colnames(sranames) <- "Run"

### convert sra accessions to sample names and add meta data
sranames_meta <- read.csv("../Input/SraRunTable_dpgp2.txt")
sranames_meta_all <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_all) <- c("Run","Sample")

sranames_meta <- read.csv("../Input/SraRunTable_dpgp3.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("../Input/SraRunTable_dgrp_sra2.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$strain)))
sranames_meta_sub$V2 <- gsub("DGRP-","RAL-",sranames_meta_sub$V2)
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("../Input/SraRunTable_bergman.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Library.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("../Input/SraRunTable_pool.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("../Input/SraRunTable_nuzhdin1.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("../Input/SraRunTable_clark.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("../Input/SraRunTable_ages.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)
sranames_meta_all <- na.omit(sranames_meta_all)

sranames_meta <- left_join(sranames,sranames_meta_all, by =c("Run"="Run")) ### sra names with sample names and meta data

### modify sample names for compatibility
sranames_meta$Sample <- gsub("-HE","",sranames_meta$Sample)
sranames_meta$Sample <- gsub("_HE","",sranames_meta$Sample)
sranames_meta$Sample <- gsub("diTAG-","",sranames_meta$Sample)
sranames_meta$Sample <- gsub("_new","",sranames_meta$Sample)
sranames_meta$Sample <- gsub("CO8-3","CO8N",sranames_meta$Sample)
sranames_meta$Sample <- gsub("CO10-3","CO10N",sranames_meta$Sample)
sranames_meta$Sample <- gsub("CO13-3","CO13N",sranames_meta$Sample)
individuals <- read.csv(file = "../Input/TableS1_individuals.csv")
individuals <- individuals[,c(1,4,10,11)]
datades <- read.csv(file = "../Input/TableS2_populations.csv")
individuals_des <- left_join(individuals,datades, by =c("Population"="Population.ID") )
sranames_meta_des <- inner_join(sranames_meta,individuals_des, by = c("Sample"="Stock.ID"))

### add sample names and meta data with population descriptions to genotype matrix
full_snps <- inner_join(sranames_meta_des,full_snps, by =c("Run"))


#### 6. Combine SNPs and upstream indels  ####

dpgpv <- full_snps ### combine with SNP dataset
dpgpv[,16:ncol(dpgpv)] <- apply(dpgpv[,16:ncol(dpgpv)],2,as.numeric)

#### 7. Combine duplicated samples  ####
### some samples have multiple sra accessions as they were seqeunced multiple times
strains_unique <- unique(dpgpv$Sample) ### unique sample names
fullsnp <- ""
### loop over unique samples
for (s in 1:length(strains_unique)){
  strains_unique_data <- dpgpv[dpgpv$Sample %in% strains_unique[s],]
  if (nrow(strains_unique_data) == 1){
    ### no change is each sample only has one sra accession
    fullsnp <- rbind(fullsnp,strains_unique_data)
  } else{
    strains_unique_data_mod <- strains_unique_data[1,1:16]
    ### if multiple sra accession exist, check if variant calls match among them
    ### looping over each site
    for(c in 17:ncol(strains_unique_data)){
      ### if one sra is NA and the other is called, use called from non-NA sra accession
      if (length(table(strains_unique_data[,c])) < 2){
        if(length(unique(strains_unique_data[,c])) > 1){
          geno <- unique(na.omit(strains_unique_data[,c]))
        } else {
          ### if variants call is same in all sra samples, use first call
          geno <- unique(strains_unique_data[,c])
        }
        strains_unique_data_mod <- cbind(strains_unique_data_mod,geno)
      } else {
        ### if variants call is different in  sra samples, code as NA
        strains_unique_data_mod <- cbind(strains_unique_data_mod,"NA")
      }
      
    }
    colnames(strains_unique_data_mod)[17:ncol(strains_unique_data_mod)] <- colnames(strains_unique_data)[17:ncol(strains_unique_data)]
    fullsnp <- rbind(fullsnp,strains_unique_data_mod)
  }
}
dpgpv <- fullsnp[-c(1),] ### remove empty row

dpgpv_new <- dpgpv

##only retain SNPs with allele frequency  > 0.05 in at least two populations
tmpnew <- dpgpv_new[1:16]
for (i in 17:ncol(dpgpv_new)){
  ### only working with SNPs with maximum of three alleles, including reference call
  if (ncol(table(dpgpv_new$Population,dpgpv_new[,i])) > 1 & ncol(table(dpgpv_new$Population,dpgpv_new[,i])) < 4){
    if (length(table(dpgpv_new$Population)) - sum(table(dpgpv_new$Population,dpgpv_new[,i])[,2]/rowSums(table(dpgpv_new$Population,dpgpv_new[,i])) < 0.05 | table(dpgpv_new$Population,dpgpv_new[,i])[,2]/rowSums(table(dpgpv_new$Population,dpgpv_new[,i])) > 0.95, na.rm=T) > 1){
      tmpnew <- cbind(tmpnew,dpgpv_new[i])
    }
  }
}
dpgpv_new <- cbind(tmpnew,dpgpv_new[ncol(dpgpv_new)])
####process dataframe to add back chr and pos
dpgpvt <- as.data.frame(t(dpgpv_new[,17:ncol(dpgpv_new)]))
colnames(dpgpvt) <- dpgpv_new$Sample
dpgpvt$Chr <- colsplit(rownames(dpgpvt),"_",names=c("Chr","Pos","c"))[1]$Chr
dpgpvt$Chr <- gsub(".*ndel.*","2L",dpgpvt$Chr)
dpgpvt$Pos <- as.numeric(gsub("X","",colsplit(rownames(dpgpvt),"_",names=c("Chr","Pos","c"))[2]$Pos))
dpgpvt$Pos[which(is.na(dpgpvt$Pos))] <- 3718040 ### remove 8bp indel that was called using unmodified reference

### modify variant names having filtered SNPs
### modify positions of upstream SNPs occurring before the 21 BP indel to account of the 21 BP indel
### not necessary for 8BP since insertion is derived state and 7BP is the most upstream variant
ann_snpname <- ann_snpname[names(ann_snpname) %in% colnames(dpgpv_new)] ### modfify variant names to only include filtered SNPs
ann_snpname_num <- 3717728-as.numeric(gsub("_.*","",gsub("2L_","",names(ann_snpname))))
ann_snpname2 <- ann_snpname[ann_snpname_num>(-171)] ### those downstream of 21BP indel 
ann_snpname3 <- ann_snpname[ann_snpname_num<(-171)] ### those upstream of 21BP indel 
ann_snpname_num <- ann_snpname_num[ann_snpname_num<(-171)]
ann_snpname_num <- ann_snpname_num-21 ### change numbering of variants upstream of the 21 Bp indel
ann_snpname4 <- paste("c.",ann_snpname_num,gsub(".*[0-9]","",ann_snpname3),sep="")
names(ann_snpname4) <- names(ann_snpname3)

### add indel names to naming list, 3 upstream indels and 171 bp coding indel
ann_snpname <- c(ann_snpname2,ann_snpname4)
extranames <- c("c.-439_-433del","c.-334_-333insACATTCAT","c.-171_-151del","p.Phe217_Glu273del*")
names(extranames) <- c("Indel7bp_X3718140","Indel8bp_X3718040","Indel21bp_X3717878","Indel_X3716932")
ann_snpname <- c(ann_snpname,extranames) ### snp names by id
ann_snpname_pos <- ann_snpname
names(ann_snpname_pos) <- as.numeric(gsub("2L_||.*_X||_.*","",names(ann_snpname)))
ann_snpname_pos <- rev(ann_snpname_pos[order(as.numeric(names(ann_snpname_pos)))])  ### snp names by position


### this is a new part. It check if SNPs are called in 6 BP region on either side of presumed indel. If they are not, then do not call deletion. Wasn't necessary in main analyses due to 50% lectin genotyped requirement.

dpgpvt_deleteion2 <- dpgpvt[dpgpvt$Pos < 3716910,]
missingdata <- ""
for (i in 1:(ncol(dpgpvt_deleteion2)-2)){
  tab <- table(dpgpvt_deleteion2[,i],useNA = "ifany")
  if(length(tab[names(tab) %in% NA])>0){
    missingdata <- rbind(missingdata,cbind(colnames(dpgpvt_deleteion2)[i],tab[names(tab) %in% NA]/sum(tab)))
  }
}
missingdata <- as.data.frame(missingdata)
missingdata$V2 <- as.numeric(missingdata$V2)

dpgpvt_deleteion3 <- dpgpvt[dpgpvt$Pos < 3717271 & dpgpvt$Pos > 3717081,]
missingdata2 <- ""
for (i in 1:(ncol(dpgpvt_deleteion2)-2)){
  tab <- table(dpgpvt_deleteion2[,i],useNA = "ifany")
  if(length(tab[names(tab) %in% NA])>0){
    missingdata2 <- rbind(missingdata2,cbind(colnames(dpgpvt_deleteion2)[i],tab[names(tab) %in% NA]/sum(tab)))
  }
}
missingdata2 <- as.data.frame(missingdata2)
missingdata2$V2 <- as.numeric(missingdata2$V2)


#### 10. Call 171 BP coding deletion in DPGP samples  ####

### samples designed as having deletion if NA for region between 3717102 and 3716904
### 93 samples designed as having deletion this way
### pindel calls deletions in all of these 93 samples
### however pindel also calls deletions is an extra ~25 samples which have non-NA values here with high confidence.
### manual inspection of BAM files does not indicate evidence for a deletion for this extra 25 samples. So only using these 97 which have been confirmed by pindel and manual inspection of bam files

pindel_calls <- read.csv(file="../Input/pindel_calls.csv") ## read in file containing samples where pindel has called a deletion
pindel_calls$Sample <- gsub("-HE","",pindel_calls$Sample) ### modify file names for compatability
pindel_calls$Sample <- gsub("_HE","",pindel_calls$Sample)
pindel_calls$GT <- 1 ### since file only contains samples where deletion was called
pindel_calls <- pindel_calls[!duplicated(pindel_calls$Sample),] ### remove duplicated samples, samples with multiple SRA accessions

dpgpvt_deleteion <- dpgpvt[dpgpvt$Pos < 3717081 & dpgpvt$Pos > 3716910,] ### extract SNPs within the region where deletion occurs, 9 SNPs total
dpgpvt_deleteion[is.na(dpgpvt_deleteion)] <- "NA"
Indel_X3716932 <- ""
for (n in 1:(ncol(dpgpvt_deleteion)-2)){
  ### if all 9 SNPs are NA, assign as deletion
  if (sum(dpgpvt_deleteion[,n]=="NA")==nrow(dpgpvt_deleteion)){
    Indel_X3716932[n] <- 1
    ### if even one the 9 SNPs are called, assign as no deletion
  } else if (sum(dpgpvt_deleteion[,n]=="NA")<nrow(dpgpvt_deleteion)){
    Indel_X3716932[n] <- 0
  } else {
    Indel_X3716932[n] <- NA
  }
}
Indel_X3716932 <- c(Indel_X3716932,"2L","3716932") ### the position where the deletion begins is actually 3716904, this is corrected below
dpgpvt <- rbind(dpgpvt,Indel_X3716932)
rownames(dpgpvt)[nrow(dpgpvt)] <- "Indel_X3716932"
names(Indel_X3716932) <- colnames(dpgpvt)

Indel_X3716932[names(Indel_X3716932) %in% unique(c(missingdata[missingdata$V2 > 0.5,]$V1,missingdata2[missingdata2$V2 > 0.5,]$V1))] <- NA ## this part makes samples with over 50% missing data on either side of presumed deletion as NA


### create lectin coding deletion dataframe for comparison
Indel_X3716932_df <- as.data.frame(Indel_X3716932)
Indel_X3716932_df$Sample <- rownames(Indel_X3716932_df)
Indel_X3716932_df <- left_join(Indel_X3716932_df,pindel_calls,by="Sample") ### combine pindel calls with approach used here
table(Indel_X3716932_df$Indel_X3716932,Indel_X3716932_df$GT) ###  90 samples called as having coding deletions here also called as having deletions by pindel. Pindel calls deletions is additonal 22 samples, but manual BAM inspection does not support it and CRISP has called variants in that region. So ignoring those calls. Also, pindel was not run on another 6 samples, unclear why. But since the CRISP calls have already been tested, going with that approach.

codingdelcall <- dpgpvt[nrow(dpgpvt),]
codingdelcall <- codingdelcall[-c(ncol(codingdelcall))]
codingdelcall <- codingdelcall[-c(ncol(codingdelcall))]

codingdelcall2 <- t(as.data.frame(c("2L","3716932","stop_gained","T","A","220",".",".","GT",paste(as.numeric(codingdelcall),"/",as.numeric(codingdelcall),sep=""))))
colnames(codingdelcall2) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",names(codingdelcall))
rownames(codingdelcall2) <- NA

write.table(codingdelcall2,"coding_del.vcf",quote = F, row.names = F,sep="\t")

europe_north_africa <- read.table(file="../Input/europe_north_africa.txt")
america <- read.table(file="../Input/america.txt")
southern_africa <- read.table(file="../Input/southern_africa.txt")

table(as.numeric(codingdelcall[colnames(codingdelcall) %in% europe_north_africa$V1]))[2]/length(as.numeric(codingdelcall[colnames(codingdelcall) %in% europe_north_africa$V1]))
table(as.numeric(codingdelcall[colnames(codingdelcall) %in% america$V1]))[2]/length(as.numeric(codingdelcall[colnames(codingdelcall) %in% america$V1]))
table(as.numeric(codingdelcall[colnames(codingdelcall) %in% southern_africa$V1]))[2]/length(as.numeric(codingdelcall[colnames(codingdelcall) %in% southern_africa$V1]))



