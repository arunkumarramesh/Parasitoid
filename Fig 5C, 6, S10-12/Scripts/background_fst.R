setwd("/scratch/Arun/Projects/flydata/2all")
library(tidyverse)
library(dplyr)

### Aim of script is to get significance threshold for autosomal Fsts and plot histogram of Fst in DPGP
### Note that script cannot be run alone, require larges files found on private server (Multivac)

ann <- read.table(file = "CRISP_dpgp_Annotations", header = F, sep=" ", stringsAsFactors = F) ### read in annotations for SNPs
colnames(ann) <- c("Chr","Pos","Allele","Ann")
ann <- ann[!(ann$Pos < 3718141 & ann$Pos > 3716799),] ### remove lectin region variants, will be added in later from vcf file
ann$ID <- paste(ann$Chr,ann$Pos,sep="_") ### get id by combining chr, position and annotation
ann <- ann[!(duplicated(ann$ID) | duplicated(ann$ID, fromLast = TRUE)), ] ### remove variants with duplicated ids
ann <- ann[-c(1:3)]


lectin_ann <- read.table(file="lectin_ann",header=T,stringsAsFactors = F) ### read in lectin region annotations
lectin_ann <- lectin_ann[order(2:1)]
lectin_ann$POS <- paste("2L",lectin_ann$POS,sep="_") ## generate positions from ID
colnames(lectin_ann) <- colnames(ann)

lectin_ann$Ann <- gsub("upstream_7BP_indel", "7BP",lectin_ann$Ann)
lectin_ann$Ann <- gsub("upstream_8BP_indel", "8BP",lectin_ann$Ann)
lectin_ann$Ann <- gsub("upstream_21BP_indel", "21BP",lectin_ann$Ann)
lectin_ann$Ann <- gsub("coding_indel", "Coding Deletion",lectin_ann$Ann)
lectin_ann$ID <- gsub("2L_3716932", "2L_3717081",lectin_ann$ID)

### custom make indel annotations
#annindel <- as.data.frame(cbind(c("7BP","8BP","21BP","Coding Deletion"),c("2L_3718140","2L_3718040","2L_3717878","2L_3717081")))
#colnames(annindel) <- colnames(ann)

ann <- rbind(ann,lectin_ann) ### add lectin snp and indel annotations to background annotations file
ann <- ann[order(ann$ID),] ## order annotations by ID


perlocusfstall <- read.table(file = "dpgp_vcftools_POP.weir.fst", stringsAsFactors = F, header = T) ### read in file with annotations for all SNPs
perlocusfstall <- perlocusfstall[!is.na(perlocusfstall$WEIR_AND_COCKERHAM_FST),] ### remove sites without Fst calls
perlocusfstall <- perlocusfstall[!perlocusfstall$CHR %in% "X",] ### only analyse autosomal sites

## this section identified SNPs in short introns (i.e. between XX and XX)
introns <- read.table(file = "dmel-all-r5.13_intron.gff",sep="\t") ## open file containing all intron positions
introns <- introns[introns$V5 - introns$V4 < 66,] ## short introns defined as introns of 65bp of less (https://doi.org/10.1093/molbev/msq046)
introns <- introns[introns$V3 %in% "intron",]
introns <- introns[introns$V5 - introns$V4 > 49,] ## but ensure it is at least 50bp since hard to get 8-30 bp for things shorter
## this part gets based 8-30 bp in short introns, the sites that are most free to evolve (https://doi.org/10.1093/molbev/msq046)
introns[introns$V7 %in% "+",]$V4 <- introns[introns$V7 %in% "+",]$V4+7
introns[introns$V7 %in% "+",]$V5 <- introns[introns$V7 %in% "+",]$V4+22
introns[introns$V7 %in% "-",]$V5 <- introns[introns$V7 %in% "-",]$V5-7
introns[introns$V7 %in% "-",]$V4 <- introns[introns$V7 %in% "-",]$V5-22
introns <- introns[c(1,4,5)]
introns <- introns[introns$V1 %in% c("2L","2R","3L","3R"),] ## only retain introns in autosomes

## note this loop takes very long to run (~12 hrs)! Goal is to identify all snps used for calculating FST within short introns
newlist <- ""
for (i in 1:nrow(perlocusfstall)){ ## loop over each site where fst was calculated
  print(i)
  introns_part <- introns[introns$V1 %in% perlocusfstall$CHR[i],] ## select all introns that occur in the chromsome that matches the site being studied
  ## the next few lines checks if postion in within the range of the intron bed. Basically, if a site is within any of those ranges, the  start position column value will be negative and the end position value will be positive
  introns_part$V4 <- introns_part$V4 - perlocusfstall$POS[i] 
  introns_part$V5 <- introns_part$V5 - perlocusfstall$POS[i]
  ## if the site occurs at the start of end of introns, will be coded as 0 or artifically make then - or +
  if(length(table(introns_part$V4 %in% 0))>1){
    introns_part[introns_part$V4 %in% 0,]$V4 <- -1
  }
  if(length(table(introns_part$V5 %in% 0))>1){
    introns_part[introns_part$V5 %in% 0,]$V5 <- 1
  }
  if(nrow(introns_part[ introns_part$V4 * introns_part$V5 < 0, ])>0){ ## to get rows where first column value is negative and second value is positive. Then if site in within the range of even a single short intron, include it for analyses
    newlist <- rbind(newlist, cbind(perlocusfstall$CHR[i],perlocusfstall$POS[i]))
  }
}
newlist <- newlist[-c(1),]
newlist <- as.data.frame(newlist)
newlist$V2 <- as.numeric(as.character(newlist$V2))
write.csv(newlist, file="newlist.csv") ## save list of short intron snps for sites being analysed for fst

lectin_masked_fst <- read.table(file="lectin_dpgp_stop_gain_masked_vcftools_POP.weir.fst",header=T,stringsAsFactors = F) ### read in file with annotations for lectin region variants
lectin_masked_fst_indel <- lectin_masked_fst[lectin_masked_fst$POS %in% c(3718140,3718040,3717878,3716932),] ### extract fst from indels
lectin_masked_fst_indel[lectin_masked_fst_indel$POS == 3716932,]$POS <- 3717081  ### change to correct position for 171 bp coding indel
colnames(lectin_masked_fst_indel) <- colnames(perlocusfstall) 
perlocusfstall <- rbind(perlocusfstall,lectin_masked_fst_indel) ### combine lectin indel fst to SNP fst
perlocusfstall <- perlocusfstall[order(perlocusfstall$CHR,perlocusfstall$POS),] ### order fst df by chr and positon
perlocusfstall$WEIR_AND_COCKERHAM_FST <- pmax(perlocusfstall$WEIR_AND_COCKERHAM_FST,0) ## make negative fst values  to 0
perlocusfstall$ID <- paste(perlocusfstall$CHR,perlocusfstall$POS,sep="_") ### make id by combining chr and position
snpsOfInterest <- perlocusfstall[perlocusfstall$CHR %in% "2L" & perlocusfstall$POS < 3718050 & perlocusfstall$POS > 3716812,]$ID ### get lectin region snp ids
perlocusfstall <- perlocusfstall  %>% mutate( is_highlight=ifelse(ID %in% snpsOfInterest, "yes", "no")) ### highlight column if lectin region snps
perlocusfstall <- left_join(perlocusfstall,ann,by="ID") ### add annotations by ID column
perlocusfstall <- perlocusfstall  %>% mutate( shortintron=ifelse(ID %in% newlist$ID, "yes", "no")) ### highlight column if short intron snps
lectinregion <- perlocusfstall[perlocusfstall$CHR %in% "2L" & perlocusfstall$POS < 3717728 & perlocusfstall$POS > 3716880,] ### extract lectin region fst
perlocusfstall2 <- perlocusfstall[perlocusfstall$shortintron %in% c("yes"),] ## only considering silent variants here

### calculate number of mutations and 1% Fst percentile (99% - 0.412153,  99.9%, 0.6445234, 8440512 SNPs)
quantile(perlocusfstall$WEIR_AND_COCKERHAM_FST,0.99,na.rm=T)
quantile(perlocusfstall$WEIR_AND_COCKERHAM_FST,0.999,na.rm=T)
dim(perlocusfstall)
### calculate number of mutations and 1% Fst percentile for silent mutations (99% - 0.4070053,  99.9%, 0.614813, 20783 SNPs)
quantile(perlocusfstall2$WEIR_AND_COCKERHAM_FST,0.99,na.rm=T)
quantile(perlocusfstall2$WEIR_AND_COCKERHAM_FST,0.999,na.rm=T)
dim(perlocusfstall2)

### extract lectin region fst with some surrounding region
lectin_1kb_arms <- perlocusfstall[perlocusfstall$CHR %in% "2L" & perlocusfstall$POS < 3719050 & perlocusfstall$POS > 3715513,] 
dim(lectin_1kb_arms)
write.csv(lectin_1kb_arms,file="lectin_1kb_arms.csv",row.names=F) 

### plot histogram of fst with lectin coding region variants highlighed by lines and 1% threshold highlighted as red region
fsthist <- ggplot(perlocusfstall2, aes( WEIR_AND_COCKERHAM_FST)) +
  geom_rect(aes(xmin=quantile(perlocusfstall2$WEIR_AND_COCKERHAM_FST,0.99,na.rm=T),xmax=max(perlocusfstall2$WEIR_AND_COCKERHAM_FST,na.rm=T),ymin=0,ymax=8),fill = "pink", alpha = 0.03)+
  geom_histogram(bins = 1000,aes(y = ..density..))+
  geom_vline(data=lectinregion[!is.na(lectinregion$Ann),],aes(xintercept=WEIR_AND_COCKERHAM_FST, color=Ann))+
  ylim(0,8)+
  xlab("Fst")+
  ylab("Percentage of sites")

png("fsthist.png",width = 700)
fsthist
dev.off()



