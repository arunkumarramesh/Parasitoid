library(methods)
library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)
library(data.table)

### Aim of script is to process and filter vcf files to make it compatible with vcftools for estimating Fst per site
### The original vcf file was split by chromosome and into five further chucks to aid computation, a script is associated with each chunk
### This is an example script for one of those chunks.
### Note that script cannot be run alone, require larges files found on private server (Multivac)

### convert sra accessions to sample names and add meta data

sranames_meta <- read.csv("SraRunTable_dpgp2.txt")
sranames_meta_all <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_all) <- c("Run","Sample")

sranames_meta <- read.csv("SraRunTable_dpgp3.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("SraRunTable_dgrp_sra2.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$strain)))
sranames_meta_sub$V2 <- gsub("DGRP-","RAL-",sranames_meta_sub$V2)
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("SraRunTable_bergman.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Library.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("SraRunTable_pool.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("SraRunTable_nuzhdin1.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("SraRunTable_clark.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)

sranames_meta <- read.csv("SraRunTable_ages.txt")
sranames_meta_sub <- as.data.frame(cbind(as.character(sranames_meta$Run),as.character(sranames_meta$Sample.Name)))
colnames(sranames_meta_sub) <- c("Run","Sample")
sranames_meta_all <- rbind(sranames_meta_all,sranames_meta_sub)
sranames_meta_all <- na.omit(sranames_meta_all)

individuals <- read.csv(file = "TableS1_individuals.csv")
individuals <- individuals[,c(1,4,10,11)]
datades <- read.csv(file = "TableS2_populations.csv")
individuals_des <- left_join(individuals,datades, by =c("Population"="Population.ID") )

### get sample names
names_files <- read.table("file.list",stringsAsFactors = FALSE)
for (n in 1:nrow(names_files)) {
  names_files[n,] <- gsub("_paired_wasp_mapped_reorder_readgroup_mq20_sort_marked_realign.bam","",names_files[n,])
}
names_files <- as.character(names_files$V1)


#### Mask regions with ibd, admixture and heterozygosity, from DPGP analyses
### regions for masking obtained from DPGP paper
ibd <- read.table(file = "ibd_filter_tracts.txt",sep="\t")
admix <- read.table(file = "admixture_filter_tracts.txt",sep="\t")
het <- read.table(file = "het_filter_tracts.txt",sep="\t")
fil <- rbind(ibd,admix,het)
fil$V2 <- gsub("Chr","",fil$V2)

transactFile <- 'CRISP_dpgp_2R_1-5500000_ann.vcf'
index <- 0
chunkSize <- 1000
con <- file(description=transactFile,open="r")
dataChunk <- read.table(con, nrows=chunkSize,stringsAsFactors = FALSE)
full <- dataChunk

### repeat over 1000 chunks of variants
repeat {
  index <- index + 1
  print(paste('Processing rows:', index * chunkSize))
  
  full <- full[full$V8 %like% "VT=SNV;",] ### only keep SNPs
  full <- full[nchar(full$V4) < 2 & nchar(full$V5) < 2,] ### only keep biallelic SNPs
  
  ##separate genotype and meta data 
  gt <- as.data.frame(full[,10:ncol(full)]) ### genotype data
  fix <- as.data.frame(full[,1:8]) ### meta data
  colnames(gt) <- names_files ## add sra sample names
  
  ### process vcf to extract genotype for each sample
  new2 <- ""
  for (i in 1:ncol(gt)){
    gt_sub <- gt[i] ### per sample unprocessed genotype
    new <- gt_sub %>% separate(colnames(gt_sub), into =c("MLAC","GQ","DP","ADf","ADr","ADb"), sep = ":") ### MLAC contains genotypes called by CRISP
    new <- new %>% separate(MLAC, into = c("ref","gt"),sep="/") ### after / is genotype of sample, expecting single call in homozygous genomes
    if(nrow(new[new$ref != 0,]) > 0){
      new[new$ref != 0,]$gt <- NA
    }
    genosample <- as.data.frame(new$gt)
    colnames(genosample) <- names_files[i] ### add SRA names to samples
    genosample[as.numeric(new$GQ ) < 10 | as.numeric(new$GQ ) %in% NA,] <- NA  ### if genotype quality less than 10 or NA, make genotype call NA
    new2 <- cbind(new2,genosample)
  }
  new2 <- new2[,-1] ### remove empty row
  full <- cbind(fix,new2) ### add genotype and meta data
  
  ### convert gt calls to numeric
  full[,9:ncol(full)] <- apply(full[,9:ncol(full)],2,as.character)
  full[,9:ncol(full)] <- apply(full[,9:ncol(full)],2,as.numeric)

  #### Get variants SNPEff annotations from info column in vcf and added to ann_allsnps df  
  ann=""
  ann_allsnps=""
  for (i in 1:nrow(full)){
    altalleles <- str_split(full[i,5],",")[[1]]
    ann[i] <- strsplit(full$V8[i], split ="\\|")[[1]][2]
    ann_allsnps <- rbind(ann_allsnps,c(full$V2[i],altalleles,ann[i]))
  }
  ann_allsnps <- ann_allsnps[-c(1),] ## remove blank line from annotations
  
  full_snps <- full[9:ncol(full)] ### get genotype calls only
  rownames(full_snps) <- paste(full$V1,full$V2,ann,sep="_") ## rownames combination of chr and pos
  full_snps$ann <- ann ## annotation column for genotype matrix
  full_snps$Pos <- full$V2 ## positions column for genotype matrix
  full_snps <- full_snps[!is.na(full_snps$ann),] ### remove sites with no or multiple annotations
  full_snps <- as.data.frame(t(full_snps[1:(ncol(full_snps)-2)])) #transpose genotype matrix
  full_snps$Run <- rownames(full_snps) ### add run column using rownames
  sranames <- full_snps$Run ### separate vector for run names
  sranames <- as.data.frame(sranames)
  colnames(sranames) <- "Run"
  
  ### process sample names for comparability
  sranames_meta <- left_join(sranames,sranames_meta_all, by =c("Run"="Run"))
  sranames_meta$Sample <- gsub("-HE","",sranames_meta$Sample)
  sranames_meta$Sample <- gsub("_HE","",sranames_meta$Sample)
  sranames_meta$Sample <- gsub("diTAG-","",sranames_meta$Sample)
  sranames_meta$Sample <- gsub("_new","",sranames_meta$Sample)
  sranames_meta$Sample <- gsub("CO8-3","CO8N",sranames_meta$Sample)
  sranames_meta$Sample <- gsub("CO10-3","CO10N",sranames_meta$Sample)
  sranames_meta$Sample <- gsub("CO13-3","CO13N",sranames_meta$Sample)
  
  ### combine genotype and meta data
  sranames_meta_des <- inner_join(sranames_meta,individuals_des, by = c("Sample"="Stock.ID"))
  dpgpv <- inner_join(sranames_meta_des,full_snps, by =c("Run"))
  
  dpgpv[,17:ncol(dpgpv)] <- apply(dpgpv[,17:ncol(dpgpv)],2,as.numeric) ### make genotype calls numeric
  dpgpv[,1:16] <- apply(dpgpv[,1:16],2,as.character) ### make meta data character vector
  
  ### some samples have multiple sra accessions as they were seqeunced multiple times so combine them here
  
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
  dpgpv <- fullsnp[-c(1),]
  dpgpv <- dpgpv[!dpgpv$Sample %in% "ZW184",] ### first remove sample designated as outlier in global diversity lines (GDL)
  
  Popover3 <- read.table(file="Popover3") ### samples analysed for lectin-24A molecular evolution
  dpgpv <- dpgpv[dpgpv$Sample %in% Popover3$V1,]  ### background analyses restricted to samples those analysed in lectin-24a analyses!
  dpgpv_new <- dpgpv
  rownames(dpgpv_new)<- dpgpv_new$Sample
  
  rownames(dpgpv_new)<- dpgpv_new$Sample ### sample names now rownames
  dpgpvt <- as.data.frame(t(dpgpv_new[,17:ncol(dpgpv_new)])) ### extract just genotype data
  dpgpvt$Chr <- colsplit(rownames(dpgpvt),"_",names=c("Chr","Pos","c"))[1]$Chr ### add chr column
  dpgpvt$Pos <- as.numeric(colsplit(rownames(dpgpvt),"_",names=c("Chr","Pos","c"))[2]$Pos) ### add Pos column
  
  ### basically, loop over samples and replaces masked region with NA
  fulldata <-  dpgpvt[1:(ncol(dpgpvt)-2)]
  ### loop over sample
  for (c in 1:(ncol(dpgpvt)-2)){
    newdat <- cbind(dpgpvt[c],dpgpvt["Chr"],dpgpvt["Pos"]) ### isolate sample data, chr and position from genotype matrix
    mask <- fil[fil$V1 %in% colnames(dpgpvt[c]),] ### get all regions to be masked for a given sample
    newdat_all <- ""
    newdat_all <- as.data.frame(newdat_all)
    colnames(newdat_all) <- colnames(dpgpvt[c])
    ###if there are regions to mask, proceed
    if(nrow(mask) > 0){
      for(ch in 1:length(unique(mask$V2))){
        ### loop over each chromosome
        mask_ch <- mask[mask$V2 %in% unique(mask$V2)[ch],] ### get regions to be masked for each chr
        newdat_ch <- newdat[newdat$Chr %in% unique(mask$V2)[ch],] ### isolate chr for each samples
        ###if there are regions to mask for a given chromosome, proceed
        if (nrow(mask_ch) > 0 & nrow(newdat_ch) > 0){
          ### basically mask by start and stop positions
          for(i in 1:nrow(newdat_ch)){
            for(j in 1:nrow(mask_ch)){
              if((as.numeric(newdat_ch[i,]$Pos) >= as.numeric(mask_ch[j,]$V3)) & (as.numeric(newdat_ch[i,]$Pos) <= as.numeric(mask_ch[j,]$V4))){
                newdat_ch[i,1] <- NA
              }
            }
          }
        }
        colnames(newdat_ch)[1] <- colnames(dpgpvt[c])
        newdat_all <- rbind(newdat_all,newdat_ch[1])
      }
      newdat_all <- newdat_all[-c(1),]
      fulldata[c] <- newdat_all
    } else {
      ###for samples there are no region to mask, copy into new df
      fulldata[c] <- dpgpvt[c]
    }
  }
  
  fulldata <- as.data.frame(t(fulldata))
  fulldata <- cbind(dpgpv_new[1:16],fulldata)
  
  ##only keep SNPs with AF > 0.05 in at least two pop
  tmpnew <- fulldata[1:16] ### create tmp tmrix containing meta data
  for (i in 17:ncol(fulldata)){
    if (ncol(table(fulldata$Population,fulldata[,i])) > 1 & ncol(table(fulldata$Population,fulldata[,i])) < 4){
      if (length(table(fulldata$Population)) - sum(table(fulldata$Population,fulldata[,i])[,2]/rowSums(table(fulldata$Population,fulldata[,i])) < 0.05 | table(fulldata$Population,fulldata[,i])[,2]/rowSums(table(fulldata$Population,dpgpv_new[,i])) > 0.95, na.rm=T) > 1){
        tmpnew <- cbind(tmpnew,fulldata[i]) ### add filtered variants to tmp matrix
      }
    }
  }
  fulldata <- tmpnew ### make tmp matrix as original
  fulldata <- fulldata[-c(1:16)] ### remove meta data
  fulldata <- t(fulldata) ### transpose matrix

  
  write.table(fulldata, paste('header',gsub("ann.vcf","",transactFile),sep='_'), quote=F)
  write.table(fulldata, paste(gsub("ann.vcf","parse",transactFile),'geno',index,sep='_'), col.names = F, quote=F) ### annotations file
  write.table(ann_allsnps, paste(gsub("ann.vcf","parse",transactFile),'ann',index,sep='_'),row.names = F,quote=F) ### genotype file with metadata
  
  if (nrow(dataChunk) != chunkSize){
    print('Processed all files!')
    break}
  ### process next chunk of variants
  dataChunk <- read.table(con, nrows=chunkSize,stringsAsFactors = FALSE)
  full <- dataChunk
  if (index > 628) break ### stop when exceed number of lines in file

}
close(con)

