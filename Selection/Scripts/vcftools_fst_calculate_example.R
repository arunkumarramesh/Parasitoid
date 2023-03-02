library(dplyr)
library(methods)

### Aim of script is to do per site individual filtering and calculate Fst using vcftools
### This is an example script. The variants were split into 15 files and processed simultaneously for faster processing.

transactFile <- 'xaa_h'  ### variant matrix file, see DGN_analysis.sh for how they were created
index <- 0
chunkSize <- 1000  ### processing 1000 lines at a time to avoid holding too much memory
header <- read.table(file="dpgp_masked2_header",stringsAsFactors = FALSE) ### contains sample names
header <- as.character(header[1,])
con <- file(description=transactFile,open="r")
dataChunk <- read.table(con, nrows=chunkSize,stringsAsFactors = FALSE)
dpgp_masked2 <- dataChunk
pop_sample <- read.table("pop_sample",header = T)  ### file contains region and population designations for samples
fst_vals <- "" 
### repeat over 1000 chunks of variants
repeat {
  index <- index + 1
  print(paste('Processing rows:', index * chunkSize))
  colnames(dpgp_masked2) <- header

  ### now looping over each site
  for (i in 1:nrow(dpgp_masked2)){
    
    gt <- as.data.frame(t(dpgp_masked2[i,10:ncol(dpgp_masked2)]))  ### just gets the genotypes for all samples
    gt <- as.data.frame(gt[!gt == "./.",])  ## remove samples where genotypes are not called
    gt$Sample <- rownames(gt)
    gt <- left_join(gt,pop_sample,by="Sample") ### add population information for samples
    sample_count_per_pop <- table(as.character(gt$Population))  ### count number of samples per population
    sample_count_per_pop <- sample_count_per_pop[sample_count_per_pop > 3] ### only keep populations with four or more samples
    if(length(names(sample_count_per_pop)) > 24){ ### only calculate fst if there are at least 25 populations
      line <- dpgp_masked2[i,]
      
      ###create a mock per site vcf file
      write.table(line,file="xaa.vcf",quote=F, row.names = F) 
      system("sed -i 's/CHROM/#CHROM/' xaa.vcf")
      system("cat tester_wh_header2 xaa.vcf >xaa2.vcf")
      system("sed -i 's/ /\t/g'  xaa2.vcf")
      ### run vcf tools on mock vcf file and capture per site fst output 
      a <- system(paste("vcftools --vcf xaa2.vcf --weir-fst-pop ",gsub(",","_POP --weir-fst-pop ", toString(names(sample_count_per_pop))), "_POP --stdout",sep=""),intern = T)
      fst_vals <- rbind(fst_vals,a)  ### add fsts to a common dataframe
    }
  }
  
  if (nrow(dataChunk) != chunkSize){
    print('Processed all files!')
    break}
  ### process next chunk of variants
  dataChunk <- read.table(con, nrows=chunkSize,stringsAsFactors = FALSE)
  dpgp_masked2 <- dataChunk
  if (index > 4340) break ### stop when exceed number of lines in file
  
}
#close(con)

fst_vals <- as.data.frame(fst_vals[-c(1),])
rownames(fst_vals) <- seq(1,nrow(fst_vals))  ### give arbitrary rownames
library(tidyr)

fst_vals <- separate(fst_vals,V2,sep="\t",into=c("CHR","POS","WEIR_AND_COCKERHAM_FST"))[2:4] ### separate fst values into columns
write.table(fst_vals,file="dpgp_vcftools_POP.weir.fst_xaa",row.names = F,col.names=F,quote = F) ### write output



