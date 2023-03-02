setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(lmerTest)
library(ggpubr)
library(reshape2)
library(car)

## script is to analyse the encapsulation rates in adult DGRP lines and their association with Lectin-24A expression

## read in file containing encapsulation rates in adult DGRPs
DGRP_screen <- read.table(file="DGRP_Lboulardi_Screen_Results_Filtered.csv",sep = ";",header = T)
DGRP_screen <- DGRP_screen[DGRP_screen$Wasp_Strain %in% "G486",] ## only keep rates for G486 infection 
DGRP_screen$DGRP <- as.character(DGRP_screen$DGRP)
DGRP_screen$highlight <- "others"
DGRP_screen[DGRP_screen$DGRP %in% "437",]$highlight <- "437"
DGRP_screen[DGRP_screen$DGRP %in% "892",]$highlight <- "892"
## get  mean encapsulation rate per line, used for ordering DGRP lines
DGRP_screen_mean <- DGRP_screen %>%
  dplyr::group_by(DGRP) %>%
  dplyr::summarise(Encapsulation_Proportion=mean(Encapsulation_Proportion))
DGRP_screen_mean <- DGRP_screen_mean[order(DGRP_screen_mean$Encapsulation_Proportion),]

## getting means for online DGRP association analysis
DGRP_screen_sum_mean <- DGRP_screen
DGRP_screen_sum_mean$sus <- (DGRP_screen_sum_mean$Capsules/DGRP_screen_sum_mean$Encapsulation_Proportion)-DGRP_screen_sum_mean$Capsules
DGRP_screen_sum_mean[is.na(DGRP_screen_sum_mean$sus),]$sus <- 0
DGRP_screen_sum_mean <- DGRP_screen_sum_mean %>%
  dplyr::group_by(DGRP) %>%
  dplyr::summarise(Capsules=sum(Capsules),sus=sum(sus))
DGRP_screen_sum_mean$Encapsulation_Proportion <- DGRP_screen_sum_mean$Capsules/(DGRP_screen_sum_mean$Capsules + DGRP_screen_sum_mean$sus)
DGRP_screen_sum_mean <- DGRP_screen_sum_mean[!is.na(DGRP_screen_sum_mean$Encapsulation_Proportion),]
write.table(DGRP_screen_sum_mean[c(1,4)],file="DGRP_screen_mean.csv",row.names = F, col.names = F, quote = F,sep = ",")

## plot ordered boxplot of encapsulation rates in the DGRPs
dgrpscreenplot <- ggplot(data=DGRP_screen,aes(x=factor(DGRP,levels=DGRP_screen_mean$DGRP),y=Encapsulation_Proportion,fill=highlight))+
  geom_boxplot()+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1))+
  ylab("Proportion wasps encapsulated") +
  xlab("DGRP line") +
  scale_fill_manual(values=c("red", "blue", "white")) +
  theme(legend.title= element_blank()) +
  theme(legend.position = c(0.1, 0.7))

pdf(file="dgrpscreenplot.pdf",height=3.5,width=16)
dgrpscreenplot
dev.off()

## plot of infection rate with encapsulation rate
ggscatter(DGRP_screen, x = "Infection_Proportion", y = "Encapsulation_Proportion",add = "reg.line") + 
  stat_cor(method = "pearson", label.x = 0.7, label.y = 1) +
  ylab("Proportion wasps encapsulated") +
  xlab(expression('Infection Proportion'))

## file contains the lectin-24a upstream indel haplotype calls for the DGRPs, obtained from dgn_lectin.R
upstreamindelsgt <- read.csv(file="upstreamindelsgt2.txt")
colnames(upstreamindelsgt) <- c("DGRP","indel7bp","indel8bp","indel21bp")
upstreamindelsgt$DGRP <- gsub("RAL-","",upstreamindelsgt$DGRP)
## only keep DGRP lines filly genotypes from all 3 indels
upstreamindelsgt <- upstreamindelsgt[which(upstreamindelsgt$indel21bp %in% c(0,1)),]
upstreamindelsgt <- upstreamindelsgt[which(upstreamindelsgt$indel8bp %in% c(0,1)),]
upstreamindelsgt <- upstreamindelsgt[which(upstreamindelsgt$indel7bp %in% c(0,1)),]
## convert 0,1 calls to insertion and deletion
upstreamindelsgt$indel21bp <- gsub("1","D",upstreamindelsgt$indel21bp)
upstreamindelsgt$indel21bp <- gsub("0","I",upstreamindelsgt$indel21bp)
upstreamindelsgt$indel8bp <- gsub("1","D",upstreamindelsgt$indel8bp)
upstreamindelsgt$indel8bp <- gsub("0","I",upstreamindelsgt$indel8bp)
upstreamindelsgt$indel7bp <- gsub("1","D",upstreamindelsgt$indel7bp)
upstreamindelsgt$indel7bp <- gsub("0","I",upstreamindelsgt$indel7bp)
## generate haplotypes for upstream indels
upstreamindelsgt$haplotype2 <- paste(upstreamindelsgt$indel7bp,upstreamindelsgt$indel8bp,upstreamindelsgt$indel21bp,sep=" ")
## grouping haplotypes based whether they can express lectin-24a, based on ASE results
upstreamindelsgt$haplotype  <- ""
upstreamindelsgt[upstreamindelsgt$haplotype2 %in% c("I I I","I D I","D I I"),]$haplotype <- "Expresses Lectin-24A"
upstreamindelsgt[upstreamindelsgt$haplotype2 %in% c("D D I","D D D"),]$haplotype <- "Does not express Lectin-24A"

table(upstreamindelsgt$haplotype)

## merge data containing encapsulation rates and upstream indel haplotypes
DGRP_screen_h <- left_join(DGRP_screen,upstreamindelsgt,by="DGRP")
DGRP_screen_h <- DGRP_screen_h[!is.na(DGRP_screen_h$haplotype),] ## only keep lines where indel haplotypes are known
t.test(data=DGRP_screen_h,Encapsulation_Proportion~haplotype)

## to see how many lines of each lectin-24a expression haplotype exist
DGRP_screen_h %>%
  dplyr::group_by(DGRP) %>%
  dplyr::summarise(haplotype=haplotype) %>% 
  dplyr::distinct() %>%
  dplyr::group_by(haplotype) %>%
  dplyr::summarise(table(haplotype))

## plot containing violin plot of encapsualtion rates by upstream indel haplotype
pdf(file="indel_haplotype_encapsulation.pdf",height=3.5,width=3.8)
ggplot(data=DGRP_screen_h,aes(x=haplotype,y=Encapsulation_Proportion))+
  geom_violin() +
  geom_jitter(width=0.2,height=0,size=0.9,shape=1,color="grey") +
  ylab("Proportion wasps encapsulated") +
  xlab("Upstream indel haplotype") +
  annotate("text",x=1.5,y=1.06,label=expression(italic('t')*'=-3.07, '*italic('df')*'=230, '*italic('p')*'=0.002'), size = 5) +
  theme_classic2() +
  geom_segment(aes(x = 0.9, y = 1, xend = 2.1, yend = 1)) +
  scale_x_discrete(labels=c("Does not express","Express")) +
  theme(axis.text=element_text(size=14)) +
  theme(axis.text.x=element_text(color=NULL)) +
  theme(axis.title = element_text(size = 15))   
dev.off()

## gets the female and male lines means for the DGRP transcriptomes from 10.1101/gr.257592.119
#female.dgrp.exp <- read.table(file="female.genvar.line.mean.txt",header=T)
#male.dgrp.exp <- read.table(file="male.genvar.line.mean.txt",header=T)
## extract lectin expression values only
#female.dgrp.exp <- female.dgrp.exp[female.dgrp.exp$GENE %in% "FBgn0040104",]
#male.dgrp.exp <- male.dgrp.exp[male.dgrp.exp$GENE %in% "FBgn0040104",]

## combine male and female values
all <- cbind(t(female.dgrp.exp),t(male.dgrp.exp))
all <- all[-c(1),]
colnames(all) <- c("female","male")
all <- as.data.frame(all)
all$DGRP <- gsub("_F","",rownames(all))
all$DGRP <- gsub("X","",all$DGRP)
all <- melt(all,id.var="DGRP") ## make into 2 columns with sex as a seperate column
colnames(all)[2:3] <- c("sex","log2FPKM")
all <- left_join(all,DGRP_screen,by="DGRP") ## combine with encapsulation rates
all$log2FPKM <- as.numeric(all$log2FPKM)
all <- left_join(all,upstreamindelsgt,by="DGRP")  ## combine with upstream indel haplotype calls

all_expression <- all[c(1:3,17)] %>% distinct() ## only keep unique expression estimates, duplicated caused by previous merging step
## first test for differences in females with posthoc test to identify differences between haplotypes
model <- lm(data=all_expression[all_expression$sex %in% "female",],log2FPKM ~ as.factor(haplotype2))
tukey.test <- TukeyHSD(aov(model),  conf.level=0.95)
tukey.test$`as.factor(haplotype2)`
plot(tukey.test , las=1 , col="brown") # Tukey test representation

## next test for differences in males with posthoc test to identify differences between haplotypes
model <- lm(data=all_expression[all_expression$sex %in% "male",],log2FPKM ~ as.factor(haplotype2))
tukey.test <- TukeyHSD(aov(model),  conf.level=0.95)
tukey.test$`as.factor(haplotype2)`
plot(tukey.test , las=1 , col="brown") # Tukey test representation

## plotting upstream indel haplotype by adult expression (log2FPKM)
pdf(file="haplotype2_adult_FPKM.pdf",height=4,width=6)
ggplot(data=all_expression[!is.na(all_expression$haplotype2),],aes(x=haplotype2,y=log2FPKM)) +
  facet_grid(~sex) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Upstream indel haplotype") +
  ylab(expression('Average log'[2] *italic(' Lectin-24A')*' FPKM')) +
  geom_jitter(width=0.35,height=0,size=0.9, shape=1) +
  geom_segment(aes(x=1,xend=2,y=3,yend=3)) +
  geom_segment(aes(x=3,xend=5,y=5,yend=5)) +
  geom_text(aes(x=1.5,y=3.5,label="a")) +
  geom_text(aes(x=4,y=5.5,label="b"))
dev.off()

## plotting adult expression (log2FPKM) in 437 and 892
pdf(file="adult_FPKM_encapsulation_437_892.pdf",height=4,width=4)
ggplot(data=all[all$DGRP %in% c("437","892"),],aes(x=DGRP,y=log2FPKM,fill=DGRP))+
  geom_col() +
  facet_grid(~sex) +
  scale_fill_manual(values=c("red","blue")) +
  ylab(expression('Average log'[2] *italic(' Lectin-24A')*' FPKM')) +
  xlab("DGRP line") +
  theme(legend.position="none")
dev.off()



