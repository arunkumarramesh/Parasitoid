setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(reshape2)
library(ggplot2)
library(dplyr)
library(car)

## script is to plot qPCR values for the 8 candidate genes in chr2 qtl

candidates <- read.csv(file="candidate 8 genes 437 vs 892 011122.csv") ## read in qPCR values
candidates <- candidates[-c(4:14)] ## remove ct columns as not needed, only using delta ct against rpl32
candidates <- melt(candidates,id.vars=c(1:3))  ## make genes as a column
## change names
candidates$variable <- gsub("_Dct","",candidates$variable)
candidates$variable <- gsub("Bark","bark",candidates$variable)
candidates$variable <- gsub("Capu","capu",candidates$variable)
candidates$variable <- gsub("SrCIV3","Sr-CIV",candidates$variable)
candidates$variable <- gsub("Lectin","Lectin-24A",candidates$variable)
candidates$Treatment <- gsub("con","Uninfected",candidates$Treatment)
candidates$Treatment <- gsub("G486","Infected",candidates$Treatment)
candidates$Hours <- gsub("18h","18 hpi",candidates$Hours)
candidates$Hours <- gsub("6h","6 hpi",candidates$Hours)
candidates$Line <- as.character(candidates$Line)

## get lectin specific values
lectin <- candidates[candidates$variable %in% "Lectin-24A",]

deltacttimecourse <- ggplot(lectin,aes(y=value,x=as.factor(Line)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color=Line),width=0.35,height=0,size=1.5)+
  facet_nested(factor(Hours,levels = c("6 hpi","18 hpi"))~factor(Treatment,levels = c("Uninfected","Infected")))+
  ylab(expression(paste(Delta,"CT")))+
  ylab(expression('Log'[2]*' relative '*italic('Lectin-24A')*' expression'))+
  scale_color_manual(values=c("red","blue"))+
  scale_shape_manual(values=c(19,1))+
  xlab("Drosophila Line")+
  labs(fill = "DGRP line")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+
  labs(color = "")+
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 13)) +
  theme(strip.text = element_text(size = 13))

deltacttimecourse

pdf(file="expression_time_course.pdf",height=3.5,width=3)
deltacttimecourse
dev.off()

Anova(lm(data=lectin,formula=value ~Line + Hours+ Treatment))

## now do a test for interaction between DGRP line and treatment on expression
pvals <- ""
for (g in 1:length(table(candidates$variable))){
  for (h in 1:length(table(candidates$Hours))){
    pvals <- rbind(pvals,cbind(names(table(candidates$variable))[g],names(table(candidates$Hours))[h],Anova(lm(data=candidates[candidates$variable %in%  names(table(candidates$variable))[g] & candidates$Hours %in% names(table(candidates$Hours))[h],],formula=value ~Line * Treatment))[3,2:4]))
  }
}
pvals <- pvals[-c(1),]
colnames(pvals)  <- c("variable","Hours","df","F","Pvalue")
pvals$Pvalue <- as.numeric(pvals$Pvalue)
pvals$F <- as.numeric(pvals$F)
colnames(pvals)[1:2] <- c("Gene","Time post infection") 
write.csv(pvals,file="intearction_pval.csv",row.names = F) ## write interaction pvalues to file

## calculate mean and standard error for detla ct
candidates2 <- candidates %>% 
  dplyr::group_by(Line,Treatment,Hours,variable) %>% 
  dplyr::summarise(mean=mean(value),se=(sd(value)/sqrt(4)))

## generate plot
pdf(file="candidate_genes.pdf",height=4,width=8)
ggplot(data=candidates2, aes(x=factor(Treatment,levels=c("Uninfected","Infected")), y=mean, color=Line,group=Line)) + 
  geom_point()+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(0.05)) +
  facet_grid(factor(Hours,levels=c("6 hpi","18 hpi"))~variable) +
  xlab("Treatment") +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1)) +
  ylab(expression('Average log'[2]*' relative expression +/- standard error')) +
  scale_color_discrete(name = "DGRP line") + 
  theme(strip.text.x = element_text(face = "italic"))
dev.off()



