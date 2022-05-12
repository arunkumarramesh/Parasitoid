library(car)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(lme4)
library(binom)
#set wd

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dataset = read.csv("G486_Encapsulation_Assays_Raw_Data.csv", header=T, fileEncoding="UTF-8-BOM")
head(dataset)


##remove C2xC3
dataset <- dataset[dataset$DGRP %in% c("437","892"),]

#calculate Resistance as a proportion of resistant larvae per sample
datasummary <- dataset %>% 
  dplyr::group_by(DGRP)%>% 
  dplyr::summarise(Capsules=sum(Capsules),Non_Capsules =sum(Non_Capsules ))
  
#Define binomial CIs
binomial_CIs <- binom.confint(x = datasummary$Capsules, n = as.numeric(rowSums(datasummary[2:3])), methods = "prop.test")
datasummary <- cbind(datasummary,binomial_CIs)

#plot bars with mean resistance per stock overlayed with individual samples per stock - color by combination and shape by researcher                     

p1 = ggplot(datasummary,aes(x=as.factor(DGRP),y=mean,fill=as.factor(DGRP)))+
  geom_bar(stat="identity",width=0.8)+
  scale_fill_manual(values=c("red","blue"))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.4)+ 
  geom_text(aes(y= upper+0.1, label=paste("n=",n,sep="")), size=3, angle=0)+
  geom_jitter(data=dataset,aes(x=as.factor(DGRP),y=Encapsulation_Ratio),width=0.25,height=0,size=1, alpha=0.3)+
  theme_bw()+
  annotate("text", x = 1.5, y = 0.6, label = "****")+
  ylim(0,0.6)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(size=10), legend.position="none") +
  ylab('Proportion Wasps Encapsulated') +
  xlab('Genotype') +
  theme_classic()+
  theme(legend.position = "none")

p1

glmer_model_comb<-glmer(data=dataset, cbind(Capsules,Non_Capsules)~DGRP+(1|Replica),family="binomial")
summary(glmer_model_comb)

pdf(file="dgrp_encapsulation.pdf",height=2.67,width=1.6)
p1 
dev.off()
