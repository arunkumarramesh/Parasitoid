setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(plyr)
library(ggh4x)
library(dplyr)
library(reshape2)
library(cowplot)

### now looking at fold change in lectin expression at 6h and 18h
injection_infection_qPCR <- read.csv(file = "injection_infection_qPCR.csv")
injection_infection_qPCR$Line <- as.factor(injection_infection_qPCR$Line)
injection_infection_qPCR$delta_ct_lectinnort <- rowMeans(injection_infection_qPCR[8:9]) -  injection_infection_qPCR$Lectin_no_RT
injection_infection_qPCR$Time <- gsub("6h","6",injection_infection_qPCR$Time)
injection_infection_qPCR$Time <- gsub("18h","18",injection_infection_qPCR$Time)
injection_infection_qPCR2 <- cbind(injection_infection_qPCR[2:5],injection_infection_qPCR[10],injection_infection_qPCR[14])
injection_infection_qPCR2$Time <- gsub("6","6 hpi",injection_infection_qPCR2$Time)
injection_infection_qPCR2$Time <- gsub("18","18 hpi",injection_infection_qPCR2$Time)
injection_infection_qPCR2 <- melt(injection_infection_qPCR2, id.vars = c(1:4))
injection_infection_qPCR2$variable <- mapvalues(injection_infection_qPCR2$variable, from = c("delta_ct_lectin","delta_ct_lectinnort"),to=c("Lectin-24A","NRT"))

anova(lm(data = injection_infection_qPCR2[injection_infection_qPCR2$variable %in% "NRT",], value ~ Line + Treatment + Time)) ### NRT values differ by time and line, cannot directly compare 6 and 18h, line may still be OK since fold difference is very large
anova(lm(data = injection_infection_qPCR, Lectin_no_RT ~ Line + Treatment + Time))
anova(lm(data = injection_infection_qPCR[injection_infection_qPCR$Treatment %in% "Control",], delta_ct_lectin ~ Line  + Time)) ### slight difference in delta ct in control samples. So 437 may have slightly higher basal expression.
anova(lm(data = injection_infection_qPCR[!injection_infection_qPCR$Treatment %in% "Control",], fold_induction_lectin ~ Line  + Treatment + Time)) ### line clearly important for fold induction of lectin.

injection_infection_qPCR2$Treatment <- gsub("Control","Uninfected",injection_infection_qPCR2$Treatment)

nrt <- injection_infection_qPCR2[injection_infection_qPCR2$variable %in% "NRT",] %>% 
  dplyr::group_by(Line,Treatment,Time)%>% 
  dplyr::summarise(nrt=max(value))
nrt$id <- paste(nrt$Line,nrt$Treatment,nrt$Time,sep="_")
nrt <- nrt[-c(1:3)]

lectin <- injection_infection_qPCR2[injection_infection_qPCR2$variable %in% "Lectin-24A",]
lectin$id <- paste(lectin$Line,lectin$Treatment,lectin$Time,sep="_")
lectin <- left_join(lectin,nrt,by="id")
lectin$detection <- lectin$value > lectin$nrt
lectin$detection <- gsub(T,"ADT",lectin$detection)
lectin$detection <- gsub(F,"BDT",lectin$detection)

labdat <- data.frame(Treatment=c("Uninfected","Infection","Injection"), Time=c("18 hpi","18 hpi","18 hpi"),lab=c("ADT","",""))
labdat2 <- data.frame(Treatment=c("Uninfected","Infection","Injection"), Time=c("18 hpi","18 hpi","18 hpi"),lab=c("BDT","",""))

deltacttimecourse <- ggplot(lectin,aes(y=value,x=as.factor(Line)))+
  geom_boxplot(outlier.shape = NA)+
  #  geom_jitter(aes(color=detection),width=0.35,height=0,size=1.5)+
  geom_jitter(aes(color=Line,shape=detection),width=0.35,height=0,size=1.5)+
  facet_nested(factor(Time,levels = c("6 hpi","18 hpi"))~factor(Treatment,levels = c("Uninfected","Infection","Injection")))+
  ylab(expression(paste(Delta,"CT")))+
  ylab(expression('Log'[2]*' relative '*italic('Lectin-24A')*' expression'))+
  scale_color_manual(values=c("red","blue"))+
  scale_shape_manual(values=c(19,1))+
  xlab("Drosophila Line")+
  labs(fill = "DGRP line")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")+
  labs(color = "")
 # geom_text(data=labdat,aes(x=1.5,y= -3.6,label=lab),color=c("orange"))+
 # geom_text(data=labdat2,aes(x=1.5,y= -4.6,label=lab),color=c("purple"))

deltacttimecourse

pdf(file="expression_time_course.pdf",height=3.5,width=3.5)
deltacttimecourse
dev.off()


### plot of fold change in lectin expression at 6 and 18h
foldchange <- ggplot(data = injection_infection_qPCR[!injection_infection_qPCR$Treatment %in% "Control",], aes(x = factor(Time,levels=c("6","18")), y = fold_induction_lectin, fill = Line))+
  geom_boxplot()+
  geom_jitter(width=0.25,height=0,size=1, alpha=0.3)+
  facet_nested(. ~ Treatment, scales = "free_y")+
  #facet_nested(. ~ Treatment + Line, scales = "free_y")+
  xlab("Hours Post Infection")+
  labs(fill = "DGRP line")+
  ylab(expression('Log'[2]*'FC '*italic('Lectin-24A')))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
foldchange

pdf(file="foldchange_time_course.pdf",height=2.6,width=6)
foldchange
dev.off()
