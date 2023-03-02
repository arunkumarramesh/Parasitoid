library(ggplot2)
library(plyr)
library(Hmisc)
library(lme4)
library(multcomp)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data = read.csv('DGRP_Dominance_Results_Combined.csv', header=T, fileEncoding="UTF-8-BOM")
head(data)
data$ID=as.factor(paste(data$Female,'x',data$Male))

### all the different DGRP crosses that were analysed
data$ID=factor(data$ID, levels=c("437 x 437", "892 x 892", "437 x 892","892 x 437", 
                                 "566 x 566", "566 x 892", "892 x 566", 
                                 "589 x 589", "748 x 748", "589 x 748", "748 x 589"))
#remove unparasitised larvae
data = data[which(!data$Resistance_Phenotype=="<NA>"),]

#calculate Resistance as a proportion of resistant larvae per sample
data_summary = ddply(data,.(Cross_Type,ID,Sex),  summarise,
                     Proportion_Resistance = sum(Resistance_Phenotype=='Resistant')/(sum(Resistance_Phenotype=='Resistant')+sum(Resistance_Phenotype=='Susceptible')),
                     Lower_Confidence_Interval_Value = binconf(x=sum(Resistance_Phenotype=='Resistant'), n=(sum(Resistance_Phenotype=='Resistant')+sum(Resistance_Phenotype=='Susceptible')), alpha=0.05, method='wilson')[2],
                     Upper_Confidence_Interval_Value = binconf(x=sum(Resistance_Phenotype=='Resistant'), n=(sum(Resistance_Phenotype=='Resistant')+sum(Resistance_Phenotype=='Susceptible')), alpha=0.05, method='wilson')[3],
                     total=(sum(Resistance_Phenotype=='Resistant')+sum(Resistance_Phenotype=='Susceptible')))




######### To plot data
data2=data_summary[data_summary$ID %in% c("437 x 437", "892 x 892", "437 x 892",  "892 x 437"),]

p1 = ggplot(data2, aes(x=ID, y=Proportion_Resistance,ymin=Lower_Confidence_Interval_Value, ymax=Upper_Confidence_Interval_Value,fill=Sex))+
  geom_bar(stat='identity', 
                    width=0.8,position=position_dodge(width=0.85))+
  geom_errorbar( width=.4,position=position_dodge(width=0.85))+ 
  
  ylim(0,1)+
  guides(fill=FALSE)+
  theme_bw()+
  scale_fill_manual(values=c("Male"="lightblue3","Female"="lightsalmon"))+
  theme(panel.grid = element_blank())+
  theme(panel.border =element_blank())+
  theme(axis.line = element_line(colour="black"))+
    theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  ylab('Proportion Wasps Encapsulated\n')+
  xlab('\n\nDrosophila Genotype')+
  annotate('text', x = 1.2, y = 0.93, label = "\u2640 larvae",colour="lightsalmon",size=5)+
  annotate('text', x = 1.2, y = 1, label = "\u2642 larvae",colour="lightblue3",size=5)+
  geom_text(aes(x=ID,y=Upper_Confidence_Interval_Value+0.035, label=total), size=4.5,colour = "grey30",position=position_dodge(width=0.85))+
  theme(axis.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 17))  + 
  theme(axis.text.x = element_text(colour ="white")) +
  theme(axis.ticks = element_blank())
  


#use this to embed male and female symbols in pdf. Think only works on mac (?)
p1
quartz.save(type = "pdf", file="autosomal_dominant.pdf",height=5,width=3.6)



#cartoon of chromosomes

space=0.2
s=0.28

ax=c(1,2+s,4,5+s,7,8+s,10,11+s)
Fem_ax=ax[c(1,3,5,7)]

x=c(4,4,4,
    8,8,8,
    4,8,4,
    4,8,8)
aut=c(4,4,4,4,
      8,8,8,8,
      4,8,4,8,
      4,8,4,8)


pdf(file='chromosome cartoon2.pdf',height=4,width=7.3)
par(bg = "white") 
plot(NULL, xlim=c(0,13), ylim=c(-1,3), ylab="y label", xlab="x lablel")


#X chromsomes
ax2=sort(c(ax,Fem_ax+space))
for (i in 1:length(ax2)){
  colour=ifelse(x[i]==4,"red","blue")
  lines(x=c(ax2[i],ax2[i]),y=c(2,2.5),col=colour,lwd=5)
}
#2 chromsomes
ax2=sort(c(ax,ax+space))
for (i in 1:length(ax2)){
  colour=ifelse(aut[i]==4,"red","blue")
  lines(x=c(ax2[i],ax2[i]),y=c(1,1.7),col=colour,lwd=5)
}
#3 chromsomes
ax2=sort(c(ax,ax+space))
for (i in 1:length(ax2)){
  colour=ifelse(aut[i]==4,"red","blue")
  lines(x=c(ax2[i],ax2[i]),y=c(0,0.7),col=colour,lwd=5)
}



text(0,2.3,'X',cex=3)
text(0,1.32,'2',cex=3)
text(0,0.35,'3',cex=3)
text(12.7,2,'437',col='red',cex=2.5)
text(12.7,0.5,'892',col='blue',cex=2.5)

dev.off()
### Stats to test effect of cross on encapsulation
data_comb1 <- data[data$ID %in% c("437 x 437", "892 x 892", "437 x 892", "892 x 437"),]
data_comb1 <- droplevels(data_comb1)
data_comb1$Replica_ID = paste(data_comb1$Replica,data_comb1$Dissection_Date, sep=',' )
glmer_model_comb<-glmer(data=data_comb1, as.factor(Resistance_Phenotype)~ID+Sex+(1|Replica_ID)+(1|Dissection_Date),family="binomial")
summary(glmer_model_comb)
#dissection date explains no variance- remove from model
glmer_model_comb1<-glmer(data=data_comb1, as.factor(Resistance_Phenotype)~ID+Sex+(1|Replica_ID),family="binomial")

glmer_model_comb1_2 <-glmer(data=data_comb1, as.factor(Resistance_Phenotype)~ID+(1|Replica_ID),family="binomial")

anova(glmer_model_comb1,glmer_model_comb1_2)

summary(glmer_model_comb1)
summary(glht(glmer_model_comb1, linfct = mcp(ID = "Tukey")), test = adjusted("holm"))

