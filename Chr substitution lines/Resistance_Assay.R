library(car)
library(ggplot2)
library(plyr)
library(Hmisc)
library(lme4)
library(scales)
#set wd

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dataset = read.csv("2018_08_09_Resistance_Assay_Combo1_2.csv", header=T, fileEncoding="UTF-8-BOM")
head(dataset)


##remove C2xC3
dataset <- dataset[!dataset$Stock %in% "C2xC3",]

#change order of factor levels in Stock
dataset$Stock = factor(dataset$Stock, levels=c("DGRP437", "DGRP892", 
                                               "C3_X(DGRP892)_II(DGRP892)_III(DGRP437)",
                                               "C5_X(DGRP437)_II(DGRP892)_III(DGRP437)",
                                               "C2_X(DGRP892)_II(DGRP437)_III(DGRP892)",
                                               "C6_X(DGRP437)_II(DGRP437)_III(DGRP892)",
                                               "C4_X(DGRP892)_II(DGRP437)_III(DGRP437)"))



#create a new column with a unique ID for each sample
dataset$Sample_ID = paste(dataset$Stock, dataset$Treatment, dataset$Replica, dataset$Dissection_Date, sep= "_")

#count the number of observation in each unique Sample ID and create a data frame with it
data_summary_counts = plyr::count(dataset, c("Sample_ID"))

#filter the data keeping only the ones with 10 observation and remove NA
dataset <- dataset[!is.na(dataset$Resistance_Phenotype),]

dataset = dataset[dataset$Sample_ID %in% data_summary_counts[data_summary_counts$freq > 10,]$Sample_ID,]

#calculate Resistance as a proportion of resistant larvae per sample
data_summary = ddply(dataset, .(Researcher, Replica, Stock, Dissection_Date), summarise,
                     Proportion_Resistance= sum(Resistance_Phenotype== "Resistant")/(sum(Resistance_Phenotype== "Resistant")+sum(Resistance_Phenotype== "Susceptible")),
                     tot = sum(Resistance_Phenotype== "Resistant")+sum(Resistance_Phenotype== "Susceptible")) 

#Define standard error function
se <- function(x) sd(x)/sqrt(length(x))

#calculate Resistance as a proportion of resistant larvae per stock
data_summary_2 = ddply(data_summary, .(Stock), summarise, 
                       Mean_Proportion_Resistance = mean(Proportion_Resistance),
                       se = se(Proportion_Resistance),
                       Total = sum(tot))
                       
#plot bars with mean resistance per stock overlayed with individual samples per stock - color by combination and shape by researcher                     

palette = c('green0', 'green10', 'green20', 'green30', 'green40', 'green50','green60','green70')

p1 = ggplot(data_summary_2,aes(x=Stock,y=Mean_Proportion_Resistance))+
  geom_bar(stat="identity",fill='grey',width=0.8)+
  geom_errorbar(data=data_summary_2, aes(ymin=Mean_Proportion_Resistance-se, ymax=Mean_Proportion_Resistance+se), width=.4)+ 
  geom_text(aes(y= Mean_Proportion_Resistance+se+0.05, label=Total), size=4.5, angle=0)+
  geom_jitter(data=data_summary,aes(x=Stock,y=Proportion_Resistance),width=0.25,height=0,size=0.8, alpha=0.3)+
  theme_bw()+
  annotate("text", x = 1, y = 1, label = "a")+
  annotate("text", x = 3, y = 1, label = "b")+
  annotate("text", x = 5.5, y = 1, label = "c")+
  annotate("text", x = 7, y = 1, label = "d")+
  annotate("segment", x = 0.7, xend = 1.3, y = 0.98, yend = 0.98)+
  annotate("segment", x = 1.7, xend = 4.3, y = 0.98, yend = 0.98)+
  annotate("segment", x = 4.7, xend = 6.3, y = 0.98, yend = 0.98)+
  annotate("segment", x = 6.7, xend = 7.3, y = 0.98, yend = 0.98)+
  ylim(0,1)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(size=10), legend.position="none") +  
  ylab('Proportion Wasps Encapsulated\n') +
  xlab('\n\n\nDrosophila Genotype') +
  scale_x_discrete(labels = c("437", "892",
                              expression(X[892]*II[892]*III[437]),
                              expression(X[437]*II[892]*III[437]),
                              expression(X[892]*II[437]*III[892]),
                              expression(X[437]*II[437]*III[892]),
                              expression(X[892]*II[437]*III[437]),
                              expression(X[892]*II[437/892]*III[437/892])))+ 
  theme_classic() +
  theme(axis.text=element_text(size=15)) +
  theme(axis.title = element_text(size = 16))  + 
  theme(axis.text.x = element_text(colour ="white"))



p1

pdf(file="chromosome_2.pdf",height=5.5,width=3.5)

  p1 

dev.off()

svg(file="chromosome_2.svg",height=5.5,width=3.5)

p1 

dev.off()

#cartoon of chromosomes
x=c(4,8,8,4,8,4,8)
two=c(4,8,8,8,4,4,4)
three=c(4,8,4,4,8,8,4)


pdf(file='chromosome cartoon2.pdf',height=4,width=7)
par(bg = "white")
plot(NULL, xlim=c(0,9), ylim=c(-1,3), ylab="y label", xlab="x lablel")

for (i in 1:7){
  colour=ifelse(x[i]==4,"red","blue")
  lines(x=c(i,i),y=c(2,2.4),col=colour,lwd=6)
}

for (i in 1:7){
  colour=ifelse(x[i]==4,"red","blue")
  lines(x=c(i+0.2,i+0.2),y=c(2,2.4),col=alpha(colour,0.3),lwd=6)
}

for (i in 1:7){
  colour=ifelse(two[i]==4,"red","blue")
  lines(x=c(i,i),y=c(1,1.6),col=colour,lwd=6)
}

for (i in 1:7){
  colour=ifelse(two[i]==4,"red","blue")
  lines(x=c(i+0.2,i+0.2),y=c(1,1.6),col=colour,lwd=6)
}

for (i in 1:7){
  colour=ifelse(three[i]==4,"red","blue")
  lines(x=c(i,i),y=c(0,0.6),col=colour,lwd=6)
}

for (i in 1:7){
  colour=ifelse(three[i]==4,"red","blue")
  lines(x=c(i+0.2,i+0.2),y=c(0,0.6),col=colour,lwd=6)
}

text(0.2,2.1,'X',cex=3)
text(0.2,1.2,'2',cex=3)
text(0.2,0.3,'3',cex=3)
text(8.5,2,'437',col='red',cex=2.5)
text(8.5,1,'892',col='blue',cex=2.5)

dev.off()

######################################################


### Stats to test effect of chromosome on encapsulation
library(multcomp)
data_comb1 = dataset
data_comb1$Replica_ID = paste(data_comb1$Replica,data_comb1$Dissection_Date, sep=',' )
glmer_model_comb<-glmer(data=data_comb1, as.factor(Resistance_Phenotype)~Stock+(1|Replica_ID)+(1|Dissection_Date),family="binomial")
summary(glmer_model_comb)
#dissection date explains no variance- remove from model,NOT TRUE ANY LONGER!!!
#glmer_model_comb1<-glmer(data=data_comb1, as.factor(Resistance_Phenotype)~Stock+(1|Replica_ID),family="binomial")
#summary(glmer_model_comb1)

summary(glht(glmer_model_comb, linfct = mcp(Stock = "Tukey")), test = adjusted("holm"))


#test for epistasis

dataset$chrII=vapply(strsplit(as.character(dataset$Stock),"_"), `[`, 3, FUN.VALUE=character(1))
dataset$chrIII=vapply(strsplit(as.character(dataset$Stock),"_"), `[`, 4, FUN.VALUE=character(1))

dataset$chrII[dataset$Stock=="DGRP437"]="II(DGRP437)"
dataset$chrIII[dataset$Stock=="DGRP437"]="III(DGRP437)"
dataset$chrII[dataset$Stock=="DGRP892"]="II(DGRP892)"
dataset$chrIII[dataset$Stock=="DGRP892"]="III(DGRP892)"

data_summary_chromosome = ddply(dataset, .(Sample_ID,chrII,chrIII), summarise,
                     resist= sum(Resistance_Phenotype== "Resistant"),
                     suscept=sum(Resistance_Phenotype== "Susceptible")) 
data_summary_chromosome$prop_resist=data_summary_chromosome$resist/(data_summary_chromosome$resist+data_summary_chromosome$suscept)



data_summary_chromosome$replica=as.factor(1:nrow(data_summary_chromosome))
model2<-glmer(data=data_summary_chromosome, cbind(resist,suscept)~chrII*chrIII+(1|replica),family="binomial")
Anova(model2)


#sanity check: from single chromosome effects we get the combined effect
#note this is the multiplicative model of epistasis
#https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347(01)02213-3#GLOSS12

data_summary_chromosome$II_III=paste(data_summary_chromosome$chrII,data_summary_chromosome$chrIII)
table(data_summary_chromosome$II_III)
x=tapply(data_summary_chromosome$prop_resist,data_summary_chromosome$II_III,mean)
x=x/x[1]
x
#predict the effect of double susceptible allele from single chromosomes
x[2]*x[3]
#uncannily like the observed value!: --NOT TRUE ANY LONGER
x[4]


##compare C2xC3 (het) with C4_X(DGRP892)_II(DGRP437)_III(DGRP437)
dataset = read.csv("2018_08_09_Resistance_Assay_Combo1_2.csv", header=T, fileEncoding="UTF-8-BOM")
dataset2 <- dataset[dataset$Stock %in% c("C4_X(DGRP892)_II(DGRP437)_III(DGRP437)","C2xC3"),]
fisher.test(table(dataset2$Stock,dataset2$Resistance_Phenotype))