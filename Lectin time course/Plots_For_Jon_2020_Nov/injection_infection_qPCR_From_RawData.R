#load required packages
library(ggplot2)
library(plyr)
library(multcomp)
library(Hmisc)
library(reshape2)
library(stringr)

#load dataset 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
list.files()
data = read.table('injection_infection_qPCR.csv', sep=',',dec='.', head=T, check.names = F)
head(data)

data$Time = factor(data$Time, levels=c('6h', '18h'))


#Calculate DeltaCt
data$DeltaCtLectin = data$lectin_CT - data$rpl_CT
data$DeltaCtIM1 = data$IM1_CT - data$rpl_CT
data$DeltaCtCG33462 = data$CG33462_CT - data$rpl_CT


#Summarise data from control treatment: Calculate average of Ct and calculate Delta Ct. Test data will be normalised to this value
data_control = data[data$Treatment=='Control',]
data_control_summary = ddply(data_control, .(Line, Treatment, Time), summarise, 
                            mean_Rpl = mean(rpl_CT),
                            mean_Lectin = mean(lectin_CT), 
                            mean_IM1 = mean(IM1_CT),
                            mean_CG33462 = mean(CG33462_CT))

data_control_summary$mean_DeltaCtLectin = data_control_summary$mean_Lectin - data_control_summary$mean_Rpl
data_control_summary$mean_DeltaCtIM1 = data_control_summary$mean_IM1 - data_control_summary$mean_Rpl
data_control_summary$mean_DeltaCtCG33462 = data_control_summary$mean_CG33462 - data_control_summary$mean_Rpl


#Create dataset with test data and value of buffer injections for normalisation, per day of injecion
data_infection = data[!data$Treatment=='Control',]

mdata= merge(data_infection, data_control_summary, by=c('Line', 'Time'))

mdata$Fold_Change_Lectin = 2^(-(mdata$DeltaCtLectin - mdata$mean_DeltaCtLectin))
mdata$Fold_Change_IM1 = 2^(-(mdata$DeltaCtIM1 - mdata$mean_DeltaCtIM1))
mdata$Fold_Change_CG33462 = 2^(-(mdata$DeltaCtCG33462 - mdata$mean_DeltaCtCG33462))

subset= mdata[,c(1,2,3,19,20,21)]
msubset = melt(subset, id.vars=c('Line', 'Time', 'Treatment.x'))


###Plot boxplots
stage.labs <- c(
  `Fold_Change_Lectin` = "lectin-24A",
  `Fold_Change_IM1` = "Bomanin Short 1",
  `Fold_Change_CG33462` = "CG33462", 
  `437` = "DGRP-437", 
  `892` = "DGRP-892")

p=(ggplot(msubset, aes(x=Time, y=log2(value))))+
  geom_boxplot(aes(color=Treatment.x, fill=Treatment.x))+
  scale_color_manual(values=c('grey0', 'grey40'))+
  scale_fill_manual(values=c('grey0', 'grey40'))+
  ylab('log2(Fold change in gene expression)')+
  labs(color = "Treatment", fill= "Treatment")+
  theme_bw()+
  facet_grid(Line~variable, labeller = as_labeller(stage.labs))
p

pdf(file="Gene_Expression_Boxplots.pdf",height=5,width=5)
p
dev.off()




mdata_summary = ddply(mdata, .(Line,Selection_Regime), summarise, 
                      mean_Fold_Change_Atilla = mean(Fold_Change_Atilla),
                      mean_Fold_Change_PPO3 = mean(Fold_Change_PPO3),
                      mean_Fold_Change_Vkg = mean(Fold_Change_Vkg),
                      mean_Fold_Change_CG8157 = mean(Fold_Change_CG8157))

mdata_summary_2 = ddply(mdata_summary, .(Selection_Regime), summarise, 
                      mean_Fold_Change_Atilla = mean(mean_Fold_Change_Atilla),
                      mean_Fold_Change_PPO3 = mean(mean_Fold_Change_PPO3),
                      mean_Fold_Change_Vkg = mean(mean_Fold_Change_Vkg),
                      mean_Fold_Change_CG8157 = mean(mean_Fold_Change_CG8157))

#Plot fold change
p_Atilla = ggplot(mdata_summary_2, aes(x=Selection_Regime, y=log2(mean_Fold_Change_Atilla)))+
  geom_bar(stat='identity', fill='grey80')+
  geom_point(data=mdata_summary, aes(x=Selection_Regime, y=log2(mean_Fold_Change_Atilla)))
p_Atilla

p_PPO3 = ggplot(mdata_summary_2, aes(x=Selection_Regime, y=log2(mean_Fold_Change_PPO3)))+
  geom_bar(stat='identity', fill='grey80')+
  geom_point(data=mdata_summary, aes(x=Selection_Regime, y=log2(mean_Fold_Change_PPO3)))
p_PPO3

pdf(file="Atilla.pdf",height=4,width=6)
p_Atilla
dev.off()

pdf(file="PPO3.pdf",height=4,width=6)
p_PPO3
dev.off()
