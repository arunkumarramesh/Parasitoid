#Load required libraries
library(plyr)
library(ggplot2)
library(Hmisc)

#Set up wd and load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()
data = read.table('2019_03_21_Fluorescent_stock_cross.csv', sep=',',dec='.', head=T)
head(data)

data_summary = ddply(data,.(Cross),summarise,
                     Resistant = sum(Phenotype=='R'),
                     Susceptible = sum(Phenotype=='S'))

data_summary$Total=data_summary$Resistant+data_summary$Susceptible

data_summary_2=ddply(data_summary,.(Cross),summarise,
                     Proportion_Resistant = binconf(x=Resistant, n=Total, alpha=0.05, method='wilson')[1],
                     Lower_Confidence_Interval_Value = binconf(x=Resistant, n=Total, alpha=0.05, method='wilson')[2],
                     Upper_Confidence_Interval_Value = binconf(x=Resistant, n=Total, alpha=0.05, method='wilson')[3],
                     Total = sum(Total))


##Plot
p = ggplot(data_summary_2, aes(x=Cross, y=Proportion_Resistant))+
  geom_bar(stat='identity', aes(fill=Cross))+
  scale_fill_manual(values=c('red1', 'dodgerblue3', 'red1', 'dodgerblue3'))+ guides(fill="none")+
  geom_errorbar(aes(ymin=Lower_Confidence_Interval_Value, ymax=Upper_Confidence_Interval_Value), width=.4)+ 
  geom_text(aes(y=Upper_Confidence_Interval_Value+0.05, label=Total), size=3, angle=0)+
  ylab(expression('Proportion of Melanised Wasp Embryo'))+
  xlab(NULL)+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_x_discrete(labels=c("DGRP-437",
                            "DGRP-892",
                            "Hml>GFP\n\ Msn-mCherry\n\ DGRP-437",
                            "Hml>GFP\n\ Msn-mCherry\n\ DGRP-892"))+
  theme_classic2()
p

pdf(file="Resistance_Assay_Fluorescent_Stocks.pdf",height=3.5,width=5)
p
dev.off()

#Stats
data_summary_new = ddply(data,.(Cross,Date),summarise,
                     Resistant = sum(Phenotype=='R'),
                     Susceptible = sum(Phenotype=='S'))
m1 <- glm(data=data_summary_new, cbind(Resistant,Susceptible) ~ Cross+Date,family=binomial())
summary(m1)
anova(m1,test="Chisq")

data_summary_new2 = ddply(data,.(Cross),summarise,
                         Resistant = sum(Phenotype=='R'),
                         Susceptible = sum(Phenotype=='S'))
fisher.test(data_summary_new2[3:4,2:3])
