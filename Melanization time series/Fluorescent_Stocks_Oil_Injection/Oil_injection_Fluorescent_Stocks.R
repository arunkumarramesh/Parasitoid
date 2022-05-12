#load required packages
library(ggplot2)
library(plyr)
library(Hmisc)
library(tidyverse)

#Set up 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #This changes the path director to the folder where the file is saved
list.files()
data = as.data.frame(read_csv('Fluorescent_Stocks_Oil_Injections.csv'))
head(data)

data$Cross = paste(data$Female, data$Male, sep='x')

#Create data summaries
data_summary = ddply(data,.(Cross),summarise,
                     Melanised = sum(Melanised),
                     Non_Melanised = sum(Non_Melanised),
                     Total = sum(Total))

data_summary_2=ddply(data_summary,.(Cross),summarise,
                     Proportion_Melanised = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[1],
                     Lower_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[2],
                     Upper_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[3],
                     Total = sum(Total))



##Plot
p = ggplot(data_summary_2, aes(x=Cross, y=Proportion_Melanised))+
  geom_bar(stat='identity', aes(fill=Cross))+ guides(fill="none")+
  geom_errorbar(aes(ymin=Lower_Confidence_Interval_Value, ymax=Upper_Confidence_Interval_Value), width=0.3)+
  geom_text(aes(y=Upper_Confidence_Interval_Value+0.05, label=Total), size=3, angle=0)+
  scale_fill_manual(values=c('red1', 'dodgerblue3'))+
  scale_y_continuous(limits = c(0, 0.9),expand = c(0, 0),breaks=seq(0,1,.1))+
  scale_x_discrete(labels=c("Hml>GFP\n\ Msn-mCherry\n\ DGRP-437",
                            "Hml>GFP\n\ Msn-mCherry\n\ DGRP-892"))+

  xlab(NULL)+
  ylab('Proportion Oil Droplets Melanized')+
  theme_classic2()
p



## Save plot
pdf(file="Oil_Melanization_Fluorescent_Stocks.pdf",height=3.5,width=3)
p
dev.off()

## fishers test to compare melanization in the two lines
fisher.test(cbind(data[5:6]))
