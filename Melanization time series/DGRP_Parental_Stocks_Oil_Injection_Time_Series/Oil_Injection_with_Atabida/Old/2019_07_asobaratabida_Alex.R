#Load required libraries
library(plyr)
library(ggplot2)
library(Hmisc)

#Set up wd and load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()
data = read.table('2019_07_asobaratabida.csv', sep=',',dec='.', head=T)
head(data)
data$DGRP_Line = as.factor(data$DGRP_Line)

data_summary=ddply(data,.(DGRP_Line, Time),summarise,
                     Proportion_Melanised = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[1],
                     Lower_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[2],
                     Upper_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[3],
                     Total = sum(Total))


##Plot
p = ggplot(data_summary, aes(x=DGRP_Line, y=Proportion_Melanised))+
  geom_bar(stat='identity', aes(fill=DGRP_Line))+
  scale_fill_manual(values=c('red1', 'dodgerblue3', 'red1', 'dodgerblue3'))+ guides(fill=FALSE)+
  geom_errorbar(aes(ymin=Lower_Confidence_Interval_Value, ymax=Upper_Confidence_Interval_Value), width=.4)+ 
  geom_text(aes(y=Lower_Confidence_Interval_Value+0.05, label=Total), size=3, angle=0)+
  ylab(expression('Proportion of Melanised Oil Droplets'))+
  xlab(NULL)+
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, color='grey30', size=10),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(colour = "grey30", face='bold', size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"))+
  facet_grid(.~Time)
p

pdf(file="Oil_Melanization_Assay_with_Atabida_24h_PI.pdf",height=4.5,width=4)
p
dev.off()

#Stats - to correct
m1 <- glm(data=data, cbind(data$Melanised, data$Non_Melanised)~DGRP_Line*Duration+Dissection_Date,family=binomial())
summary(m1)


