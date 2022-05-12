library(ggplot2)
library(plyr)
library(Hmisc)
library(stringr)
library(car)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #This changes the path director to the folder where the file is saved
list.files()
data = read.table('2019_07_asobaratabida.csv', sep=',',dec='.', head=T)
head(data)

data$DGRP_Line = as.factor(data$DGRP_Line)
str(data)

#Create data summaries
data_summary = ddply(data,.(DGRP_Line, Time),summarise,
                     Melanised = sum(Melanised),
                     Non_Melanised = sum(Non_Melanised),
                     Total = sum(Total))

data_summary_2=ddply(data_summary,.(DGRP_Line, Time),summarise,
                     Proportion_Melanised = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[1],
                     Lower_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[2],
                     Upper_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[3],
                     Total = sum(Total))


##Plot
p1 = ggplot(data_summary_2, aes(x=Time, y=Proportion_Melanised))+
  geom_point(aes(color=DGRP_Line))+
  geom_errorbar(aes(ymin=Lower_Confidence_Interval_Value, ymax=Upper_Confidence_Interval_Value, color=DGRP_Line), width=0)+ 
  geom_path(aes(group=DGRP_Line, color=DGRP_Line))+
  scale_fill_manual(values=c('red1', 'dodgerblue3'))+
  scale_color_manual(values=c('red1', 'dodgerblue3'))+
  scale_y_continuous(limits = c(0, 1),expand = c(0, 0),breaks=seq(0,1,.1))+
  scale_x_continuous(breaks=unique(data_summary_2$Time))+
  theme_bw()+
  xlab(expression('Hours post injection'))+
  ylab(expression('Proportion of Melanised Oil Droplets'))+
  theme_classic2()+
  theme(legend.position = c(0.85, 0.20))+ 
  guides(color=guide_legend("DGRP line"))

p1

## Save plot
pdf(file="Oil_Melanization_Atabida.pdf",height=3.4,width=3)
p1
dev.off()

#Stats - to correct
m1 <- glm(data=data, cbind(data$Melanised, data$Non_Melanised)~DGRP_Line + Duration,family=binomial())
summary(m1)
Anova(m1)

