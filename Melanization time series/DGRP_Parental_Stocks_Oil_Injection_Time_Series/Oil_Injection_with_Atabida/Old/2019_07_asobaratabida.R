library(ggplot2)
library(plyr)
library(Hmisc)
library(stringr)
list.files()
data = read.table('2019_07_asobaratabida`.csv', sep=',',dec='.', head=T)
head(data)
list.files()
data$DGRP_Line = as.factor(data$DGRP_Line)
str(data)
data$Time = as.factor(data$Time)
data_summary = ddply(data,.(DGRP_Line, Time),summarise,
                       +                      Melanised = sum(Melanised),
                       +                      Non.melanised = sum(Non.melanised),
                       +                      Total = sum(Total))
data_summary_2=ddply(data_summary,.(DGRP_Line, Time),summarise,
                       +                      Proportion_Melanised = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[1],
                       +                      Lower_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[2],
                       +                      Upper_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[3],
                       +                      Total = sum(Total))
p = ggplot(data_summary_2, aes(x=DGRP_Line, y=Proportion_Melanised))+
    geom_bar(stat='identity')+
    #  scale_fill_manual(values=c('grey60','grey20'))+
    geom_errorbar(aes(ymin=Lower_Confidence_Interval_Value, ymax=Upper_Confidence_Interval_Value), width=.4)+ 
  facet_wrap(~Time)
   geom_text(aes(y=Upper_Confidence_Interval_Value+0.05, label=Total), size=3, angle=90)+
    ylab(expression('Proportion of meanised oil droplets'))+
     xlab(NULL)+
     theme_bw()+
   scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0))+
  theme(panel.grid = element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 0, color='grey30', size=8),
                    axis.text.y = element_text(size=15),
                    axis.title.y = element_text(colour = "grey30", face='bold', size=20),
                   panel.border = element_blank(),
                  axis.line = element_line(colour = "grey30"))
p
