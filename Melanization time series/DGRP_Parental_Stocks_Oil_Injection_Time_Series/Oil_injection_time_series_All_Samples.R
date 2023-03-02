library(ggplot2)
library(plyr)
library(Hmisc)
library(stringr)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #This changes the path director to the folder where the file is saved
list.files()
data = read.table('Oil_injection_time_series_All_Samples.csv', sep=',',dec='.', head=T)
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
  xlab(expression('Time Post Injection'))+
  ylab(expression('Proportion of Melanised Oil Droplets'))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, face='bold', color='grey30'),
        axis.title.x = element_text(colour = "grey30", face='bold'),
        axis.title.y = element_text(colour = "grey30", face='bold'),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"),
        legend.position = c(0.15, 0.85))

p1

## Save plot
pdf(file="Oil_Melanization_Rec_Lines.pdf",height=3,width=4)
p1
dev.off()
