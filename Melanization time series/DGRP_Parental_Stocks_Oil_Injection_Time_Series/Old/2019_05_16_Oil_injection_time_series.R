library(ggplot2)
library(plyr)
library(Hmisc)
library(stringr)
list.files()
data = read.table('2019_05_16_Oil_injection_time_series.csv', sep=',',dec='.', head=T)
head(data)
data$DGRP_Line = as.factor(data$DGRP_Line)

data_summary=ddply(data[1:6,],.(DGRP_Line,Time),summarise,
                   Proportion_Melanised = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[1],
                   Lower_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[2],
                   Upper_Confidence_Interval_Value = binconf(x=Melanised, n=Total, alpha=0.05, method='wilson')[3],
                   Total = sum(Total))
p = ggplot(data, aes(x=Time, y=Proportion_Melanised))+
  geom_point(aes(color=DGRP_Line))+
  geom_path(aes(group=DGRP_Line, color=DGRP_Line))+
  scale_y_continuous(limits = c(0, 1),expand = c(0, 0),breaks=seq(0,1,.1))+
  theme_bw()+
  xlab('Time Post Injection')+
  ylab('Proportion Melanized')+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, face='bold', color='grey30'),
        axis.title.x = element_text(colour = "grey30", face='bold'),
        axis.title.y = element_text(colour = "grey30", face='bold'),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"))
p
