library(ggplot2)
library(plyr)
library(reshape2)

#set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()
data = read.table("injection_infection_qPCR.csv", sep = ",", dec = ".", header = TRUE)
head(data)
str(data)

#Filter data
data = data[!data$Treatment=='Control',]
data = data[,c(2,3,4,16,19,22)]
head(data)

mdata = melt(data, id.vars = c('Line', 'Treatment', 'Time'))
head(mdata)

#Create data summaries
data_summary = ddply(mdata,.(Line, Treatment, Time, variable),summarise,
                     value = mean(value, na.rm = T))

data_summary

###Plot data
stage.labs <- c(
  `Female` = "Female larvae",
  `Male` = "Male larvae")

p = ggplot(data_summary, aes(x=Time, y=log2(value)))+
  geom_boxplot(aes(color=Treatment), position=position_dodge(width=1))+
  geom_point(data=mdata, aes(color=Treatment), position=position_dodge(width=1))+
#  ylim(0,100)+
  facet_grid(Line~variable)
p

Crystal_Cells = ggplot(data_summary_2, aes(x=Selection_Regime, y=mean_Crystal_Cells))+
  geom_bar(stat='identity',  width=1,  aes(fill=Selection_Regime))+
  geom_point(data=data_summary, position = position_dodge(width = .9), aes(y=mean_Crystal_Cells, color=Selection_Regime), size=2.5)+
  scale_color_manual(values = c(rep("grey0",3)))+
  scale_fill_manual(name='Selection regime', values = c('#0072b2','#d55e00'))+
  guides(color=FALSE)+
  theme_bw()+
  xlab(expression('Selection regime'))+
  ylab(expression(paste('Crystal Cells')))+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour = "grey30", face='bold',size=15),
        axis.text.y = element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"), 
        legend.position = 'none')+
  facet_grid(.~Sex, labeller = as_labeller(stage.labs))
Crystal_Cells


##Save plots
pdf(file="Crystal_Cells.pdf",height=4,width=3)
Crystal_Cells
dev.off()


###Stats
library(car)
library(lme4)

m = lmer(Crystal_Cells ~ Selection_Regime + Sex + (1|Line), data=data)
Anova(m, type='II')

