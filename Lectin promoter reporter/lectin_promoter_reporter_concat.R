## Author: Olivia Zhou
# This is code for visualizing the lectin-24A promoter reporter assay data

library(ggplot2)
library(plyr)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()

#import data
data = read.table("Promoter_reporter_concat.csv", 
                  sep = ",", dec = ".", header = TRUE)
head(data)
str(data)

data$Treatment = factor(data$Treatment, levels = c("no_infection", "infection"))

data$Label = factor(data$Label, 
                       levels = c('LP437_6xSNP_8A', 'LP437_7A', 'LP437_8A', 
                                  "LP437_2xSNP","LP437_4xSNP",
                                  'LP437_21A', 'LP437_7A_8A',   
                                  "DGRP437",  "DGRP892"))

str(data)

###Stats
library(lme4)
library(car)
library(multcomp)

LP437_noninf_all = data$Value[data$Line == "LP437" & data$Treatment == "no_infection"]/
  data$Protein_Quantity_ug[data$Line == "LP437" & data$Treatment == "no_infection"]

LP437_inf_all = data$Value[data$Line == "LP437" & data$Treatment == "infection"]/
  data$Protein_Quantity_ug[data$Line == "LP437" & data$Treatment == "infection"]


m1= lmer(data=data, log10(Value/Protein_Quantity_ug)~Label + (1|Extr.day), subset= Treatment=='no_infection')
summary(m1)
Anova(m1)
m1.glm <- summary(glht(m1, mcp(Label="Tukey")))
labels_noinf <- cld(m1.glm)$mcletters$Letters
labels_noinf

m2= lmer(data=data, Value/Protein_Quantity_ug~Label + (1|Extr.day), subset= Treatment=='infection')
summary(m2)
Anova(m2)
m2.glm <- summary(glht(m2, mcp(Label="Tukey")))
labels_inf <- cld(m2.glm)$mcletters$Letters
labels_inf


#plot data on log scale
p = ggplot(data, aes(x=Label, y=log10(Value/Protein_Quantity_ug)))+
  theme_bw() +
  theme(legend.position = c(0.15, 0.1),
        legend.title = element_blank(),
        axis.text = element_text(size = 12.5),
        legend.text = element_text(size = 12.5)) +
  coord_cartesian(ylim = c(2.5, 5.5)) +
  geom_boxplot(aes(color = Treatment), position = "identity", outlier.shape = NA)+
  geom_point(aes(color = Treatment), position = position_jitter(width = 0.3, height = 0), 
             size = 1.0)+
  ylab(expression(log[10]*"(Relative Fluorescence Intensity / Total Protein (ug))")) + 
  theme(text = element_text(size = 12.5)) + 
  scale_color_hue(labels = c("no infection", "infection")) +
  xlab("Genotype")
p

#add Tukey HSD labels

for(i in 1:9) {
  p <- p + annotate("text", label = labels_noinf[[i]],
                    x = i, y = 5.25, size = 5, color = "#F8766D", fontface = "bold")
}

for(i in 1:9) {
  p <- p + annotate("text", label = labels_inf[[i]],
                    x = i, y = 5.5, size = 5, color = "#00BFC4", fontface = "bold")
}

p

#flip x and y axes
p + coord_flip()


ggsave("lectin-24A promoter reporter assay Tukey.pdf", width = 7, height = 6)
