## Author: Olivia Zhou
# This is code for analyzing the encapsulation assay data for lectin-24A overexpression  
# and plotting the figure 

library(multcomp)
library(lme4)
library(car)
library(ggplot2)
library(binom)
library(emmeans)
library(tidyverse)
library(tidyr)
library(stringr)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()

#import data as encapsulation per total larvae per genotype
data = read.table("UAS-lectin round 2 all.csv", 
                  sep = ",", dec = ".", header = TRUE)

str(data)
data

data$genotype = factor(data$genotype, 
                       levels=c("2Ar5", "UAS-lectin-1", "UAS-lectin-2"))
data$date = factor(data$date, levels = c("2703", "3003", "3103", "2704", "2804"))
data

#fit binomial glm accounting for random variation between vials
model = glmer(cbind(data$resistant, data$susceptible)  ~ genotype + (1|vial), data = data, 
              family = "binomial")

summary(model)

anova(model)

#TUKEY post hoc test on binomial glm
TUKEY.glm <- glht(model, mcp(genotype = "Tukey"))
summary(TUKEY.glm)


#generate TUKEY post hoc test labels
LABELS <- cld(TUKEY.glm, decreasing = TRUE)
LABELS
LABELS.df <- data.frame(LABELS$mcletters$Letters)
LABELS.df$genotype <- rownames(LABELS.df)
LABELS.df

#aggregate data by genotype
data <- subset(data, select = -c(vial, date))
data
data_sum <- aggregate(. ~ genotype, data, sum)
data_sum

#generate binomial proportion confidence intervals
binomial_CIs <- binom.confint(x = data_sum$resistant, n = data_sum$n, methods = "prop.test")
data_sum$lower <- binomial_CIs$lower
data_sum$upper <- binomial_CIs$upper
data_sum

data_sum$resistance = data_sum$resistant / data_sum$n * 100

data_sum$n <- paste("n=", data_sum$n, sep = "")
data_sum


#add TUKEY labels
data_sum = merge(data_sum, LABELS.df, by.x = "genotype", by.y = "genotype")

data_sum

data_sum$genotype = paste(data_sum$genotype, "x C7-GAL4", sep = " ")

data_sum

# remove n from labels, explain this in legend
data_sum$n <- gsub("n=","",data_sum$n)

#bar plot 
p <- ggplot(data = data_sum, aes(x = factor(genotype, levels = c("UAS-lectin-1 x C7-GAL4",
                                            "UAS-lectin-2 x C7-GAL4",
                                            "2Ar5 x C7-GAL4")),
                                            y = resistance, label = n)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_bar(stat = "identity", fill = "gray") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = .2) +
  geom_text(aes(y = upper*100), nudge_y = 5, size = 4) +
  scale_x_discrete(labels = function(x)
    stringr::str_wrap(x, width = 15)) +
  scale_y_continuous(expand = c(0, 1)) +
  xlab("Genotype") +
  ylab("Parasitoid\nMelanization Rate (%)") + 
  coord_cartesian(ylim=c(0,55), clip="off") +
  labs(fill = "Tukey's test \n p < 0.05") +
  geom_line(data = tibble(x=c(1, 2), y=c(47, 47)), 
            aes(x=x, y=y), 
            inherit.aes=FALSE) + 
  geom_text(data = tibble(x=1.5, y = 51), 
            aes(x = x, y = y, label = "a"), size = 4, 
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=3, y = 51), 
            aes(x = x, y = y, label = "b"), size = 4, 
            inherit.aes = FALSE)
p

ggsave("UAS-lectin overexpression resistance assay binomial glm.pdf", width = 2.4, height = 2.8) 

