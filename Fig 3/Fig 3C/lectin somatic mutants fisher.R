## Author: Shuyu Olivia Zhou
# This is code for analyzing the encapsulation assay data for somatic mutants of lectin-24A
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

#import data
data = read.table("Lectin somatic mutants.csv", 
                  sep = ",", dec = ".", header = TRUE)

str(data)
data

data$male = factor(data$male, 
                       levels=c("2Ar5", "gRNAs (X)", "68A4", "gRNAs (III)"))
data$chromosome = factor(data$chromosome, levels=c("III", "X"))
str(data)

#separate X and III chromosome somatic mutagenesis data 
X_data <- subset(data, chromosome == "X")
X_data
III_data <- subset(data, chromosome == "III")
III_data

####### X chromosome guides #######
#Fishers exact test
X_data_fisher <- subset(X_data, select = c(male, resistant, susceptible))
X_data_fisher
X_fisher <- aggregate(.~ male, X_data_fisher, sum)
rownames(X_fisher) <- X_fisher[,1]
X_fisher[,1] <- NULL
X_fisher
X_fisher_test <- fisher.test(X_fisher)
X_fisher_test$p.value


########### III chromosome guides ############
#Fishers exact test
III_data_fisher <- subset(III_data, select = c(male, resistant, susceptible))
III_data_fisher
III_fisher <- aggregate(.~ male, III_data_fisher, sum)
rownames(III_fisher) <- III_fisher[,1]
III_fisher[,1] <- NULL
III_fisher
III_fisher_test <- fisher.test(III_fisher)
III_fisher_test$p.value



######## ALL ########
#binomial confidence intervals
data_subset <- subset(data, select = c(male, resistant, susceptible, n))
data_sum <- aggregate(. ~ male, data_subset, sum)
binomial_CIs <- binom.confint(x = data_sum$resistant, n = data_sum$n, methods = "prop.test")
data_sum$lower <- binomial_CIs$lower
data_sum$upper <- binomial_CIs$upper
data_sum
data_sum$resistance = data_sum$resistant / data_sum$n * 100
data_sum$n <- paste("n=", data_sum$n, sep = "")
data_sum
data_sum$chromosome = NA
data_sum[data_sum$male == "68A4" | data_sum$male == "gRNAs (III)", ]$chromosome <- "III"
data_sum[data_sum$male == "2Ar5" | data_sum$male == "gRNAs (X)", ]$chromosome <- "X"
data_sum
# remove n from labels, explain this in legend
data_sum$n <- gsub("n=","",data_sum$n)
data_sum$male <- gsub("gRNAs.*","gRNA",data_sum$male)

p <- ggplot(data = data_sum, aes(x = male, y = resistance, label = n)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", fill = "gray") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = .2) +
  geom_text(aes(y = upper*100), nudge_y = 5, size = 4) +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  xlab("Genotype") +
  ylab("Parasitoid Melanization Rate (%)") + 
  facet_wrap( ~ chromosome, scales = "free_x") 
p

ggsave("Lectin somatic mutants encapsulation rates.pdf", width = 2.3, height = 2.8)
