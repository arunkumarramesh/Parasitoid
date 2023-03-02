## Author: Olivia Zhou
# This is code for analyzing the lectin-24A mutant resistance assay 
# and plotting the figure 

library(multcompView)
library(multcomp)
library(ggplot2)
library(binom)
library(emmeans)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list.files()

#import data as encapsulation per total larvae per genotype
data = read.table("lectin mutant resistance assay.csv", 
                  sep = ",", dec = ".", header = TRUE)

data
data$genotype = factor(data$genotype, 
                    levels=c('C4/C4', 'C3/C4', 'C3/C3', 
                             'C3/L69', 'L69/L69'))
#fit binomial glm
model = glm(cbind(data$resistant, data$susceptible)  ~ genotype, data = data, family = "binomial")
summary(model)

#TUKEY post hoc test on binomial glm
TUKEY.glm <- glht(model, mcp(genotype = "Tukey"))
summary(TUKEY.glm)

#generate labels
LABELS <- cld(TUKEY.glm, decreasing = TRUE)
LABELS
LABELS.df <- data.frame(LABELS$mcletters$Letters)
LABELS.df$genotype <- rownames(LABELS.df)


#generate binomial proportion confidence intervals
binomial_CIs <- binom.confint(x = data$resistant, n = data$n, methods = "prop.test")
data$lower <- binomial_CIs$lower
data$upper <- binomial_CIs$upper
data
data$genotype = factor(data$genotype, 
                       levels=c('C4/C4', 'C3/C4', 'C3/C3', 
                                'C3/L69', 'L69/L69'))

data$n <- paste("n=", data$n, sep = "")
data

#add labels for genotypes
genotype_labels <- c("resistant / resistant", 
                     "susceptible / resistant", 
                     "susceptible / susceptible",
                     "susceptible / mutant",
                     "mutant / mutant")

data = merge(data, LABELS.df)
data

#x-axis title for figure
xtitle = expression(paste(italic("Lectin-24A"), " allele"))

#bar plot 
p <- ggplot(data = data, aes(x = genotype, y = resistance, label = n)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), axis.line = element_line(),panel.border = element_blank(),panel.grid.minor = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 25)), 
        axis.text.x=element_blank()) +
  geom_bar(stat = "identity", fill = "gray") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = .2) +
  geom_text(aes(y = upper*100), nudge_y = 5, size = 3) +
  scale_x_discrete(labels = genotype_labels, expand = c(0,1)) + 
  scale_y_continuous(expand = c(0, 1)) +
  ylab("Parasitoid Encapsulation Rate (%)") + 
  xlab(xtitle) + coord_cartesian(ylim=c(0,110), clip="off") +
  labs(fill = "Tukey's test \n p < 0.05") + 
  geom_line(data = tibble(x=c(1, 2), y=c(100, 100)), 
            aes(x=x, y=y), 
            inherit.aes=FALSE) + 
  geom_text(data = tibble(x=1.5, y = 103), 
            aes(x = x, y = y, label = "a"), size = 3, 
            inherit.aes = FALSE) + 
  geom_text(data = tibble(x=3, y = 103), 
            aes(x = x, y = y, label = "bc"), size = 3, 
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=4, y = 103), 
            aes(x = x, y = y, label = "b"), size = 3, 
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=5, y = 103), 
            aes(x = x, y = y, label = "c"), size = 3, 
            inherit.aes = FALSE)


re_positions <- c(0.9, 1.1, 2.1, 4.1, 4.9, 5.1)
for(i in re_positions) {
  p <- p + geom_segment(x = i, y = -5, xend = i, yend = -10, 
                        lineend = "round", size = 1.25, col = "red")
}

sus_positions <- c(1.9, 2.9, 3.1, 3.9)
for(i in sus_positions) {
  p <- p + geom_segment(x = i, y = -5, xend = i, yend = -10, 
                        lineend = "round", size = 1.25, col = "blue")
}

mut_positions <- c(4.1, 4.9, 5.1)
for(i in mut_positions) {
  p <- p + geom_point(x = i, y = -7, shape = 25, 
                      fill = "white", color = "black", size = 1.5)
}

p <- p + annotate(geom = "text", x = 0, y = -8, label = "resistant", 
                  color = "red", size = 3)
p <- p + annotate(geom = "text", x = 0, y = -13, label = "susceptible", 
                  color = "blue", size = 3)
p
ggsave("lectin resistance assay binomial glm.pdf", width = 4, height = 4) 
