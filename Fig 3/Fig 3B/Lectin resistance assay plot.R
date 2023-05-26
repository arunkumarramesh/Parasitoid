## Author: Olivia Zhou
# This is code for analyzing the lectin-24A mutant resistance assay 
# and plotting the figure 

library(multcompView)
library(ggplot2)
library(binom)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#import data
data = read.table("lectin mutant resistance assay full.csv", 
                  sep = ",", dec = ".", header = TRUE)
str(data)
data$cross = factor(data$cross, 
                    levels=c('C4/C4', 'C3/C4', 'C3/C3', 
                             'C3/L69', 'L69/L69'))

#fit linear model to data
model = lm(data$resistance ~ data$cross)
#ANOVA
ANOVA = aov(model)
summary(ANOVA)

#Tukey test for each pair of treatment
TUKEY <- TukeyHSD(x = ANOVA, 'data$cross', conf.level = 0.95)
TUKEY
plot(TUKEY, las = 1)

#group the treatments that are not significantly different by Tukey
generate_label_df <- function(TUKEY, variable){
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  Tukey.labels$cross = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$cross), ]
  return(Tukey.labels)
}

#generate Tukey labels 
LABELS <- generate_label_df(TUKEY, "data$cross")
colnames(LABELS) <- c("TUKEY", "cross")

#import data as encapsulation per total larvae per genotype
data = read.table("lectin mutant resistance assay.csv", 
                  sep = ",", dec = ".", header = TRUE)

data

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

data = merge(data, LABELS, by.x = "genotype", by.y = "cross")
data

#x-axis title for figure
xtitle = expression(paste(italic("Lectin-24A"), " allele"))

# remove n from labels, explain this in legend
data$n <- gsub("n=","",data$n)

#bar plot 

p <- ggplot(data = data, aes(x = genotype, y = resistance, label = n)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 25)), 
        axis.text.x=element_blank()) +
  geom_bar(stat = "identity", fill = "gray") +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = .2) +
  geom_text(aes(y = upper*100), nudge_y = 5, size = 4) +
  scale_x_discrete(labels = genotype_labels, expand = c(0,1)) + 
  scale_y_continuous(expand = c(0, 1)) +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 13)) +
  ylab("Parasitoid\nMelanization Rate (%)") + 
  xlab(xtitle) + coord_cartesian(ylim=c(0,110), clip="off") +
  labs(fill = "Tukey's test \n p < 0.05") + 
  geom_line(data = tibble(x=c(1, 2), y=c(100, 100)), 
            aes(x=x, y=y), 
            inherit.aes=FALSE) + 
  geom_line(data = tibble(x=c(3, 5), y=c(100, 100)), 
            aes(x=x, y=y), 
            inherit.aes=FALSE) +
  geom_text(data = tibble(x=1.5, y = 103), 
            aes(x = x, y = y+2, label = "a"), size = 4.5, 
            inherit.aes = FALSE) + 
  geom_text(data = tibble(x=4, y = 103+2), 
            aes(x = x, y = y, label = "b"), size = 4.5, 
            inherit.aes = FALSE)

re_positions <- c(0.85, 1.15, 2.15, 4.15, 4.85, 5.15)
for(i in re_positions) {
  p <- p + geom_segment(x = i, y = -4, xend = i, yend = -10, 
                        lineend = "round", size = 1.25, col = "red")
}

sus_positions <- c(1.85, 2.85, 3.15, 3.85)
for(i in sus_positions) {
  p <- p + geom_segment(x = i, y = -4, xend = i, yend = -10, 
                        lineend = "round", size = 1.25, col = "blue")
}

mut_positions <- c(4.15, 4.85, 5.15)
for(i in mut_positions) {
  p <- p + geom_point(x = i, y = -7, shape = 25, 
                      fill = "white", color = "black", size = 1.5)
}

p <- p + annotate(geom = "text", x = -1, y = -7, label = "resistant", 
                  color = "red", size = 4)
p <- p + annotate(geom = "text", x = -1, y = -14, label = "susceptible", 
                  color = "blue", size = 4)
p
ggsave("lectin resistance assay.pdf", width = 2.8, height = 3.3) 
