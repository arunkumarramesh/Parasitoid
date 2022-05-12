###Script to test if experimental strategy or cage had an effect on chr II risk ratio
library(dplyr)
library(ggplot2)
library(psych)
library(polycor)
library(zoo)
library(lmerTest)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##file with all meta data for experimental strategies and genotypes of adult flies with capsules
all_samples <- read.csv(file = "all_samples.csv", header = T, stringsAsFactors = F)
all_samples<- all_samples[-c(1:6)]
all_samples<- all_samples[-c(2)]
all_samples<- all_samples[-c(3)]
all_samples<- all_samples[-c(18:25)]
all_samples<- all_samples[-c(12:16)]

##create unique experimental strategy id
all_samples$ID <- paste(all_samples$Egg.lay.period,all_samples$Infection.start.date,all_samples$Infection.start.date...Egg.lay.date,sep="_")

##calculate risk ratio by experimental strategy
all_samples2 <- cbind(all_samples["ID"],all_samples["Genotype..3.and.27.cM."])
all_samples2 <- all_samples2 %>% group_by(ID) %>% count(Genotype..3.and.27.cM.)
all_samplesHH <- all_samples2[all_samples2$Genotype..3.and.27.cM. %in% "HH",]
all_samplesSS <- all_samples2[all_samples2$Genotype..3.and.27.cM. %in% "SS",]
all_samplesrr <- inner_join(all_samplesHH,all_samplesSS,by="ID")
all_samplesrr$rr <- all_samplesrr$n.x/all_samplesrr$n.y
all_samplesrr <- cbind(all_samplesrr["ID"],all_samplesrr["rr"])
all_samples <- left_join(all_samples,all_samplesrr,by="ID")

##calculate number of recombinants by experimental strategy
all_samples2_nrec <- all_samples2[all_samples2$Genotype..3.and.27.cM. %in% c("HH","SS"),]
all_samples2_nrec <- all_samples2_nrec %>% group_by(ID) %>% tally(n)
all_samples2_rec <- all_samples2[all_samples2$Genotype..3.and.27.cM. %in% c("HS","SH"),]
all_samples2_rec <- all_samples2_rec %>% group_by(ID) %>% tally(n)
all_samplesrec <- inner_join(all_samples2_nrec,all_samples2_rec,by="ID")
all_samplesrec$n <- all_samplesrec$n.x + all_samplesrec$n.y
all_samplesrec$rec <- all_samplesrec$n.y / all_samplesrec$n
all_samplesrec <- cbind(all_samplesrec["ID"],all_samplesrec["rec"])
all_samples <- left_join(all_samples,all_samplesrec,by="ID")

##calculate proportion of recombinant samples amplifying wasp DNA test by experimental strategy
all_samples_wasp <- all_samples[all_samples$Wasp.DNA.amplified %in% c("Y","N"),]
all_samples_wasp <- cbind(all_samples_wasp["ID"],all_samples_wasp["Wasp.DNA.amplified"])
all_samples_wasp <- all_samples_wasp %>% group_by(ID) %>% count(Wasp.DNA.amplified)
all_samples_wasp_y <- all_samples_wasp[all_samples_wasp$Wasp.DNA.amplified %in% "Y",]
all_samples_wasp_n <- all_samples_wasp[all_samples_wasp$Wasp.DNA.amplified %in% "N",]
all_samples_wasp <- inner_join(all_samples_wasp_y,all_samples_wasp_n,by = "ID")
all_samples_wasp$waspstatus <- all_samples_wasp$n.x/(all_samples_wasp$n.x+all_samples_wasp$n.y)
all_samples_wasp <-  cbind(all_samples_wasp["ID"],all_samples_wasp["waspstatus"])
all_samples <- left_join(all_samples,all_samples_wasp,by="ID")

##remove duplicated rows
all_samples <- all_samples[-c(9:11)]
all_samples <- all_samples[-c(9)]
all_samples <- all_samples[!duplicated(all_samples),]
all_samples$Scaled.total.collected <- all_samples$Total.flies.collected.on.that.date/all_samples$No..vials.made.on.that.date

##correlation matrix for factors
cor.plot(hetcor(all_samples[complete.cases(all_samples),]))

##test effect of each variable in experimental strategy on risk ratio
anova(lmerTest::lmer(data=all_samples,rr~ Infection.period*Infection.start.date...Egg.lay.date*Egg.lay.period+(1|Cage)))

ggplot(data=all_samples,aes(x=as.factor(Infection.period),y=waspstatus))+
  geom_boxplot()+
  geom_point()

##plot qtl map by cage

par(mfrow=c(1,3))
HSSH <- read.csv("HS-SH,15-5.csv")
HSSH <- HSSH[HSSH$Cage %in% c(1),]
#create matrix of marker genotypes
genotypes1=HSSH[,grep("cM",names(HSSH))]
names(genotypes1)=gsub(".cM", "", names(genotypes1))
names(genotypes1)=gsub("X", "", names(genotypes1))
ind=order(as.numeric(names(genotypes1)))
genotypes1=data.matrix(genotypes1[,ind])
##only samples with wasp 1 amplification
genotypes1 <- genotypes1[HSSH$Wasp.DNA.amplified=="Y",]
#impute missing genotypes when flanking values same
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
genotypes2=t(apply(genotypes1,1,na.approx))
ind=rowSums(!t(apply(genotypes2,1,is.wholenumber)))==0
genotypes=genotypes2[ind,]
#test if genotypes differ from 50:50
genotypes[genotypes==2]=0
chi2=vector()
for(i in 1:ncol(genotypes)){
  x=prop.test(sum(genotypes[,i]),nrow(genotypes),p=0.5)$statistic
  #next line changes sign depnding on whether in correct direction
  chi2[i]=ifelse((sum(genotypes[,i])/nrow(genotypes)>0.5),
                 x,-x)
}
locations=as.numeric(colnames(genotypes1))
#pdf(file="qtl mapping in adults_removegreateronedayofinfection.pdf",height=3.8,width=3.8)
par(mar=c(5.1,4.3,2.1,1.2))
plot(chi2~locations,type="l",
     xlab="",
     ylab=expression(chi^2),
     las=1,bty="l",
     main=paste("Cage 1, n=",nrow(genotypes1),sep=""))
points(chi2~locations,type="p",pch=20)
#dev.off()


HSSH <- read.csv("HS-SH,15-5.csv")
HSSH <- HSSH[HSSH$Cage %in% c(2),]
#create matrix of marker genotypes
genotypes1=HSSH[,grep("cM",names(HSSH))]
names(genotypes1)=gsub(".cM", "", names(genotypes1))
names(genotypes1)=gsub("X", "", names(genotypes1))
ind=order(as.numeric(names(genotypes1)))
genotypes1=data.matrix(genotypes1[,ind])
##only samples with wasp 1 amplification
genotypes1 <- genotypes1[HSSH$Wasp.DNA.amplified=="Y",]
#impute missing genotypes when flanking values same
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
genotypes2=t(apply(genotypes1,1,na.approx))
ind=rowSums(!t(apply(genotypes2,1,is.wholenumber)))==0
genotypes=genotypes2[ind,]
#test if genotypes differ from 50:50
genotypes[genotypes==2]=0
chi2=vector()
for(i in 1:ncol(genotypes)){
  x=prop.test(sum(genotypes[,i]),nrow(genotypes),p=0.5)$statistic
  #next line changes sign depnding on whether in correct direction
  chi2[i]=ifelse((sum(genotypes[,i])/nrow(genotypes)>0.5),
                 x,-x)
}
locations=as.numeric(colnames(genotypes1))
#pdf(file="qtl mapping in adults_removegreateronedayofinfection.pdf",height=3.8,width=3.8)
par(mar=c(5.1,4.3,2.1,1.2))
plot(chi2~locations,type="l",
     xlab="",
     ylab=expression(chi^2),
     las=1,bty="l" ,
     main=paste("Cage 2, n=",nrow(genotypes1),sep=""))
points(chi2~locations,type="p",pch=20)
#dev.off()


HSSH <- read.csv("HS-SH,15-5.csv")
HSSH <- HSSH[HSSH$Cage %in% c(3),]
#create matrix of marker genotypes
genotypes1=HSSH[,grep("cM",names(HSSH))]
names(genotypes1)=gsub(".cM", "", names(genotypes1))
names(genotypes1)=gsub("X", "", names(genotypes1))
ind=order(as.numeric(names(genotypes1)))
genotypes1=data.matrix(genotypes1[,ind])
##only samples with wasp 1 amplification
genotypes1 <- genotypes1[HSSH$Wasp.DNA.amplified=="Y",]
#impute missing genotypes when flanking values same
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
genotypes2=t(apply(genotypes1,1,na.approx))
ind=rowSums(!t(apply(genotypes2,1,is.wholenumber)))==0
genotypes=genotypes2[ind,]
#test if genotypes differ from 50:50
genotypes[genotypes==2]=0
chi2=vector()
for(i in 1:ncol(genotypes)){
  x=prop.test(sum(genotypes[,i]),nrow(genotypes),p=0.5)$statistic
  #next line changes sign depnding on whether in correct direction
  chi2[i]=ifelse((sum(genotypes[,i])/nrow(genotypes)>0.5),
                 x,-x)
}
locations=as.numeric(colnames(genotypes1))
#pdf(file="qtl mapping in adults_removegreateronedayofinfection.pdf",height=3.8,width=3.8)
par(mar=c(5.1,4.3,2.1,1.2))
plot(chi2~locations,type="l",
     xlab="",
     ylab=expression(chi^2),
     las=1,bty="l" ,
     main=paste("Cage 3, n=",nrow(genotypes1),sep=""))
points(chi2~locations,type="p",pch=20)

mtext("Chromosome 2 Map Position (cM)",
      side=1,outer=T,line=-2.2,cex=1.1, at = 0.5)
#dev.off()