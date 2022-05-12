setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#this is the drop in chi2 to get confidence intervals on location. Derived empirically from simulation script.
chi_drop=4.6

library(zoo)
HSSH <- read.csv("HS-SH,15-5.csv")

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

#get location CIs

locations=as.numeric(colnames(genotypes1))

outside=which(chi2<max(chi2-4.6))
peak=which(chi2==max(chi2))
low=locations[outside[max(which(outside<peak))]]
up=locations[outside[min(which(outside>peak))]]

low
up


pdf(file="qtl mapping in adults.pdf",height=3.8,width=3.8)
par(mar=c(5.1,4.3,2.1,1.2))

plot(chi2~locations,type="l",
     xlab="Chromosome 2 Map Position (cM)",
     ylab=expression(chi^2),
     las=1,bty="l" )
points(chi2~locations,type="p",pch=20)
rect(low,-15,up,40, col= rgb(0,0,1.0,alpha=0.15), border = NA)
#remove this line: use for guessing where to put markers
#abline(h=max(chi2)-chi_drop,col="red",lwd=3)
#beware this is not interval mapping so pay little attention to line between markers

dev.off()


svg(file="qtl mapping in adults.svg",height=3.8,width=3.8)
par(mar=c(5.1,4.3,2.1,1.2))

plot(chi2~locations,type="l",
     xlab="Chromosome 2 Map Position (cM)",
     ylab=expression(chi^2),
     las=1,bty="l" )
points(chi2~locations,type="p",pch=20)
rect(low,-15,up,40, col= rgb(0,0,1.0,alpha=0.15), border = NA)
#remove this line: use for guessing where to put markers
#abline(h=max(chi2)-chi_drop,col="red",lwd=3)
#beware this is not interval mapping so pay little attention to line between markers

dev.off()

