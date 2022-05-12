setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#this is the oberved risk ratio
#from non-recombinant flies (adults with capsules) (same  marker genotype at 3 and 27cM)
#thre were 738 flies RR genotype and 260 SS flies
observedRR=906/282


##function to simulate recombinant genotypes
genosim=function(nrecombs=298, risk=observedRR){

#number genetically resistant recombinants
nR=rbinom(1, nrecombs, risk/(risk+1))

#recomninants at 3 and 27 cM, so 240 0.1cM intervals
#1,0 recombinants
mat=matrix(1,nrow=1000,ncol=240)
breakpoint=sample(238, 1000,replace=T)+1
for (i in 1:1000){
  mat[i,breakpoint[i]:240]=0
}

#0,1 recombinants
mat2=matrix(0,nrow=1000,ncol=240)
breakpoint=sample(238, 1000,replace=T)+1
for (i in 1:1000){
  mat2[i,breakpoint[i]:240]=1
}

#shuffle the two
mat3=rbind(mat,mat2)[sample(1:2000),]

#split into resistant and susceptible at position 70 (ie 7cm from start, 10cM) 
#assign R as genotype 1
matR=mat3[mat3[,70]==1,]
matS=mat3[mat3[,70]!=1,]

#sample resistant and susceptible lines following risk ratio
genotypes=rbind(matR[1:nR,],matS[1:(nrecombs-nR),])
genotypes
}



genosim2=function(N=251, R=observedRR,ldrop=1.5,chidrop=4.5){
      
    QTL=vector()
    CIsize=vector()
    CIgene=vector()
    lodCIsize=vector()
    lodCIgene=vector()
    chiCIsize=vector()
    chiCIgene=vector()
    
    
    for(i in 1:10000){
      genos=genosim(nrecombs=N, risk=R)
      x=which(colSums(genos)==max(colSums(genos)))
      #note this is to randomly choose marker when ties
      QTL[i]=x[sample(length(x))][1]
    
       #get chi2 scores
            chi2=vector()
            for(j in 1:ncol(genos)){
              x=prop.test(sum(genos[,j]),nrow(genos),p=0.5)$statistic
              #next line changes sign depnding on whether in correct direction
              chi2[j]=ifelse((sum(genos[,j])/nrow(genos)>0.5),
                             x,-x)
            }
            
            ci_chi=which(chi2>(max(chi2)-chidrop))
            chiCIsize[i]=(max(ci_chi)-min(ci_chi))/10
            chiCIgene[i]=ifelse(min(ci_chi)<=70&max(ci_chi)>=70,"gene in CI","gene not in CI")
            
            #convert to lod 
#            lod=chi2/(2*log(10))
##            ci=which(lod>(max(lod)-ldrop))
#            lodCIsize[i]=(max(ci)-min(ci))/10
#            lodCIgene[i]=ifelse(min(ci)<=70&max(ci)>=70,"gene in CI","gene not in CI")
            
    }
    #use this for lod drop conf intervals
#    c(table(CIgene)[2]/length(CIgene),
#    mean(CIsize),
#    table(lodCIgene)[2]/length(lodCIgene),
#    mean(lodCIsize)
#    )
    #use this for chi-squared ci
    c(table(chiCIgene)[2]/length(chiCIgene),
#      mean(chiCIsize),
#      table(chiCIgene)[2]/length(chiCIgene),
      mean(chiCIsize)
    )
}
###############run simulations

#########################
##test different chi-squared drops. called lods as old script!




lods=40:60/10
effect_L=matrix(nrow=length(lods),ncol=2)



for(k in 1:length(lods)){
  effect_L[k,]=genosim2(N=251, R=observedRR,chidrop=lods[k])
  print(k)
}
effect_L2=cbind(lods,effect_L[,1:2])
colnames(effect_L2)[2:3]=c("proportion_qtl_containing_gene","mean_qtl_size")
#write.csv(effect_L2,file="QTL simulations different chi drops.csv")

effect_L2=read.csv(file="QTL simulations different chi drops.csv")

pdf(file="QTL simulations different chi drops.pdf",height=4,width=5)
par(mfrow=c(1,1))


plot(100*effect_L2[,3]~effect_L2[,2],
     ylim=c(0,9),
     xlab="chi2 drop",
     ylab="% times gene outside 95% CI interval",
     main='risk ratio=3.21,N=298')
abline(h=5,col="red",lwd=3)

dev.off()

