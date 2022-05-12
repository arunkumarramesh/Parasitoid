library(qtl)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## read genotypes from cross into rqtl
mapthis <- read.cross("csv", file="qtl_chr3rest_with5cm.csv", estimate.map=FALSE)

### initally creating a genetic map using marker segreation ratios, to see how linked markers are in current cross
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>9)) ### only keep individuals typed for more than nine markers
mapthis <- est.rf(mapthis) ## estimates recombination fraction (rf)
lg <- formLinkageGroups(mapthis, max.rf=0.5, min.lod=6) ### assign markers into linkage groups
table(lg[,2])
rf <- pull.rf(mapthis) ### get rf
lod <- pull.rf(mapthis, what="lod") ### get lod scores
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score") ## plot rf v lod scores
mapthis <- calc.errorlod(mapthis, error.prob=0.01) ### calculate error probabilities
newmap <- est.map(mapthis, error.prob=0.01) ### estimate map with new error probabilities, to see how genetic map changes with greater genotyping error
plotMap(newmap,show.marker.names = TRUE)
mapthis <- replace.map(mapthis, newmap) ### replace old map with new one

#run this code to use standard drosophila genetic map
markers=colnames(read.csv("qtl_chr3rest_with5cm.csv")) ## column names contain recombination map position
markers2=markers[2:length(markers)]
markers2=gsub('Chr2L_','',markers2)
markers2=gsub('Chr2R_','',markers2)
markers2=gsub('cM','',markers2)
markers2=as.numeric(markers2)
flymap=as.list(markers2) 
flymap=newmap
flymap[[1]]<-markers2
names(flymap[[1]])=markers[2:length(markers)]
flymap
mapthis <- replace.map(mapthis, flymap) ### use recombination estimates from drosophila genetic map rather than calculating here using marker ratios. used further
##############################

plotGeno(mapthis, chr=1, ind=c(1:115))

mapthis <- calc.genoprob(mapthis, step=1, error.prob=0.05) ### calculate genotype probabilities along map

out.hk <- scanone(mapthis, method="hk") ## genome scan with single qtl model
summary(out.hk)

operm.hk <- scanone(mapthis, method="hk", n.perm=1000) ### here carrying out a permutation test

summary(operm.hk, alpha=0.05)
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE)

#get confidence intervals
cM_in_qtl=out.hk$pos[out.hk$lod>(max(out.hk$lod)-1.5)]
low=min(cM_in_qtl)
up=max(cM_in_qtl)
low
up

#composite interval mapping, see if more than one qtl exists
out.cim <- cim(mapthis,window=10)

#pdf(file="qtl.pdf",height=5,width=5)
plot(out.hk,xlab="Chromosome 2 Map Position (cM)")
plot(out.cim,add=T,col="red")
rect(low,-5,up,20, col= rgb(0,0,1.0,alpha=0.15), border = NA)
#dev.off()
quartz.save(type = "pdf", file="qtl.pdf",height=4,width=4)

svg(file="qtl.svg",height=4,width=4)
plot(out.hk,xlab="Chromosome 2 Map Position (cM)", frame.plot = FALSE)
plot(out.cim,add=T,col="red")
rect(low,-5,up,20, col= rgb(0,0,1.0,alpha=0.15), border = NA)
box(bty="l")
dev.off()

##to check of phenotype vs marker
comb1 <- read.csv(file="qtl_chr3rest_with5cm.csv")
table(comb1[c(1,2)])
