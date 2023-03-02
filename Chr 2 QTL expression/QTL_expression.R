setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(LSD)
library(dplyr)
library(org.Dm.eg.db)

### plot expression of QTL genes in wasp, oil injection data
###  a modification of the original fatbody and hemocyte DE plots to highlight QTL genes

all.top.table.hemocyte.waspvcontrol <- read.csv(file="all.top.table.hemocyte.waspvcontrol.csv",row.names = 1) ### read in DE genes for hemocytes
all.top.table.fatbody.waspvcontrol <- read.csv(file="all.top.table.fatbody.waspvcontrol.csv",row.names = 1) ### read in DE genes for fat body
all.top.table.hemocyte.oilvcontrol <- read.csv(file="all.top.table.hemocyte.oilvcontrol.csv",row.names = 1) ### read in DE genes for hemocytes
all.top.table.fatbody.oilvcontrol <- read.csv(file="all.top.table.fatbody.oilvcontrol.csv",row.names = 1) ### read in DE genes for fat body
qtl_genes <- read.table(file ="qtl_genes.txt")

nrow(all.top.table.hemocyte.waspvcontrol[all.top.table.hemocyte.waspvcontrol$logFC > 2,])
nrow(all.top.table.hemocyte.waspvcontrol[all.top.table.hemocyte.waspvcontrol$logFC < -2,])

nrow(all.top.table.fatbody.waspvcontrol[all.top.table.fatbody.waspvcontrol$logFC > 2,])
nrow(all.top.table.fatbody.waspvcontrol[all.top.table.fatbody.waspvcontrol$logFC < -2,])


all.top.table.fatbody.waspvcontrol[rownames(all.top.table.fatbody.waspvcontrol) %in% qtl_genes$V1,]$logFC > 2
all.top.table.fatbody.waspvcontrol[rownames(all.top.table.fatbody.waspvcontrol) %in% qtl_genes$V1,]$logFC < -2
fatqtlgenes <- all.top.table.fatbody.waspvcontrol[rownames(all.top.table.fatbody.waspvcontrol) %in% qtl_genes$V1,]
nrow(fatqtlgenes)
mapIds(org.Dm.eg.db, as.character(rownames(fatqtlgenes[fatqtlgenes$logFC > 2,])), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
#mapIds(org.Dm.eg.db, as.character(rownames(fatqtlgenes[fatqtlgenes$logFC < -2,])), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

all.top.table.hemocyte.waspvcontrol[rownames(all.top.table.hemocyte.waspvcontrol) %in% qtl_genes$V1,]$logFC > 2
all.top.table.hemocyte.waspvcontrol[rownames(all.top.table.hemocyte.waspvcontrol) %in% qtl_genes$V1,]$logFC < -2
hemoqtlgenes <- all.top.table.hemocyte.waspvcontrol[rownames(all.top.table.hemocyte.waspvcontrol) %in% qtl_genes$V1,]
nrow(hemoqtlgenes)
mapIds(org.Dm.eg.db, as.character(rownames(hemoqtlgenes[hemoqtlgenes$logFC > 2,])), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
mapIds(org.Dm.eg.db, as.character(rownames(hemoqtlgenes[hemoqtlgenes$logFC < -2,])), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

### in this plot, genes within the chromosome 2 QTL are highlighted against all other DE genes in fat body and hemocytes
cairo_pdf("qtl_genes.pdf", height = 5, width = 10,family = "Arial")

par(mfrow=c(1,2))

### this part basically highlights log2FC between wasp homogenATE or oil against control on separate axes for hemocytes. The line indicating wasp and injury effect is indicated.
heatscatter(all.top.table.hemocyte.waspvcontrol$logFC,all.top.table.hemocyte.oilvcontrol$logFC,ylim=c(-6,6),xlim=c(-6,6),ylab =expression('Oil v. Control log'[2]*'(FC)'),xlab =expression('Wasp homogenate v. Control log'[2]*'(FC)'), main = "Haemocytes",nrcol=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(0,1,lty=2)
text(4.5,-0.6,paste("Wasp",intToUtf8(8593)),cex=0.9)
text(4,4.7,paste("Injury",intToUtf8(8593)),cex=0.9,srt = 45)
text(-4.5,0.6,paste("Wasp",intToUtf8(8595)),cex=0.9)
text(-4,-4.7,paste("Injury",intToUtf8(8595)),cex=0.9,srt = 45)
text(0,7.2,"Haemocyte",cex=1.2)

### here highlight the QTL genes
all.top.table.hemocyte.waspvcontrol_qtl <- all.top.table.hemocyte.waspvcontrol
all.top.table.hemocyte.waspvcontrol_qtl$gene <- rownames(all.top.table.hemocyte.waspvcontrol_qtl)
all.top.table.hemocyte.oilvcontrol_qtl <- all.top.table.hemocyte.oilvcontrol
all.top.table.hemocyte.oilvcontrol_qtl$gene <- rownames(all.top.table.hemocyte.oilvcontrol_qtl)
all.top.table.hemocyte.waspvcontrol_qtl <- left_join(all.top.table.hemocyte.waspvcontrol_qtl,qtl_genes,by= c("gene"="V1"))
all.top.table.hemocyte.oilvcontrol_qtl <- left_join(all.top.table.hemocyte.oilvcontrol_qtl,qtl_genes,by= c("gene"="V1"))
siggenesall <- c(rownames(all.top.table.hemocyte.waspvcontrol[all.top.table.hemocyte.waspvcontrol$adj.P.Val < 0.05,]),rownames(all.top.table.hemocyte.oilvcontrol[all.top.table.hemocyte.oilvcontrol$adj.P.Val < 0.05,]))
siggenesall <- siggenesall[!duplicated(siggenesall)]
points(x = all.top.table.hemocyte.waspvcontrol[rownames(all.top.table.hemocyte.waspvcontrol) %in% siggenesall,]$logFC,
       y = all.top.table.hemocyte.oilvcontrol[rownames(all.top.table.hemocyte.oilvcontrol) %in% siggenesall,]$logFC,
       pch = 19, col = "red", cex = 0.6)
points(x = all.top.table.hemocyte.waspvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2),]$logFC,
       y = all.top.table.hemocyte.oilvcontrol_qtl[!is.na(all.top.table.hemocyte.oilvcontrol_qtl$V2),]$logFC,
       pch = 19, col = "black", cex = 1)
text(x = all.top.table.hemocyte.waspvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC > 2,]$logFC,
     y = 0.7+all.top.table.hemocyte.oilvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC > 2,]$logFC,
     all.top.table.hemocyte.waspvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC > 2,]$V2, cex =1, font = 3)
text(x = -1.5+all.top.table.hemocyte.waspvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC < (-2) & !all.top.table.hemocyte.waspvcontrol_qtl$gene %in% "FBgn0031558",]$logFC,
     y = all.top.table.hemocyte.oilvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC < (-2) & !all.top.table.hemocyte.waspvcontrol_qtl$gene %in% "FBgn0031558",]$logFC,
     all.top.table.hemocyte.waspvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC < (-2) & !all.top.table.hemocyte.waspvcontrol_qtl$gene %in% "FBgn0031558",]$V2, cex = 1, font = 3)
text(x = -0.7+all.top.table.hemocyte.waspvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC < (-2) & all.top.table.hemocyte.waspvcontrol_qtl$gene %in% "FBgn0031558",]$logFC,
     y = -0.7+all.top.table.hemocyte.oilvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC < (-2) & all.top.table.hemocyte.waspvcontrol_qtl$gene %in% "FBgn0031558",]$logFC,
     all.top.table.hemocyte.waspvcontrol_qtl[!is.na(all.top.table.hemocyte.waspvcontrol_qtl$V2) & all.top.table.hemocyte.waspvcontrol_qtl$logFC < (-2) & all.top.table.hemocyte.waspvcontrol_qtl$gene %in% "FBgn0031558",]$V2, cex = 1, font = 3)

### this part basically highlights log2FC between wasp + oil or oil against control on separate axes for fat body. The line indicating wasp and injury effect is indicated.
heatscatter(all.top.table.fatbody.waspvcontrol$logFC,all.top.table.fatbody.oilvcontrol$logFC,ylim=c(-6,6),xlim=c(-6,6),ylab =expression('Oil v. Control log'[2]*'(FC)'),xlab =expression('Wasp homogenate v. Control log'[2]*'(FC)'), main = "Fat body",nrcol=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(0,1,lty=2)
text(4.5,-0.6,paste("Wasp",intToUtf8(8593)),cex=0.9)
text(4,4.7,paste("Injury",intToUtf8(8593)),cex=0.9,srt = 45)
text(-4.5,0.6,paste("Wasp",intToUtf8(8595)),cex=0.9)
text(-4,-4.7,paste("Injury",intToUtf8(8595)),cex=0.9,srt = 45)
text(0,7.2,"Fat body",cex=1.1)

### here highlight the QTL genes
all.top.table.fatbody.waspvcontrol_qtl <- all.top.table.fatbody.waspvcontrol
all.top.table.fatbody.waspvcontrol_qtl$gene <- rownames(all.top.table.fatbody.waspvcontrol_qtl)
all.top.table.fatbody.oilvcontrol_qtl <- all.top.table.fatbody.oilvcontrol
all.top.table.fatbody.oilvcontrol_qtl$gene <- rownames(all.top.table.fatbody.oilvcontrol_qtl)
all.top.table.fatbody.waspvcontrol_qtl <- left_join(all.top.table.fatbody.waspvcontrol_qtl,qtl_genes,by= c("gene"="V1"))
all.top.table.fatbody.oilvcontrol_qtl <- left_join(all.top.table.fatbody.oilvcontrol_qtl,qtl_genes,by= c("gene"="V1"))
points(x = all.top.table.fatbody.waspvcontrol[rownames(all.top.table.fatbody.waspvcontrol) %in% rownames(all.top.table.fatbody.waspvcontrol[all.top.table.fatbody.waspvcontrol$adj.P.Val < 0.05,]),]$logFC,
       y = all.top.table.fatbody.oilvcontrol[rownames(all.top.table.fatbody.oilvcontrol) %in% rownames(all.top.table.fatbody.waspvcontrol[all.top.table.fatbody.waspvcontrol$adj.P.Val < 0.05,]),]$logFC,
       pch = 19, col = "red", cex = 0.6)
points(x = all.top.table.fatbody.waspvcontrol_qtl[!is.na(all.top.table.fatbody.waspvcontrol_qtl$V2),]$logFC,
       y = all.top.table.fatbody.oilvcontrol_qtl[!is.na(all.top.table.fatbody.oilvcontrol_qtl$V2),]$logFC,
       pch = 19, col = "black", cex = 1)
text(x = all.top.table.fatbody.waspvcontrol_qtl[!is.na(all.top.table.fatbody.waspvcontrol_qtl$V2) & all.top.table.fatbody.waspvcontrol_qtl$logFC > 2,]$logFC,
     y = 0.7+all.top.table.fatbody.oilvcontrol_qtl[!is.na(all.top.table.fatbody.waspvcontrol_qtl$V2) & all.top.table.fatbody.waspvcontrol_qtl$logFC > 2,]$logFC,
     all.top.table.fatbody.waspvcontrol_qtl[!is.na(all.top.table.fatbody.waspvcontrol_qtl$V2) & all.top.table.fatbody.waspvcontrol_qtl$logFC > 2,]$V2, cex = 1, font = 3)

dev.off()

