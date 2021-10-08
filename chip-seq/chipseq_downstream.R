d9<-readPeakFile('/Users/zhulab/Desktop/download/rxx_bed/bed.txt')
d9ifa<-readPeakFile('/Users/zhulab/Desktop/download/rxx_bed/4D9IFA_BKDL210030914-1a_summits.bed')
stat1<-readPeakFile('/Users/zhulab/Desktop/download/rxx_bed/stat1.txt')
covplot(d9)
covplot(d9ifa)
covplot(stat1)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene

ann_d9ifa = annotatePeak(d9ifa, annoDb="org.Mm.eg.db",TxDb=txdb)
plotAnnoPie(ann_d9ifa)

ann_d9 = annotatePeak(d9, annoDb="org.Mm.eg.db",TxDb=txdb)
plotAnnoPie(ann_d9)



library(Vennerable)
ann_stat1<-annotatePeak(stat1, annoDb="org.Mm.eg.db",TxDb=txdb)
plotAnnoPie(ann_stat1, pie3D = F,col=NA)

gene_d9<-ann_d9@anno$SYMBOL
gene_stat1<-ann_stat1@anno$SYMBOL
gene_d9ifa<-ann_d9ifa@anno$SYMBOL

file<-list(gene_d9,gene_stat1,gene_d9ifa)

library(Vennerable)
vennplot(file, by='Vennerable')






peakAnnoList <- lapply(c(a), annotatePeak, 
                       TxDb=txdb,tssRegion=c(-3000, 3000))
plotAnnoBar(peakAnnoList)



genes <- lapply(peakAnnoList, function(i) 
  as.data.frame(i)$geneId)
