library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
pheno_data<-read.table('~/test_r.txt',header = T)
AS<-ballgown(dataDir = "/home/zhushu/zhouty19990625/download/rxx_data/ballgown",samplePattern = "S",pData = pheno_data)
save.image('/home/zhushu/zhouty19990625/download/rxx_data/output_test/ballgown.Rdata')

# 这里滤掉了样本间差异少于一个转录本的数据
AS_filt = subset(AS,"rowVars(texpr(AS)) >1",genomesubset=TRUE)
# 比较男和女的基因表达差异
results_transcripts = stattest(AS_filt, feature="transcript",covariate="feature",getFC=TRUE, meas="FPKM")
# 确认组间有差异的基因
results_genes = stattest(AS_filt, feature="gene", covariate="feature", getFC=TRUE, meas="FPKM")
# 对结果增加基因名和基因ID
results_transcripts = data.frame(geneNames=ballgown::geneNames(AS_filt), geneIDs=ballgown::geneIDs(AS_filt), results_transcripts)
results_genes = data.frame(geneNames=ballgown::geneNames(AS_filt), geneIDs=ballgown::geneIDs(AS_filt), results_genes)
# 按照P值排序（从小到大）
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

plotTranscripts(ballgown::geneIDs(AS)[115993],AS, main=c('Gene Nlrp9b in sample S215WT_L2_337X37'), sample=c('S215WT_L2_337X37'))
plotTranscripts(ballgown::geneIDs(AS)[115993],AS, main=c('Gene Nlrp9b in sample S316WT_L2_338X38'), sample=c('S316WT_L2_338X38'))
plotTranscripts(ballgown::geneIDs(AS)[115993],AS, main=c('Gene Nlrp9b in sample S418KO_L2_339X39'), sample=c('S418KO_L2_339X39'))
plotTranscripts(ballgown::geneIDs(AS)[115993],AS, main=c('Gene Nlrp9b in sample S519KO_L2_340X40'), sample=c('S519KO_L2_340X40'))
plotTranscripts(ballgown::geneIDs(AS)[115993],AS, main=c('Gene Nlrp9b in sample S620KO_L2_341X41'), sample=c('S620KO_L2_341X41'))
plotTranscripts(ballgown::geneIDs(AS)[115993],AS, main=c('Gene Nlrp9b in sample S723WT_L2_342X42'), sample=c('S723WT_L2_342X42'))

plotMeans('stringtie_merged.16167', AS_filt,groupvar="feature",legend=T)
test_AS_cut<-results_transcripts[which(results_transcripts$geneIDs%in%unchange_gene$id),]
write.table(fdft1,'/Users/tingyue/Desktop/upload/SE.MATS.JC.fdft1.txt',quote = F,row.names = F,sep = '\t')
#先把results_transcripts和all_sig的交集找出来：

diff_trans_SE<-results_transcripts[which(results_transcripts$geneNames%in%all_sig$geneSymbol),]
plotMeans('stringtie_merged.16167', AS_filt,groupvar="feature",legend=T)
#########################
# 把结果写到csv文件
write.csv(results_transcripts, "chrX_transcript_results.csv", row.names=FALSE)
write.csv(results_genes, "chrX_gene_results.csv", row.names=FALSE)
# 赋予调色板五个指定颜色
tropical= c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
# 提取FPKM值
fpkm = texpr(bg_chrX,meas="FPKM")

#方便作图将其log转换，+1是为了避免出现log2(0)的情况
fpkm = log2(fpkm+1)
# 作图
boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')
# 就单个转录本的查看在样品中的分布
ballgown::transcriptNames(bg_chrX)[12]
ballgown::geneNames(bg_chrX)[12]
plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2), main=paste(ballgown::geneNames(bg_chrX)[12],' : ', ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab="Sex", ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))
#查看某一基因位置上所有的转录本
plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))
# 以性别为区分查看基因表达情况
plotMeans('MSTRG.575', bg_chrX_filt,groupvar="sex",legend=FALSE)