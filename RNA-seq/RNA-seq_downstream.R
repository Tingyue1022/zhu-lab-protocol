###########合并数据#############
KO1<-read.table('/Users/tingyue/Desktop/download/output2/S418KO_L2_339X39.count',sep = "\t",col.names = c("gene_id","KO1"))
KO2<-read.table('/Users/tingyue/Desktop/download/output2/S519KO_L2_340X40.count',sep = "\t",col.names = c("gene_id","KO2"))
KO3<-read.table('/Users/tingyue/Desktop/download/output2/S620KO_L2_341X41.count',sep = "\t",col.names = c("gene_id","KO3"))
WT1<-read.table('/Users/tingyue/Desktop/download/output2/S215WT_L2_337X37.count',sep = "\t",col.names = c("gene_id","WT1"))
WT2<-read.table('/Users/tingyue/Desktop/download/output2/S316WT_L2_338X38.count',sep = "\t",col.names = c("gene_id","WT2"))
WT3<-read.table('/Users/tingyue/Desktop/download/output2/S723WT_L2_342X42.count',sep = "\t",col.names = c("gene_id","WT3"))
#/Users/tingyue/Desktop/download/output2/S215WT_L2_337X37.count
#/Users/tingyue/Desktop/download/output2/S215WT_L2_337X37.count
#/Users/tingyue/Desktop/download/output2/S418KO_L2_339X39.count
#/Users/tingyue/Desktop/download/output2/S519KO_L2_340X40.count
#Users/tingyue/Desktop/download/output2/S620KO_L2_341X41.count
#/Users/tingyue/Desktop/download/output2/S723WT_L2_342X42.count


KO12<-merge(KO1,KO2,by='gene_id')
KO123<-merge(KO12,KO3,by='gene_id')
WT12<-merge(WT1,WT2,by='gene_id')
WT123<-merge(WT12,WT3,by='gene_id')
wdc<-merge(KO123,WT123,by='gene_id')
wdc<-wdc[-c(1,2,3,4,5),]
rm (KO1,KO2,KO3,WT1,WT2,WT3,KO12,KO123,WT12,WT123)
rownames(wdc)<-wdc$gene_id
wdc<-wdc[,-1]
class(wdc$KO1)



##################################input matrix######
non<-read.csv('/Users/zhulab/Desktop/download/matrix_htd.count',header = F,sep='')
colnames(non)=c('sample','gene','count')
library(reshape2)
non_counts=dcast(non,formula=gene~sample)
non_counts<-non_counts[-c(55488:55491  ),]
non_counts<-merge(ann3,non_counts,by='gene')
non_counts<-as.matrix(non_counts)
rownames(non_counts)<-non_counts[,2]
non_counts<-non_counts[,-c(1,2)]
non<-apply(non_counts,2,as.numeric)
rownames(non)<-rownames(non_counts)

















##################给DESeq2建立一个condition table###3
condition_table<-matrix(c('KO1','KO2','KO3',"WT1","WT2","WT3",'KO1','KO2','KO3',"WT1","WT2","WT3","CR2_neg","CR2_neg","CR2_neg","CR2_pos","CR2_pos","CR2_pos","IL17A_neg","IL17A_neg","IL17A_neg","IL17A_pos",'IL17A_pos','IL17A_pos'),ncol = 2,nrow = 12)
condition_table[,1]<-colnames(zwq)
View(condition_table)
rownames(condition_table)<-condition_table[,1]
condition_table<-condition_table[,-1]
condition_table<-as.data.frame(condition_table)
colnames(condition_table)
#[1] "condition_table"
colnames(condition_table)<-'feature'
condition_table$feature<-as.factor(condition_table$feature)
condition_table$feature<-relevel(condition_table$feature,ref = c("IL17A"))
mRNA_exprSet<-wdc[,c(1,2,3,4,5,6)]
library(DESeq2)

mRNA_exprSet<-non_counts[,c(4:9)]
mRNA_exprSet<-as.data.frame(mRNA_exprSet)
condition_table<-c[c(4:9),]
condition_table<-as.data.frame(condition_table)
rownames(condition_table)<-rownames(c)[c(4:9)]
colnames(condition_table)<-'feature'
condition_table$feature<-as.factor(condition_table$feature)
condition_table$feature<-relevel(condition_table$feature,ref = 'cov1')
#过滤表达量较低的基因，DESeq2分析之前要去除表达比较低的基因。
#在80%的样本中read counts的数量小于等于8则被过滤。
qualified_genes<-c()
for (gene_in_sheet in rownames(mRNA_exprSet)) {
  qualification<-mRNA_exprSet[gene_in_sheet,]<=8
  if(sum(qualification)<0.8*length(mRNA_exprSet)){
    qualified_genes<-append(qualified_genes,gene_in_sheet)
  }
  
}
mRNA_expr_for_DESeq<-mRNA_exprSet[qualified_genes,]

#use DESeq2 to calculate different expression gene. 对于几百个样本来说，这个过程可能是非常耗时的。
dds<-DESeqDataSetFromMatrix(mRNA_expr_for_DESeq,colData = condition_table,design = ~feature)
dds_DE<-DESeq(dds)
resultsNames(dds_DE)

res_DE<-results(dds_DE,alpha = 0.05,contrast = c('feature','CR2','IL17A'),name = 'feature_CR2_vs_IL17A')
#############################！！！！运行到这里，请接着到DESeq2.R这个文件的第81行开始运行。



for_volcano<-data.frame('log2FoldChange'=res_DE$log2FoldChange,
                        'padj'=res_DE$padj,
                        'pvalue'=res_DE$pvalue,
                        'descrip'<-rep('no',length(res_DE$log2FoldChange)),
                        'gene_name'<-rownames(res_DE))

rownames(for_volcano)<-rownames(res_DE)
names(for_volcano)[4]<-'descrip'
up_sig<-intersect(which(for_volcano$log2FoldChange>1),which(for_volcano$padj<0.05)) #满足这些要求的行名
down_sig<-intersect(which(for_volcano$log2FoldChange<(-1)),which(for_volcano$padj<0.05)) #
for_volcano$descrip<-as.character(for_volcano$descrip)
for_volcano[up_sig,'descrip']<-'Up' 
for_volcano[down_sig,'descrip']<-'Down' 
for_volcano$descrip<-as.factor(for_volcano$descrip)
colnames(for_volcano)[5]<-'gene'
for_volcano<-merge(for_volcano,ann3,by='gene_id')
EnhancedVolcano(for_volcano,
                lab = for_volcano$gene_name,
                x = 'log2FoldChange',
                FCcutoff = 1.0,
                y = 'padj',
                pCutoff=0.05,
                drawConnectors = F,
                labSize=3,
                title = "CoV2 vs. gfp",
                subtitle = " ",
                boxedLabels = F)
















library(org.Mm.eg.db)
library(clusterProfiler)

up_gene<-for_volcano[which(for_volcano$descrip%in%c('Down')),]
up_gene<-up_gene$gene
up_gene<-as.character(up_gene)
up_gene<-bitr(up_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')

ego<-enrichGO(gene = up_gene$ENTREZID,
              OrgDb = 'org.Mm.eg.db',
              ont = 'BP',
              pAdjustMethod = 'BH',
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)
ggo<-gseGO(geneList = glist,
           OrgDb = 'org.Mm.eg.db',
           ont = 'BP',
           minGSSize = 10,
           maxGSSize = 500,
           pvalueCutoff = 0.5)

dim(ego)
simplify(ego)
barplot(ego,showCategory=20,col="qvalue")
ego_hr<-data.frame(ego)
#kk<-enrichMKEGG(up_gene$ENTREZID,
#               organism = 'mmu',
#                keyType = 'kegg',
#                pAdjustMethod = 'BH',
#                pvalueCutoff = 0.05,
#                qvalueCutoff = 0.1)
#barplot(kk,showCategory = 10,title = 'KEGG_UP')
up_ego<-data.frame(ego)
barplot(ego,showCategory=11,col="qvalue")
go_up<-ego[c(1:9),c(1,2,8,6)]
####!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
go_down<-ego[c(1:8),c(1,2,8,6)]
GO<-rbind(go_up,go_down)
GO$geneID<-str_replace_all(GO$geneID,"/",",")
names(GO)<-c('ID',"Term","Genes",'adj_pval')
GO$Category<-'BP'

up_gene<-for_volcano[which(for_volcano$descrip==c('Down')|for_volcano$descrip==c('Up')),]
names(up_gene)[c(2,3,4,6)]<-c('logFC','adj.P.Val','P.Value','ID')
library(GOplot)
circ<-circle_dat(GO,up_gene)
GOBar(circ,display = 'multiple',)
GOBubble(circ, labels = 3)
barplot(ego,showCategory=11,col="qvalue")

##############PCA
res<-vst(dds_DE,blind = F)
exp<-assay(res)
par(cex = 0.7)
par(mar=c(1,3,1,1))
n.sample=ncol(exp)
cols <- rainbow(n.sample*1.2)
boxplot(exp,las=2,col=cols,)
PCA_data<-plotPCA(res,intgroup='feature',returnData=T)
plotPCA(res,intgroup='feature')
ggplot(PCA_data,aes(PCA_data$PC1,PCA_data$PC2,color=PCA_data$feature))+stat_ellipse(aes(fill=PCA_data$feature),type="norm",geom="polygon",alpha=0,color=NA)+guides(fill=F)+geom_point(size=4)+labs(x='PC1:57% variance',y='PC2:26% variance',color="")+geom_text(label=paste(rownames(PCA_data)),colour="black",size=0)

exp<-as.data.frame(exp)
exp$gene<-rownames(exp)
colnames(for_volcano)[5]<-'gene'
result_table2<-merge(for_volcano,exp,by="gene")

for_volcano2[which(for_volcano2$descrip%in%c('Down')),][which(for_volcano2[which(for_volcano2$descrip%in%c('Down')),]$gene%in%c(for_volcano[which(for_volcano$descrip%in%c('Down')),]$gene)),]$descrip

down_pos<-result_table[which(result_table$descrip%in%c('Up')),]
down_all<-result_table2[which(result_table2$descrip%in%c('Up')),]
down_unsig<-down_pos[which(down_pos$gene%in%c(down_all$gene)),]
table<-result_table
table$descrip<-as.character(table$descrip)
table[which(table$gene%in%c(down_unsig$gene)),]$descrip<-rep('Up_unsig',38)

#########GSEA


go_gsea<-read.gmt('/Users/tingyue/Downloads/c5.go.bp.v7.3.entrez.gmt')
gene<-for_volcano
gene <- dplyr::distinct(gene,gene$X.gene_name.....rownames.res_DE.,.keep_all=TRUE)
gene=bitr(gene,fromType="ENTREZID",toType="SYMBOL",OrgDb="org.Hs.eg.db")
genelist2<-gene[,c(6,2)]
colnames(genelist2)[1]<-'SYMBOL'
genelist<-merge(genelist,genelist2,by='SYMBOL')
genelist<-bitr(genelist$gene_symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
#genelist=sort(genelist$log2FoldChange,decreasing = T) 
genelist<-genelist[order(genelist$log2FoldChange,decreasing = T
                         ),]
genelist<-for_volcano[,c(6,2)]
genelist2<-for_volcano[,c(6,2)]
genelist<-bitr(genelist$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(genelist2)[1]<-'SYMBOL'
genelist<-merge(genelist,genelist2,by='SYMBOL')
glist<-genelist$log2FoldChange
names(glist)<-genelist$ENTREZID
glist <- sort(glist,decreasing = T)
ggo<-gseGO(geneList = glist,
           OrgDb = org.Hs.eg.db,
           ont = 'BP',
           minGSSize = 10,
           maxGSSize = 500,
           pvalueCutoff = 0.5,
           eps=0)


gseaplot2(ggo,2,color="red",pvalue_table = T,base_size = 14)
gseaplot2(ggo,1:10,base_size = 18)
gseaplot2(ggo,geneSetID=c('GO:0071396'),color="red",pvalue_table = T,base_size = 15)
ggo_hr<-data.frame(ggo)

ggo_hr$core_enrichment<-str_replace_all(ggo_hr$core_enrichment,"/","','")
ggo_hr2$core_enrichment<-str_replace_all(ggo_hr2$core_enrichment,"/","','")
wamup$geneID<-str_replace_all(wamup$geneID,"/","','")

wdc_exp[which(wdc_exp$ENTREZID%in%c('257632','20202','12775','14537','12986','245195','20310','11501','11689','12768','58218','20302','20305','17395','19204','19261','16409','16175','12977','56792','16414','106512','74748','16176','14131','20343','217303','57765','12771','67133','20299','16952','16170','246177','20963','17698','110168','12721','74734','12490','20868','378425','19354','16149','105855','20344','18578','27226','17329','12458','277360','24088','20308','13616','15139','20303','21838','14345','13051','18613','20345','19260','20737','20750','100952','14127','21936','13601','19277','16408','216892','12766','14254','193385','56212','17101','15894','19271','57781','15945','17874','64095','321019','71660','109700')),]



########################3
gsea_rxx<-GSEA(geneList = glist,TERM2GENE = go_gsea,pvalueCutoff = 0.05)
glist<-genelist[,3]
names(glist) <- as.character(genelist[,2])
glist <- sort(glist,decreasing = T)





