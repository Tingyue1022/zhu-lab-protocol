library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

# Load the PBMC dataset
data_1 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/1/filtered_feature_bc_matrix/")
data_3 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/3/filtered_feature_bc_matrix/")
data_2 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/2/filtered_feature_bc_matrix/")
data_4 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/4/filtered_feature_bc_matrix/")
data_5 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/5/filtered_feature_bc_matrix/")
data_8 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/8/filtered_feature_bc_matrix/")
data_7 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/7/filtered_feature_bc_matrix/")
data_9 <- Read10X(data.dir = "~/Desktop/download/sc_hkx/9/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data). min.cells这个的设置是怎么确定的？我看有的教程设置的是60啊？
sc1 <- CreateSeuratObject(counts = data_1, project = "1", min.cells = 3, min.features = 200,umi.assay = "RNA")
sc3 <- CreateSeuratObject(counts = data_3, project = "3", min.cells = 3, min.features = 200,umi.assay = "RNA")
sc2 <- CreateSeuratObject(counts = data_2, project = "2", min.cells = 3, min.features = 200,umi.assay = "RNA")
sc4 <- CreateSeuratObject(counts = data_4, project = "4", min.cells = 3, min.features = 200,umi.assay = "RNA")
sc5 <- CreateSeuratObject(counts = data_5, project = "5", min.cells = 3, min.features = 200,umi.assay = "RNA")
sc8 <- CreateSeuratObject(counts = data_8, project = "8", min.cells = 3, min.features = 200,umi.assay = "RNA")
sc7 <- CreateSeuratObject(counts = data_7, project = "7", min.cells = 3, min.features = 200,umi.assay = "RNA")
sc9 <- CreateSeuratObject(counts = data_9, project = "9", min.cells = 3, min.features = 200,umi.assay = "RNA")
#The values they picked here are somewhat arbitrary, but min.cells helps limit the number of genes used by removing those 
#unlikely to play any part in differentiating groups of cells due to being expressed in very few cells. In general, most 
#genes removed will be those with zero counts across all cells. min.features removes dead cells cells and empty droplets 
#where few genes are detected.
#总之感觉这个值不用太管，后面质控会进行更严格的筛选。
#1     2     3     4     5     7     8     9 
#10692  9563  9742  9528 11738 10497 11820 12878 
alldata@meta.data$group<-rep('non',53916)
alldata@meta.data[which(alldata@meta.data$orig.ident%in%c("1","3")),12]<-"a"
alldata@meta.data[which(alldata@meta.data$orig.ident%in%c("2","4")),12]<-"b"
alldata@meta.data[which(alldata@meta.data$orig.ident%in%c("5","8")),12]<-"c"
alldata@meta.data[which(alldata@meta.data$orig.ident%in%c("7","9")),12]<-"d"

alldata.list$aa@meta.data$group2<-rep('non',39723)
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("1","3")),21]<-"a"
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("2","4")),21]<-"b"
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("8")),21]<-"c"
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("9")),21]<-"d"
WT
KO

WT$type<-'WT'
KO$type<-'KO'

alldata<- merge(sc1, c(sc3,sc2,sc4,sc5,sc8,sc7,sc9), add.cell.ids = c("1","3","2","4","5","8","7","9"))#如果是有更多的组，在c()里面增加就完事了
alldata<- merge(sc1, c(sc3,sc2,sc4), add.cell.ids = c("1","3","2","4"))
alldata<- merge(sc5, c(sc8,sc7,sc9), add.cell.ids = c("5","8","7","9"))
data.filt<-merge(a_clean,c(b_clean),add.cell.ids = c("a","b"))

alldata <- PercentageFeatureSet(alldata, "^mt-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^Rp[sl]", col.name = "percent_ribo")
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.还有血小板
alldata <- PercentageFeatureSet(alldata, "^Hb[^(p)]", col.name = "percent_hb")
alldata <- PercentageFeatureSet(alldata, "Pecam1|Pf4", col.name = "percent_plat")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats,在WT这个Seurat文件中和的meta.data加了一列

grep ("^mt-", rownames(alldata[["RNA"]]),value = T)
grep ("^Rp[sl]", rownames(alldata[["RNA"]]),value = T)
grep ("^Hb[^(p)]", rownames(alldata[["RNA"]]),value = T)
grep ("Pecam1|Pf4", rownames(alldata[["RNA"]]),value = T)

#使用这个命令查看‘MT’是大写还是小写，核糖体的gene_symbol可能和网上的教程不太一样，不然线粒体出来的比例太低（是0）
#plotQC:画图检测测序的质量
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb","percent_plat")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.0000001, ncol = 3) + 
  NoLegend()
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)#看基因数量和读数总数的相关性
#开始过滤基因和细胞
selected_c <- WhichCells(alldata, expression = nFeature_RNA > 200)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]

data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt) #这里会滤掉一部分细胞和基因，比如刚开始WT7490个细胞，KO9111细胞，一共是16601个，基因是17279，这里过滤之后基因16839，细胞16545，但是我在CreateSeuratObject那一步设置的参数也是200和3啊？为什么这一步还是会过滤细胞呢？
table(data.filt$orig.ident)
# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
library(Matrix)
par(mar = c(4, 8, 2, 1))
C <- data.filt@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", 
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
#Malt1高属于正常现象
#过滤线粒体高的和没有核糖体表达的
selected_mito <- WhichCells(data.filt, expression = percent_mito < 20)
selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 0.05)

# and subset the object to only keep those cells
data.filt <- subset(data.filt, cells = selected_mito)
data.filt <- subset(data.filt, cells = selected_ribo)

dim(data.filt)

table(data.filt$orig.ident)#这部弄完之后KO剩下1484,WT剩下1382

VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()
##过滤掉线粒体基因，Malat1，线粒体基因等等
dim(data.filt)

# Filter MALAT1
data.filt <- data.filt[!grepl("Malat1", rownames(data.filt)), ]
dim(data.filt)
# Filter Mitocondrial
data.filt <- data.filt[!grepl("^mt-", rownames(data.filt)), ]
dim(data.filt)
# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
data.filt<- data.filt[ ! grepl('^Rp[sl]', rownames(data.filt)), ]
dim(data.filt)
# Filter Hemoglobin gene (optional if that is a problem on your data)
data.filt <- data.filt[!grepl("^Hb[^(p)]", rownames(data.filt)), ]
dim(data.filt)

#看眼KO和WT是不是同一个性别的，我这都没有Xist，应该都是雄的吧，其实应该看下Y染色体上基因的表达
VlnPlot(data.filt, features = c("Xist"))


data.filt <- NormalizeData(data.filt, normalization.method = "LogNormalize", scale.factor = 10000)

###############for mouse gene only###################

# Basic function to convert mouse to human gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

#There may be issues with timeout for the Biomart server, just keep trying,我他妈一直连不上啊，这图就没做

m.s.genes <- convertHumanGeneList(cc.genes$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)



data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes, 
                              s.features = cc.genes$s.genes)
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
        ncol = 4, pt.size = 0.1)

#Here, we will use DoubletFinder to predict doublet cells. But before doing doublet detection we need to run
#scaling, variable gene selection and pca, as well as UMAP for visualization. These steps will be explored in
#more detail in coming exercises.
suppressMessages(require(DoubletFinder))

#找到差异基因
data.filt = FindVariableFeatures(data.filt, verbose = T,selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data.filt), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#归一化 But before doing doublet detection we need to run scaling, variable gene selection and pca, as well as UMAP for visualization. 
#These steps will be explored in more detail in coming exercises.


data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"), 
                      verbose = T)
data.filt = RunPCA(data.filt, verbose = T, npcs = 20)
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = T)


library(DoubletFinder)
sweep.res <- paramSweep_v3(data.filt) 
sweep.stats <- summarizeSweep(sweep.res,GT = FALSE) 
bcmvn <- find.pK(sweep.stats) 
mpk<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
#这个推测出来的pk值在0.0005～0.1之间都是可以的


nExp <- round(ncol(data.filt) * 0.06)  # expect 6% doublets
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]



cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())
#We should expect that two cells have more detected genes than a single cell, lets check if our predicted 
#doublets also have more detected genes in general.
VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
dim(data.filt)
table(data.filt$orig.ident) #这样处理完KO剩下1400，WT剩下1323。。。第二次跑线粒体剩下了10，上回是5，这回KO2356，WT1946

#############################QC结束，开始降维####################################
#############################QC结束，开始降维####################################

suppressWarnings(suppressMessages(data.filt <- FindVariableFeatures(data.filt, selection.method = "vst", 
                                                                  nfeatures = 2000, verbose = FALSE, assay = "RNA")))
top20 <- head(VariableFeatures(data.filt), 20)

LabelPoints(plot = VariableFeaturePlot(data.filt), points = top20, repel = TRUE)

data.filt <- ScaleData(data.filt, vars.to.regress = c("percent_mito", "nFeature_RNA"), 
                     assay = "RNA")

data.filt <- RunPCA(data.filt, npcs = 50, verbose = T)

plot_grid(ncol = 3, DimPlot(data.filt, reduction = "pca", group.by = "orig.ident", 
                            dims = 1:2), DimPlot(data.filt, reduction = "pca", group.by = "orig.ident", dims = 3:4), 
          DimPlot(data.filt, reduction = "pca", group.by = "orig.ident", dims = 5:6))

#To identify which genes (Seurat) or metadata paramters (Scater/Scran) contribute the most to each PC
VizDimLoadings(data.filt, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)

#We can also plot the amount of variance explained by each PC.
data.filt <- JackStraw(data.filt, num.replicate = 100)
data.filt<- ScoreJackStraw(data.filt, dims = 1:20)
JackStrawPlot(data.filt, dims = 1:20)

ElbowPlot(data.filt, reduction = "pca", ndims = 50)

#tSNE
data.filt <- RunTSNE(data.filt, reduction = "pca", dims = 1:12, perplexity = 50, max_iter = 1000, 
                   theta = 0.5, eta = 200, num_threads = 0)
plot_grid(ncol = 1, DimPlot(data.filt, reduction = "tsne", group.by = "c"))

#UMAP
data.filt <- RunUMAP(data.filt, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30, 
                   n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
DimPlot(data.filt, reduction = "umap", group.by = "orig.ident",dims = 1:2,split.by = 'group')
DimPlot(data.filt, reduction = "umap", group.by = "orig.ident",dims = 1:2)
data.filt <- RunUMAP(data.filt,reduction.name = "UMAP10_on_PCA", reduction = "pca", dims = 1:30, n.components = 10, n.neighbors = 25, 
                     n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
data.filt <- RunUMAP(data.filt,reduction.name = "UMAP14_on_PCA", reduction = "pca", dims = 1:30, n.components = 14, n.neighbors = 25, 
                     n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
plot_grid(ncol = 3, DimPlot(data.filt, reduction = "umap", group.by = "orig.ident") + 
            ggplot2::ggtitle(label = "UMAP_on_PCA"), 
          DimPlot(data.filt, reduction = "UMAP10_on_PCA",group.by = "orig.ident", dims = 1:2) + ggplot2::ggtitle(label = "UMAP10_on_PCA"), 
          DimPlot(data.filt, reduction = "UMAP14_on_PCA", group.by = "orig.ident", dims = 1:2) + 
            ggplot2::ggtitle(label = "UMAP14_on_PCA"))
#dims:Which dimensions to use as input features, used only if features is NULL.
#n.components:The dimension of the space to embed into.就是用了几个PC

#以上的降维手段是基于PCA的结果，但是我发现WT和KO之间好像有一群细胞没有重叠啊，试一下基于scaled data和graph的降维
data.filt <- RunUMAP(data.filt, reduction.name = "UMAP_on_ScaleData", features = data.filt@assays$RNA@var.features, 
                   assay = "RNA", n.components = 2, n.neighbors = 30, n.epochs = 200, min.dist = 0.3, 
                   learning.rate = 1, spread = 1)

# Build Graph
data.filt <- FindNeighbors(data.filt, reduction = "pca", graph.name = "SNN", assay = "RNA", 
                         k.param = 20, features = data.filt@assays$RNA@var.features)
# Run UMAP on a graph,我靠还要下载uamp-learn。。算了算了回头再说，下载好了重新启动一下试试
data.filt <- RunUMAP(data.filt, reduction.name = "UMAP_on_Graph", graph = "SNN", assay = "RNA")

p1 <- DimPlot(data.filt, reduction = "umap", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_PCA")
p2 <- DimPlot(data.filt, reduction = "UMAP_on_ScaleData", group.by = "orig.ident") + 
  ggplot2::ggtitle(label = "UMAP_on_ScaleData")
p3 <- DimPlot(alldata, reduction = "UMAP_on_Graph", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_Graph")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), p2 + NoLegend() + 
                                                 NoAxes(),  leg, nrow = 2), ncol = 1, widths = c(1))

myfeatures <- c("Cd3e", "Cd4", "Cd8a", "Nkg7", "Gnly", "Ms4a1", "Cd14", "Lyz", "Ms4a7", 
                "Fcgr3a", "Cd326", "Fcgr1a","Cd19","Epcam","Lgr5")
FeaturePlot(data.filt, reduction = "umap", dims = 1:2, features =c("Cd4","Cd8a"), ncol = 3, 
            order = T)

###################################use CCA for removing batch effects######################
###########################################################################################
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_03_integration.html#Session_Info

##################################graph clustering###############################################
##################################graph clustering###############################################
#library(rafalib)
#library(clustree)
data.filt <- FindNeighbors(data.filt, dims = 1:30, k.param = 60, prune.SNN = 1/15)
#这个FindNeighbor就是搞了个k-nn图出来，好像是unweighted graph

#take a look at the kNN graph
library(pheatmap)
pheatmap(data.filt@graphs$RNA_nn[1:200, 1:200], col = c("white", "black"), border_color = "grey90", 
         legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2)
#Once the graph is built, we can now perform graph clustering. The clustering is done respective to a resolution 
#which can be interpreted as how coarse you want your cluster to be. Higher resolution means higher number of clusters.

# Clustering with louvain (algorithm 1)
for (res in c(1,1.25,1.5)) {
  data.filt <- FindClusters(data.filt, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}

# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.
library(enrichplot)
plot_grid(ncol = 3, DimPlot(data.filt, reduction = "umap", group.by = "RNA_snn_res.1.5",label = T) + 
            ggtitle("louvain_1.5"), DimPlot(data.filt, reduction = "umap", group.by = "RNA_snn_res.1.25",label = T) + 
            ggtitle("louvain_1.25"), DimPlot(data.filt, reduction = "umap", group.by = "RNA_snn_res.1",label = T) + 
            ggtitle("louvain_1"))
DimPlot(data.filt, reduction = "umap", group.by = "RNA_snn_res.1.5",label = T)
#check clusters
library('clustree')
clustree(data.filt@meta.data, prefix = "RNA_snn_res.")

#QC看下不同组之间的bias（这组就是clustering之后得到的组
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb","percent_plat")
VlnPlot(data.filt, group.by = "RNA_snn_res.1", features = feats, pt.size = 0.01, ncol = 3) + 
  NoLegend()

############################DEG####################################
############################DEG####################################
#In single cell, differential expresison can have multiple functionalities such as of identifying marker genes for cell populations,
#as well as differentially regulated genes across conditions (healthy vs control).

# Set the identity as louvain with resolution 0.5
sel.clust = "RNA_snn_res.1.5_rerange"

alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)

# plot this clustering
plot_grid(ncol = 2, DimPlot(data.filt,group.by = "RNA_snn_res.0.5",label = T) + NoAxes(), DimPlot(data.filt, group.by = "group") + 
            NoAxes())

markers_genes <- FindAllMarkers(alldata, logfc.threshold = 0.2, test.use = "wilcox", 
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, 
                                assay = "RNA")

#top 25 DEG across different clusters and plot them
top50 <- markers_genes %>% group_by(cluster) %>% top_n(-50, p_val_adj)
top50
library(rafalib)
mypar(2, 5, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F), horiz = T, 
          las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
  abline(v = c(0, 0.25), lty = c(1, 2))
}
#这画图不work啊，用热图看吧。。

#画个热图，每个cluster找五个padj最小的，第14cluster排名一样的比较多，所以基因的数量就比较多
top5 <- markers_genes %>% group_by(cluster) %>% top_n(-5, p_val_adj)
# create a scale.data slot for the selected genes
data.filt <- ScaleData(data.filt, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(data.filt, features = as.character(unique(top5$gene)), group.by = sel.clust, 
          assay = "RNA")
DotPlot(data.filt, features = rev(as.character(unique(top5$gene))), group.by = sel.clust, 
        assay = "RNA") + coord_flip()
#marker gene小提琴图
top1 <- top5 %>% group_by(cluster) %>% top_n(-1, p_val)
VlnPlot(data.filt, features = as.character(unique(top1$gene)), ncol = 8, group.by = sel.clust, 
        assay = "RNA", pt.size = 0)

#我看了下热图刚画的，有的cluster之间有点像啊，找找cluster之间的差异看下能不能合并
#DEG&GO analysis

sel.clust = "RNA_snn_res.1.5_rerange"

alldata <- SetIdent(alldata, value = sel.clust)
markers_diff <- FindMarkers(alldata, ident.1 = c('04'), ident.2 = c('02','03','05'), min.pct=0.25)
markers_diff_sig <- markers_diff[ markers_diff$p_val_adj < 0.05, ]

FeaturePlot(data.filt, reduction = "umap", dims = 1:2, features = c('Ppia'), ncol = 3, 
            order = T)
#GSEA is better for a single pathway

markers_diff_sig<-markers_diff

markers_diff_sig$gene_name<-rownames(markers_diff_sig)
genelist<-markers_diff_sig[,c(6,2)]
genelist2<-markers_diff_sig[,c(6,2)]
genelist<-bitr(genelist$gene_name,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
colnames(genelist2)[1]<-'SYMBOL'
genelist<-merge(genelist,genelist2,by='SYMBOL')
glist<-genelist$avg_log2FC
names(glist)<-genelist$ENTREZID
glist <- sort(glist,decreasing = T)
ggo<-gseGO(geneList = glist,
           OrgDb = org.Mm.eg.db,
           ont = 'BP',
           minGSSize = 10,
           maxGSSize = 500,
           pvalueCutoff = 0.999999,
           eps=0)


gseaplot2(ggo,430,color="red",pvalue_table = T,base_size = 14)
gseaplot2(ggo,1:10,base_size = 18)
gseaplot2(ggo,geneSetID=c('GO:0071396'),color="red",pvalue_table = T,base_size = 15)
ggo_hr<-data.frame(ggo)

which(ggo_hr$ID=="GO:0006915")
which(ggo_hr$ID=="GO:0008632")



#GO analysis
library(org.Mm.eg.db)
library(clusterProfiler)
up_gene<-rownames(markers_diff_sig)
up_gene<-as.character(up_gene)
up_gene<-bitr(up_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')
ego<-enrichGO(gene = up_gene$ENTREZID,
              OrgDb = 'org.Mm.eg.db',
              ont = 'BP',
              pAdjustMethod = 'BH',
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              readable = T)
dim(ego)
simplify(ego)
barplot(ego,showCategory=20,col="qvalue")
markers_diff
#不信的话画个图看下
FeaturePlot(data.filt, reduction = "umap", dims = 1:2, features = c("Nlrp6" , "Nlrp9b" , "Cgas" , "Gsdmc"  , "Cd80","Cd86","Casp3","Casp7" ), ncol = 3, 
            order = T,split.by = "orig.ident")
diff_3_8<-markers_diff

#下面开始找组间的差异基因，这个组间指的就是WT和KO组之间的

#对于组间的差异表达基因，你得一个cluster一个cluster找，我先吧cluster0给他提出来
cell_selection <- subset(data.filt, cells = colnames(data.filt)[data.filt@meta.data[, sel.clust] ==0])
cell_selection <- SetIdent(cell_selection, value = "type")
# Compute differentiall expression
DGE_cell_selection2 <- FindAllMarkers(cell_selection, logfc.threshold = 0.2, test.use = "DESeq2", 
                                     min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, 
                                     assay = "RNA")


####################################predict cell types#############
# Load the human marker table
markers <- read.delim("~/Desktop/singlecell_project/Mouse_cell_markers.txt")
markers <- markers[markers$speciesType == "Mouse", ]
markers <- markers[markers$cancerType == "Normal", ]

# Filter by tissue (to reduce computational time and have tissue-specific
# classification) 
sort(unique(markers$tissueType))
#grep('blood',unique(markers$tissueType),value = T) 
markers <- markers [markers$tissueType %in% c('Blood','blood vessel',"Basilar membrane","Colon","Colon epithelium",'Lymph node',"Epithelium",
                                              "Gastrointestinal tract","Ileum","Intestinal crypt","Intestine","Lymphoid tissue","Mesenteric lymph node",
                                              "Peyer patch","Small intestine"), ]

celltype_list <- lapply(unique(markers$cellName), function(x) {
  x <- paste(markers$geneSymbol[markers$cellName == x], sep = ",")
  x <- gsub("[[]|[]]| |-", ",", x)
  x <- unlist(strsplit(x, split = ","))
  x <- unique(x[!x %in% c("", "NA", "family")])
})
names(celltype_list) <- unique(markers$cellName)
# celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]} )
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) < 100]
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) > 1]

# run differential expression in our dataset, using clustering at resolution 2
data.filt <- SetIdent(data.filt, value = sel.clust)
DGE_table <- FindAllMarkers(data.filt, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, 
                            assay = "RNA")

# split into a list
DGE_list <- split(DGE_table, DGE_table$cluster)


# run fgsea for each of the clusters in the list
library(fgsea)
res <- lapply(DGE_list, function(x) {
  gene_rank <- setNames(x$avg_log2FC, x$gene)
  fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank, nperm = 10000)
  return(fgseaRes)
})
names(res) <- names(DGE_list)

# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
  x[x$pval < 0.01, ]
})
res <- lapply(res, function(x) {
  x[x$size > 5, ]
})
res <- lapply(res, function(x) {
  x[order(x$NES, decreasing = T), ]
})

# show top 3 for each cluster.
lapply(res, head, 3)

##手动注释

FeaturePlot(data.filt, reduction = "umap", dims = 1:2, features = c('Thy1','Il2ra','Cd4','Cd8a','Cd160','Dapk2','Cd200r4','Nkg7','Cd244a'), ncol = 3, 
            order = T)
gene<-c("Ciita","Cd74","H2-Aa","Il10")
gene<-"Il10"
FeaturePlot(alldata.list$aa, reduction = "umap", dims = 1:2, features = gene, ncol = 1,order = T,split.by = 'group2')
FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = c("Cd74","H2-Aa","H2-Ab1","H2-Eb1"), ncol = 1,order = T,cols = c("lightgrey", "red"))
FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = c("Cd74","H2-Aa","H2-Ab1","H2-Eb1"), ncol = 2,order = T)
FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = c( 'Gsdma', 'Gsdma2','Gsdma3','Gsdmb','Gsdmc1','Gsdmc2','Gsdmc3','Gsdmc4',"Gsdmd","Gsdme",'Nlrp6','Nlrp9b','Casp3','Casp7','Casp1','Irf7'), ncol = 4,order = T)

alldata.list2 <- SplitObject(alldata.list$bb, split.by = "orig.ident")
alldata.list3 <- SplitObject(alldata.list$aa, split.by = "orig.ident")

FeaturePlot(alldata.list2$`8`, reduction = "umap", dims = 1:2, features = c( 'Epcam','Ptprc','Gsdmd','Casp3','Casp7','Casp1'), ncol = 3,order = T,keep.scale = "feature")

FeaturePlot(alldata.list2$`8`, reduction = "umap", dims = 1:2, features = c( 'Epcam','Ptprc','Gsdmd','Casp3','Casp7','Casp1'), ncol = 3,order = T,keep.scale = "feature",min.cutoff = 2)

FeaturePlot(alldata.list2$`8`, reduction = "umap", dims = 1:2, features = c("Cd74","H2-Aa","H2-Ab1","H2-Eb1"), ncol = 2,order = T)
FeaturePlot(alldata.list2$`7`, reduction = "umap", dims = 1:2, features = c("H2-Aa","Cd74","H2-Eb1","Ciita"), ncol = 2,order = T)
FeaturePlot(alldata.list3$`4`, reduction = "umap", dims = 1:2, features = c("H2-Aa","Cd74","H2-Eb1","Ciita"), ncol = 2,order = T)
FeaturePlot(alldata.list3$`3`, reduction = "umap", dims = 1:2, features = c("H2-Aa","Cd74","H2-Eb1","Ciita"), ncol = 2,order = T)

FeaturePlot(alldata.iec, reduction = "umap", dims = 1:2, features = c("Cd80","Cd86","Cd40","Icosl",'Cd274'), ncol = 3,order = T)

FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = c("Epcam"), ncol = 1,order = T, cols = c("lightgrey", "red"),keep.scale = NULL)

FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = c( 'Gsdma', 'Gsdma2','Gsdma3'), ncol = 3,order = T)
FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = c('Lag3'), ncol = 2,order = T,split.by = "orig.ident")
FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = c("Ido1"), ncol = 1,order = T)
FeaturePlot(alldata.list$aa, reduction = "umap", dims = 1:2, features = c('Tgfb1'), ncol = 1,order = T,split.by = 'group2')
new.cluster.ids <- c("0 Enterocyte (SI distal)", "1 Intestinal epithelial cell", "2 Enterocyte", "3 Programmed death?", 
                     "4 Acinar cell","5 Cholangiocyte", "6 Paneth cell", "7 Intestinal stem cell", "8 Goblet cell", 
                     "9 Goblet cell_2", "10 Tuft cell", "11 Tuft cell_2", "12 Enterochromaffin cell", 
                     "13 CD8+ T cell", "14 CD4+ T cell", "15 CD8+ T_eff", "16 Follicular B cell (late stage)","17 Plasma cell","18 Plasma cell_2","19 Antigen activated B cell","20 Plasma cell_3",
                     "21 Follicular B cell","22 Macrophage","23 Plasmacytoid dendritic cell","24 Fibroblast")
#new.cluster.ids <- c("0 Enterocyte_Cyp2b10 high", "1 Enterocyte", "2 Stem cell", "3 B cell", 
#                     "4 T cell","5 Enterocyte", "6 Enterocyte", "7 T cell", "8 Epithelial cell", 
#                     "9 Goblet cell", "10 B cell", "11 Macrophage", "12", 
#                     "13 Fibroblast", "14", "15 B cell", "16","17","18 Enteroendocrine","19","20",
#                     "21","22 T cell","23","24 Tuft cell")
names(new.cluster.ids) <- levels(data.filt)
data.filt.plot <- RenameIdents(data.filt, new.cluster.ids)
DimPlot(alldata, label = T,cols =c("grey95","grey", "darkred", "red", "black","blue","green","pink","orange", "yellow","grey20", "yellow green", "green", "skyblue", "blue", "purple","green", "skyblue", "blue", "purple","green", "skyblue", "blue", "purple","blue", "purple"))
VlnPlot(data.filt, features = "Il10",split.by = "group", 
        assay = "RNA", pt.size = 0.01,)
VlnPlot(alldata, features = c("Casp3", "Casp7","Gsdmd","Casp1"), 
        assay = "RNA", pt.size = 0)
#casp3 casp7 gsdmd casp1

DimPlot(alldata,group.by = "RNA_snn_res.1.5_rerange",label = T,cols ='alphabet')
DimPlot(alldata,group.by = "RNA_snn_res.1.5_rerange",label = T,cols ='alphabet2')
DimPlot(alldata,group.by = "RNA_snn_res.1.5_rerange",label = T,cols ='glasbey')
DimPlot(alldata,group.by = "RNA_snn_res.1.5_rerange",label = T,cols ='polychrome')
DimPlot(alldata,group.by = "RNA_snn_res.1.5_rerange",label = T,cols ='stepped')
#
alldata.list$aa@meta.data$group2<-rep('non',39723)
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("1","3")),21]<-"a"
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("2","4")),21]<-"b"
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("8")),21]<-"c"
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("9")),21]<-"d"
#
data.filt@meta.data$label<-rep('non',53916)
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("0")),20]<-"Enterocyte (SI distal)"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("1")),20]<-"Intestinal epithelial cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("6")),20]<-"Enterocyte"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("5")),20]<-"Programmed death?"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("8")),20]<-"Acinar cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("12")),20]<-"Cholangiocyte"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("14")),20]<-"Paneth cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("2")),20]<-"Intestinal stem cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("9")),20]<-"Goblet cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("21")),20]<-"Goblet cell_2"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("24")),20]<-"Tuft cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("23")),20]<-"Tuft cell_2"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("18")),20]<-"Enterochromaffin cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("4")),20]<-"CD8+ T cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("7")),20]<-"CD4+ T cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("22")),20]<-"CD8+ T_eff"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("3")),20]<-"Follicular B cell (late stage)"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("10")),20]<-"Plasma cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("15")),20]<-"Plasma cell_2"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("17")),20]<-"Antigen activated B cel"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("19")),20]<-"Plasma cell_3"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("16")),20]<-"Follicular B cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("11")),20]<-"Macrophage"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("20")),20]<-"Plasmacytoid dendritic cell"
data.filt@meta.data[which(data.filt@meta.data$RNA_snn_res.0.5%in%c("13")),20]<-"Fibroblast"


new.cluster.ids <- c("1 Intestinal stem cell", "2 Enterocyte", "3 Intestinal epithelial barrier cell (Reg1+, lower crypt)", "4 Enterocyte", 
                     "5 Enterocyte","6 Goblet cell", "7 Goblet cell (non-canonical, Hes1+, Muc2-)", "8 Paneth cell", "9 Tuft cell", 
                     "10 Enterochromaffin cell", "11 CD8aa T cell", "12 CD8+ T cell", "13 Activated CD8+ T cell", 
                     "14 CD4+ T cell", "15 NK cell", "16 NK1.1+ CD8+ T cell", "17 Plasma cell_1","18 Germinal center B cell/early plasmablast (Ki67+)",
                     "19 Plasma cell_2","20 Plasma cell_3","21 Follicular B cell",
                     "22 Macrophage","23 Plasmacytoid dendritic cell","24 Fibroblaast","25 Pancreas cell","26 Apoptotic cell")
names(new.cluster.ids) <- levels(alldata)
data.filt.plot <- RenameIdents(alldata, new.cluster.ids)
DimPlot(data.filt.plot, label = T,cols ='polychrome')
alldata@meta.data$RNA_snn_res.1.5_rerange<-rep("non",25899)
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("7")),22]<-"01"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("8")),22]<-"01"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("6")),22]<-"02"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("10")),22]<-"02"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("4")),22]<-"03"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("2")),22]<-"04"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("5")),22]<-"04"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("24")),22]<-"05"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("15")),22]<-"06"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("29")),22]<-"07"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("20")),22]<-"08"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("30")),22]<-"09"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("26")),22]<-"10"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("0")),22]<-"11"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("14")),22]<-"12"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("13")),22]<-"13"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("16")),22]<-"14"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("19")),22]<-"15"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("28")),22]<-"16"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("3")),22]<-"17"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("9")),22]<-"17"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("11")),22]<-"17"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("23")),22]<-"18"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("18")),22]<-"19"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("22")),22]<-"20"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("25")),22]<-"21"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("12")),22]<-"22"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("27")),22]<-"23"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("17")),22]<-"24"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("21")),22]<-"25"
alldata@meta.data[which(alldata@meta.data$RNA_snn_res.1.5%in%c("1")),22]<-"26"
DimPlot(alldata, label = T,cols ='polychrome')



data.filt@meta.data$split<-rep('non',53916)
data.filt@meta.data[which(data.filt@meta.data$orig.ident%in%c("3")),20]<-"aa"
data.filt@meta.data[which(data.filt@meta.data$orig.ident%in%c("4")),20]<-"aa"
data.filt@meta.data[which(data.filt@meta.data$orig.ident%in%c("7")),20]<-"aa"
data.filt@meta.data[which(data.filt@meta.data$orig.ident%in%c("8")),20]<-"aa"

#############DEG_ann##################################shiji
DGE_table <- FindAllMarkers(data.filt, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, 
                            assay = "RNA")
DGE_list <- split(DGE_table, DGE_table$cluster)
unlist(lapply(DGE_list, nrow))

test_list<-split(markers_genes,markers_genes$cluster)
unlist(lapply(test_list, nrow))


##don't run!!! check the difference between two "FindAllMarkers" paramrters!
markers_genes <- FindAllMarkers(data.filt, logfc.threshold = 0.2, test.use = "wilcox", 
                               min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, 
                                assay = "RNA")
ref<-read.csv('/Users/tingyue/Library/Containers/com.tencent.xinWeChat/Data/Library/Application\ Support/com.tencent.xinWeChat/2.0b4.0.9/3fdd4a306d1ecc6b382cdff9ed879555/Message/MessageTemp/eeb76e3b1fda5224e0094e823e8a29af/File/Combined.combined.integrated_resolution0.8.csv')
reference_markers<-ref
reference_markers <- reference_markers[order(reference_markers$avg_log2FC, decreasing = T),]
top50_cell_selection <- reference_markers %>% group_by(cluster) %>% top_n(-100, p_val) %>% 
  top_n(50, avg_log2FC)

# Transform the markers into a list
ref_list = split(top50_cell_selection$gene, top50_cell_selection$cluster)

unlist(lapply(ref_list, length))

#GSEA
suppressPackageStartupMessages(library(fgsea))

# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
  gene_rank <- setNames(x$avg_log2FC, x$gene)
  fgseaRes <- fgsea(pathways = ref_list, stats = gene_rank, nperm = 10000)
  return(fgseaRes)
})
names(res) <- names(DGE_list)

# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
  x[x$pval < 0.1, ]
})
res <- lapply(res, function(x) {
  x[x$size > 2, ]
})
res <- lapply(res, function(x) {
  x[order(x$NES, decreasing = T), ]
})
res

new.cluster.ids <- c("0", "1", "2", "3", 
                     "CD4+ T cell","5", "6", "7", "8", 
                     "9", "DC", "11", "Paneth cell", 
                     "B cell", "Intestinal stem cell", "15", "16")
names(new.cluster.ids) <- levels(data.filt)
data.filt.plot <- RenameIdents(data.filt, new.cluster.ids)
DimPlot(data.filt.plot, label = T,)
VlnPlot(data.filt.plot$group%in%c("c","d"), features = "Ciita",split.by = "group", 
        assay = "RNA", pt.size = 0)

###################################################################




WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^mt-")
KO[["percent.mt"]] <- PercentageFeatureSet(KO, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(KO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(KO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(KO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

FeaturePlot(data.filt, reduction = "umap", dims = 1:2, 
            features = c('Il10'), ncol = 4,order = T,split.by = "orig.ident")



group (a/b/c/d)
total)percent (cluster4+cluster7+cluster22)/(total cell number)
4)percent cluster4/(cluster4+cluster7+cluster22)
7)percent cluster7/(cluster4+cluster7+cluster22)
22)percent cluster22/(cluster4+cluster7+cluster22)

a
total)18.16% 2164/11918
4)76.85% 1663/2164
7)21.44% 464/2164
22)1.71% 37/2164

b
total)16.80% 1846/10991
4)71.18% 1314/1846
7)26.44% 488/1846
22)2.38% 44/1846

c
total)17.30% 2650/15314
4)60.00% 1590/2650
7)38.15% 1011/2650
22)1.849% 49/2650

d
total)11.48% 1802/15693
4)54.94% 990/1802
7)44.06% 794/1802
22)0.999% 18/1802



dim(data.filt@meta.data[which(data.filt@meta.data$group=="d"),])
dim(data.filt@meta.data[which(data.filt@meta.data$group=="d"&data.filt@meta.data$RNA_snn_res.0.5=="4"),])
dim(data.filt@meta.data[which(data.filt@meta.data$group=="d"&data.filt@meta.data$RNA_snn_res.0.5=="7"),])
dim(data.filt@meta.data[which(data.filt@meta.data$group=="d"&data.filt@meta.data$RNA_snn_res.0.5=="22"),])
dim(data.filt@meta.data[which(data.filt@meta.data$group=="d"&data.filt@meta.data$RNA_snn_res.0.5%in%c('4','7','22')),])

#split and plot
data.filt.plot<-alldata
data.filt.plot@meta.data$split<-rep("non",25899)

data.filt.plot@meta.data[which(data.filt.plot@meta.data$orig.ident%in%c("3","4")),20]<-"aa"
data.filt.plot@meta.data[which(data.filt.plot@meta.data$orig.ident%in%c("7","8")),20]<-"bb"
alldata.list <- SplitObject(data.filt.plot, split.by = "split")

#11474
alldata.list$aa@meta.data$split_c<-rep("non",11474)
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("3")),23]<-"2"
alldata.list$aa@meta.data[which(alldata.list$aa@meta.data$orig.ident%in%c("4")),23]<-"1"

alldata.list$bb@meta.data$split_c<-rep("non",14425)
alldata.list$bb@meta.data[which(alldata.list$bb@meta.data$orig.ident%in%c("7")),23]<-"6"
alldata.list$bb@meta.data[which(alldata.list$bb@meta.data$orig.ident%in%c("8")),23]<-"5"

VlnPlot(alldata.list$aa, features = "H2-Eb1",split.by = "split_c",combine  = T,split.plot = T,
        assay = "RNA", pt.size = 0,cols = c("skyblue","coral"),idents = c("01","02","03","04","05","17","18","19","20","21","22","23"))
VlnPlot(alldata.list$bb, features = "H2-Eb1",split.by = "split_c",combine  = T,split.plot = T,
        assay = "RNA", pt.size = 0,cols = c("#6CA6CD","#CD5B45"),idents = c("01","02","03","04","05","17","18","19","20","21","22","23"))




FeaturePlot(alldata.list$aa,features = "Il10",label = F,order = T)
FeaturePlot(alldata.list$aa,split.by = "orig.ident",features = c("Pdcd1","Cd274","Ifng"),label = F,order = T)
VlnPlot(alldata.list$bb, features = "Il10",split.by = "group",combine  = T,split.plot = T,
        assay = "RNA", pt.size = 0.0001)
FeaturePlot(alldata.list$bb,features = "Ciita",label = F,order = T)
FeaturePlot(alldata.list$bb,split.by = "group",features = "Il10",label = F,order = T)

alldata.list$bb@meta.data$split2<-rep("non",31007)
alldata.list$bb@meta.data[which(alldata.list$bb@meta.data$orig.ident%in%c("5","9")),21]<-"aaa"
alldata.list$bb@meta.data[which(alldata.list$bb@meta.data$orig.ident%in%c("8","7")),21]<-"bbb"
alldata.list2 <- SplitObject(alldata.list$bb, split.by = "split2")
VlnPlot(alldata, features = "Cd74",split.by = "orig.ident",combine  = T,split.plot = T,
        assay = "RNA", pt.size = 0.0001)
alldata.list$aa@meta.data$split2<-rep("non",22909)


#calculate positive cell number(ratio)
gene<-c("Cd74", "H2-Aa")
p1<-FeaturePlot(alldata.list$bb,features = gene,label = F,order = T)
p2<-FeaturePlot(alldata.list$bb,split.by = "group",features = gene,label = F,order = T)
FeaturePlot(data.filt,split.by = "group",features = "H2-Eb1",label = F,order = T)
exprs <- data.frame(FetchData(object = alldata, vars=c("Il10","Foxp3","Ciita","Cd80","Cd86","Cd40","Cd274","Icosl","Cd70","Lag3")))
exprs$id<-rownames(exprs)
meta<-alldata@meta.data
meta$id<-rownames(meta)
exprs<-merge(exprs,meta[,c(1,12,22,23)],by="id")


ab<-c("3","4","7","8")
cd<-c("c","d","5","8","7","9")
gene<-"Foxp3"
gr<-ab
par(mfrow=c(2,2))
#par(cex=1.5,font=2)
#layout(matrix(c(1,1, 2, 2, byrow = TRUE)))
#1
a<-dim(exprs[which(exprs$orig.ident==paste0(gr[1])&exprs$RNA_snn_res.1.5_rerange=="14"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[1])&exprs$RNA_snn_res.1.5_rerange=="14"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
a1<-a[1]
b1<-b[1]
r1<-(a1-b1)
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p1<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("#8B475D","#668B8B"),init.angle = 90,main = "Gsdmd f/f vil-cre")

a<-dim(exprs[which(exprs$orig.ident==paste0(gr[2])&exprs$RNA_snn_res.1.5_rerange=="14"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[2])&exprs$RNA_snn_res.1.5_rerange=="14"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
a1<-a[1]
b1<-b[1]
r1<-(a1-b1)
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p1<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("#8B475D","#668B8B"),init.angle = 90,main = "Gsdmd f/f")

a<-dim(exprs[which(exprs$orig.ident==paste0(gr[3])&exprs$RNA_snn_res.1.5_rerange=="14"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[3])&exprs$RNA_snn_res.1.5_rerange=="14"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
a1<-a[1]
b1<-b[1]
r1<-(a1-b1)
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p1<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("#8B475D","#668B8B"),init.angle = 90,main = "Gsdmd D88A")

a<-dim(exprs[which(exprs$orig.ident==paste0(gr[4])&exprs$RNA_snn_res.1.5_rerange=="14"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[4])&exprs$RNA_snn_res.1.5_rerange=="14"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
a1<-a[1]
b1<-b[1]
r1<-(a1-b1)
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p1<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("#8B475D","#668B8B"),init.angle = 90,main = "Gsdmd +/+")

chisq.test(matrix(c(4,14,62,114),nrow = 2,ncol = 2),correct=F)
chisq.test(matrix(c(17,10,183,128),nrow = 2,ncol = 2),correct=F)
####################################################

#2
a<-dim(exprs[which(exprs$group==paste0(gr[2])&exprs$RNA_snn_res.0.5=="0"),])
b<-dim(exprs[which(exprs$group==paste0(gr[2])&exprs$RNA_snn_res.0.5=="0"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
a2<-a[1]
b2<-b[1]
r2<-(a2-b2)
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p2<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("coral","gray85"),init.angle = 90,main = paste0(gr[2]))

#3

a<-dim(exprs[which(exprs$orig.ident==paste0(gr[3])&exprs$RNA_snn_res.0.5=="0"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[3])&exprs$RNA_snn_res.0.5=="0"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p3<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("coral","gray90"),init.angle = 90,main = paste0(gr[3]))


#4
a<-dim(exprs[which(exprs$orig.ident==paste0(gr[4])&exprs$RNA_snn_res.0.5=="0"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[4])&exprs$RNA_snn_res.0.5=="0"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p4<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("black","gray85"),init.angle = 90,main = paste0(gr[4]))

#5

a<-dim(exprs[which(exprs$orig.ident==paste0(gr[5])&exprs$RNA_snn_res.0.5=="0"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[5])&exprs$RNA_snn_res.0.5=="0"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p5<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("black","gray85"),init.angle = 90,main = paste0(gr[5]))


#6
a<-dim(exprs[which(exprs$orig.ident==paste0(gr[6])&exprs$RNA_snn_res.0.5=="0"),])
b<-dim(exprs[which(exprs$orig.ident==paste0(gr[6])&exprs$RNA_snn_res.0.5=="0"&exprs[,paste0(gene)]>0),])
b[1]/a[1]
n<-round(b[1]/a[1],5)*100
m<-(100-n)
p6<-pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("black","gray85"),init.angle = 90,main = paste0(gr[6]))

s<-chisq.test(matrix(c(b1,b2,r1,r2),nrow = 2,ncol = 2))
round(s$p.value,5)















































a1
b1
r1
a2
b2
r2




































s<-chisq.test(matrix(c(12,21,999,773),nrow = 2,ncol = 2))
round(s$p.value,5)
text(50,50,"Pearson's Chi-squared Test P=paste0(num)")


a<-dim(exprs[which(exprs$group=="d"&exprs$RNA_snn_res.0.5=="7"),])
b<-dim(exprs[which(exprs$group=="d"&exprs$RNA_snn_res.0.5=="7"&exprs[,paste0(gene)]>0),])
21/794

a<-dim(exprs[which(exprs$orig.ident=="5"&exprs$RNA_snn_res.0.5=="7"),])
b<-dim(exprs[which(exprs$orig.ident=="5"&exprs$RNA_snn_res.0.5=="7"&exprs$Foxp3>0),])
b[1]/a[1]
78/1195
0.06527197


a<-dim(exprs[which(exprs$orig.ident=="7"&exprs$RNA_snn_res.0.5=="7"),])
b<-dim(exprs[which(exprs$orig.ident=="7"&exprs$RNA_snn_res.0.5=="7"&exprs$Foxp3>0),])
b[1]/a[1]
236/1237
0.1907842
a<-dim(exprs[which(exprs$orig.ident=="8"&exprs$RNA_snn_res.0.5=="7"),])
b<-dim(exprs[which(exprs$orig.ident=="8"&exprs$RNA_snn_res.0.5=="7"&exprs$Foxp3>0),])
b[1]/a[1]
79/1307
0.06044376
a<-dim(exprs[which(exprs$orig.ident=="9"&exprs$RNA_snn_res.0.5=="7"),])
b<-dim(exprs[which(exprs$orig.ident=="9"&exprs$RNA_snn_res.0.5=="7"&exprs$Foxp3>0),])
b[1]/a[1]
142/849
0.1672556

#pie chart
par(cex=1.3,font=2)
n<-2.016
m<-(100-n)
pie(c(n,m),labels = c(paste0(n,"%"),NA),col = c("black","gray85"),init.angle = 90,main = "9")
legend("right",legend = c("Il10 pos","Il10 neg"),fill = c("black","gray85"))

#chi-sq test
s<-chisq.test(matrix(c(12,21,999,773),nrow = 2,ncol = 2))
round(s$p.value,5)

ggpie(data=c(n,m),x=c(n,m),labels = c(paste0(n,"%"),NA))












#chipseq code, run for now
chip15<-readPeakFile('/Users/tingyue/Desktop/download/chip15_r_peaks.narrowPeak')
chip16<-readPeakFile('/Users/tingyue/Desktop/download/chip16_r_peaks.narrowPeak')
library('TxDb.Hsapiens.UCSC.hg38.knownGene')
txdb <- TxDb.Mmusculus.UCSC.hg38.knownGene
covplot(chip15,weightCol = 5)


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(chip15, windows=promoter)

data<-alldata@meta.data[,c(1,22)]
data$RNA_snn_res.1.5_rerange<-factor(data$RNA_snn_res.1.5_rerange)
data$number <- 1
data <- ddply(data,'orig.ident',transform,percent = 1/sum(number)*100)
colName = c("black","grey", "darkred", "red", "black","blue","green","black","pink","grey","orange", "yellow","grey20", "yellow green", "green", "skyblue", "blue", "purple")#无限制，随便添
colors = colorRampPalette(colName)(26)
ggplot(data,aes(orig.ident,percent,fill=cell))+
  geom_bar(stat="identity",position="stack")+
  ggtitle("")+
  theme_classic()+
  theme(axis.ticks = element_line(size=0.5),
        axis.text.y = element_text(size=10,hjust = 0.5),
        axis.text.x = element_text(size=10,hjust = 0.9,angle = 45,vjust = 0.9),
        axis.title.y = element_text(size=12,hjust = 0.5),
        axis.title.x = element_text(size=12,hjust = 0.5))+
  guides(fill=guide_legend(title=NULL))+scale_fill_manual(values = colors)


colName=c("skyblue","red","blue","coral","skyblue3","darkred","coral3","blueviolet",'brown1','deepskyblue')
#cell type
data$cell<-rep("non",25899)
data[which(data$RNA_snn_res.1.5%in%c("01")),3]<-"Intestinal stem cell"
data[which(data$RNA_snn_res.1.5%in%c("02","03","04","05")),3]<-"Enterocyte"
data[which(data$RNA_snn_res.1.5%in%c("06","07","08")),3]<-"Goblet cell & Paneth cell"
data[which(data$RNA_snn_res.1.5%in%c("09")),3]<-"Tuft cell"
data[which(data$RNA_snn_res.1.5%in%c("10")),3]<-"Enterochromaffin cell"
data[which(data$RNA_snn_res.1.5%in%c("11","12","13","14","16")),3]<-"T cell"
data[which(data$RNA_snn_res.1.5%in%c("15")),3]<-"NK cell"
data[which(data$RNA_snn_res.1.5%in%c("17","18","19","20","21")),3]<-"B cell"
data[which(data$RNA_snn_res.1.5%in%c("22")),3]<-"Macrophage"
data[which(data$RNA_snn_res.1.5%in%c("23")),3]<-"Plasmacytoid dendritic cell"
data[which(data$RNA_snn_res.1.5%in%c("24")),3]<-"Fibroblaast"
data[which(data$RNA_snn_res.1.5%in%c("25")),3]<-"Pancreas cell"
data[which(data$RNA_snn_res.1.5%in%c("26")),3]<-"Apoptotic cell"
data$cell<-factor(data$cell,levels = c("Intestinal stem cell","Enterocyte","Goblet cell & Paneth cell","Tuft cell","Enterochromaffin cell",
                                       "T cell","NK cell","B cell","Macrophage","Plasmacytoid dendritic cell","Fibroblaast","Pancreas cell",
                                       "Apoptotic cell"))

#Plasmacytoid dendritic cell","24 Fibroblaast","25 Pancreas cell","26 Apoptotic cell"
#change the order of the clusters 

#Seurat 
transfer.anchors <- FindTransferAnchors(reference = data.filt, query = alldata, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = data.filt$label, 
                            dims = 1:30)
alldata2 <- AddMetaData(object = alldata, metadata = predictions)








##a script to cauculae positve cell percentage of each cluster
num<-c("01","02","03","04","05","17","18","19","20","21","22","23")
#para
gene<-"Ciita"
gr1<-"7"
gr2<-"8"
per_1<-c()
per_2<-c()
pval<-c()

for (i in c(1:12)) {print(paste0(num[i]));
  a<-dim(exprs[which(exprs$orig.ident==paste0(gr1)&exprs$RNA_snn_res.1.5_rerange==paste0(num[i])),]);
  b<-dim(exprs[which(exprs$orig.ident==paste0(gr1)&exprs$RNA_snn_res.1.5_rerange==paste0(num[i])&exprs[,paste0(gene)]>0),]);
  a1<-a[1];b1<-b[1];r1<-(a1-b1);n<-round(b[1]/a[1],5)*100;
  per_1<-append(per_1,n);
  a<-dim(exprs[which(exprs$orig.ident==paste0(gr2)&exprs$RNA_snn_res.1.5_rerange==paste0(num[i])),]);
  b<-dim(exprs[which(exprs$orig.ident==paste0(gr2)&exprs$RNA_snn_res.1.5_rerange==paste0(num[i])&exprs[,paste0(gene)]>0),]);
  a2<-a[1];b2<-b[1];r2<-(a2-b2);n<-round(b[1]/a[1],5)*100;
  per_2<-append(per_2,n);
  s<-chisq.test(matrix(c(a1,a2,r1,r2),nrow = 2,ncol = 2),correct = F);
  p<-round(s$p.value,5);
  pval<-append(pval,p)
  
  }

tab<-as.data.frame(num)
tab<-rbind(tab,tab)

#!something change
tab$orig.idents<-c(rep("2",length(num)),rep("1",length(num)))
tab$percent<-c(per_1,per_2)
ggplot(tab,aes(x=num,y=percent)) + geom_bar(stat="identity", aes(fill=factor(orig.idents)),position="dodge")+scale_fill_manual(values = c("skyblue3","coral3"))+
  theme_classic()+labs(x="",y='Percent')+theme(axis.ticks = element_line(size=0.5),
                             axis.text.y = element_text(size=10,hjust = 0.5),
                             axis.text.x = element_text(size=10,hjust = 0.9,angle = 45,vjust = 0.9),
                             axis.title.y = element_text(size=12,hjust = 0.5),
                             axis.title.x = element_text(size=12,hjust = 0.5))
#####################################################positve cell number:CD4,CD8,il10+cd4+,cd4+helios+foxp3+
gr<-c("3","4","7","8")
cd4<-c()
cd8<-c()
cd4_il10<-c()
cd4_foxp3<-c()
#for (i in c(1:4)) {print(paste0(gr[i]));
#  a<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs$RNA_snn_res.1.5_rerange==c("13")),]);#CD4
#  b<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs$RNA_snn_res.1.5_rerange%in%c("10","11","12","15")),]);#CD8
#  c<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs$RNA_snn_res.1.5_rerange==c("13")&exprs[,"Il10"]>0),]);#CD4+Il10+
#  d<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs$RNA_snn_res.1.5_rerange==c("13")&exprs[,"Foxp3"]>0&exprs[,"Ikzf2"]>0),]);#CD4+Foxp3+helios+
#  cd4<-append(cd4,a[1]);
# cd8<-append(cd8,b[1]);
# cd4_il10<-append(cd4_il10,c[1]);
# cd4_foxp3<-append(cd4_foxp3,d[1])
  
  
#}
for (i in c(1:4)) {print(paste0(gr[i]));
  a<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs[,"Cd4"]>0),]);#CD4
  b<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs[,"Cd8a"]>0),]);#CD8
  c<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs[,"Cd4"]>0&exprs[,"Il10"]>0),]);#CD4+Il10+
  d<-dim(exprs[which(exprs$orig.ident==paste0(gr[i])&exprs[,"Cd4"]>0&exprs[,"Foxp3"]>0&exprs[,"Ikzf2"]>0),]);#CD4+Foxp3+helios+
  cd4<-append(cd4,a[1]);
  cd8<-append(cd8,b[1]);
  cd4_il10<-append(cd4_il10,c[1]);
  cd4_foxp3<-append(cd4_foxp3,d[1])
  
  
}
cd4<-as.data.frame(cd4)
cd8<-as.data.frame(cd8)
cd4_il10<-as.data.frame(cd4_il10)
cd4_foxp3<-as.data.frame(cd4_foxp3)
cell_num<-cbind(cd4,cd8,cd4_il10,cd4_foxp3)
###############################monocle3##############
cds <- as.cell_data_set(alldata)
cds <- cluster_cells(cds,k=20,reduction_method = "UMAP",partition_qval = 0.05,resolution = 1)
p1 <- plot_cells(cds, show_trajectory_graph = F)
p2 <- plot_cells(cds, color_cells_by = "RNA_snn_res.1.5_rerange", show_trajectory_graph = FALSE,labels_per_group = 1)
p2
wrap_plots(p1, p2)

alldata.cd4t<-subset(data.filt,RNA_snn_res.1.5_rerange%in%c("14"))
alldata.t<-subset(alldata,RNA_snn_res.1.5_rerange%in%c("14"))
#Seurat process
data.filt<-alldata.t
data.filt <- NormalizeData(data.filt, normalization.method = "LogNormalize", scale.factor = 10000)
suppressWarnings(suppressMessages(data.filt <- FindVariableFeatures(data.filt, selection.method = "vst", 
                                                                    nfeatures = 2000, verbose = FALSE, assay = "RNA")))
data.filt <- ScaleData(data.filt, vars.to.regress = c("percent_mito", "nFeature_RNA"), 
                       assay = "RNA")
data.filt <- RunPCA(data.filt, npcs = 50, verbose = T)
data.filt <- RunUMAP(data.filt, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30, 
                     n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)


DimPlot(data.filt, reduction = "umap", split.by = "orig.ident",dims = 1:2)
DimPlot(data.filt, reduction = "umap", group.by = "RNA_snn_res.1.5_rerange",dims = 1:2,cols = "polychrome")
DimPlot(alldata.cd4t, reduction = "umap", group.by = "RNA_snn_res.1.5_rerange",dims = 1:2,cols = "polychrome")
DimPlot(data.filt, reduction = "umap", group.by = "RNA_snn_res.1.5_rerange",dims = 1:2,cols = "polychrome",split.by = "orig.ident")

FeaturePlot(data.filt, reduction = "umap", dims = 1:2, features = c("Foxp3","Il10"), ncol = 2,order = T,split.by = "orig.ident")
FeaturePlot(alldata.cd4t, reduction = "umap", dims = 1:2, features = c("Il10"), ncol = 2,order = T,split.by = "orig.ident")
###############
alldata.list <- SplitObject(data.filt, split.by = "orig.ident")
###############
cds<-as.cell_data_set(alldata.list$`8`)
#
#cds <- preprocess_cds(cds, num_dim = 50)
#cds <- reduce_dimension(cds)
#plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "RNA_snn_res.1.5")
#
plot_cells(cds, show_trajectory_graph = F,color_cells_by = "RNA_snn_res.1.5_rerange")
cds <- cluster_cells(cds,k=20,reduction_method = "UMAP")
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "RNA_snn_res.1.5_rerange",
           label_groups_by_cluster=FALSE,
           label_leaves=T,
           label_branch_points=T)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)



num<-c("01","02","03","04","05")
#para
gene<-"Lag3"
gr1<-"8"
gr2<-"4"
per_1<-c()
per_2<-c()
pval<-c()

for (i in c(1:5)) {print(paste0(num[i]));
  a<-dim(exprs[which(exprs$orig.ident==paste0(gr1)&exprs$RNA_snn_res.1.5_rerange==paste0(num[i])),]);
  b<-dim(exprs[which(exprs$orig.ident==paste0(gr1)&exprs$RNA_snn_res.1.5_rerange==paste0(num[i])&exprs[,paste0(gene)]>0),]);
  a1<-a[1];b1<-b[1];r1<-(a1-b1);n<-round(b[1]/a[1],5)*100;
  per_1<-append(per_1,n);
  
}

lag3<-per_1
tab<-as.data.frame(num)
tab<-rbind(tab,tab)

#!something change
tab$orig.idents<-c(rep("2",length(num)),rep("1",length(num)))
tab$percent<-c(per_1,per_2)
ggplot(num,aes(x=num,y=percent)) + geom_bar(stat="identity", aes(fill=factor(feature,levels = c("Cd80","Cd86","Cd40","Lag3","Cd274"))),position="dodge")+scale_fill_manual(values = c("red","coral","skyblue","skyblue3","coral3"))+
  theme_classic()+labs(x="",y='Percent')+theme(axis.ticks = element_line(size=0.5),
                                               axis.text.y = element_text(size=10,hjust = 0.5),
                                               axis.text.x = element_text(size=10,hjust = 0.9,angle = 45,vjust = 0.9),
                                               axis.title.y = element_text(size=12,hjust = 0.5),
                                               axis.title.x = element_text(size=12,hjust = 0.5))



cd70<-as.data.frame(cd70)
lag3<-as.data.frame(lag3)
eee<-cbind(cd80,cd86,cd40,cd274,lag3)
num<-c(rep("01",5),rep("02",5),rep("03",5),rep("04",5),rep("05",5))
num<-as.data.frame(num)
num$feature<-rep('non',25)
num$feature[1:5]<-c('Cd80',"Cd86","Cd40","Cd274","Lag3")
num$feature[6:10]<-c('Cd80',"Cd86","Cd40","Cd274","Lag3")
num$feature[11:15]<-c('Cd80',"Cd86","Cd40","Cd274","Lag3")
num$feature[16:20]<-c('Cd80',"Cd86","Cd40","Cd274","Lag3")
num$feature[21:25]<-c('Cd80',"Cd86","Cd40","Cd274","Lag3")

num$percent<-rep(1,25)
num$percent[1:5]<-as.numeric(eee[1,])
num$percent[6:10]<-as.numeric(eee[2,])
num$percent[11:15]<-as.numeric(eee[3,])
num$percent[16:20]<-as.numeric(eee[4,])
num$percent[21:25]<-as.numeric(eee[5,])


#transfer
transfer.anchors <- FindTransferAnchors(reference = reference, query = ctrl, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = reference$RNA_snn_res.1.5_rerange, 
                            dims = 1:30)
ctrl <- AddMetaData(object = ctrl, metadata = predictions)

DimPlot(ctrl, group.by = "predicted.id", label = T, repel = T,cols = 'polychrome') 




#ref
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Session_Info
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_02_dim_reduction.html#Session_Info
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_03_integration.html#Session_Info
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_04_clustering.html#Session_Info
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_05_dge.html#Session_Info
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_06_celltype.html#Session_Info
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Session_Info







