#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("singscore")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#install.packages("devtools")
#library(devtools)
#devtools::install_github('dviraran/SingleR')
#对细胞的注释

setwd("/")
#设置工作空间

library(limma)
library(Seurat)
library(dplyr)
library(magrittr)

rt <- read.table("GSE176509_APC_TPM.txt",sep="\t",header=T,check.names=F)
#读取不同格式文件时注意参数修改
rt <- as.matrix(rt)
rownames(rt) = rt[,1] 
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))

data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data <- avereps(data)
#如果基因存在多行中会取均值

pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)
#min.cells为基因存在样本最小数，需要根据实际情况选择，min.features = 50基因最小存在细胞数

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#使用PercentageFeatureSet函数计算线粒体基因的百分比

pdf(file="04.featureViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#保存对于基因特征的小提琴图
#第一个图是在样品中对所有细胞检测到的基因数量小提琴图，横坐标为样品名，纵坐标为每个细胞中包含的基因数量；
#第二个图为样品所有细胞的中基因序列数量小提琴图，横坐标为样品名，纵坐标为每个细胞基因的数量。
#第三个图为样品所有细胞的线粒体比例小提琴图，横坐标为样品名，纵坐标为每个细胞线粒体比例。
#小提琴图展示了任意位置的密度，可以知道哪些位置的密度较高。


pbmc <- subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 7) 
#对数据进行过滤，线粒体比例大于7%进行删掉，因为线粒体基因比例过高的细胞，会干扰细胞分群
pdf(file="04.featureCor.pdf",width=10,height=6) 
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()
#保存基因特征相关性图
#随着测序数量的增加，单细胞所检测到的基因数量也在增加
#这两个是有一定关联性的，在作图之前需要对含有线粒体基因的细胞进行筛选，因为过高会干扰细胞分群，我所选择的是高于7%的进行筛选，
#第一个图是对于线粒体比例于测序数据的关系的散点图即测序深度和线粒体基因含量的关系，
#第二个是对基因于测序数据序列的关系的散点图，即测序深度于基因数量的关系。


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#对数据进行标准化
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#提取那些在细胞间变异系数较大的基因
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)#top10是对基因
CombinePlots(plots = list(plot1, plot2))
dev.off()
#保存基因特征方差图
#第一个图，是进行筛选，是筛选基因表达最高的前1500个（表达量要根据实际情况进行更改）,在这个图中横坐标是对于基因在所有细胞的平均表达值，纵坐标是对于基因的离差值，离差值值越大表示基因的可变性越大。
#第二个图是对表达最显著的前十个基因进行标注。