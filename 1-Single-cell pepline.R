# ===== 加载R包 =====
library(data.table)    # 高效数据读写
library(tidyverse)     # 数据处理套件
library(Seurat)        # 单细胞分析核心框架
library(reshape2)      # 数据重塑
library(Matrix)        # 稀疏矩阵操作
library(dplyr)         # 数据操作
library(umap)          # UMAP降维
library(parallel)      # 并行计算
library(ggplot2)       # 可视化

#########################################################################
#                                                                       #
#                        Single cell RNA-seq                            #
#                  水稻幼苗单细胞RNA测序分析流程                           #
#                                                                       #
#########################################################################

# ===== 第1部分：通过bulk RNA-seq鉴定原生质体制备敏感基因 =====
# 原因：单细胞制备过程中需要酶解细胞壁，这会影响基因表达
# 目的：识别并排除受原生质体制备影响的基因，减少技术偏差

# 读取基因长度文件，用于计算TPM
# TPM = Transcripts Per Million，标准化基因表达量的指标
genelength <- read.table("Gene.length",sep="\t",header = T)[,c(1,6)]
head(genelength)  # 查看数据结构
str(genelength)   # 查看数据类型
genelength$Geneid <- as.character(genelength$Geneid)

# ----- 处理重复1的原生质体后(AP)数据 -----
# AP = After Protoplasting，原生质体制备后
# BP = Before Protoplasting，原生质体制备前
rep1_AP <- read.table("HTSeq_out_rep1_AP.txt")  # 读取HTSeq计数结果
rep1_AP <- rep1_AP[grep("AT",rep1_AP$V1),]       # 过滤，只保留拟南芥基因(AT开头)
rep1_AP$V1 <- as.character(rep1_AP$V1)
colnames(rep1_AP)[2] <- "AP.rep1.count"

# 合并基因长度信息
mer <- merge(rep1_AP,genelength,by.x="V1",by.y="Geneid")
head(mer)
colnames(mer)[1] <- "gene"

# 计算TPM的步骤：
# 1. RPK = Reads Per Kilobase = 计数 / 基因长度(kb)
# 2. TPM = RPK / 总RPK * 10^6
mer$rpk_AP <- mer[,2]/mer[,3]                    # 计算RPK
mer$AP.rep1 <- mer$rpk_AP*1e6/sum(mer[,4])       # 计算TPM
bulk_rice <- mer[,c("gene","AP.rep1")]
# ----- 处理重复1的原生质体前(BP)数据 -----
rep1_BP <- read.table("HTSeq_out_rep1_BP.txt")
rep1_BP <- rep1_BP[grep("AT",rep1_BP$V1),]
rep1_BP$V1 <- as.character(rep1_BP$V1)
colnames(rep1_BP)[2] <- "BP.rep1.count"
mer <- merge(rep1_BP,genelength,by.x="V1",by.y="Geneid")
colnames(mer)[1] <- "gene"
mer$rpk_BP <- mer[,2]/mer[,3]
mer$BP.rep1 <- mer$rpk_BP*1e6/sum(mer[,4])
bulk_rice <- merge(bulk_rice,mer[,c("gene","BP.rep1")],by="gene")

# ----- 处理重复2的原生质体后(AP)数据 -----
rep2_AP <- read.table("HTSeq_out_rep2_AP.txt")
rep2_AP <- rep2_AP[grep("AT",rep2_AP$V1),]
rep2_AP$V1 <- as.character(rep2_AP$V1)
colnames(rep2_AP)[2] <- "AP.rep2.count"
mer <- merge(rep2_AP,genelength,by.x="V1",by.y="Geneid")
colnames(mer)[1] <- "gene"
mer$rpk_AP <- mer[,2]/mer[,3]
mer$AP.rep2 <- mer$rpk_AP*1e6/sum(mer[,4])
bulk_rice <- merge(bulk_rice,mer[,c("gene","AP.rep2")],by="gene")

# ----- 处理重复2的原生质体前(BP)数据 -----
rep2_BP <- read.table("HTSeq_out_rep2_BP.txt")
rep2_BP <- rep2_BP[grep("AT",rep2_BP$V1),]
rep2_BP$V1 <- as.character(rep2_BP$V1)
colnames(rep2_BP)[2] <- "BP.rep2.count"
mer <- merge(rep2_BP,genelength,by.x="V1",by.y="Geneid")
colnames(mer)[1] <- "gene"
mer$rpk_BP <- mer[,2]/mer[,3]
mer$BP.rep2 <- mer$rpk_BP*1e6/sum(mer[,4])
bulk_rice <- merge(bulk_rice,mer[,c("gene","BP.rep2")],by="gene")

# 计算两个重复的平均TPM
bulk_rice$mean_AP <- rowMeans(bulk_rice[,c("AP.rep1","AP.rep2")])
bulk_rice$mean_BP <- rowMeans(bulk_rice[,c("BP.rep1","BP.rep2")])

# ===== 鉴定原生质体制备敏感基因 =====
# 计算fold change = AP/BP，比较原生质体制备前后的表达变化
bulk_rice$fc_rep1 <- bulk_rice$AP.rep1/bulk_rice$BP.rep1
head(bulk_rice)
bulk_rice$fc_rep2 <- bulk_rice$AP.rep2/bulk_rice$BP.rep2
head(bulk_rice)

# 绘制火山图可视化两个重复的一致性
ggplot(bulk_rice,aes(log2(fc_rep1),log2(fc_rep2)))+
  geom_point(size=.5)+my_theme

# 设定阈值，筛选敏感基因
# cutoff = 3 表示 |log2FC| > 3，即表达变化超过8倍
cutoff <- 3
bulk_rice$group<- "N"  # 默认不敏感
# 在两个重复中都表现出显著变化的基因被标记为敏感基因
bulk_rice$group[(log2(bulk_rice$fc_rep1)>cutoff&log2(bulk_rice$fc_rep2)>cutoff)|
                  (log2(bulk_rice$fc_rep1)< -cutoff&log2(bulk_rice$fc_rep2)< -cutoff)]<- "Y"

# ===== 第2部分：Seurat单细胞分析 =====

# ----- 读取单细胞表达矩阵 -----
# fread比read.table更快，适合大文件
all_integ <- as.data.frame(fread("all.csv"))
rownames(all) <- all$V1
all <- all[,-1]

# ----- 排除原生质体制备敏感基因 -----
# 这一步很关键，移除受技术处理影响的基因，保留生物学真实信号
all_integ <- all_integ[!(rownames(all_integ) %in% bulk_rice[bulk_rice$group=="Y",]$gene),]

# ----- 创建Seurat对象 -----
# 从列名中提取样本信息 (格式: barcode-sample)
sample <- separate(data = data.frame(cell=colnames(all_integ)), col = "cell",
                   into = c("barcode", "sample"), sep = "-")$sample
mydata <- CreateSeuratObject(counts = all_integ, project = "mydata_scRNAseq")
mydata@meta.data$sample <- sample

# 设置并行计算参数
future::plan("multiprocess", workers = 10)           # 使用10个核心
options(future.globals.maxSize = 80 * 1024^3)        # 内存限制80GB

# ----- SCTransform标准化 -----
# SCTransform是Seurat推荐的标准化方法，比传统方法更有效
# 它可以消除技术变异，同时保留生物学变异
rice.list <- SplitObject(mydata, split.by = "sample")  # 按样本分割
for (i in 1:length(rice.list)) {
  rice.list[[i]] <- SCTransform(rice.list[[i]], verbose = T)
}

# ----- 样本整合 -----
# 使用锚点整合不同样本，消除批次效应
rice.features <- SelectIntegrationFeatures(object.list = rice.list,nfeatures = 3000)
rice.list <- PrepSCTIntegration(object.list = rice.list, anchor.features = rice.features,
                                verbose = T)
reference_dataset <- which(names(rice.list) == "Ctrl1")  # 以Ctrl1为参考
rice.anchors <- FindIntegrationAnchors(object.list = rice.list, normalization.method = "SCT",
                                       anchor.features = rice.features, reference = reference_dataset)
rice.integrated <- IntegrateData(anchorset = rice.anchors, normalization.method = "SCT")
# ----- 降维分析 -----
# PCA (主成分分析) - 将高维数据投影到主要变异方向
rice.integrated <- RunPCA(object = rice.integrated, verbose = FALSE, npcs = 100)

# UMAP (统一流形近似与投影) - 非线性降维，用于可视化
# 使用Seurat内置UMAP
rice.integrated <- RunUMAP(object = rice.integrated, dims = 1:10,
                           min.dist = 0.05, n.neighbors = 5, seed.use = 100)

# 使用umap包进行更精细的UMAP分析
# metric="pearson2" 使用皮尔逊相关性作为距离度量
data <- as.data.frame(rice.integrated@reductions[["pca"]]@cell.embeddings)
umap <- umap::umap(data[,1:100],
                   n_neighbors=10,metric="pearson2",
                   min_dist=0.01,random_state=39)

# 提取UMAP坐标用于可视化
plt <- as.data.frame(umap$layout)
plt$sample <- sample
plt$cell <- rownames(plt)
data <- as.matrix(plt[colnames(rice.integrated),c("V1","V2")])
rownames(data) <- plt[colnames(rice.integrated),]$cell
colnames(data) <- c("UMAP_1","UMAP_2")
rice.integrated@reductions$umap@cell.embeddings <- data

# ----- 细胞聚类 -----
# 使用图聚类方法对细胞进行分群
DefaultAssay(rice.integrated) <- "integrated"
rice.integrated <- FindNeighbors(rice.integrated, reduction = "pca", dims = 1:100)  # 构建KNN图
rice.integrated <- FindClusters(rice.integrated, resolution = 0.75, n.start = 10)   # 聚类

# 重新排列聚类编号以匹配文献中的顺序
plt$cluster <- rice.integrated@meta.data$seurat_clusters
plt$cell <- rownames(plt)
# 映射原始聚类ID到新的编号顺序
plt <- merge(plt,data.frame(c2=c(0:28),stringsAsFactors = F,
                            cluster=c(24,13,3,23,7,21,17,1,10,9,8,5,6,14,27,20,19,22,4,12,2,0,15,26,11,18,25,16,28)),by="cluster")
plt$cluster <- plt$c2

# 绘制UMAP图，展示聚类结果
ggplot(plt,aes(x=V1,y=V2,color=cluster))+
  geom_point(alpha=0.5,size=.1)
# ===== 第3部分：Log标准化和数据准备 =====
# 对原始RNA数据进行log标准化，用于后续差异表达分析
# SCTransform适合整合，LogNormalize适合差异表达

# 创建新的Seurat对象用于log标准化
mydata_allgene <- CreateSeuratObject(counts = all_integ, project = "mydata_scRNAseq")
DefaultAssay(mydata_allgene) <- "RNA"

# Log标准化：将每个细胞的总计数归一化到10^6，然后取log
mydata_allgene <- NormalizeData(mydata_allgene, normalization.method = "LogNormalize",scale.factor = 1e6)

# 添加样本和聚类信息
mydata_allgene$sample <- plt[colnames(mydata_allgene),]$sample

# 识别高变基因 (HVGs)
mydata_allgene <- FindVariableFeatures(mydata_allgene, selection.method = "vst", nfeatures = 2000)

# 数据缩放：对所有基因进行标准化，使均值为0，方差为1
all.genes <- rownames(mydata_allgene)
mydata_allgene <- ScaleData(mydata_allgene, features = all.genes)

# 运行PCA和UMAP
mydata_allgene <- RunPCA(object = mydata_allgene, verbose = FALSE, npcs = 10)
mydata_allgene <- RunUMAP(object = mydata_allgene, dims = 1:10,
                          min.dist = 0.05, n.neighbors = 5, seed.use = 100)

# 将之前计算的UMAP坐标和聚类信息转移到新对象
data <- as.matrix(plt[colnames(mydata_allgene),c("V1","V2")])
rownames(data) <- plt[colnames(mydata_allgene),]$cell
colnames(data) <- c("UMAP_1","UMAP_2")
mydata_allgene@reductions$umap@cell.embeddings <- data
mydata_allgene@meta.data$seurat_clusters <- plt[colnames(mydata_allgene),]$cluster

# ===== 脚本完成 =====
# 输出对象：
# - rice.integrated: 整合后的Seurat对象 (SCTransform标准化)
# - mydata_allgene: Log标准化的Seurat对象 (用于差异表达)
# - plt: 包含UMAP坐标、样本和聚类信息的数据框