# ===== 加载R包 =====
library(data.table)
library(tidyverse)
library(Seurat)
library(reshape2)
library(Matrix)
library(dplyr)
library(parallel)

###########################################################################
#                                                                         #
#        Determination of the organ origin in the Seedling samples        #
#                  鉴定幼苗样本中细胞的器官来源                               #
#                                                                         #
###########################################################################

# ===== 第1部分：筛选叶偏向基因和根偏向基因 =====
# 目的：找出在叶或根中特异性高表达的基因，用于判断细胞来源

# ----- 批次1分析 -----
# 设置细胞标识为样本
Idents(mydata_allgene) <- "sample"

# 提取批次1的样本：S_A1(叶) 和 S_R1(根)
ae.ro1@meta.data$Sample <- substr(ae.ro1@meta.data$sample,1,3)  # 提取样本类型(A或R)

# 创建聚类-样本组合标识 (例如: "0_S_A", "0_S_R")
ae.ro1$cluster.stress <- paste(ae.ro1$seurat_clusters, ae.ro1$Sample, sep = "_")
Idents(ae.ro1) <- "cluster.stress"
unique(ae.ro1$cluster.stress)

DefaultAssay(ae.ro1) <- "RNA"  # 使用RNA assay进行差异表达
ae.ro1$cluster.stress <- as.character(ae.ro1$cluster.stress)

# 定义差异表达分析函数
# 使用Wilcoxon检验比较同一聚类中叶和根细胞的表达差异
fun <- function(x){
  response.ae.ro1 <- FindMarkers(ae.ro1, slot="data",test.use = "wilcox",
                                 ident.1 = as.character(paste0(x,"_S_A")),   # 叶细胞
                                 ident.2 = as.character(paste0(x,"_S_R")),   # 根细胞
                                 verbose = T, logfc.threshold = 0, pseudocount.use=1)
  response.ae.ro1$cluster <- x
  response.ae.ro1$gene <- rownames(response.ae.ro1)
  rownames(response.ae.ro1) <- NULL
  return(response.ae.ro1)
}

# 并行计算所有聚类的差异表达 (使用28核)
clu <- c(0:22,26:28)  # 排除23, 24, 25聚类
deg.ae.ro1 <- do.call('rbind',parallel::mclapply(clu,function(x){fun(x)},mc.cores = 28))
deg.ae.ro1$batch <- "batch1"
# ----- 批次2分析 -----
# 同样的流程处理第二批样本
Idents(mydata_allgene) <- "sample"
ae.ro2 <- subset(mydata_allgene,idents = c("S_A2","S_R2"))
ae.ro2@meta.data$Sample <- substr(ae.ro2@meta.data$sample,1,3)
ae.ro2$cluster.stress <- paste(ae.ro2$seurat_clusters, ae.ro2$Sample, sep = "_")
Idents(ae.ro2) <- "cluster.stress"
DefaultAssay(ae.ro2) <- "RNA"
ae.ro2$cluster.stress <- as.character(ae.ro2$cluster.stress)

# 定义差异表达函数
fun <- function(x){
  response.ae.ro2 <- FindMarkers(ae.ro2, slot="data",test.use = "wilcox",
                                 ident.1 = as.character(paste0(x,"_S_A")),
                                 ident.2 = as.character(paste0(x,"_S_R")),
                                 verbose = T, logfc.threshold = 0, pseudocount.use=1)
  response.ae.ro2$cluster <- x
  response.ae.ro2$gene <- rownames(response.ae.ro2)
  rownames(response.ae.ro2) <- NULL
  return(response.ae.ro2)
}

# 批次2包含聚类25
clu <- c(0:22,25:28)
deg.ae.ro2 <- do.call('rbind',parallel::mclapply(clu,function(x){fun(x)},mc.cores = 28))
table(deg.ae.ro2$cluster)
deg.ae.ro2$batch <- "batch2"

# 将批次2的聚类25结果补充到批次1
deg.ae.ro1 <- rbind(deg.ae.ro1,deg.ae.ro2[deg.ae.ro2$cluster==25,])
deg.ae.ro1[deg.ae.ro1$cluster==25,]$batch <- "batch1"

# 合并两个批次的结果
deg.ae.ro <- rbind(deg.ae.ro1,deg.ae.ro2)
colnames(deg.ae.ro)[3:4] <- c("pct.aerial","pct.root")  # 重命名列

# 标记基因表达方向：ae = aerial(叶), ro = root(根)
deg.ae.ro$direction <- "ae"
deg.ae.ro[deg.ae.ro$avg_logFC <0,]$direction <- "ro"

# 创建基因-聚类唯一标识
deg.ae.ro$gene_cluster <- paste0(deg.ae.ro$gene,"-",deg.ae.ro$cluster)

# 计算log2 fold change (基于平均表达量)
deg.ae.ro$log2fc_avg <- log2(exp(deg.ae.ro$avg_logFC))

# 处理低表达比例：将低于0.1的比例设为0.1，避免log计算问题
prop <- 0.1
deg.ae.ro[deg.ae.ro$avg_logFC>0&deg.ae.ro$pct.root<0.1,]$pct.root <- prop
deg.ae.ro[deg.ae.ro$avg_logFC<0&deg.ae.ro$pct.aerial<0.1,]$pct.aerial <- prop

# 计算基于表达比例的log2FC
deg.ae.ro$log2fc_pct <- log2(deg.ae.ro$pct.aerial/deg.ae.ro$pct.root)

# 综合评分 = 平均表达FC + 表达比例FC
# 分数越高，基因越偏向叶；分数越低，越偏向根
deg.ae.ro$score <- deg.ae.ro$log2fc_avg+deg.ae.ro$log2fc_pct
head(deg.ae.ro)

# ===== 第2部分：筛选在胁迫条件下不变的基因 =====
# 目的：排除受胁迫影响的基因，保留器官特异性基因
# 胁迫条件：HS(热胁迫), ID(缺铁), LN(低氮), Ct(对照)

clu <- c(0:22,25:28)
deg.ae.ro$Ct <- 0
for (i in clu){
  for(j in c("batch1","batch2")){
    deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$Ct <- log2(rowMeans(expm1(mydata_allgene@assays$RNA@data[deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$gene,plt[plt$cluster==i&plt$Sample=="Ct",]$cell]))+1)
  }
}
deg.ae.ro$HS <- 0
for (i in clu){
  for(j in c("batch1","batch2")){
    deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$HS <- log2(rowMeans(expm1(mydata_allgene@assays$RNA@data[deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$gene,plt[plt$cluster==i&plt$Sample=="HS",]$cell]))+1)
  }
}
deg.ae.ro$ID <- 0
for (i in clu){
  for(j in c("batch1","batch2")){
    deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$ID <- log2(rowMeans(expm1(mydata_allgene@assays$RNA@data[deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$gene,plt[plt$cluster==i&plt$Sample=="ID",]$cell]))+1)
  }
}
# 整合两个批次的结果，筛选显著差异基因 (p < 0.01)
deg.twobatch <- reshape2::dcast(deg.ae.ro[deg.ae.ro$p_val<0.01,],gene_cluster+cluster+gene+Ct+HS+ID+LN+direction~batch,value.var = "score")
deg.twobatch <- na.omit(deg.twobatch)

# 计算两个批次的平均分数
deg.twobatch$mean_score <- rowMeans(deg.twobatch[,c("batch1","batch2")])

# 鉴定胁迫响应基因
# 条件：任一胁迫条件与对照差异 > 1 (log2FC)
deg.twobatch$is_stressresponse <-"Y"
deg.twobatch[abs(deg.twobatch$HS-deg.twobatch$Ct)<1&
               abs(deg.twobatch$ID-deg.twobatch$Ct)<1&
               abs(deg.twobatch$LN-deg.twobatch$Ct)<1, ]$is_stressresponse <-"N"
table(deg.twobatch$is_stressresponse)

# ===== 第3步：选择器官特异性标记基因 =====
# 每个聚类选择top 3个叶偏向或根偏向基因 (排除胁迫响应基因)

num <- 3

# 叶偏向标记基因
marker.ae <- deg.twobatch[deg.twobatch$direction=="ae"&deg.twobatch$is_stressresponse=="N",] %>%
  group_by(cluster) %>% top_n(num,mean_score) %>% as.data.frame()
length(marker.ae$gene)
length(unique(marker.ae$gene))

# 根偏向标记基因
marker.ro <- deg.twobatch[deg.twobatch$direction=="ro"&deg.twobatch$is_stressresponse=="N",] %>%
  group_by(cluster) %>% top_n(num,abs(mean_score)) %>% as.data.frame()
length(marker.ro$gene)
length(unique(marker.ro$gene))

# 标记选中的差异表达基因
deg.twobatch$is_deg <- "N"
deg.twobatch[deg.twobatch$direction=="ae"&deg.twobatch$gene_cluster %in% marker.ae$gene_cluster,]$is_deg <- "Y"
deg.twobatch[deg.twobatch$direction=="ro"&deg.twobatch$gene_cluster %in% marker.ro$gene_cluster,]$is_deg <- "Y"
table(deg.twobatch$is_deg,deg.twobatch$cluster,deg.twobatch$direction)

deg.twobatch <- deg.twobatch[order(deg.twobatch$is_deg),]

# ===== 第3部分：混合样本的Wilcoxon检验 =====
# 将两个批次的样本混合进行差异表达分析

# 设置细胞标识
Idents(mydata_allgene) <- "sample"

# 提取所有幼苗样本
ae.ro <- subset(mydata_allgene,idents = c("S_A1","S_A2","S_R1","S_R2"))
ae.ro

# 创建样本类型标识
ae.ro@meta.data$Sample <- substr(ae.ro@meta.data$sample,1,3)

# 创建聚类-样本组合标识
ae.ro$cluster.stress <- paste(ae.ro$seurat_clusters, ae.ro$Sample, sep = "_")
Idents(ae.ro) <- "cluster.stress"
unique(ae.ro$cluster.stress)

# 可视化UMAP图
DimPlot(ae.ro,reduction = 'umap',split.by = 'sample',label = T,group.by = 'seurat_clusters')
DefaultAssay(ae.ro) <- "RNA"
ae.ro$cluster.stress <- as.character(ae.ro$cluster.stress)

# 定义差异表达函数
fun <- function(x){
  response.ae.ro <- FindMarkers(ae.ro, slot="data",test.use = "wilcox",
                                ident.1 = as.character(paste0(x,"_S_A")),   # 叶细胞
                                ident.2 = as.character(paste0(x,"_S_R")),   # 根细胞
                                verbose = T, logfc.threshold = 0, pseudocount.use=1)
  response.ae.ro$cluster <- x
  response.ae.ro$gene <- rownames(response.ae.ro)
  rownames(response.ae.ro) <- NULL
  return(response.ae.ro)
}

# 并行计算差异表达
clu <- c(0:25,27:28)
deg.ae.ro.mix2rep <- do.call('rbind',parallel::mclapply(clu,function(x){fun(x)},mc.cores = 28))
table(deg.ae.ro.mix2rep$cluster)
head(deg.ae.ro.mix2rep)
dim(deg.ae.ro.mix2rep)

# 计算log2FC和方向
deg.ae.ro.mix2rep$log2fc_avg <- log2(exp(deg.ae.ro.mix2rep$avg_logFC))
deg.ae.ro.mix2rep$direction <- "ae"
deg.ae.ro.mix2rep[deg.ae.ro.mix2rep$avg_logFC <0,]$direction <- "ro"

# 创建基因-聚类唯一标识
deg.ae.ro.mix2rep$gene_cluster <- paste0(deg.ae.ro.mix2rep$gene,"-",deg.ae.ro.mix2rep$cluster)

# 标记选中的差异表达基因
deg.ae.ro.mix2rep$is_deg <- "N"
deg.ae.ro.mix2rep[deg.ae.ro.mix2rep$direction=="ae"&deg.ae.ro.mix2rep$gene_cluster %in% marker.ae$gene_cluster,]$is_deg <- "Y"
deg.ae.ro.mix2rep[deg.ae.ro.mix2rep$direction=="ro"&deg.ae.ro.mix2rep$gene_cluster %in% marker.ro$gene_cluster,]$is_deg <- "Y"
deg.ae.ro.mix2rep <- deg.ae.ro.mix2rep[order(deg.ae.ro.mix2rep$is_deg),]

# 创建样本分组标识
plt$sample2 <- plt$sample
plt[plt$sample %in% c("S_A1","S_A2","S_R1","S_R2"),]$sample2 <- "Small"

# 定义分数计算函数
# 对于每个聚类和样本，计算叶和根的平均Z-score表达量
plt.ae.ro <- data.frame(stringsAsFactors = F)
deg <- deg.ae.ro.mix2rep

score <- function(i) {
  tmp <- data.frame(stringsAsFactors = F)
  for (j in sample2){
    plt.tmp <- plt[plt$cluster==i&plt$sample2==j,]

    # 提取叶偏向基因表达
    gene.ae <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC>0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ae,plt.tmp$cell]))

    # 计算叶分数 (Z-score标准化后的平均表达)
    if (length(gene.ae)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))  # 行标准化(Z-score)
        dat[dat=="NaN"] <- min(dat[dat!="NaN"])  # 处理NaN
      } else {
        dat[dat==0] <- -5  # 无表达设为-5
      }
      plt.tmp$aerial <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$aerial <- scale(dat)
    } else {
      plt.tmp$aerial <- -5
    }
    }

    # 提取根偏向基因表达
    gene.ro <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC<0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ro,plt.tmp$cell]))

    # 计算根分数
    if (length(gene.ro)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))
        dat[dat=="NaN"] <- min(dat[dat!="NaN"])
      } else {
        dat[dat==0] <- -5
      }
      plt.tmp$root <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$root <- scale(dat)
    } else {
      plt.tmp$root <- -5
    }
    }

    # 归一化分数到非负范围
    plt.tmp$aerial <- plt.tmp$aerial-min(plt.tmp$aerial)
    plt.tmp$root <- plt.tmp$root-min(plt.tmp$root)

    tmp <- rbind(tmp,plt.tmp)
  }
  return(tmp)
}

# 并行计算主要聚类的分数
clu1 <- c(0:13,15:22,25,27:28)
plt.ae.ro <- do.call('rbind',parallel::mclapply(clu1,function(i){score(i)},mc.cores = 28))

# 标注样本来源
plt.ae.ro[plt.ae.ro$sample %in% c("S_A1","S_A2"),]$Sample <- "Aerial"
plt.ae.ro[plt.ae.ro$sample %in% c("S_R1","S_R2"),]$Sample <- "Root"

# 根据分数判断组织来源
plt.score.1 <- plt.ae.ro
plt.score.1$tissue <- "N"
plt.score.1[plt.score.1$root< plt.score.1$aerial,]$tissue <- "Aerial"  # 叶分数更高
# 对特殊聚类(14, 26)进行单独处理
# 这些聚类可能需要不同的样本组合
plt.ae.ro <- data.frame(stringsAsFactors = F)
deg <- deg.ae.ro.mix2rep

score <- function(i) {
  tmp <- data.frame(stringsAsFactors = F)
  # 特殊聚类使用所有胁迫条件样本
  sample3 <- c("Ctrl1","Ctrl2","Ctrl3","Ctrl4","ID1","ID2","ID3","ID4","LN1","LN2","LN3","LN4","Small")
  for (j in sample3){
    plt.tmp <- plt[plt$cluster==i&plt$sample2==j,]
    gene.ae <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC>0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ae,plt.tmp$cell]))
    if (length(gene.ae)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))
        dat[dat=="NaN"] <- min(dat[dat!="NaN"])
      } else {
        dat[dat==0] <- -5
      }
      plt.tmp$aerial <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$aerial <- scale(dat)
    } else {
      plt.tmp$aerial <- -5
    }
    }

    gene.ro <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC<0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ro,plt.tmp$cell]))
    if (length(gene.ro)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))
        dat[dat=="NaN"] <- min(dat[dat!="NaN"])
      } else {
        dat[dat==0] <- -5
      }
      plt.tmp$root <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$root <- scale(dat)
    } else {
      plt.tmp$root <- -5
    }
    }
    plt.tmp$aerial <- plt.tmp$aerial-min(plt.tmp$aerial)
    plt.tmp$root <- plt.tmp$root-min(plt.tmp$root)
    tmp <- rbind(tmp,plt.tmp)
  }
  return(tmp)
}

# 并行计算特殊聚类
clu2 <- c(14,26)
plt.ae.ro <- do.call('rbind',parallel::mclapply(clu2,function(i){score(i)},mc.cores = 8))
plt.ae.ro[plt.ae.ro$sample %in% c("S_A1","S_A2"),]$Sample <- "Aerial"
plt.ae.ro[plt.ae.ro$sample %in% c("S_R1","S_R2"),]$Sample <- "Root"

# 判断组织来源
plt.score.2 <- plt.ae.ro
plt.score.2$tissue <- "N"
plt.score.2[plt.score.2$root< plt.score.2$aerial,]$tissue <- "Aerial"
plt.score.2[plt.score.2$root> plt.score.2$aerial,]$tissue <- "Root"

# ===== 合并所有结果 =====
# 将主要聚类和特殊聚类的结果合并
plt.score <- rbind(plt.score.1,plt.score.2)

# 最终判断每个细胞的器官来源
plt.score$tissue <- "N"
plt.score[plt.score$root< plt.score$aerial,]$tissue <- "Aerial"
plt.score[plt.score$root> plt.score$aerial,]$tissue <- "Root"

# ===== 脚本完成 =====
# 输出：plt.score - 包含每个细胞的器官来源信息
# tissue列: "Aerial" = 叶, "Root" = 根, "N" = 无法判断