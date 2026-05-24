# ===== 加载R包 =====
library(data.table)
library(tidyverse)
library(Seurat)
library(reshape2)
library(Matrix)
library(dplyr)
library(parallel)
library(patchwork)

#########################################################
#                                                       #
#                          MICI                         #
#     Marker-based Identification of Cell type Identity #
#            基于标记基因的细胞类型鉴定方法                 #
#                                                       #
#########################################################

# ===== 第1部分：叶组织(Lef)细胞类型鉴定 =====

# 读取标记基因表 (Table S5)
marker <- as.data.frame(fread('Table S5',header = T))
marker$`Gene ID` <- gsub("_","-",marker$`Gene ID`)  # 统一基因ID格式
rownames(marker) <- 1:dim(marker)[1]
colnames(marker)[1] <- "gene"

# 筛选叶组织相关标记基因 (Location为Both或Leaf)
marker1 <-  marker[marker$Location %in% c("Both","Leaf"),]
dim(marker1)

# 将Leaf列中的多个组织拆分为多行
# 例如：一个基因标记多个组织时，会变成多行
marker2 <- marker1 %>% separate_rows(Leaf,sep = ",")
dim(marker2)

# 选择叶组织样本
sam <- c("S_A1","S_A2")

# 提取表达数据：只保留标记基因在叶样本中的表达
expr_data <- as.data.frame(mydata_allgene@assays$RNA@data[
  intersect(rownames(mydata_allgene@assays$RNA@data),unique(marker1$gene)),plt[plt$sample %in% sam,]$cell])
# ===== 第2步：计算每个标记基因的权重 =====
# 权重 = 基因在各聚类中表达的方差
# 方差越大，说明该基因在不同细胞类型间差异越大，区分能力越强

weight_data <- data.frame(stringsAsFactors = F)
for (i in 0:28){
  for (j in 1:length(marker1$gene)) {
    # 计算基因在每个聚类中的平均表达量 (log2转换)
    weight_data[j,i+1] <- log2(rowMeans(expm1(expr_data[marker1$gene[j],plt[plt$cluster==i&plt$sample %in% sam,]$cell]))+1)
  }
}
rownames(weight_data) <- marker1$gene
# 计算方差作为权重
weight_data$weight <- apply(weight_data,1,var)
weight_data$gene <- rownames(weight_data)
weight_data <- merge(weight_data,marker1[,c("gene","Tissue","Leaf")],by="gene")

# ===== 第3步：计算每个细胞的MICI分数 =====
# MICI分数 = 加权后的标记基因表达之和
# 用于判断每个细胞最可能属于哪种细胞类型

# 对表达数据进行Z-score标准化 (行-wise，即按基因)
expr_scale_data <- as.data.frame(t(scale(t(expr_data))))

# 合并标记基因信息
expr_scale_data <- merge(marker2,expr_scale_data,by.x='gene',by.y="row.names")
expr_scale_data <- merge(weight_data[,c("gene","weight")],expr_scale_data,by="gene")

# 用权重加权表达值
expr_scale_data[,plt[plt$sample %in% sam,]$cell] <- expr_scale_data[,plt[plt$sample %in% sam,]$cell]*expr_scale_data$weight

# 定义函数：计算每个细胞类型的MICI分数
fun <- function(x) {
  tmp <- expr_scale_data[expr_scale_data$Leaf==unique(expr_scale_data$Leaf)[x],plt[plt$sample %in% sam,]$cell]
  tmp1 <- as.data.frame(t(as.data.frame(colSums(tmp))))
  rownames(tmp1) <- unique(dat$Leaf)[x]
  return(tmp1)
}

# 并行计算所有细胞类型的MICI分数
MICI_out <- do.call('rbind',parallel::mclapply(1:length(unique(expr_scale_data$Leaf)),
                                               function(x){fun(x)},mc.cores = length(unique(expr_scale_data$Leaf))))

# 为每个细胞分配细胞类型 (选择MICI分数最高的类型)
MICI_out$cell_type <- as.character(rownames(MICI_out))
MICI_result <-apply(MICI_out[,-dim(MICI_out)[2]], 2, function(x){MICI_out[,dim(MICI_out)[2]][which(x==max(x))]}) %>% as.data.frame()
colnames(MICI_result) <- c("identified_ct")

# 合并聚类信息
MICI_result <- merge(plt,MICI_result,by.x="cell",by.y="row.names")
MICI_result$identified_ct <- as.character(MICI_result$identified_ct)

# 将MICI分数 < 2 的细胞标记为"Unknown" (低置信度)
tmp <- as.data.frame(t(MICI_out[,-dim(MICI_out)[2]]))
unct <- rownames(tmp[rowMax(as.matrix(tmp))<2,])
MICI_result[MICI_result$cell %in% unct,]$identified_ct <- "Unknown"

# ===== 第4步：将聚类分配给细胞类型 (叶组织) =====
# 根据MICI结果，将每个聚类分配到最主要的细胞类型

# 统计每个聚类中各细胞类型的细胞数量
dat1 <- as.data.frame(table(MICI_result$identified_ct,MICI_result$cluster))
colnames(dat1) <- c("cell_type","cluster","number_of_cell")

# 统计每个聚类的总细胞数
dat2 <- as.data.frame(table(MICI_result$cluster))
colnames(dat2) <- c("cluster","number_of_cluster")

# 合并并计算比例
dat3 <- merge(dat1,dat2,by="cluster")
dat3$ratio <- signif(dat3$number_of_cell/dat3$number_of_cluster,3)

# 选择每个聚类中比例最高的细胞类型
dat4 <- dat3 %>% group_by(cluster) %>% top_n(1,ratio) %>% as.data.frame()

# 创建聚类-细胞类型映射矩阵
dat5 <- reshape2::dcast(dat3,cluster~cell_type,value.var = "ratio")

# 最终结果：聚类到细胞类型的分配
cluster_assign <- merge(dat4[,c("cluster","cell_type")],dat5,by="cluster")

# ===== 第2部分：根组织(Root)细胞类型鉴定 =====
# 与叶组织分析流程相同，但使用根组织标记基因

# 重新读取标记基因表
marker$`Gene ID` <- gsub("_","-",marker$`Gene ID`)
rownames(marker) <- 1:dim(marker)[1]
colnames(marker)[1] <- "gene"

# 筛选根组织相关标记基因 (Location为Both或Root)
marker1 <-  marker[marker$Location %in% c("Both","Root"),]
dim(marker1)

# 拆分多组织标记基因
marker2 <- marker1 %>% separate_rows(Leaf,sep = ",")
dim(marker2)

# 选择根组织样本
sam <- c("S_R1","S_R2")

# 指定根组织相关的聚类 (排除23, 24聚类，这些是叶组织特有的)
clu <- c(0:22,25,26:28)

# 提取表达数据
# ===== 根组织：计算标记基因权重 =====
weight_data <- data.frame(stringsAsFactors = F)
for (i in 1:length(clu)){
  for (j in 1:length(marker1$gene)) {
    # 计算基因在每个聚类中的平均表达量
    tmp[j,i] <- log2(rowMeans(expm1(dat[marker1$gene[j],plt[plt$cluster==clu[i]&plt$sample %in% sam,]$cell]))+1)
  }
}
rownames(weight_data) <- marker1$gene
weight_data$weight <- apply(weight_data,1,var)  # 方差作为权重
weight_data$gene <- rownames(weight_data)
weight_data <- merge(weight_data,marker1[,c("gene","Tissue","Leaf")],by="gene")

# ===== 根组织：计算MICI分数 =====
# Z-score标准化
expr_scale_data <- as.data.frame(t(scale(t(expr_data))))

# 合并标记基因信息
expr_scale_data <- merge(marker2,expr_scale_data,by.x='gene',by.y="row.names")
expr_scale_data <- merge(weight_data[,c("gene","weight")],expr_scale_data,by="gene")

# 用权重加权表达值
expr_scale_data[,plt[plt$sample %in% sam,]$cell] <- expr_scale_data[,plt[plt$sample %in% sam,]$cell]*expr_scale_data$weight

# 计算每个细胞类型的MICI分数
fun <- function(x) {
  tmp <- expr_scale_data[expr_scale_data$Leaf==unique(expr_scale_data$Leaf)[x],plt[plt$sample %in% sam,]$cell]
  tmp1 <- as.data.frame(t(as.data.frame(colSums(tmp))))
  rownames(tmp1) <- unique(dat$Leaf)[x]
  return(tmp1)
}

# 并行计算
MICI_out <- do.call('rbind',parallel::mclapply(1:length(unique(expr_scale_data$Leaf)),
                                               function(x){fun(x)},mc.cores = length(unique(expr_scale_data$Leaf))))

# 为每个细胞分配细胞类型
MICI_out$cell_type <- as.character(rownames(MICI_out))
MICI_result <-apply(MICI_out[,-dim(MICI_out)[2]], 2, function(x){MICI_out[,dim(MICI_out)[2]][which(x==max(x))]}) %>% as.data.frame()
colnames(MICI_result) <- c("identified_ct")
MICI_result <- merge(plt,MICI_result,by.x="cell",by.y="row.names")
MICI_result$identified_ct <- as.character(MICI_result$identified_ct)

# 低置信度细胞标记为Unknown
tmp <- as.data.frame(t(MICI_out[,-dim(MICI_out)[2]]))
unct <- rownames(tmp[rowMax(as.matrix(tmp))<2,])
MICI_result[MICI_result$cell %in% unct,]$identified_ct <- "Unknown"

# ===== 根组织：将聚类分配给细胞类型 =====
dat1 <- as.data.frame(table(MICI_result$identified_ct,MICI_result$cluster))
colnames(dat1) <- c("cell_type","cluster","number_of_cell")
dat2 <- as.data.frame(table(MICI_result$cluster))
colnames(dat2) <- c("cluster","number_of_cluster")
dat3 <- merge(dat1,dat2,by="cluster")
dat3$ratio <- signif(dat3$number_of_cell/dat3$number_of_cluster,3)
dat4 <- dat3 %>% group_by(cluster) %>% top_n(1,ratio) %>% as.data.frame()
dat5 <- reshape2::dcast(dat3,cluster~cell_type,value.var = "ratio")
cluster_assign <- merge(dat4[,c("cluster","cell_type")],dat5,by="cluster")

# ===== 脚本完成 =====
# 输出：cluster_assign - 聚类到细胞类型的映射表