# ===== 加载R包 =====
library(data.table)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(monocle)      # 发育轨迹推断核心包
library(parallel)
library(ggpubr)

#################################################################################
#                                                                               #
#          Reconstruction of developmental trajectory by monocle                #
#                    使用Monocle重建发育轨迹                                    #
#                                                                               #
#################################################################################

# ===== 数据准备 =====
# 目的：从叶组织细胞中重建发育轨迹
# 选择聚类0, 1, 17-22的叶组织细胞，排除叶和根的原生质体制备样本

# 提取整合后的标准化数据
data_input <- as.matrix(rice.integrated@assays[["integrated"]]@scale.data)

# 筛选目标细胞
# 条件：特定聚类 + 叶组织 + 排除原生质体制备样本
cell_input <- plt.score[plt.score$cluster %in% c(0,1,17:22)&plt.score$tissue=="Aerial"&
                plt.score$Sample!="Aerial"&plt.score$Sample!="Root",]

# 随机抽样10,000个细胞 (Monocle对大数据集较慢)
cell_input$random <- sample(1:dim(cell_input)[1],dim(cell_input)[1])
cell_input <- cell_input[order(cell_input$random),]
data_input <- data_input[,colnames(data_input) %in% cell_input$cell[1:10000]]

# ===== 创建Monocle对象 =====
# Monocle需要三个组件：
# 1. 表达矩阵 (genes x cells)
# 2. 基因注释 (gene_annotation)
# 3. 细胞信息 (sample_sheet)

# 基因注释
gene_annotation <- data.frame(row.names = rownames(data_input),gene_short_name =rownames(data_input))
gene_annotation$gene_short_name <- gsub("LOC-","",gene_annotation$gene_short_name)

# 转换为稀疏矩阵格式
raw_select <- as.matrix(data_input)

# 细胞信息 (这里只有一个样本组)
sample_sheet <- data.frame(row.names = colnames(raw_select),Library= rep("mes",dim(raw_select)[2]))

# 创建AnnotatedDataFrame对象
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = sample_sheet)

# 创建CellDataSet对象
mono <- newCellDataSet(as(raw_select, "sparseMatrix"),phenoData = pd,featureData = fd)

# 估计大小因子 (用于标准化)
mono <- estimateSizeFactors(mono)

# 检测表达基因 (最低表达阈值0.1)
mono <- detectGenes(mono, min_expr = 0.1)

# ===== 降维分析 =====
# 设置用于轨迹推断的基因
ordering_genes <- rownames(df)
mono <- setOrderingFilter(mono, ordering_genes)

# 使用DDRTree进行降维 (Discriminative Dimensionality Reduction via Learning a Tree)
# 这是Monocle 2的核心算法，可以学习细胞的树状发育结构
mono <- reduceDimension(mono, reduction_method = "DDRTree",norm_method="none",pseudo_expr=0)

# ===== 计算拟时序 =====
# 拟时序(pseudotime)代表细胞在发育过程中的位置
# 值越小越接近发育起点，越大越接近终点

mono <- orderCells(mono)

# 可视化轨迹
plot_cell_trajectory(mono, color_by="Pseudotime")  # 按拟时序着色
plot_cell_trajectory(mono, color_by="State")        # 按状态着色

# 设置根状态 (发育起点)
mono <- orderCells(mono,root_state = 2)
plot_cell_trajectory(mono, color_by="Pseudotime")   # 重新可视化

# ===== 脚本完成 =====
# 输出：mono - 包含发育轨迹信息的Monocle对象
# 可用于分析基因表达随发育时间的变化
