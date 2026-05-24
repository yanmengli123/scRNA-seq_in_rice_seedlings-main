##############################################################################
#                                                                            #
#   水稻幼苗单细胞转录组完整分析流程                                           #
#   基于 Wang et al. 2021 JGG 论文复现                                       #
#   数据来源: GSA CRA004082                                                   #
#   GitHub: https://github.com/Yuwang-art/scRNA-seq_in_rice_seedlings        #
#                                                                            #
##############################################################################

# ============================================================================
# 0. 安装和加载 R 包
# ============================================================================

# 首次运行请安装所有包
# install.packages(c("data.table", "tidyverse", "reshape2", "dplyr",
#                     "ggplot2", "patchwork", "ggpubr", "parallel", "umap"))
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("Seurat")  # 注意: 使用 Seurat v4
# install.packages("Seurat")      # 或直接从 CRAN 安装

library(data.table)
library(tidyverse)
library(Seurat)
library(reshape2)
library(Matrix)
library(dplyr)
library(umap)
library(parallel)
library(ggplot2)

# 设置全局参数
options(future.globals.maxSize = 80 * 1024^3)  # 80GB 内存上限

# ============================================================================
# 1. 读取 Cell Ranger 输出数据
# ============================================================================

# 假设你已经用 Cell Ranger 处理了每个样本
# 每个样本的输出在: cellranger_output/<sample>/outs/filtered_feature_bc_matrix/

# --- 修改这里的路径 ---
data_dir <- "./cellranger_output"  # Cell Ranger 输出的根目录

# 先读取一个对照样本测试
# 如果你只有部分样本，按需修改这个列表
sample_names <- c("Leaf1", "Leaf2", "Root1", "Root2",
                   "Seedling1", "Seedling2", "Seedling3", "Seedling4")

# 读取每个样本的 10x 矩阵
data_list <- list()
for (sample in sample_names) {
  matrix_path <- file.path(data_dir, sample, "outs", "filtered_feature_bc_matrix")
  if (file.exists(matrix_path)) {
    data_list[[sample]] <- Read10X(data.dir = matrix_path)
    cat("Loaded:", sample, "- Genes:", nrow(data_list[[sample]]),
        "Cells:", ncol(data_list[[sample]]), "\n")
  } else {
    cat("WARNING: Path not found for", sample, ":", matrix_path, "\n")
  }
}

# ============================================================================
# 2. 创建 Seurat 对象 + 质控
# ============================================================================

# 为每个样本创建 Seurat 对象
seurat_list <- list()
for (sample in names(data_list)) {
  seu <- CreateSeuratObject(counts = data_list[[sample]],
                            project = sample,
                            min.cells = 3,        # 基因至少在3个细胞中表达
                            min.features = 200)    # 细胞至少检测到200个基因
  seu$sample <- sample
  # 添加线粒体基因比例（水稻线粒体基因以 "Os" 开头的不算，需要用 MT 前缀）
  # 水稻中通常用叶绿体/线粒体基因比例来评估细胞质量
  # 水稻线粒体基因在 RGAP7 注释中以 LOC_Os12g... 或特定编号开头
  # 这里先不做线粒体过滤，原论文也没有这步
  seurat_list[[sample]] <- seu
}

# 合并所有样本（如果只需要简单合并）
# mydata <- merge(seurat_list[[1]], y = seurat_list[-1],
#                 add.cell.ids = names(seurat_list))

# 或者分别处理每个样本（推荐，后续做整合）
cat("Total samples loaded:", length(seurat_list), "\n")

# ============================================================================
# 3. SCTransform 标准化（每个样本单独做）
# ============================================================================

for (i in 1:length(seurat_list)) {
  cat("Running SCTransform on:", names(seurat_list)[i], "\n")
  seurat_list[[i]] <- SCTransform(seurat_list[[i]], verbose = TRUE)
}

# ============================================================================
# 4. 多样本整合 (Integration)
# ============================================================================

# 选择整合特征
rice.features <- SelectIntegrationFeatures(object.list = seurat_list,
                                            nfeatures = 3000)

# 准备 SCT 整合
seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                   anchor.features = rice.features,
                                   verbose = TRUE)

# 以第一个样本（如 Seedling-1）作为参考
reference_dataset <- 1

# 找锚点
rice.anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                        normalization.method = "SCT",
                                        anchor.features = rice.features,
                                        reference = reference_dataset,
                                        verbose = TRUE)

# 整合数据
rice.integrated <- IntegrateData(anchorset = rice.anchors,
                                  normalization.method = "SCT",
                                  verbose = TRUE)

cat("Integrated object:", nrow(rice.integrated), "genes x",
    ncol(rice.integrated), "cells\n")

# ============================================================================
# 5. PCA 降维
# ============================================================================

rice.integrated <- RunPCA(object = rice.integrated, verbose = FALSE, npcs = 100)

# 可视化 PCA（查看前几个 PC）
# DimPlot(rice.integrated, reduction = "pca")
# ElbowPlot(rice.integrated, ndims = 50)

# ============================================================================
# 6. UMAP 降维
# ============================================================================

# 方法1: 使用 Seurat 内置 UMAP
rice.integrated <- RunUMAP(object = rice.integrated,
                            dims = 1:10,
                            min.dist = 0.05,
                            n.neighbors = 5,
                            seed.use = 100)

# 方法2: 使用 umap 包（与原论文一致，可选）
# data <- as.data.frame(rice.integrated@reductions[["pca"]]@cell.embeddings)
# umap_result <- umap::umap(data[, 1:100],
#                            n_neighbors = 10,
#                            metric = "pearson2",
#                            min_dist = 0.01,
#                            random_state = 39)
# plt <- as.data.frame(umap_result$layout)
# colnames(plt) <- c("UMAP_1", "UMAP_2")
# plt$cell <- rownames(plt)
# 将结果放回 Seurat 对象...

# ============================================================================
# 7. 聚类
# ============================================================================

DefaultAssay(rice.integrated) <- "integrated"
rice.integrated <- FindNeighbors(rice.integrated, reduction = "pca", dims = 1:100)
rice.integrated <- FindClusters(rice.integrated, resolution = 0.75, n.start = 10)

# 查看聚类结果
cat("Number of clusters:", length(unique(rice.integrated$seurat_clusters)), "\n")
table(rice.integrated$seurat_clusters)

# ============================================================================
# 8. 可视化 UMAP 聚类图
# ============================================================================

# 8a. 按聚类着色
p1 <- DimPlot(rice.integrated, reduction = "umap",
              group.by = "seurat_clusters",
              label = TRUE, pt.size = 0.1) +
  ggtitle("Cell Clusters") +
  theme(legend.text = element_text(size = 8))

# 8b. 按样本着色
p2 <- DimPlot(rice.integrated, reduction = "umap",
              group.by = "sample",
              pt.size = 0.1) +
  ggtitle("Samples")

# 保存图片
ggsave("UMAP_clusters.pdf", p1, width = 10, height = 8)
ggsave("UMAP_samples.pdf", p2, width = 12, height = 8)

# 组合图
library(patchwork)
combined <- p1 + p2
ggsave("UMAP_combined.pdf", combined, width = 20, height = 8)

# ============================================================================
# 9. Log-Normalization（用于后续差异分析）
# ============================================================================

# 原论文在聚类后还做了一次 log-normalization
# 这是因为 SCTransform 的数据用于整合和聚类
# 而 log-normalized 数据用于差异表达分析

mydata_allgene <- CreateSeuratObject(counts = rice.integrated@assays$RNA@counts,
                                      project = "rice_allgene")
DefaultAssay(mydata_allgene) <- "RNA"
mydata_allgene <- NormalizeData(mydata_allgene,
                                 normalization.method = "LogNormalize",
                                 scale.factor = 1e6)
mydata_allgene <- FindVariableFeatures(mydata_allgene,
                                        selection.method = "vst",
                                        nfeatures = 2000)
all.genes <- rownames(mydata_allgene)
mydata_allgene <- ScaleData(mydata_allgene, features = all.genes)

# 复制聚类和 UMAP 信息
mydata_allgene$sample <- rice.integrated$sample
mydata_allgene$seurat_clusters <- rice.integrated$seurat_clusters
mydata_allgene@reductions$umap <- rice.integrated@reductions$umap
mydata_allgene@reductions$pca <- rice.integrated@reductions$pca

# ============================================================================
# 10. 鉴定 Cluster-Specific 基因（marker 基因）
# ============================================================================

DefaultAssay(mydata_allgene) <- "RNA"
cluster_markers <- FindAllMarkers(mydata_allgene,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25)

# 保存 marker 基因结果
write.csv(cluster_markers, "cluster_specific_markers.csv", row.names = FALSE)

# 每个 cluster 的 top5 marker
top5 <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# 可视化 top5 marker 的热图
pdf("Top5_markers_heatmap.pdf", width = 16, height = 12)
DoHeatmap(mydata_allgene, features = top5$gene) + NoLegend()
dev.off()

# ============================================================================
# 11. 细胞类型注释 (MICI 方法 - 简化版)
# ============================================================================

# 原论文使用 Table S5 中的 marker 基因做 MICI 注释
# 这里提供一个简化版：直接用已知的水稻细胞类型 marker 基因

# 水稻已知的 marker 基因（来自论文 Table S6）
# 叶片细胞类型 markers
leaf_markers <- list(
  mesophyll = c("LOC_Os02g05830", "LOC_Os12g17600",  # RBCS1, RBCS2
                "LOC_Os01g41710"),                       # CAB2R
  epidermis = c("LOC_Os04g48530"),                      # SLAC1
  procambium = c("LOC_Os02g08100"),                     # 4CL3
  vascular_initial = c("LOC_Os04g55590"),               # WOX4
  fiber = c("LOC_Os10g06000")
)

# 根部细胞类型 markers
root_markers <- list(
  cortex = c("LOC_Os02g41904"),
  xylem = c("LOC_Os01g73980"),
  pericycle = c("LOC_Os08g37300"),
  root_hair = c("LOC_Os07g35860"),
  endodermis = c("LOC_Os08g03450", "LOC_Os01g16890")
)

# 可视化 marker 基因表达
all_markers <- unlist(c(leaf_markers, root_markers))
# 过滤出在数据中存在的基因
available_markers <- all_markers[all_markers %in% rownames(mydata_allgene)]

if (length(available_markers) > 0) {
  pdf("Cell_type_markers_featureplot.pdf", width = 16, height = 20)
  print(FeaturePlot(mydata_allgene, features = available_markers,
                    ncol = 3, pt.size = 0.1))
  dev.off()
} else {
  cat("NOTE: Marker genes use MSU LOC format (LOC_Os...).\n")
  cat("If your data uses different gene IDs, you need to convert them.\n")
  cat("Check rownames(mydata_allgene) to see the format.\n")
}

# ============================================================================
# 12. 按样本统计细胞比例
# ============================================================================

# 统计每个样本中各 cluster 的细胞数
cluster_sample_table <- table(mydata_allgene$seurat_clusters,
                                mydata_allgene$sample)
print(cluster_sample_table)

# 保存
write.csv(as.data.frame(cluster_sample_table),
          "cluster_sample_distribution.csv", row.names = FALSE)

# ============================================================================
# 13. 保存 Seurat 对象
# ============================================================================

saveRDS(rice.integrated, "rice_integrated_seurat.rds")
saveRDS(mydata_allgene, "rice_allgene_seurat.rds")
cat("Seurat objects saved!\n")

# ============================================================================
# 14. 差异表达分析（以胁迫响应为例）
# ============================================================================

# 如果你有 LN/HS/ID 样本，可以做胁迫响应分析
# 以低氮(LN) vs 对照(Ctrl) 为例

# 假设对照样本为 Seedling-1~4，低氮样本为 LN-1~4
# ctrl_cells <- WhichCells(mydata_allgene, expression = sample %in%
#                            c("Seedling1", "Seedling2", "Seedling3", "Seedling4"))
# ln_cells <- WhichCells(mydata_allgene, expression = sample %in%
#                           c("LN1", "LN2", "LN3", "LN4"))

# 对每个 cluster 做差异分析
# for (clu in levels(Idents(mydata_allgene))) {
#   markers <- FindMarkers(mydata_allgene,
#                          ident.1 = "LN", ident.2 = "Ctrl",
#                          group.by = "condition",
#                          subset.ident = clu,
#                          logfc.threshold = 0, min.pct = 0)
#   write.csv(markers, paste0("DEG_cluster", clu, "_LN_vs_Ctrl.csv"))
# }

# ============================================================================
# 15. 拟时序分析（Monocle2，可选）
# ============================================================================

# 如果想做发育轨迹分析，需要安装 monocle2
# BiocManager::install("monocle")
# library(monocle)

# 基本流程（参考原论文脚本4）：
# 1. 选择特定 cluster 的细胞（如 mesophyll 相关 clusters: 0,1,17-22）
# 2. 从 Seurat 对象提取数据
# 3. 创建 CellDataSet
# 4. 降维 (DDRTree)
# 5. 计算 pseudotime
# 6. 可视化轨迹

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("Output files:\n")
cat("  - UMAP_clusters.pdf\n")
cat("  - UMAP_samples.pdf\n")
cat("  - cluster_specific_markers.csv\n")
cat("  - Top5_markers_heatmap.pdf\n")
cat("  - Cell_type_markers_featureplot.pdf\n")
cat("  - rice_integrated_seurat.rds\n")
cat("  - rice_allgene_seurat.rds\n")
cat("========================================\n")
