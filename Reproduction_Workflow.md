# 水稻幼苗单细胞RNA测序论文复现完整流程

> **论文**: Single-cell transcriptome atlas of the leaf and root of rice seedlings
> **数据编号**: GSA CRA004082 / BioProject PRJCA004855
> **物种**: 水稻 (*Oryza sativa*)

---

## 目录

1. [数据概述](#1-数据概述)
2. [计算环境要求](#2-计算环境要求)
3. [第一阶段：原始数据下载](#3-第一阶段原始数据下载)
4. [第二阶段：Cell Ranger 处理](#4-第二阶段cell-ranger-处理)
5. [第三阶段：Seurat 单细胞分析](#5-第三阶段seurat-单细胞分析)
6. [第四阶段：细胞类型鉴定 (MICI)](#6-第四阶段细胞类型鉴定-mici)
7. [第五阶段：器官来源鉴定](#7-第五阶段器官来源鉴定)
8. [第六阶段：发育轨迹重建](#8-第六阶段发育轨迹重建)
9. [第七阶段：高级分析与可视化](#9-第七阶段高级分析与可视化)
10. [常见问题与解决方案](#10-常见问题与解决方案)

---

## 1. 数据概述

### 1.1 样本清单

| 样本名 | 样本类型 | 重复 | 处理条件 | 数据类型 |
|--------|---------|------|---------|---------|
| Leaf-1, Leaf-2 | 叶组织 | 2 | 对照 | 10x scRNA-seq |
| Root-1, Root-2 | 根组织 | 2 | 对照 | 10x scRNA-seq |
| Seedling-1 ~ Seedling-4 | 幼苗 | 4 | 对照 | 10x scRNA-seq |
| HS-1 ~ HS-4 | 幼苗 | 4 | 热胁迫 (Heat Stress) | 10x scRNA-seq |
| ID-1 ~ ID-4 | 幼苗 | 4 | 缺铁 (Iron Deficiency) | 10x scRNA-seq |
| LN-1 ~ LN-4 | 幼苗 | 4 | 低氮 (Low Nitrogen) | 10x scRNA-seq |
| QC replicate 1/2 | 质控 | 2 | 水稻+拟南芥混合 | 10x scRNA-seq |
| Pre-protoplasting 1/2 | 原生质体前 | 2 | bulk | FASTQ |
| Post-protoplasting 1/2 | 原生质体后 | 2 | bulk | FASTQ |
| Mesophyll cell bulk | 叶肉细胞 | 2 | Ctrl/HS | FASTQ |

### 1.2 数据存储结构

```
CRA004082/
├── CRA004082.csv          # 元数据（下载链接）
├── CRA004082.txt          # 实验accession列表
├── 1-s2.0-...-mmc2.xlsx  # Table S1-S14（标记基因等）
└── 1-s2.0-...-mmc1.pdf   # 补充材料说明
```

### 1.3 关键编号对照

| 编号类型 | 编号 | 说明 |
|---------|------|------|
| GSA项目 | CRA004082 | 数据提交号 |
| BioProject | PRJCA004855 | NCBI项目号 |
| 样本ID | SAMC350942 ~ SAMC350969 | BioSample编号 |
| 实验ID | CRX235962 ~ CRX235989 | 实验编号 |
| 运行ID | CRR279079 ~ CRR279180 | 测序运行编号 |

---

## 2. 计算环境要求

### 2.1 硬件要求

| 阶段 | CPU | 内存 | 磁盘 | 时间估算 |
|------|-----|------|------|---------|
| Cell Ranger | 16核+ | 64GB+ | 2TB+ | 3-7天 |
| Seurat分析 | 8核+ | 32GB+ | 500GB | 1-2天 |
| MICI鉴定 | 8核 | 16GB+ | 100GB | 数小时 |
| 轨迹分析 | 4核 | 16GB | 50GB | 数小时 |

### 2.2 软件依赖

```bash
# Cell Ranger (需要10x Genomics官方下载)
# 版本: Cell Ranger 7.0+
# 下载: https://www.10xgenomics.com/support/software/cell-ranger/downloads

# R 包
install.packages(c("Seurat", "monocle", "ggplot2", "dplyr", "tidyverse",
                   "data.table", "reshape2", "Matrix", "patchwork", "ggpubr",
                   "umap", "parallel"))

# 参考基因组
# 水稻参考基因组: IRGSP-1.0 (RAP-DB)
# 下载: https://rapdb.dna.affrc.go.jp/
```

### 2.3 方案选择

| 方案 | 适用条件 | 优点 | 缺点 |
|------|---------|------|------|
| **A: 从头处理** | 有高性能服务器 | 完全掌控，学习全流程 | 耗时长，资源需求高 |
| **B: 下载已处理数据** | 个人电脑 | 快速，资源需求低 | 依赖他人处理质量 |
| **C: 混合方案** | 中等配置 | 平衡学习与效率 | 需要一定经验 |

---

## 3. 第一阶段：原始数据下载

### 3.1 创建目录结构

```bash
mkdir -p rice_scRNAseq/{raw_data,cellranger_output,seurat_analysis,results,figures}
cd rice_scRNAseq/raw_data
```

### 3.2 下载单细胞数据（10x Genomics）

**方案A：下载全部样本（推荐用于完整复现）**

```bash
# 从CRA004082.csv解析下载链接
# 每个样本4个lane，共20个样本 = 80个tar文件

# 示例：下载 Leaf-1 (4个lane)
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279087/CRR279087.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279088/CRR279088.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279089/CRR279089.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279090/CRR279090.tar

# 示例：下载 Root-1 (4个lane)
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279095/CRR279095.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279096/CRR279096.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279097/CRR279097.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279098/CRR279098.tar

# ... 以此类推下载所有样本
```

**方案B：批量下载脚本**

```bash
#!/bin/bash
# download_all.sh - 从CSV批量下载

CSV_FILE="CRA004082.csv"
mkdir -p ../raw_data

# 提取10x数据的下载链接
grep "10x genomics" $CSV_FILE | cut -d',' -f7 | while read url; do
    filename=$(basename $url)
    echo "Downloading $filename..."
    wget -P ../raw_data/ $url
done

# 提取fastq数据的下载链接
grep "fastq" $CSV_FILE | cut -d',' -f7 | tr '|' '\n' | while read url; do
    filename=$(basename $url)
    echo "Downloading $filename..."
    wget -P ../raw_data/ $url
done
```

**方案C：下载已处理数据（推荐用于快速学习）**

```bash
# 从GEO下载 Liu et al. 2021 的水稻根部数据
# GSE146035 包含处理好的 Seurat 对象
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146035/suppl/

# 或者从 figshare 下载本论文作者可能提供的处理后数据
```

### 3.3 解压数据

```bash
cd ../raw_data

# 解压所有tar文件
for f in CRR*.tar; do
    echo "Extracting $f..."
    tar xf $f
done

# 查看解压后的文件结构
ls -la
# 每个tar解压后会得到：
# - CRR279xxx_1.fastq.gz (Read 1，包含barcode和UMI)
# - CRR279xxx_2.fastq.gz (Read 2，包含插入序列)
# - CRR279xxx_3.fastq.gz (Index read，可选)
```

### 3.4 下载Bulk RNA-seq数据

```bash
# 下载原生质体制备前后的bulk数据
# Pre-protoplasting replicate 1
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279175/CRR279175_f1.fastq.gz
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279175/CRR279175_r2.fastq.gz

# Post-protoplasting replicate 1
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279177/CRR279177_f1.fastq.gz
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279177/CRR279177_r2.fastq.gz

# Mesophyll cell bulk RNA-seq (Ctrl)
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279179/CRR279179_f1.fastq.gz
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279179/CRR279179_r2.fastq.gz
```

---

## 4. 第二阶段：Cell Ranger 处理

### 4.1 下载水稻参考基因组

```bash
cd ../cellranger_output

# 从RAP-DB下载水稻参考基因组 (IRGSP-1.0)
wget https://rapdb.dna.affrc.go.jp/data/irgsp1/IRGSP-1.0_genome.fasta.gz
gunzip IRGSP-1.0_genome.fasta.gz

# 下载基因注释文件
wget https://rapdb.dna.affrc.go.jp/data/irgsp1/IRGSP-1.0_genes.gff3.gz
gunzip IRGSP-1.0_genes.gff3.gz
```

### 4.2 构建Cell Ranger参考

```bash
# 构建参考基因组索引
cellranger mkref \
    --genome=rice_genome \
    --fasta=IRGSP-1.0_genome.fasta \
    --genes=IRGSP-1.0_genes.gff3 \
    --nthreads=16

# 这会生成 rice_genome/ 目录，包含：
# - fasta/ (基因组序列)
# - genes/ (基因注释)
# - star/ (STAR索引)
```

### 4.3 合并同一样本的多个Lane

```bash
# 创建样本目录
mkdir -p fastqs/{Leaf1,Leaf2,Root1,Root2,Seedling1,Seedling2,Seedling3,Seedling4}
mkdir -p fastqs/{HS1,HS2,HS3,HS4,ID1,ID2,ID3,ID4,LN1,LN2,LN3,LN4}

# 将同一样本的4个lane合并到同一目录
# Leaf-1 示例
mv raw_data/CRR279087*.fastq.gz fastqs/Leaf1/
mv raw_data/CRR279088*.fastq.gz fastqs/Leaf1/
mv raw_data/CRR279089*.fastq.gz fastqs/Leaf1/
mv raw_data/CRR279090*.fastq.gz fastqs/Leaf1/

# ... 以此类推
```

### 4.4 运行Cell Ranger

```bash
# 对每个样本运行Cell Ranger count
# 示例：处理 Leaf-1

cellranger count \
    --id=Leaf1_result \
    --transcriptome=./rice_genome \
    --fastqs=./fastqs/Leaf1 \
    --sample=Leaf1 \
    --expect-cells=10000 \
    --nthreads=16 \
    --memgb=64

# 重要参数说明：
# --id: 输出目录名
# --transcriptome: 参考基因组路径
# --fastqs: FASTQ文件目录
# --sample: 样本名前缀
# --expect-cells: 预期细胞数 (10x通常10,000-100,000)
# --nthreads: 线程数
# --memgb: 内存限制
```

### 4.5 批量处理脚本

```bash
#!/bin/bash
# run_cellranger_all.sh

samples=("Leaf1" "Leaf2" "Root1" "Root2" "Seedling1" "Seedling2" "Seedling3" "Seedling4")
samples+=("HS1" "HS2" "HS3" "HS4" "ID1" "ID2" "ID3" "ID4" "LN1" "LN2" "LN3" "LN4")

for sample in "${samples[@]}"; do
    echo "Processing $sample..."
    cellranger count \
        --id=${sample}_result \
        --transcriptome=./rice_genome \
        --fastqs=./fastqs/$sample \
        --sample=$sample \
        --expect-cells=10000 \
        --nthreads=16 \
        --memgb=64
    echo "$sample completed."
done
```

### 4.6 Cell Ranger输出

每个样本会在 `outs/` 目录生成：

```
${sample}_result/
└── outs/
    ├── filtered_feature_bc_matrix/  # 过滤后的表达矩阵
    │   ├── barcodes.tsv.gz          # 细胞条形码
    │   ├── features.tsv.gz          # 基因信息
    │   └── matrix.mtx.gz            # 表达矩阵
    ├── raw_feature_bc_matrix/       # 未过滤的原始矩阵
    ├── metrics_summary.csv          # 质控指标
    ├── web_summary.html             # 质控报告 (在浏览器打开)
    └── possorted_genome_bam.bam     # 比对文件
```

### 4.7 质控指标检查

```bash
# 查看质控报告
cat ${sample}_result/outs/metrics_summary.csv

# 关键指标：
# - Estimated Number of Cells: 估计的细胞数
# - Mean Reads per Cell: 平均每个细胞的reads数
# - Median Genes per Cell: 中位基因数/细胞
# - Total Genes Detected: 检测到的总基因数
# - Fraction Reads in Cells: 细胞内reads比例 (>70%为好)
```

---

## 5. 第三阶段：Seurat 单细胞分析

### 5.1 合并所有样本的表达矩阵

```R
# 准备合并所有Cell Ranger输出
library(Seurat)
library(Matrix)

# 读取所有样本的10x数据
samples <- c("Leaf1", "Leaf2", "Root1", "Root2",
             "Seedling1", "Seedling2", "Seedling3", "Seedling4",
             "HS1", "HS2", "HS3", "HS4",
             "ID1", "ID2", "ID3", "ID4",
             "LN1", "LN2", "LN3", "LN4")

# 创建合并的表达矩阵
all_data <- list()
for (sample in samples) {
    path <- paste0("./cellranger_output/", sample, "_result/outs/filtered_feature_bc_matrix/")
    data <- Read10X(path)
    colnames(data) <- paste0(colnames(data), "-", sample)
    all_data[[sample]] <- data
}

# 合并为一个大矩阵
all_matrix <- do.call(cbind, all_data)

# 保存为CSV供后续分析
write.csv(as.matrix(all_matrix), "all.csv")
```

### 5.2 原生质体制备敏感基因鉴定

```R
library(data.table)
library(ggplot2)

# 读取bulk RNA-seq数据
genelength <- read.table("Gene.length", sep="\t", header=TRUE)[, c(1,6)]

# 处理Post-protoplasting数据
rep1_AP <- read.table("HTSeq_out_rep1_AP.txt")
rep1_AP <- rep1_AP[grep("AT", rep1_AP$V1),]
# ... (完整代码见 1-Single-cell pepline.R)

# 计算TPM
# TPM = (reads / gene_length_kb) / (total_rpk / 10^6)

# 鉴定敏感基因：|log2FC| > 3
cutoff <- 3
bulk_rice$group <- "N"
bulk_rice$group[(log2(bulk_rice$fc_rep1) > cutoff & log2(bulk_rice$fc_rep2) > cutoff) |
                (log2(bulk_rice$fc_rep1) < -cutoff & log2(bulk_rice$fc_rep2) < -cutoff)] <- "Y"

# 排除敏感基因
sensitive_genes <- bulk_rice[bulk_rice$group == "Y", "gene"]
all_integ <- all_integ[!rownames(all_integ) %in% sensitive_genes,]
```

### 5.3 Seurat对象创建与标准化

```R
library(Seurat)
library(future)

# 创建Seurat对象
sample <- separate(data = data.frame(cell = colnames(all_integ)),
                   col = "cell", into = c("barcode", "sample"), sep = "-")$sample
mydata <- CreateSeuratObject(counts = all_integ, project = "rice_scRNAseq")
mydata@meta.data$sample <- sample

# 设置并行计算
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 80 * 1024^3)

# SCTransform标准化
rice.list <- SplitObject(mydata, split.by = "sample")
for (i in 1:length(rice.list)) {
    rice.list[[i]] <- SCTransform(rice.list[[i]], verbose = TRUE)
}
```

### 5.4 样本整合

```R
# 选择整合特征
rice.features <- SelectIntegrationFeatures(object.list = rice.list, nfeatures = 3000)
rice.list <- PrepSCTIntegration(object.list = rice.list,
                                anchor.features = rice.features)

# 找到锚点
reference_dataset <- which(names(rice.list) == "Ctrl1")
rice.anchors <- FindIntegrationAnchors(object.list = rice.list,
                                       normalization.method = "SCT",
                                       anchor.features = rice.features,
                                       reference = reference_dataset)

# 整合数据
rice.integrated <- IntegrateData(anchorset = rice.anchors,
                                 normalization.method = "SCT")
```

### 5.5 降维与聚类

```R
# PCA降维
rice.integrated <- RunPCA(object = rice.integrated, verbose = FALSE, npcs = 100)

# UMAP降维
rice.integrated <- RunUMAP(object = rice.integrated,
                           dims = 1:10,
                           min.dist = 0.05,
                           n.neighbors = 5,
                           seed.use = 100)

# 构建KNN图并聚类
DefaultAssay(rice.integrated) <- "integrated"
rice.integrated <- FindNeighbors(rice.integrated, reduction = "pca", dims = 1:100)
rice.integrated <- FindClusters(rice.integrated, resolution = 0.75, n.start = 10)

# 可视化
DimPlot(rice.integrated, reduction = "umap", label = TRUE)
```

### 5.6 聚类重编号

```R
# 原始聚类ID到新编号的映射（根据论文Figure 1）
cluster_mapping <- data.frame(
    original = 0:28,
    new = c(24,13,3,23,7,21,17,1,10,9,8,5,6,14,27,20,19,22,4,12,2,0,15,26,11,18,25,16,28)
)

# 应用重编号
rice.integrated@meta.data$seurat_clusters_new <- cluster_mapping$new[
    match(rice.integrated@meta.data$seurat_clusters, cluster_mapping$original)
]
```

---

## 6. 第四阶段：细胞类型鉴定 (MICI)

### 6.1 MICI算法原理

**MICI = Marker-based Identification of Cell type Identity**

核心思想：
1. 使用已知的标记基因（来自 Table S5）
2. 计算每个标记基因的权重（基于聚类间方差）
3. 计算每个细胞的MICI分数（加权表达之和）
4. 将聚类分配给MICI分数最高的细胞类型

### 6.2 准备标记基因

```R
# 从mmc2.xlsx读取Table S5
library(readxl)
marker <- read_excel("1-s2.0-S1673852721001673-mmc2.xlsx", sheet = "Table S5")

# 分离叶和根的标记基因
marker_leaf <- marker[marker$Location %in% c("Both", "Leaf"),]
marker_root <- marker[marker$Location %in% c("Both", "Root"),]
```

### 6.3 计算MICI分数

```R
# 叶组织细胞类型鉴定
sam <- c("S_A1", "S_A2")  # 叶样本
expr_data <- as.data.frame(mydata_allgene@assays$RNA@data[
    intersect(rownames(mydata_allgene@assays$RNA@data), unique(marker_leaf$gene)),
    plt[plt$sample %in% sam,]$cell
])

# 计算每个标记基因的权重（方差）
weight_data <- data.frame()
for (i in 0:28) {
    for (j in 1:length(marker_leaf$gene)) {
        weight_data[j, i+1] <- log2(rowMeans(expm1(
            expr_data[marker_leaf$gene[j],
                     plt[plt$cluster == i & plt$sample %in% sam,]$cell
        ])) + 1)
    }
}
rownames(weight_data) <- marker_leaf$gene
weight_data$weight <- apply(weight_data, 1, var)

# Z-score标准化
expr_scale_data <- as.data.frame(t(scale(t(expr_data))))

# 用权重加权
expr_scale_data[, plt[plt$sample %in% sam,]$cell] <-
    expr_scale_data[, plt[plt$sample %in% sam,]$cell] * weight_data$weight

# 计算每个细胞类型的MICI分数
MICI_scores <- colSums(expr_scale_data[, plt[plt$sample %in% sam,]$cell])

# 为每个细胞分配细胞类型
MICI_result <- data.frame(cell = names(MICI_scores), score = MICI_scores)
MICI_result$cell_type <- marker_leaf$Tissue[match(MICI_result$gene_cluster, marker_leaf$gene)]
```

### 6.4 聚类分配

```R
# 统计每个聚类中各细胞类型的细胞数
cluster_celltype <- table(MICI_result$cell_type, MICI_result$cluster)

# 计算比例
cluster_celltype_ratio <- prop.table(cluster_celltype, margin = 2)

# 选择每个聚类中比例最高的细胞类型
cluster_assignment <- data.frame(
    cluster = colnames(cluster_celltype_ratio),
    cell_type = rownames(cluster_celltype_ratio)[apply(cluster_celltype_ratio, 2, which.max)]
)
```

---

## 7. 第五阶段：器官来源鉴定

### 7.1 鉴定器官偏向基因

```R
# 差异表达分析：叶 vs 根
# 使用Wilcoxon检验

# 批次1分析
Idents(mydata_allgene) <- "sample"
ae.ro1 <- subset(mydata_allgene, idents = c("S_A1", "S_R1"))

# 为每个聚类找差异基因
deg_results <- list()
for (clu in 0:28) {
    deg <- FindMarkers(ae.ro1,
                       ident.1 = paste0(clu, "_S_A"),
                       ident.2 = paste0(clu, "_S_R"),
                       test.use = "wilcox",
                       slot = "data")
    deg$cluster <- clu
    deg$gene <- rownames(deg)
    deg_results[[clu + 1]] <- deg
}

deg_all <- do.call(rbind, deg_results)
```

### 7.2 排除胁迫响应基因

```R
# 计算各胁迫条件下的表达
deg_all$Ct <- 0  # Control
deg_all$HS <- 0  # Heat Stress
deg_all$ID <- 0  # Iron Deficiency
deg_all$LN <- 0  # Low Nitrogen

for (i in 0:28) {
    for (j in c("batch1", "batch2")) {
        # 计算每个基因在各条件下的平均表达
        # ... (完整代码见 3-Determination of the organ origin...)
    }
}

# 筛选：胁迫条件与对照差异 < 1 的基因
deg_all$is_stressresponse <- "Y"
deg_all[abs(deg_all$HS - deg_all$Ct) < 1 &
        abs(deg_all$ID - deg_all$Ct) < 1 &
        abs(deg_all$LN - deg_all$Ct) < 1,]$is_stressresponse <- "N"

# 排除胁迫响应基因
deg_filtered <- deg_all[deg_all$is_stressresponse == "N",]
```

### 7.3 选择器官特异性标记基因

```R
# 每个聚类选择top 3个叶偏向和根偏向基因
marker_ae <- deg_filtered[deg_filtered$direction == "ae",] %>%
    group_by(cluster) %>%
    top_n(3, mean_score)

marker_ro <- deg_filtered[deg_filtered$direction == "ro",] %>%
    group_by(cluster) %>%
    top_n(3, abs(mean_score))
```

### 7.4 计算器官来源分数

```R
# 为每个细胞计算叶分数和根分数
score_function <- function(i) {
    # 提取叶偏向基因表达
    gene_ae <- marker_ae[marker_ae$cluster == i, "gene"]
    dat_ae <- expm1(expr_data[gene_ae, ])

    # Z-score标准化
    dat_ae_scaled <- t(scale(t(dat_ae)))

    # 计算叶分数（平均Z-score）
    aerial_score <- colMeans(dat_ae_scaled, na.rm = TRUE)

    # 同样计算根分数
    gene_ro <- marker_ro[marker_ro$cluster == i, "gene"]
    dat_ro <- expm1(expr_data[gene_ro, ])
    dat_ro_scaled <- t(scale(t(dat_ro)))
    root_score <- colMeans(dat_ro_scaled, na.rm = TRUE)

    return(data.frame(cell = colnames(dat_ae),
                      aerial = aerial_score,
                      root = root_score))
}

# 应用函数
scores <- do.call(rbind, lapply(0:28, score_function))

# 判断器官来源
scores$tissue <- ifelse(scores$root < scores$aerial, "Aerial",
                ifelse(scores$root > scores$aerial, "Root", "Unknown"))
```

---

## 8. 第六阶段：发育轨迹重建

### 8.1 准备Monocle输入

```R
library(monocle)

# 选择目标细胞
# 聚类 0, 1, 17-22 的叶组织细胞
cell_input <- plt.score[
    plt.score$cluster %in% c(0, 1, 17:22) &
    plt.score$tissue == "Aerial" &
    plt.score$Sample != "Aerial" &
    plt.score$Sample != "Root",]

# 随机抽样10,000细胞
set.seed(123)
cell_input <- cell_input[sample(nrow(cell_input), 10000),]

# 提取表达数据
data_input <- as.matrix(rice.integrated@assays[["integrated"]]@scale.data)
data_input <- data_input[, cell_input$cell]
```

### 8.2 创建Monocle对象

```R
# 基因注释
gene_annotation <- data.frame(
    gene_short_name = gsub("LOC-", "", rownames(data_input)),
    row.names = rownames(data_input)
)

# 细胞信息
sample_sheet <- data.frame(
    Library = rep("mes", ncol(data_input)),
    row.names = colnames(data_input)
)

# 创建CellDataSet
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
mono <- newCellDataSet(
    as(data_input, "sparseMatrix"),
    phenoData = pd,
    featureData = fd
)

# 估计大小因子
mono <- estimateSizeFactors(mono)
mono <- detectGenes(mono, min_expr = 0.1)
```

### 8.3 降维与轨迹推断

```R
# 设置排序基因
ordering_genes <- rownames(data_input)
mono <- setOrderingFilter(mono, ordering_genes)

# DDRTree降维
mono <- reduceDimension(mono,
                        reduction_method = "DDRTree",
                        norm_method = "none",
                        pseudo_expr = 0)

# 计算拟时序
mono <- orderCells(mono)

# 可视化轨迹
plot_cell_trajectory(mono, color_by = "Pseudotime")
plot_cell_trajectory(mono, color_by = "State")

# 设置根状态（发育起点）
mono <- orderCells(mono, root_state = 2)
plot_cell_trajectory(mono, color_by = "Pseudotime")
```

---

## 9. 第七阶段：高级分析与可视化

### 9.1 UMAP可视化

```R
library(ggplot2)
library(patchwork)

# 绘制UMAP图
p1 <- DimPlot(rice.integrated, reduction = "umap",
              group.by = "seurat_clusters",
              label = TRUE, pt.size = 0.1) +
    ggtitle("Clusters")

p2 <- DimPlot(rice.integrated, reduction = "umap",
              group.by = "sample",
              pt.size = 0.1) +
    ggtitle("Samples")

p3 <- DimPlot(rice.integrated, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.1) +
    ggtitle("Cell Types")

# 组合图
p1 + p2 + p3
ggsave("figures/UMAP_overview.pdf", width = 18, height = 6)
```

### 9.2 细胞类型比例分析

```R
# 计算每个样本中各细胞类型的比例
celltype_proportions <- table(mydata$cell_type, mydata$sample)
celltype_proportions <- prop.table(celltype_proportions, margin = 2)

# 绘制堆叠柱状图
df <- as.data.frame(celltype_proportions)
colnames(df) <- c("CellType", "Sample", "Proportion")

ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Type Proportions Across Samples")

ggsave("figures/celltype_proportions.pdf", width = 10, height = 8)
```

### 9.3 差异表达分析

```R
# 找所有细胞类型的标记基因
DefaultAssay(mydata) <- "RNA"
Idents(mydata) <- "cell_type"

all_markers <- FindAllMarkers(mydata,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# 选择top 10标记基因
top10 <- all_markers %>%
    group_by(cluster) %>%
    top_n(10, avg_logFC)

# 绘制热图
DoHeatmap(mydata, features = top10$gene) + NoLegend()
ggsave("figures/marker_heatmap.pdf", width = 15, height = 20)
```

### 9.4 基因表达可视化

```R
# 特定基因的UMAP可视化
FeaturePlot(mydata,
            features = c("LOC_Os01g01050", "LOC_Os01g01070"),
            reduction = "umap",
            pt.size = 0.1)

# 小提琴图
VlnPlot(mydata,
        features = c("LOC_Os01g01050", "LOC_Os01g01070"),
        group.by = "cell_type",
        pt.size = 0)
```

### 9.5 胁迫响应分析

```R
# 比较对照和胁迫条件下的细胞类型比例变化
conditions <- c("Ctrl", "HS", "ID", "LN")

proportion_changes <- list()
for (cond in conditions) {
    cells_cond <- colnames(mydata)[mydata$condition == cond]
    prop_cond <- table(mydata$cell_type[cells_cond]) / length(cells_cond)
    proportion_changes[[cond]] <- prop_cond
}

# 绘制比较图
df_compare <- do.call(rbind, lapply(names(proportion_changes), function(x) {
    data.frame(Condition = x,
               CellType = names(proportion_changes[[x]]),
               Proportion = as.numeric(proportion_changes[[x]]))
}))

ggplot(df_compare, aes(x = CellType, y = Proportion, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Type Proportions Under Different Conditions")
```

---

## 10. 常见问题与解决方案

### 10.1 Cell Ranger 问题

**问题**: Cell Ranger 报错 "No valid barcodes found"
**解决**:
```bash
# 检查FASTQ文件格式
zcat sample_1.fastq.gz | head -4
# 确认barcode在正确的位置

# 尝试指定chemistry
cellranger count --chemistry=threeprime ...
```

**问题**: 内存不足
**解决**:
```bash
# 减少线程数，增加内存
cellranger count --nthreads=8 --memgb=64 ...
```

### 10.2 Seurat 问题

**问题**: SCTransform 太慢
**解决**:
```R
# 使用vst方法替代
rice.list[[i]] <- SCTransform(rice.list[[i]], method = "glmGamPoi")
```

**问题**: 整合后批次效应仍然明显
**解决**:
```R
# 增加锚点数
rice.anchors <- FindIntegrationAnchors(..., k.filter = 200)
```

### 10.3 Monocle 问题

**问题**: reduceDimension 报错
**解决**:
```R
# 检查基因数量
length(ordering_genes)
# 应该 > 100

# 尝试减少基因数
ordering_genes <- ordering_genes[1:1000]
```

**问题**: 轨迹图看起来混乱
**解决**:
```R
# 调整参数
mono <- reduceDimension(mono,
                        reduction_method = "DDRTree",
                        max_components = 2,
                        ...)
```

---

## 附录

### A. 文件路径参考

| 文件 | 路径 | 说明 |
|------|------|------|
| Cell Ranger输出 | `./cellranger_output/${sample}_result/outs/` | 每个样本的结果 |
| 表达矩阵 | `./seurat_analysis/all.csv` | 合并后的表达矩阵 |
| Seurat对象 | `./seurat_analysis/rice_integrated.rds` | 整合后的Seurat对象 |
| Monocle对象 | `./seurat_analysis/monocle_object.rds` | 轨迹分析对象 |

### B. 关键参数总结

| 参数 | 值 | 说明 |
|------|-----|------|
| SCTransform特征数 | 3000 | 整合用的高变基因数 |
| PCA维度 | 100 | 使用的主成分数 |
| UMAP维度 | 10 | UMAP使用的维度数 |
| 聚类分辨率 | 0.75 | 产生29个聚类 |
| 原生质体FC阈值 | 3 (log2) | 筛选敏感基因的阈值 |
| 轨迹分析细胞数 | 10,000 | Monocle分析的抽样数 |

### C. 输出文件清单

```
results/
├── cluster_assignment.csv      # 聚类到细胞类型的映射
├── organ_origin.csv            # 细胞的器官来源
├── pseudotime.csv              # 拟时序信息
├── differentially_expressed/   # 差异表达基因列表
└── marker_genes/               # 各细胞类型的标记基因

figures/
├── UMAP_overview.pdf           # UMAP总览图
├── celltype_proportions.pdf    # 细胞类型比例图
├── marker_heatmap.pdf          # 标记基因热图
├── organ_origin.pdf            # 器官来源图
└── trajectory.pdf              # 发育轨迹图
```

---

**文档创建时间**: 2024-05-24
**最后更新**: 2024-05-24
**作者**: Claude Code Assistant
