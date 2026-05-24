#!/bin/bash
##############################################################################
#  水稻幼苗单细胞转录组 - Cell Ranger 处理流程
#  数据: GSA CRA004082 (Wang et al. 2021 JGG)
##############################################################################

# ============================================================================
# 0. 前置要求
# ============================================================================
# - Cell Ranger v3.1.0 (论文使用的版本)
#   下载: https://www.10xgenomics.com/support/software/cell-ranger/downloads
# - 至少 32GB RAM, 200GB 磁盘空间
# - 水稻参考基因组

# ============================================================================
# 1. 设置目录结构
# ============================================================================

WORK_DIR="/home/user/rice_scrna"  # 修改为你的工作目录
mkdir -p ${WORK_DIR}/{raw_data,reference,fastqs,cellranger_output}
cd ${WORK_DIR}

echo "=============================="
echo "Step 1: 准备水稻参考基因组"
echo "=============================="

# 下载水稻参考基因组 IRGSP-1.0
# 方法1: 从 RAP-DB 下载
cd ${WORK_DIR}/reference
wget https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz
gunzip IRGSP-1.0_genome.fasta.gz

wget https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_gene.gff3.gz
gunzip IRGSP-1.0_representative_gene.gff3.gz

# 用 Cell Ranger 构建参考索引
cellranger mkref --genome=rice \
  --fasta=IRGSP-1.0_genome.fasta \
  --genes=IRGSP-1.0_representative_gene.gff3 \
  --nthreads=16

echo "参考基因组构建完成: ${WORK_DIR}/reference/rice/"

# ============================================================================
# 2. 下载原始 FASTQ 数据 (从 GSA CRA004082)
# ============================================================================
# 注意: 每个样本有多个 lane，需要全部下载
# 总数据量约 2TB，建议按样本逐步下载

cd ${WORK_DIR}/raw_data

# --- 对照组: Seedling-1 (4 lanes, ~125GB) ---
echo "下载 Seedling-1..."
mkdir -p Seedling-1
cd Seedling-1
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279103/CRR279103.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279104/CRR279104.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279105/CRR279105.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279106/CRR279106.tar
for f in *.tar; do tar xf $f; done
cd ..

# --- 叶片: Leaf-1 (4 lanes, ~90GB) ---
echo "下载 Leaf-1..."
mkdir -p Leaf-1
cd Leaf-1
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279087/CRR279087.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279088/CRR279088.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279089/CRR279089.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279090/CRR279090.tar
for f in *.tar; do tar xf $f; done
cd ..

# --- 根部: Root-1 (4 lanes, ~80GB) ---
echo "下载 Root-1..."
mkdir -p Root-1
cd Root-1
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279095/CRR279095.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279096/CRR279096.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279097/CRR279097.tar
wget ftp://download.big.ac.cn/gsa/CRA004082/CRR279098/CRR279098.tar
for f in *.tar; do tar xf $f; done
cd ..

# --- 其他样本按同样模式下载 ---
# Seedling-2: CRR279123-CRR279126
# Seedling-3: CRR279143-CRR279146
# Seedling-4: CRR279159-CRR279162
# Leaf-2:     CRR279091-CRR279094
# Root-2:     CRR279099-CRR279102
# LN-1~4, HS-1~4, ID-1~4: 查看 CRA004082.csv 获取 CRR 编号

echo "数据下载完成!"

# ============================================================================
# 3. 整理 FASTQ 文件（Cell Ranger 需要特定命名格式）
# ============================================================================
# Cell Ranger 需要 FASTQ 文件命名为: <SampleName>_S1_L00<Lane>_R1/R2/I1_001.fastq.gz
# GSA 下载的文件名格式不同，需要重命名

cd ${WORK_DIR}

# 整理函数: 将 GSA 下载的 FASTQ 重命名为 Cell Ranger 格式
organize_fastq() {
  local SAMPLE_DIR=$1
  local SAMPLE_NAME=$2
  local OUTPUT_DIR="${WORK_DIR}/fastqs/${SAMPLE_NAME}"
  mkdir -p ${OUTPUT_DIR}

  local LANE=1
  for TARBALL in ${WORK_DIR}/raw_data/${SAMPLE_DIR}/*.tar; do
    tar tf ${TARBALL} | head -20  # 查看 tar 内的文件名格式
    # 解压后重命名
    tar xf ${TARBALL} -C ${OUTPUT_DIR}/
  done

  # 根据实际文件名进行重命名
  # GSA 的 10x 文件通常解压后是: <Run>_1.fastq.gz, <Run>_2.fastq.gz
  cd ${OUTPUT_DIR}
  for R1 in *_1.fastq.gz; do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    BASE=$(basename ${R1} _1.fastq.gz)
    # 重命名为 Cell Ranger 格式
    if [ -f "$R1" ] && [ -f "$R2" ]; then
      mv ${R1} "${SAMPLE_NAME}_S1_L00${LANE}_R1_001.fastq.gz"
      mv ${R2} "${SAMPLE_NAME}_S1_L00${LANE}_R2_001.fastq.gz"
      LANE=$((LANE + 1))
    fi
  done
  echo "Organized ${SAMPLE_NAME}: ${LANE} lanes"
}

# 执行整理
organize_fastq "Seedling-1" "Seedling1"
organize_fastq "Leaf-1" "Leaf1"
organize_fastq "Root-1" "Root1"
# ... 对其他样本同样操作

# ============================================================================
# 4. 运行 Cell Ranger
# ============================================================================

REF="${WORK_DIR}/reference/rice"
FASTQ_DIR="${WORK_DIR}/fastqs"
OUTPUT_DIR="${WORK_DIR}/cellranger_output"
mkdir -p ${OUTPUT_DIR}

# 运行 Cell Ranger count
run_cellranger() {
  local SAMPLE=$1
  echo "=============================="
  echo "Running Cell Ranger on: ${SAMPLE}"
  echo "=============================="

  cd ${OUTPUT_DIR}
  cellranger count \
    --id=${SAMPLE}_result \
    --transcriptome=${REF} \
    --fastqs=${FASTQ_DIR}/${SAMPLE} \
    --sample=${SAMPLE} \
    --expect-cells=10000 \
    --nthreads=16 \
    --memgb=64

  # 创建 R 代码需要的目录结构
  mkdir -p ${OUTPUT_DIR}/${SAMPLE}/outs/
  ln -sf ${OUTPUT_DIR}/${SAMPLE}_result/outs/filtered_feature_bc_matrix \
         ${OUTPUT_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix
}

# 依次处理每个样本（或用 nohup 后台运行）
run_cellranger "Seedling1"
run_cellranger "Leaf1"
run_cellranger "Root1"
# run_cellranger "Seedling2"
# run_cellranger "Seedling3"
# run_cellranger "Seedling4"
# run_cellranger "Leaf2"
# run_cellranger "Root2"
# run_cellranger "LN1"
# run_cellranger "LN2"
# run_cellranger "LN3"
# run_cellranger "LN4"
# run_cellranger "HS1"
# run_cellranger "HS2"
# run_cellranger "HS3"
# run_cellranger "HS4"
# run_cellranger "ID1"
# run_cellranger "ID2"
# run_cellranger "ID3"
# run_cellranger "ID4"

echo "=============================="
echo "Cell Ranger 处理全部完成!"
echo "下一步: 在 R 中运行 Rice_Seedling_scRNAseq_Analysis.R"
echo "=============================="

# ============================================================================
# 5. 批量运行提示（后台运行所有样本）
# ============================================================================

# 如果要一次性提交所有样本，可以用以下脚本:
# for sample in Seedling1 Seedling2 Seedling3 Seedling4 \
#                Leaf1 Leaf2 Root1 Root2 \
#                LN1 LN2 LN3 LN4 \
#                HS1 HS2 HS3 HS4 \
#                ID1 ID2 ID3 ID4; do
#   nohup bash -c "run_cellranger ${sample}" > ${sample}_cellranger.log 2>&1 &
# done
