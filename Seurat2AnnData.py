#!/usr/bin/env python3
"""
Seurat2AnnData.py - 将Seurat对象转换为AnnData格式的通用脚本

使用方法:
    python Seurat2AnnData.py <input.rds> <output.h5ad>

示例:
    python Seurat2AnnData.py 03cell_annonation.rds 03cell_annonation.h5ad
    python Seurat2AnnData.py data/seurat_obj.rds results/anndata_obj.h5ad
"""

import os
import sys
import subprocess
import tempfile
import shutil
import pandas as pd
import numpy as np
import anndata as ad
from scipy.io import mmread
from pathlib import Path

def print_usage():
    """打印使用说明"""
    print("用法: python Seurat2AnnData.py <输入RDS文件> <输出H5AD文件>")
    print()
    print("示例:")
    print("  python Seurat2AnnData.py 03cell_annonation.rds 03cell_annonation.h5ad")
    print("  python Seurat2AnnData.py data/seurat_obj.rds results/anndata_obj.h5ad")
    print()
    print("参数:")
    print("  输入RDS文件  - Seurat对象的RDS文件路径")
    print("  输出H5AD文件 - 要生成的AnnData h5ad文件路径")

def create_r_export_script(input_rds_path, temp_dir):
    """创建R脚本来导出Seurat数据"""
    r_script_content = f'''
# Seurat对象数据导出脚本
rm(list = ls())

print("开始加载R库...")
library(Seurat)
library(SeuratObject)
library(Matrix)

# 检查输入文件
input_file <- "{input_rds_path}"
if (!file.exists(input_file)) {{
  stop("错误：找不到输入文件 ", input_file)
}}

print("读取Seurat对象...")
sce <- readRDS(input_file)
print(paste("对象读取成功，包含", ncol(sce), "个细胞，", nrow(sce), "个基因"))

# 设置输出目录
output_dir <- "{temp_dir}"
if (!dir.exists(output_dir)) {{
  dir.create(output_dir, recursive = TRUE)
}}

print("导出基因表达矩阵...")
# 获取表达矩阵
expr_matrix <- GetAssayData(sce, layer = "counts")
# 保存为稀疏矩阵格式
writeMM(expr_matrix, file.path(output_dir, "matrix.mtx"))

print("导出基因信息...")
# 导出基因信息
genes <- data.frame(
  gene_id = rownames(expr_matrix),
  gene_symbol = rownames(expr_matrix)
)
write.csv(genes, file.path(output_dir, "features.csv"), row.names = FALSE)

print("导出细胞元数据...")
# 导出细胞元数据
cell_metadata <- sce@meta.data
cell_metadata$barcode <- rownames(cell_metadata)
write.csv(cell_metadata, file.path(output_dir, "metadata.csv"), row.names = FALSE)

print("导出细胞条形码...")
# 导出细胞条形码
barcodes <- data.frame(barcode = colnames(expr_matrix))
write.csv(barcodes, file.path(output_dir, "barcodes.csv"), row.names = FALSE)

# 检查是否有降维结果
if ("pca" %in% names(sce@reductions)) {{
  print("导出PCA结果...")
  pca_coords <- Embeddings(sce, reduction = "pca")
  write.csv(pca_coords, file.path(output_dir, "pca.csv"))
}}

if ("umap" %in% names(sce@reductions)) {{
  print("导出UMAP结果...")
  umap_coords <- Embeddings(sce, reduction = "umap")
  write.csv(umap_coords, file.path(output_dir, "umap.csv"))
}}

if ("tsne" %in% names(sce@reductions)) {{
  print("导出t-SNE结果...")
  tsne_coords <- Embeddings(sce, reduction = "tsne")
  write.csv(tsne_coords, file.path(output_dir, "tsne.csv"))
}}

print("R数据导出完成！")
'''
    return r_script_content

def run_r_export(input_rds_path, temp_dir):
    """运行R脚本导出Seurat数据"""
    print("第1步：使用R导出Seurat数据...")
    
    # 创建R脚本
    r_script_content = create_r_export_script(input_rds_path, temp_dir)
    r_script_path = os.path.join(temp_dir, "export_seurat.R")
    
    with open(r_script_path, 'w') as f:
        f.write(r_script_content)
    
    # 运行R脚本
    try:
        result = subprocess.run(['Rscript', r_script_path], 
                              capture_output=True, text=True, check=True)
        print(result.stdout)
        if result.stderr:
            print("R警告信息:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"R脚本执行失败: {e}")
        print("标准输出:", e.stdout)
        print("错误输出:", e.stderr)
        raise

def create_anndata_from_exports(temp_dir, output_h5ad_path, input_rds_path):
    """从导出的文件创建AnnData对象"""
    print("第2步：从导出文件创建AnnData对象...")
    
    # 检查必需文件
    required_files = ["matrix.mtx", "features.csv", "metadata.csv", "barcodes.csv"]
    for file in required_files:
        file_path = os.path.join(temp_dir, file)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"找不到必需文件: {file}")
    
    print("读取基因表达矩阵...")
    # 读取表达矩阵 (基因 x 细胞)
    X = mmread(os.path.join(temp_dir, "matrix.mtx")).T.tocsr()  # 转置为细胞 x 基因
    
    print("读取基因信息...")
    # 读取基因信息
    genes = pd.read_csv(os.path.join(temp_dir, "features.csv"))
    
    print("读取细胞元数据...")
    # 读取细胞元数据
    metadata = pd.read_csv(os.path.join(temp_dir, "metadata.csv"), index_col=0)
    
    print(f"数据维度: {X.shape[0]} 细胞 x {X.shape[1]} 基因")
    
    # 创建AnnData对象
    print("创建AnnData对象...")
    adata = ad.AnnData(
        X=X,
        obs=metadata,
        var=genes.set_index('gene_symbol')
    )
    
    # 添加降维结果（如果存在）
    embedding_files = {
        'X_pca': 'pca.csv',
        'X_umap': 'umap.csv',
        'X_tsne': 'tsne.csv'
    }
    
    for embedding_name, filename in embedding_files.items():
        file_path = os.path.join(temp_dir, filename)
        if os.path.exists(file_path):
            print(f"添加{embedding_name}结果...")
            embedding = pd.read_csv(file_path, index_col=0)
            adata.obsm[embedding_name] = embedding.values
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_h5ad_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # 保存AnnData对象
    print(f"保存AnnData对象到 {output_h5ad_path}...")
    adata.write(output_h5ad_path)
    
    # 显示转换结果
    print("\n=== 转换完成！===")
    print(f"输入文件: {input_rds_path}")
    print(f"输出文件: {output_h5ad_path}")
    print(f"数据信息:")
    print(f"  - 细胞数: {adata.n_obs:,}")
    print(f"  - 基因数: {adata.n_vars:,}")
    print(f"  - 观测元数据列: {len(adata.obs.columns)} 列")
    print(f"  - 降维结果: {list(adata.obsm.keys())}")
    
    # 检查文件大小
    file_size = os.path.getsize(output_h5ad_path) / (1024**3)  # GB
    print(f"  - 文件大小: {file_size:.2f} GB")

def main():
    """主函数"""
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("错误：参数数量不正确")
        print()
        print_usage()
        sys.exit(1)
    
    input_rds_path = sys.argv[1]
    output_h5ad_path = sys.argv[2]
    
    # 检查输入文件
    if not os.path.exists(input_rds_path):
        print(f"错误：输入文件不存在: {input_rds_path}")
        sys.exit(1)
    
    # 检查文件扩展名
    if not input_rds_path.lower().endswith('.rds'):
        print("警告：输入文件不是.rds格式，请确认文件格式正确")
    
    if not output_h5ad_path.lower().endswith('.h5ad'):
        print("警告：输出文件不是.h5ad格式，建议使用.h5ad扩展名")
    
    print("=== Seurat到AnnData转换工具 ===")
    print(f"输入文件: {input_rds_path}")
    print(f"输出文件: {output_h5ad_path}")
    print()
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # 第1步：导出Seurat数据
            run_r_export(input_rds_path, temp_dir)
            
            # 第2步：创建AnnData对象
            create_anndata_from_exports(temp_dir, output_h5ad_path, input_rds_path)
            
        except Exception as e:
            print(f"转换失败: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main() 