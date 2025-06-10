# Seurat2AnnData 转换工具

这个工具可以将Seurat对象（.rds格式）转换为AnnData对象（.h5ad格式），方便在Python的scanpy环境中进行单细胞分析。

## 功能特点

- ✅ 支持命令行参数，使用简单
- ✅ 保留完整的基因表达矩阵
- ✅ 保留所有细胞元数据
- ✅ 保留降维结果（PCA、UMAP、t-SNE）
- ✅ 自动处理临时文件清理
- ✅ 详细的转换进度显示
- ✅ 错误处理和验证

## 使用方法

### 基本语法
```bash
python Seurat2AnnData.py <输入RDS文件> <输出H5AD文件>
```

### 示例

1. **基本转换**：
```bash
python Seurat2AnnData.py 03cell_annonation.rds 03cell_annonation.h5ad
```

2. **指定路径**：
```bash
python Seurat2AnnData.py ./data/seurat_object.rds ./results/converted_data.h5ad
```

3. **完整路径**：
```bash
python Seurat2AnnData.py /path/to/input.rds /path/to/output.h5ad
```

## 输出信息

转换完成后，脚本会显示详细的数据信息：
```
=== 转换完成！===
输入文件: 03cell_annonation.rds
输出文件: 03cell_annonation.h5ad
数据信息:
  - 细胞数: 93,056
  - 基因数: 38,606
  - 观测元数据列: 9 列
  - 降维结果: ['X_pca', 'X_umap', 'X_tsne']
  - 文件大小: 3.52 GB
```

## 依赖要求

### R环境
- Seurat
- SeuratObject  
- Matrix

### Python环境
- pandas
- numpy
- anndata
- scanpy
- scipy

## 技术原理

脚本采用两步转换法：

1. **第一步（R）**：将Seurat对象导出为标准格式
   - 基因表达矩阵 → MTX格式
   - 细胞元数据 → CSV格式
   - 基因信息 → CSV格式
   - 降维结果 → CSV格式

2. **第二步（Python）**：从导出文件重建AnnData对象
   - 读取所有导出数据
   - 创建AnnData对象
   - 添加降维结果
   - 保存为h5ad格式

## 注意事项

- 确保R和Python环境都已正确配置
- 输入文件必须是有效的Seurat RDS文件
- 需要足够的磁盘空间（通常输出文件大小与输入文件相近）
- 转换过程中会创建临时文件，完成后自动清理

## 错误处理

脚本包含完善的错误处理机制：
- 检查输入文件是否存在
- 验证文件格式
- 检查必需的依赖库
- 自动清理临时文件

如果遇到问题，脚本会显示详细的错误信息。 
