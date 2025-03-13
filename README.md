# DNA自我互补结构分析工具

一个分析DNA寡核苷酸中flank区域与内部序列之间可能形成自我互补结构的工具。这些自我互补结构可能会影响PCR扩增效率，特别是当flank区域（引物结合位点）与内部序列之间形成稳定的二级结构时。

## 背景

在PCR扩增过程中，某些DNA寡核苷酸(oligo)的扩增效率显著低于其他oligo，尽管它们使用相同的引物结合位点(flank序列)。这可能是由于oligo中flank区域与内部序列之间形成稳定的自我互补结构(二级结构)，进而阻碍了PCR引物与flank区域的结合所导致的。

本工具通过计算机算法分析这种互补性，提供量化评分和可视化结果，帮助研究人员理解和预测PCR扩增效率差异。

## 安装

### 系统要求

- Python 3.6 或更高版本
- pip（Python包管理器）

### 安装步骤

1. 克隆或下载此仓库到本地

```bash
git clone https://github.com/bwin1990/highAndLowDepth.git
cd highAndLowDepth
```

2. 安装所需依赖项

```bash
pip install -r requirements.txt
```

## 使用方法

本工具提供图形用户界面(GUI)和命令行界面两种使用方式，可以进行以下分析：

1. 分析单个DNA序列的自我互补结构
2. 比较高效率和低效率oligo组的自我互补性差异

### 图形用户界面(GUI)

启动GUI界面：

```bash
python oligo_analyzer_cli.py
```

GUI界面包含两个主要标签页：

#### 单序列分析标签页

- 输入DNA序列
- 设置分析参数：
  - Flank长度：两端flank区域的长度
  - 窗口大小：内部序列滑动窗口大小
  - 最小匹配：最小互补匹配长度
  - Flank分析模式：选择分析哪一侧的flank（both=两端flank, 5p=仅5'端flank, 3p=仅3'端flank）
- 点击"分析序列"按钮进行分析
- 结果区域显示互补评分和可视化图表
- 可以保存分析结果

#### 组比较分析标签页

- 选择高效率和低效率oligo序列文件
- 设置分析参数（同单序列分析）
- 选择结果输出文件（可选）
- 点击"比较分析"按钮进行分析
- 结果区域显示两组序列的互补性差异和统计图表

### 命令行界面

#### 分析单个序列

```bash
python oligo_analyzer_cli.py single "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" --flank 10 --mode both
```

参数说明：
- `single`：表示进行单序列分析
- 序列：引号中的DNA序列
- `--flank`：指定flank区域长度，默认为20
- `--window`：指定滑动窗口大小，默认为10
- `--min-match`：指定最小匹配长度，默认为4
- `--mode`：指定flank分析模式（both, 5p, 3p），默认为both
- `--output`或`-o`：指定输出图像文件路径（可选）

#### 比较两组序列

```bash
python oligo_analyzer_cli.py compare example_data/high_efficiency_oligos.fasta example_data/low_efficiency_oligos.fasta --flank 10 --mode both --output results.xlsx --plot comparison.png
```

参数说明：
- `compare`：表示进行组间比较
- 第一个文件路径：高效率oligo序列文件
- 第二个文件路径：低效率oligo序列文件
- `--flank`：指定flank区域长度，默认为20
- `--window`：指定滑动窗口大小，默认为10
- `--min-match`：指定最小匹配长度，默认为4
- `--mode`：指定flank分析模式（both, 5p, 3p），默认为both
- `--output`或`-o`：指定结果输出文件路径（可选）
- `--plot`或`-p`：指定比较图表输出文件路径（可选）

### 支持的文件格式

- FASTA格式 (.fasta, .fa)
- 文本文件，每行一个序列 (.txt)
- CSV文件，第一列为序列 (.csv)
- TSV文件，第一列为序列 (.tsv)

### 输出格式

- 文本报告 (.txt)
- JSON格式 (.json)
- CSV格式 (.csv)
- Excel格式 (.xlsx)
- 图像格式 (.png, .jpg, .pdf, .svg)

## Flank分析模式

本工具提供三种flank分析模式：

1. **both**：同时分析5'端和3'端flank与内部序列的互补性（默认模式）
2. **5p**：只分析5'端flank与内部序列的互补性
3. **3p**：只分析3'端flank与内部序列的互补性

这使得研究人员可以更精确地分析特定端部的互补结构，特别是在只关注一侧引物结合位点的情况下。

## 示例数据

本仓库包含示例数据，位于`example_data`目录中：
- `high_efficiency_oligos.fasta`：高扩增效率的oligo序列
- `low_efficiency_oligos.fasta`：低扩增效率的oligo序列

## 结果解释

### 自我互补评分

评分越高，表示序列中flank与内部区域之间可能形成更稳定的互补结构，这可能会阻碍PCR引物的结合，从而降低扩增效率。

### 可视化图表

- **结构可视化**：显示单个序列中所有检测到的互补区域，包括位置、长度和评估的自由能
- **比较图表**：通过箱线图和小提琴图比较两组序列的自我互补性差异，包括统计信息

## 进阶用法

### 在Python脚本中使用

```python
from dna_analyzer import DNAComplementAnalyzer

# 创建分析器实例
analyzer = DNAComplementAnalyzer(flank_length=10, window_size=10, min_match=4, flank_mode='both')

# 分析单个序列
seq = "ACGTACGTACGTATCGATCGATCGATCGATCGATCGATCGATCGACGTACGTACGT"
results, score = analyzer.analyze_sequence(seq)

# 打印结果
print(f"总互补评分: {score:.2f}")
print(f"发现 {len(results)} 个潜在互补区域")

# 可视化结构
analyzer.visualize_structure(seq, results)
```

## 技术说明

该工具使用以下方法分析DNA序列：

1. **滑动窗口法**：在内部序列上滑动窗口，检测与flank区域的潜在互补
2. **评分系统**：基于匹配长度和数量计算互补性评分
3. **自由能估算**：简化版计算互补配对的预估自由能
4. **可视化**：生成直观的图形表示互补结构

## 注意事项

- 关闭GUI窗口时，程序会自动清理资源并完全退出
- 分析大量序列时可能需要较长时间，请耐心等待
- 对于非常长的序列，可能需要调整flank长度和窗口大小参数

## 开发者

如果您希望为此项目做出贡献，请参考以下步骤：

1. Fork此仓库
2. 创建您的特性分支 (`git checkout -b feature/amazing-feature`)
3. 提交您的更改 (`git commit -m 'Add some amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 打开Pull Request

## 许可证

本项目遵循MIT许可证。详情请见[LICENSE](LICENSE)文件。

## 引用

如果您在研究中使用了此工具，请引用：

```
zengcheng. (2025). DNA自我互补结构分析工具. GitHub仓库. https://github.com/bwin1990/dna-complement-analyzer
``` 