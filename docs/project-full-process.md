# 项目全流程说明

## 1. 项目定位

这是一个面向网络药理学研究的半自动分析流水线，核心目标是把“药物成分筛选 -> 靶点预测 -> 疾病靶点收集 -> 交集网络构建 -> hub 基因筛选 -> 富集分析 -> MR -> 分子对接 -> 单细胞验证 -> 汇总出图”串成一个可复现流程。

项目的执行框架是：

- `config.yaml` 负责参数配置
- `Snakefile` 负责任务编排
- `scripts/python/` 负责主流程 Python 步骤
- `scripts/R/` 负责富集、MR、单细胞等 R 步骤
- `docs/` 负责使用说明和 SOP

它不是全自动抓数系统，而是“自动分析 + 手工数据补位”的研究模板。

## 2. 适用场景

适合以下任务：

- 中药/天然药物网络药理学分析
- 单一药物对单一疾病的机制探索
- 论文方法部分的流程模板复用
- 教学演示、课题预分析、方法学复现

不适合以下任务：

- 完全无人值守的数据采集
- 大规模高通量批量药物筛选平台
- 直接作为临床结论输出系统

## 3. 项目执行总入口

### 3.1 配置入口

`config.yaml` 定义了整个项目的研究对象和阈值，包括：

- 药物名称、别名、拉丁名
- 疾病名称、缩写、MeSH、GEO accession
- ADME 阈值
- 靶点预测数据库
- 疾病靶点数据库
- STRING 参数
- hub 基因筛选参数
- 富集分析参数
- MR 参数
- 分子对接参数
- scRNA / spatial 参数
- 输出与可复现性参数

### 3.2 工作流入口

`Snakefile` 是总控文件，定义了 `step00` 到 `step10` 的依赖关系。你可以：

- 单步运行某一脚本
- 用 `snakemake --cores N` 跑全流程
- 用 `snakemake --dry-run` 预览依赖

## 4. 全流程逐步说明

### Step 00: 预检

脚本：`scripts/python/00_preflight.py`

作用：

- 检查 Python 依赖是否安装
- 检查 R 包是否安装
- 检查外部工具是否可用，如 `vina`、`obabel`、`gmx`
- 检查 `config.yaml` 是否存在占位字段
- 检查前序数据文件是否齐备
- 识别结果文件里是否包含 placeholder

这是正式运行前的入口门禁。项目里很多“看似跑完”的结果可能只是占位结果，这一步就是为了防止误用。

### Step 01: 成分收集与 ADME 筛选

脚本：`scripts/python/01_compounds.py`

输入：

- `config.yaml`
- 可选的 TCMSP 导出文件

核心逻辑：

- 统一 TCMSP 导出列名
- 支持 PubChem 名称检索
- 若缺少完整 ADME 信息，可用 RDKit 补算部分理化描述符
- 对 SMILES 做合法性校验
- 按 `OB`、`DL`、`MW`、`logP` 做筛选

输出：

- `data/processed/01_compounds_filtered.csv`
- 原始查询缓存文件

这一步的目标是拿到“可进入后续分析的候选化合物集合”。

### Step 02: 多数据库靶点预测

脚本：`scripts/python/02_targets.py`

输入：

- `01_compounds_filtered.csv`
- 手工下载的 STP / PharmMapper / CTD / SEA / TCMSP 文件

核心逻辑：

- 兼容不同数据库导出格式
- 自动识别基因列、分数字段、UniProt 字段
- 调用 `mygene` 或 MyGene.info 将标识统一成 HGNC symbol
- 尝试在线请求 SwissTargetPrediction
- 请求失败时输出手工补录说明

输出：

- `data/processed/02_targets_merged.csv`
- 原始合并缓存

这一步的目标是把药物候选成分映射成较规范的“药物潜在靶点集合”。

### Step 03: 疾病靶点收集

脚本：`scripts/python/03_disease.py`

输入：

- `config.yaml`
- 可选的 GeneCards / OMIM / DisGeNET 导出文件

核心逻辑：

- 从多个疾病数据库收集候选靶点
- 允许走 DisGeNET API
- 如果 API 不可用，转为手工导出模式
- 统一基因名并做清洗
- 按相关性或分数阈值过滤

输出：

- `data/processed/03_disease_targets.csv`

这一步构建的是“疾病相关基因集合”。

### Step 04: 交集分析与 PPI 网络

脚本：`scripts/python/04_ppi.py`

输入：

- `02_targets_merged.csv`
- `03_disease_targets.csv`

核心逻辑：

- 求药物靶点与疾病靶点交集
- 绘制 Venn 图
- 调 STRING API 拉取 PPI 网络
- 额外导出 Cytoscape 可用的 `.sif` 和节点表

输出：

- `data/processed/04_intersection_genes.csv`
- `data/processed/04_ppi_edges.csv`
- `results/figures/04_venn.pdf`

这一步把“候选靶点”变成“网络层面的核心候选集合”。

### Step 05: Hub 基因筛选

脚本：`scripts/python/05_hub_genes.py`

输入：

- `04_intersection_genes.csv`
- `04_ppi_edges.csv`
- 可选 `mcode_clusters.csv`
- 可选 GEO 表达矩阵与分组标签

核心逻辑分三层：

- 网络拓扑层：计算 degree 和 betweenness
- 模块层：读取 Cytoscape/MCODE 导出的聚类基因
- 机器学习层：对 GEO 表达矩阵运行 LASSO 和 Random Forest

最终策略：

- 对三层结果做投票
- 若有效信息不足，回退到拓扑候选

输出：

- `data/processed/05_hub_genes.csv`

这一步决定了后续富集、MR、分子对接、单细胞验证围绕哪些基因展开。

### Step 06: 富集分析

脚本：

- `scripts/python/06_enrichment.py`
- 运行期动态生成 R 脚本

输入：

- `05_hub_genes.csv`
- `04_intersection_genes.csv`

核心逻辑：

- 调用 `clusterProfiler`
- 进行 GO 和 KEGG 富集
- 可选 GSVA
- 自动生成结果图和结果表

输出：

- `data/processed/06_enrichment_go.csv`
- `data/processed/06_enrichment_kegg.csv`
- `results/figures/06_enrichment_go.pdf`
- `results/figures/06_enrichment_kegg.pdf`

这一步回答的是“这些候选靶点主要聚焦在哪些生物过程和信号通路上”。

### Step 07: MR 分析

脚本：`scripts/R/07_mr.R`

输入：

- `05_hub_genes.csv`
- eQTLGen cis-eQTL 文件
- FinnGen GWAS 文件

核心逻辑：

- 读取 hub 基因
- 以 eQTL 作为暴露
- 以 FinnGen 疾病 GWAS 作为结局
- 做格式化、clumping、harmonise、MR 分析
- 支持 IVW、Egger、Weighted Median 等方法

容错特性：

- 如果缺失 GWAS/eQTL 文件，不会直接崩溃
- 会生成 placeholder 结果和提示图

输出：

- `data/processed/07_mr_results.csv`
- `results/figures/07_mr_forest.pdf`

这一步的作用是给 hub 基因提供因果方向上的额外支持。

### Step 08: 分子对接

脚本：`scripts/python/08_docking.py`

输入：

- `01_compounds_filtered.csv`
- `05_hub_genes.csv`

核心逻辑：

- 内置一批常见靶点对应的 PDB 模板
- 自动下载 PDB 结构
- 从共晶配体估计结合中心
- 调 Open Babel 生成配体三维结构和 `pdbqt`
- 调 Vina 进行对接
- 生成 GROMACS 后续 MD 指南

输出：

- `data/processed/08_docking_scores.csv`
- `results/docking/` 下的受体、配体、对接产物及说明

这一步把网络药理学预测延伸到结构层面验证。

### Step 09: 单细胞 / 空间转录组分析

脚本：`scripts/R/09_scrna.R`

输入：

- `05_hub_genes.csv`
- `data/raw/geo_datasets/<GSE>/`

核心逻辑：

- Seurat 标准流程
- 多样本时用 Harmony 做批次校正
- 绘制 UMAP
- 统计 hub 基因在不同 cluster 的表达
- 可选 Nebulosa KDE
- 若启用 spatial，则做 Visium 空间可视化

输出：

- `results/figures/09_umap_clusters.pdf`
- `results/figures/09_hub_target_dotplot.pdf`
- `data/processed/09_cell_type_expression.csv`

这一步回答的是“关键靶点在细胞类型和空间层面是否有特异表达”。

### Step 09b: 免疫浸润补充分支

脚本：`scripts/python/09b_immune.py`

作用：

- 基于 bulk 表达矩阵做 ssGSEA 免疫浸润评分
- 输出热图和相关结果

这是一个与 Step 09 并列的补充分析模块。

### Step 10: 汇总出图

脚本：`scripts/python/10_visualization.py`

输入：

- Step 04, 05, 07, 08 的结果文件

核心逻辑：

- 汇总最终 Venn 图
- 绘制 PPI 网络图
- 绘制 hub 基因热图
- 绘制 MR forest plot
- 绘制 docking bar chart
- 生成 CETSA / SPR 模板图

输出：

- `results/figures/10_summary_figure.pdf`
- 多个子图 PDF

这是用于论文制图、汇报展示和结果交付的最后一步。

## 5. 自动化与手工介入边界

这个项目最重要的理解点是：它并不是一个完全自动采集系统。

自动化部分：

- 配置读取
- 中间数据清洗
- 统一基因符号
- 统计筛选
- STRING 网络
- 富集分析
- 对接与出图

手工部分：

- TCMSP 导出
- GeneCards / OMIM / 部分 DisGeNET 导出
- MCODE 聚类结果
- MR 所需大文件下载
- GEO 单细胞原始文件整理

因此它更准确的定位是“研究编排模板”，不是“黑盒一键出结论工具”。

## 6. 输入、输出与目录含义

### 6.1 输入层

- `config.yaml`：总参数
- `data/raw/cached/`：手工下载的数据库导出
- `data/raw/geo_datasets/`：GEO 数据

### 6.2 处理中间层

- `data/processed/`：各步骤结构化 CSV 结果

### 6.3 展示层

- `results/figures/`：论文图、展示图
- `results/tables/`：汇总表
- `results/docking/`：对接相关文件

### 6.4 管理层

- `logs/`：每步日志
- `docs/`：用户文档
- `env/`：环境定义

## 7. 可复现性设计

项目在可复现性上做了几件事：

- 关键阈值都集中在 `config.yaml`
- 支持记录查询日期
- 原始导出文件建议保存在 `data/raw/cached/`
- 各步骤输出结构相对固定
- Snakemake 明确了上下游依赖
- `00_preflight.py` 会识别 placeholder 和缺失依赖

## 8. 这个项目公开后别人如何使用

一个新用户的典型使用路径是：

1. 克隆仓库
2. 安装 Python 和 R 环境
3. 修改 `config.yaml`
4. 按文档补齐手工下载文件
5. 先运行 `00_preflight.py`
6. 再按单步或 Snakemake 跑流程
7. 到 `results/` 查看图表和汇总结果

## 9. 项目当前最适合作为哪类开源仓库

按现状看，这个仓库最适合公开成：

- 网络药理学分析模板仓库
- 教学型方法仓库
- 课题组内部标准流程仓库

如果以后要做成更强的开源项目，下一步建议是：

- 增加示例输入数据
- 增加测试数据集和 CI
- 增加 `CONTRIBUTING.md`
- 增加英文版方法说明
- 增加容器化环境，例如 Docker
