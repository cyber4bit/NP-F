# 用户操作指南

> 适用场景：需要先把项目交出去，当前先给出一份可直接给老师、同学或普通使用者阅读的操作手册。
> 本文档重点说明“用户需要做什么”，不展开算法细节。

---

## 1. 文档定位

这不是开发说明，也不是论文方法学说明，而是一份面向使用者的操作指南。

用户只需要知道四件事：

1. 先改哪里
2. 哪些数据必须手动下载
3. 下载后的文件放到哪里
4. 下一步运行什么命令

如果后续需要补充数据库截图、账号申请流程或论文级方法描述，可在此文档基础上继续完善。

---

## 2. 使用前准备

开始前请先完成以下准备：

- 已安装 Python、R，以及 `env/` 目录中的依赖
- 能正常运行 `python scripts/python/00_preflight.py --config config.yaml`
- 已打开并修改 [`config.yaml`](../config.yaml) 中的药物和疾病信息
- 已创建或确认以下目录存在：
  - `data/raw/cached/`
  - `data/raw/geo_datasets/`
  - `data/processed/`
  - `results/`

建议先完成一次环境预检：

```bash
python scripts/python/00_preflight.py --config config.yaml
```

---

## 3. 先交作业时的最低操作范围

如果当前目标是“先完成课程项目提交”，建议按下面的优先级执行：

| 优先级 | 模块 | 是否建议本轮完成 | 说明 |
|------|------|------|------|
| 必做 | Step 01 成分收集 | 是 | 至少准备 TCMSP 成分表 |
| 必做 | Step 02 药物靶点收集 | 是 | 至少准备 1 个可用来源，推荐 STP |
| 必做 | Step 03 疾病靶点收集 | 是 | 至少准备 GeneCards |
| 必做 | Step 04-06 自动分析 | 是 | 这是主流程主体 |
| 选做 | Step 05 MCODE | 否 | 没做可写明“暂未补充 Cytoscape 聚类” |
| 选做 | Step 07 MR | 否 | 缺 eQTL/GWAS 时可暂不作为正式结论 |
| 选做 | Step 08 对接 | 否 | 如未手工制备受体，只能视为流程占位 |
| 选做 | Step 09 scRNA | 否 | 可在后续版本补充 |

如果时间非常紧，最少应保证：

- `01_compounds_filtered.csv` 已生成
- `02_targets_merged.csv` 已生成
- `03_disease_targets.csv` 已生成
- `04_intersection_genes.csv` 已生成
- `05_hub_genes.csv` 已生成
- `06_enrichment_kegg.csv` 或对应图表已生成

这套结果足以说明项目主流程已经跑通。

---

## 4. 用户需要手动完成的操作

本项目不是全自动流程。以下步骤需要用户手动访问数据库、下载文件、并放到指定位置。

| 步骤 | 用户要做什么 | 文件名规范 | 放置位置 | 是否必需 |
|------|------|------|------|------|
| Step 01 | 从 TCMSP 下载药物成分表 | `tcmsp_YYYYMMDD.csv` | `data/raw/cached/` | 推荐 |
| Step 02 | 从 STP / PharmMapper / CTD / SEA 下载靶点结果 | `stp_manual.csv` 等 | `data/raw/cached/` | 推荐至少 1 个 |
| Step 03 | 从 GeneCards / OMIM / DisGeNET 下载疾病靶点 | `genecards_manual.csv` 等 | `data/raw/cached/` | 推荐至少 GeneCards |
| Step 05 | 在 Cytoscape 中运行 MCODE 并导出聚类结果 | `mcode_clusters.csv` | `data/raw/cached/` | 可选 |
| Step 07 | 下载 eQTLGen 和 FinnGen 数据 | `eqtlgen_cis_YYYYMMDD.tsv`、`finngen_r10_*.gz` | `data/raw/cached/` | MR 必需 |
| Step 09 | 下载 GEO 单细胞数据 | 以 `GSExxxxxx/` 目录保存 | `data/raw/geo_datasets/` | scRNA 必需 |

---

## 5. 各手动步骤的用户操作规范

### 5.1 Step 01：下载 TCMSP 成分表

目的：获取药物的候选化学成分及 OB、DL 等筛选指标。

用户操作：

1. 打开 `https://www.tcmsp-e.com/`
2. 搜索药材名称，例如“大血藤”或 `Sargentodoxa`
3. 进入对应药材页面
4. 打开 `Ingredients` 页面
5. 下载或复制成分表
6. 保存为 `CSV` 文件
7. 文件命名为 `tcmsp_YYYYMMDD.csv`
8. 放入 `data/raw/cached/`

完成标准：

- 文件可以正常用 Excel 或表格软件打开
- 含有成分名称
- 含有 `OB`、`DL` 等列
- 文件名带日期

运行命令：

```bash
python scripts/python/01_compounds.py --config config.yaml --tcmsp_export data/raw/cached/tcmsp_YYYYMMDD.csv
```

---

### 5.2 Step 02：下载药物靶点结果

目的：为药物成分补充潜在作用靶点。

建议最低做法：

- 至少准备 `SwissTargetPrediction` 一份结果
- 如果时间允许，再补 `PharmMapper`、`CTD`、`SEA`

用户操作：

1. 准备好 Step 01 得到的成分列表和 SMILES
2. 访问各数据库网站
3. 分别导出结果表
4. 将结果整理为 CSV
5. 保存到 `data/raw/cached/`

建议文件名：

- `stp_manual.csv`
- `pharmmapper_manual.csv`
- `ctd_manual.csv`
- `sea_manual.csv`

完成标准：

- 每个文件至少包含靶点名称或基因符号
- 文件能正常打开
- 同一数据库结果保存在单独文件中

运行命令：

```bash
python scripts/python/02_targets.py --config config.yaml
```

如需明确指定文件，可参考 [`step-by-step-sop.md`](./step-by-step-sop.md)。

---

### 5.3 Step 03：下载疾病靶点结果

目的：收集疾病相关基因，用于后续和药物靶点求交集。

建议最低做法：

- 至少准备 `GeneCards`
- 有条件时再补 `OMIM` 和 `DisGeNET`

用户操作：

1. 在 GeneCards 中搜索疾病名称，例如 `rheumatoid arthritis`
2. 导出疾病相关基因列表
3. 保存为 `genecards_manual.csv`
4. 如有时间，再分别下载 OMIM 和 DisGeNET 结果
5. 全部放入 `data/raw/cached/`

建议文件名：

- `genecards_manual.csv`
- `omim_manual.csv`
- `disgenet_manual.csv`

完成标准：

- GeneCards 文件中至少有基因符号列
- 疾病名称与 `config.yaml` 中设置一致
- 文件已放入指定目录

运行命令：

```bash
python scripts/python/03_disease.py --config config.yaml
```

---

### 5.4 Step 05：导出 MCODE 聚类结果

目的：补充 PPI 网络中的模块分析结果。

这一步是增强项，不是本轮交付的硬性前提。

用户操作：

1. 打开 Cytoscape
2. 导入 `data/processed/04_ppi_network.sif`
3. 运行 `MCODE`
4. 导出聚类节点结果
5. 保存为 `mcode_clusters.csv`
6. 放入 `data/raw/cached/`

运行命令：

```bash
python scripts/python/05_hub_genes.py --config config.yaml
```

如果暂时没有做这一步，项目仍可运行，但核心基因筛选证据会更弱。

---

### 5.5 Step 07：下载 MR 所需数据

目的：完成孟德尔随机化分析。

用户操作：

1. 从 eQTLGen 下载 cis-eQTL 数据
2. 从 FinnGen 下载对应疾病的 GWAS 数据
3. 按要求命名后放入 `data/raw/cached/`

建议文件名：

- `eqtlgen_cis_YYYYMMDD.tsv`
- `finngen_r10_RA_YYYYMMDD.gz`

运行命令：

```bash
Rscript scripts/R/07_mr.R config.yaml
```

注意：

- 如果没有真实 eQTL 或 GWAS 数据，MR 结果只能视为占位结果
- 占位结果可用于展示流程结构，不应写成正式结论

---

### 5.6 Step 09：下载 GEO 单细胞数据

目的：完成可选的单细胞分析模块。

用户操作：

1. 打开 GEO 页面
2. 下载目标数据集的原始矩阵文件
3. 以 `GSExxxxxx/` 文件夹形式保存
4. 放入 `data/raw/geo_datasets/`
5. 在 [`config.yaml`](../config.yaml) 中将 `scrna.enabled` 设为 `true`

运行命令：

```bash
Rscript scripts/R/09_scrna.R config.yaml
```

如果本轮不做单细胞分析，可以直接跳过。

---

## 6. 用户实际执行顺序

建议按下面顺序操作：

1. 修改 [`config.yaml`](../config.yaml)
2. 运行预检
3. 准备 Step 01 手工文件
4. 运行 Step 01
5. 准备 Step 02 和 Step 03 手工文件
6. 运行 Step 02 到 Step 06
7. 视时间决定是否补 Step 07 到 Step 09
8. 最后运行汇总出图

推荐命令顺序如下：

```bash
python scripts/python/00_preflight.py --config config.yaml
python scripts/python/01_compounds.py --config config.yaml --tcmsp_export data/raw/cached/tcmsp_YYYYMMDD.csv
python scripts/python/02_targets.py --config config.yaml
python scripts/python/03_disease.py --config config.yaml
python scripts/python/04_ppi.py --config config.yaml
python scripts/python/05_hub_genes.py --config config.yaml
python scripts/python/06_enrichment.py --config config.yaml
python scripts/python/10_visualization.py --config config.yaml
```

如需一键执行，可使用：

```bash
snakemake --cores 8
```

---

## 7. 交付时建议怎么表述

如果当前是课程项目、开题展示或阶段性检查，建议按以下口径表述：

- 本项目已完成网络药理学主流程的环境配置、参数配置和核心分析链路搭建
- 已支持成分筛选、药物靶点预测、疾病靶点收集、交集分析、PPI 网络、核心基因筛选和富集分析
- 部分依赖外部数据库或人工操作的模块已整理为标准化操作指南，后续可继续补充真实数据并完善结果
- MR、分子对接、scRNA 等增强模块可在后续版本继续补做

这样写更稳妥，也符合你当前“先提交、后完善”的目标。

---

## 8. 结果检查

完成后，建议重点检查以下文件是否已经生成：

- `data/processed/01_compounds_filtered.csv`
- `data/processed/02_targets_merged.csv`
- `data/processed/03_disease_targets.csv`
- `data/processed/04_intersection_genes.csv`
- `data/processed/05_hub_genes.csv`
- `data/processed/06_enrichment_kegg.csv`
- `results/figures/`

如果这些文件已经存在，说明主流程已经具备可展示性。

---

## 9. 深入文档

如果需要更细的操作说明，可继续查看：

- [`quickstart.md`](./quickstart.md)：快速开始
- [`manual-inputs.md`](./manual-inputs.md)：每个手工输入的详细下载步骤
- [`step-by-step-sop.md`](./step-by-step-sop.md)：每一步的输入、输出和报错说明
- [`research-mode.md`](./research-mode.md)：哪些结果可以用于论文
