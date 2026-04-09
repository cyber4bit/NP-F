# 手工输入数据指南 · Manual Inputs Guide

> 示例药物：血筒（大血藤）| 示例疾病：类风湿关节炎（RA）

本文档列出所有需要手工从外部数据库下载的文件。每个资源使用统一模板描述。

---

## 统一模板说明

每个手工输入都包含以下字段：

| 字段 | 说明 |
|------|------|
| 目标 | 这一步想拿到什么 |
| 来源 | 数据库/网站名 |
| 前置条件 | 账号、API key、网络要求 |
| 操作步骤 | 逐步编号 |
| 导出格式 | CSV/TSV/XLSX |
| 文件名 | 精确约定 |
| 存放路径 | 相对于项目根目录 |
| 必填列 | 下游脚本依赖的列名 |
| 自动校验 | 运行什么脚本检查 |
| 成功标志 | 什么说明下载成功 |
| 常见问题 | 编码、列名、空值等 |

---

## Step 01 — TCMSP 成分表

| 字段 | 内容 |
|------|------|
| **目标** | 获取血筒（大血藤）全部化学成分及 OB/DL 值 |
| **来源** | TCMSP (https://www.tcmsp-e.com/) |
| **前置条件** | 无需注册，需要稳定网络 |
| **操作步骤** | 1. 打开 https://www.tcmsp-e.com/ |
| | 2. 搜索框输入 `大血藤` 或 `Sargentodoxa` |
| | 3. 点击药材名进入详情页 |
| | 4. 点击 "Ingredients" 标签 |
| | 5. 点击 "Download" 或手动复制表格到 Excel |
| | 6. 保存为 CSV（UTF-8 编码） |
| **导出格式** | CSV |
| **文件名** | `tcmsp_20260409.csv`（替换为实际日期） |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `Molecule Name`, `OB(%)`, `DL`, `MW` (如有 SMILES 更好) |
| **自动校验** | `python scripts/python/00_preflight.py --step 1` |
| **成功标志** | CSV 有 >= 5 行，OB 和 DL 列有数值 |
| **常见问题** | |
| | - TCMSP 表头可能是中文，需手动改为英文 |
| | - 部分成分无 SMILES，后续 STP 预测会跳过 |
| | - 如果网站打不开，可用 TCMID 或 ETCM 替代 |

**运行命令：**
```bash
python scripts/python/01_compounds.py --config config.yaml --tcmsp_export data/raw/cached/tcmsp_20260409.csv
```

---

## Step 02 — SwissTargetPrediction（STP）

| 字段 | 内容 |
|------|------|
| **目标** | 基于成分 SMILES 预测靶点 |
| **来源** | SwissTargetPrediction (http://www.swisstargetprediction.ch/) |
| **前置条件** | 无需注册；每次只能提交 1 个 SMILES |
| **操作步骤** | 1. 打开 http://www.swisstargetprediction.ch/ |
| | 2. 粘贴成分的 Canonical SMILES |
| | 3. 选择 "Homo sapiens" |
| | 4. 点击 "Predict" |
| | 5. 等待结果页面加载 |
| | 6. 点击 "Download CSV" |
| | 7. 对每个活性成分重复 2-6 步 |
| | 8. 合并所有 CSV 为一个文件 |
| **导出格式** | CSV |
| **文件名** | `stp_manual.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `Target`(或 `Gene name`), `Uniprot ID`, `Probability` |
| **成功标志** | 有靶点且 Probability >= 0.5 的行 >= 10 |
| **常见问题** | |
| | - 每个 SMILES 提交后需等待 ~30s |
| | - 批量操作建议间隔 5 秒避免被限速 |
| | - 脚本也会尝试自动查询，但 STP 网页解析不稳定 |

---

## Step 02 — PharmMapper

| 字段 | 内容 |
|------|------|
| **目标** | 基于 3D 药效团预测靶点 |
| **来源** | PharmMapper (https://www.lilab-ecust.cn/pharmmapper/) |
| **前置条件** | 无需注册，需要 MOL2 文件 |
| **操作步骤** | 1. 将成分 SMILES 转为 MOL2（Open Babel: `obabel -:"SMILES" --gen3D -O compound.mol2`） |
| | 2. 打开 https://www.lilab-ecust.cn/pharmmapper/ |
| | 3. 上传 MOL2 文件 |
| | 4. 选择 "Human Protein Targets Only" |
| | 5. 填写邮箱，提交任务 |
| | 6. 等待邮件通知（约 2 小时） |
| | 7. 从邮件链接下载结果 CSV |
| **导出格式** | CSV |
| **文件名** | `pharmmapper_manual.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `Target Name`(或含 gene/symbol), `PDB ID`, `Fit Score` |
| **成功标志** | 有 >= 20 个靶点 |
| **常见问题** | |
| | - 网站较慢，高峰期可能 > 4 小时 |
| | - 结果靶点名可能是蛋白全名而非 gene symbol，需映射 |

---

## Step 02 — CTD

| 字段 | 内容 |
|------|------|
| **目标** | 获取化合物-基因交互关系 |
| **来源** | CTD (https://ctdbase.org/) |
| **前置条件** | 无需注册 |
| **操作步骤** | 1. 打开 https://ctdbase.org/ |
| | 2. 选择 "Chemical" 搜索 |
| | 3. 输入主要活性成分英文名（如 `sargentin`） |
| | 4. 进入化合物页面 → "Gene Interactions" 标签 |
| | 5. Filter: Organism = Homo sapiens |
| | 6. 点击 "Download" 导出 CSV |
| **导出格式** | CSV/TSV |
| **文件名** | `ctd_manual.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `Gene Symbol`, `Interaction` |
| **常见问题** | |
| | - 中药成分在 CTD 中可能无记录，这是正常的 |
| | - 如果找不到特定成分，跳过此数据库即可 |

---

## Step 02 — SEA

| 字段 | 内容 |
|------|------|
| **目标** | 基于化学相似性预测靶点 |
| **来源** | SEA (https://sea.bkslab.org/) |
| **前置条件** | 需要注册账号 |
| **操作步骤** | 1. 打开 https://sea.bkslab.org/ |
| | 2. 注册/登录 |
| | 3. 上传 SMILES 列表或 SDF 文件 |
| | 4. 提交预测 |
| | 5. 下载结果 |
| **导出格式** | CSV |
| **文件名** | `sea_manual.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `Target Name` 或 `Gene Symbol` |

---

## Step 03 — GeneCards 疾病靶点

| 字段 | 内容 |
|------|------|
| **目标** | 获取 RA 相关基因列表 |
| **来源** | GeneCards (https://www.genecards.org/) |
| **前置条件** | 批量下载需要注册（免费） |
| **操作步骤** | 1. 打开 https://www.genecards.org/ |
| | 2. 搜索: `rheumatoid arthritis` |
| | 3. 结果页面显示基因列表和 Relevance score |
| | 4. 过滤: Relevance score >= 1.0（config 中设置） |
| | 5. 点击 "Export" → CSV |
| | 6. 如果不能导出，手动复制前 500 行到 Excel 保存 |
| **导出格式** | CSV |
| **文件名** | `genecards_manual.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `Gene Symbol`（或 `Symbol`）, `Relevance score` |
| **自动校验** | `python scripts/python/00_preflight.py --step 3` |
| **成功标志** | >= 200 个基因（RA 通常有 500+） |
| **常见问题** | |
| | - 免费版有导出数量限制（~500行），付费版无限制 |
| | - 列名可能是 "Gene Symbol" 或 "Symbol"，脚本都能识别 |
| | - Relevance score 阈值可在 config.yaml 调整 |

---

## Step 03 — OMIM 疾病靶点

| 字段 | 内容 |
|------|------|
| **目标** | 获取 RA 的 OMIM 关联基因 |
| **来源** | OMIM (https://www.omim.org/) |
| **前置条件** | 批量下载需注册申请 API key |
| **操作步骤** | 1. 打开 https://www.omim.org/ |
| | 2. 搜索: `rheumatoid arthritis` |
| | 3. 筛选 "Gene" 和 "Phenotype" 类型的条目 |
| | 4. 导出基因列表 |
| | 5. 或使用 OMIM API（需申请 key） |
| **导出格式** | CSV |
| **文件名** | `omim_manual.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `Gene Symbol`（或含 gene/symbol 的列） |
| **常见问题** | |
| | - OMIM 基因数通常较少（RA ~50-100 个） |
| | - 与 GeneCards 合并后去重 |

---

## Step 03 — DisGeNET 疾病靶点

| 字段 | 内容 |
|------|------|
| **目标** | 获取 RA 的基因关联（GDA score） |
| **来源** | DisGeNET (https://www.disgenet.org/) |
| **前置条件** | 注册免费账号获取 API key |
| **操作步骤** | 方式 A（API）: |
| | 1. 注册: https://www.disgenet.org/signup/ |
| | 2. 获取 API key |
| | 3. 设置环境变量: `export DISGENET_API_KEY=<your_key>` |
| | 4. 脚本自动查询 |
| | 方式 B（手动）: |
| | 1. 登录 DisGeNET 网站 |
| | 2. 搜索: `rheumatoid arthritis` 或 MeSH ID `D001172` |
| | 3. 下载 Gene-Disease Associations |
| | 4. 保存 CSV |
| **导出格式** | CSV/TSV |
| **文件名** | `disgenet_manual.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `geneSymbol`（或 `gene_symbol`）, `score` |
| **成功标志** | >= 100 个基因（RA 通常有 300+） |
| **常见问题** | |
| | - 免费 API 有每日请求限制 |
| | - DisGeNET score >= 0.1 过滤可在 config 调整 |

---

## Step 05 — MCODE 聚类结果

| 字段 | 内容 |
|------|------|
| **目标** | 从 PPI 网络中提取蛋白质功能模块 |
| **来源** | Cytoscape + MCODE 插件 |
| **前置条件** | 安装 Cytoscape 3.8+ 和 MCODE 插件 |
| **操作步骤** | 1. 打开 Cytoscape |
| | 2. File → Import → Network from File → `data/processed/04_ppi_network.sif` |
| | 3. Apps → MCODE → Analyze Current Network |
| | 4. 参数设置（见 config.yaml `hub_selection.mcode`）: |
| |    - Degree cutoff: 2 |
| |    - Score cutoff: 0.2 |
| |    - K-core: 2 |
| |    - Max depth: 100 |
| | 5. 结果面板显示 Cluster 列表 |
| | 6. 选择 Top 3 clusters |
| | 7. 导出节点列表为 CSV |
| **导出格式** | CSV |
| **文件名** | `mcode_clusters.csv` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `gene_symbol`（或 `Node`, `Name`） |
| **成功标志** | >= 5 个基因 |
| **常见问题** | |
| | - 如果 PPI 网络太小（< 10 节点），MCODE 可能无结果 |
| | - 此步骤可选，跳过时 Step 05 会用拓扑学方法替代 |

---

## Step 07 — eQTLGen cis-eQTL 数据

| 字段 | 内容 |
|------|------|
| **目标** | 获取 hub 基因的血液 cis-eQTL 作为 MR 工具变量 |
| **来源** | eQTLGen (https://eqtlgen.org/) |
| **前置条件** | 无需注册，文件较大（~1.5 GB） |
| **操作步骤** | 1. 打开 https://eqtlgen.org/cis-eqtls.html |
| | 2. 下载 "Full cis-eQTL summary statistics" |
| | 3. 解压 .gz 文件 |
| | 4. 重命名为 `eqtlgen_cis_20260409.tsv` |
| **导出格式** | TSV |
| **文件名** | `eqtlgen_cis_YYYYMMDD.tsv`（替换日期） |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `SNP`, `beta`, `se`, `pval`, `effect_allele`, `other_allele`, `eaf`, `gene_symbol`, `N` |
| **常见问题** | |
| | - 文件非常大，建议预筛选只保留 hub 基因相关行 |
| | - 筛选命令: `awk -F'\t' 'NR==1 || $5=="TNF" || $5=="IL6"' full_file.tsv > filtered.tsv` |

---

## Step 07 — FinnGen GWAS 数据

| 字段 | 内容 |
|------|------|
| **目标** | 获取 RA 的 GWAS 汇总统计量作为 MR 结局 |
| **来源** | FinnGen (https://www.finngen.fi/en/access_results) |
| **前置条件** | 无需注册 |
| **操作步骤** | 1. 打开 https://www.finngen.fi/en/access_results |
| | 2. 搜索表型: `rheumatoid arthritis` |
| | 3. 找到 `RHEUMA_SEROPOS_RA`（血清阳性 RA）或 `M13_RHEUMA` |
| | 4. 下载 Summary statistics (.gz) |
| | 5. 保存到 cached 目录 |
| **导出格式** | TSV.GZ |
| **文件名** | `finngen_r10_RA_20260409.gz` |
| **存放路径** | `data/raw/cached/` |
| **必填列** | `rsid`, `beta`, `sebeta`, `pval`, `alt`, `ref`, `af_alt` |
| **常见问题** | |
| | - FinnGen R10 数据约 2-5 GB，需足够磁盘空间 |
| | - 也可用 IEU OpenGWAS 作为替代来源 |

---

## Step 09 — GEO scRNA-seq 数据

| 字段 | 内容 |
|------|------|
| **目标** | 获取 RA 滑膜组织单细胞数据 |
| **来源** | NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/) |
| **前置条件** | 无需注册，需要大量磁盘空间 |
| **操作步骤** | 1. 打开 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159117 |
| | 2. 找到 "Supplementary file" 部分 |
| | 3. 下载 10x matrix 文件: |
| |    - `barcodes.tsv.gz` |
| |    - `features.tsv.gz`（或 `genes.tsv.gz`） |
| |    - `matrix.mtx.gz` |
| | 4. 如提供 RDS/h5ad，下载预处理对象更方便 |
| | 5. 创建目录并放入文件 |
| **导出格式** | 10x Genomics matrix 或 RDS |
| **文件名** | 各 sample 子目录 |
| **存放路径** | `data/raw/geo_datasets/GSE159117/` |
| **成功标志** | 目录下有 3 个 .gz 文件或 1 个 .rds 文件 |
| **常见问题** | |
| | - scRNA 数据通常 > 1 GB |
| | - config.yaml 中 `scrna.enabled` 需设为 `true` |
| | - 如果只有 h5 格式，Seurat `Read10X_h5()` 可读取 |

---

## 数据完整性速查表

下载完成后，运行 preflight 检查所有手工文件：

```bash
python scripts/python/00_preflight.py --config config.yaml
```

| 文件 | 步骤 | 必需 | 替代方案 |
|------|------|------|---------|
| `tcmsp_YYYYMMDD.csv` | 01 | 推荐 | PubChem 自动搜索（成分少） |
| `stp_manual.csv` | 02 | 至少1个 | 脚本自动查询（不稳定） |
| `pharmmapper_manual.csv` | 02 | 可选 | — |
| `ctd_manual.csv` | 02 | 可选 | — |
| `sea_manual.csv` | 02 | 可选 | — |
| `genecards_manual.csv` | 03 | 推荐 | — |
| `omim_manual.csv` | 03 | 可选 | — |
| `disgenet_manual.csv` | 03 | 推荐 | API 自动查询 |
| `mcode_clusters.csv` | 05 | 可选 | 拓扑学方法替代 |
| `eqtlgen_cis_*.tsv` | 07 | MR必需 | 无（跳过MR） |
| `finngen_r10_*.gz` | 07 | MR必需 | IEU OpenGWAS |
| `geo_datasets/GSE*/` | 09 | scRNA必需 | 无（跳过scRNA） |
---

## Step 09b - MSigDB C7 immune gene sets

Optional immune infiltration analysis expects a GMT file at the path configured in
`config.yaml -> immune.gene_sets_gmt` (default: `data/raw/cached/msigdb_c7_immune.gmt`).

Related config fields:

- `immune.gene_sets_gmt`
- `immune.correlation_abs_r_min`
- `immune.p_cutoff`
