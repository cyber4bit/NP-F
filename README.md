# Network Pharmacology Pipeline · 网络药理学通用分析流程

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://python.org)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-blue.svg)](https://r-project.org)
[![Snakemake](https://img.shields.io/badge/workflow-Snakemake-green.svg)](https://snakemake.readthedocs.io)

---

## Overview · 项目简介

一个可复现的模块化**单一药物网络药理学**分析流程，设计为通用模板。

- **示例药物：** 血筒（大血藤，*Sargentodoxa cuneata*）
- **示例疾病：** 类风湿关节炎（Rheumatoid Arthritis, RA）

### 方法论来源

| Module | Reference | Key innovation |
|--------|-----------|----------------|
| scRNA-seq + SPR/CETSA | Xie et al., *Phytomedicine* 2025 | Direct target binding validation |
| MR + public DB chain | Zhao et al., *Pharmacogenomics Pers Med* 2025 | Causal gene prioritization |
| bulk + scRNA + Spatial | Zhu & He, *Front Pharmacol* 2026 | Multi-omics + spatial context |

### 重要提示

> **本流程是「半自动研究编排」，不是全自动分析系统。**
> 多个步骤需要手工从数据库下载数据。缺少数据时流程会生成 placeholder 继续执行——
> **placeholder 结果不可用于论文。** 详见 [docs/research-mode.md](docs/research-mode.md)。

---

## Workflow · 分析流程

```
00_preflight.py      → 环境与数据预检
01_compounds.py      → 成分收集 & ADME 筛选（TCMSP / PubChem）
02_targets.py        → 多数据库靶点预测（TCMSP/STP/PharmMapper/CTD/SEA）
03_disease.py        → 疾病靶点收集（GeneCards/OMIM/DisGeNET）
04_ppi.py            → 韦恩交集 + PPI 网络（STRING API）
05_hub_genes.py      → 三层筛选: 拓扑学 + MCODE + LASSO/RF
06_enrichment.py     → GO / KEGG / GSEA / GSVA 富集（R clusterProfiler）
07_mr.R              → 孟德尔随机化（TwoSampleMR, eQTLGen + FinnGen）
08_docking.py        → 分子对接（AutoDock Vina）+ MD guide（GROMACS）
09_scrna.R           → scRNA-seq / 空间转录组（Seurat + Harmony，可选）
10_visualization.py  → 出版级汇总图
```

---

## Repository Structure · 目录结构

```
np-pipeline/
├── config.yaml                  # 全局参数配置 — 首先修改这里
├── Snakefile                    # 工作流管理器
├── README.md
│
├── scripts/
│   ├── python/
│   │   ├── 00_preflight.py      # 预检脚本
│   │   ├── 01_compounds.py      # Step 01
│   │   ├── 02_targets.py        # Step 02
│   │   ├── 03_disease.py        # Step 03
│   │   ├── 04_ppi.py            # Step 04
│   │   ├── 05_hub_genes.py      # Step 05
│   │   ├── 06_enrichment.py     # Step 06 (calls R)
│   │   ├── 08_docking.py        # Step 08
│   │   ├── 10_visualization.py  # Step 10
│   │   └── utils_validate.py    # 数据验证工具
│   └── R/
│       ├── 07_mr.R              # Step 07
│       └── 09_scrna.R           # Step 09
│
├── data/
│   ├── raw/
│   │   ├── cached/              # 数据库导出缓存（带日期戳）
│   │   └── geo_datasets/        # GEO 下载数据
│   └── processed/               # 中间处理结果
│
├── results/
│   ├── figures/                 # PDF + PNG 图表
│   ├── tables/                  # CSV 结果表
│   └── docking/                 # 对接文件
│
├── docs/
│   ├── user-guide.md           # 面向用户的操作指南
│   ├── quickstart.md            # 快速开始
│   ├── manual-inputs.md         # 手工输入数据详细指南
│   ├── step-by-step-sop.md      # 逐步 SOP
│   └── research-mode.md         # 哪些结果可发表
│
├── env/
│   ├── environment.yml          # conda 环境
│   ├── requirements.txt         # pip 依赖
│   └── install_r_packages.R     # R 包安装脚本
│
└── logs/                        # 运行日志
```

---

## Quick Start · 快速开始

```bash
# 1. 安装环境
conda env create -f env/environment.yml
conda activate np-pipeline
Rscript env/install_r_packages.R

# 2. 预检
python scripts/python/00_preflight.py --config config.yaml

# 3. 修改 config.yaml（设置药物/疾病）

# 4. 按用户手册准备手工数据（见 docs/user-guide.md）

# 5. 运行
snakemake --cores 8
# 或逐步运行（见 docs/step-by-step-sop.md）
```

详细指南：[docs/user-guide.md](docs/user-guide.md)

---

## 手工数据要求

本流程**不是全自动的**。以下步骤需要手工下载外部数据：

| 步骤 | 数据 | 来源 | 详见 |
|------|------|------|------|
| 01 | TCMSP 成分表 | tcmsp-e.com | [manual-inputs.md](docs/manual-inputs.md) |
| 02 | STP/PharmMapper/CTD/SEA 靶点 | 各网站 | [manual-inputs.md](docs/manual-inputs.md) |
| 03 | GeneCards/OMIM/DisGeNET 疾病靶点 | 各网站 | [manual-inputs.md](docs/manual-inputs.md) |
| 05 | MCODE 聚类（Cytoscape） | 本地软件 | [manual-inputs.md](docs/manual-inputs.md) |
| 07 | eQTLGen + FinnGen GWAS | eqtlgen.org / finngen.fi | [manual-inputs.md](docs/manual-inputs.md) |
| 09 | GEO scRNA 数据 | ncbi.nlm.nih.gov/geo | [manual-inputs.md](docs/manual-inputs.md) |

> 自动化分工原则：能稳定 API 调用的（PubChem, STRING, DisGeNET）→ 脚本自动；
> 不稳定或需人工操作的（TCMSP, GeneCards, Cytoscape）→ **手工下载 + 自动校验 + 自动导入**。

---

## 文档导航

| 文档 | 面向人群 | 内容 |
|------|---------|------|
| [docs/user-guide.md](docs/user-guide.md) | 交付/演示/首次使用 | **面向用户的操作指南，优先看这个** |
| [docs/quickstart.md](docs/quickstart.md) | 首次使用 | 环境、配置、最小运行 |
| [docs/manual-inputs.md](docs/manual-inputs.md) | 所有用户 | **每个手工数据的下载详细步骤** |
| [docs/step-by-step-sop.md](docs/step-by-step-sop.md) | 逐步操作 | 每步的输入/输出/判据/报错 |
| [docs/research-mode.md](docs/research-mode.md) | 写论文前 | placeholder 识别、可发表判据 |

---

## Reproducibility · 可复现性清单

- [ ] 所有数据库查询日期记录在 `config.yaml → query_date`
- [ ] 原始导出文件缓存在 `data/raw/cached/`（带日期戳）
- [ ] 随机种子全局设置（`seed_global: 2024`）
- [ ] 软件版本记录在 `config.yaml → software_versions`
- [ ] GEO accession 编号已记录
- [ ] 筛选阈值（OB, DL, confidence）在 config.yaml 中可追溯

---

## Citation · 引用

If you use this pipeline, please cite the three reference papers:

```bibtex
@article{xie2025cryptotanshinone,
  title={Cryptotanshinone alleviates immunosuppression in endometriosis...},
  author={Xie, Linling and others},
  journal={Phytomedicine}, volume={136}, pages={156227}, year={2025}
}

@article{zhao2025berberine,
  title={Berberine Repairs Intestinal Mucosal Barrier by Targeting HSP90AA1 and MAPK14},
  author={Zhao, Danya and others},
  journal={Pharmacogenomics and Personalized Medicine}, volume={18}, pages={263--278}, year={2025}
}

@article{zhu2026spironolactone,
  title={A network pharmacology-guided multi-omics and spatial single-cell framework...},
  author={Zhu, Zhen and He, Xingjun},
  journal={Frontiers in Pharmacology}, volume={17}, pages={1770261}, year={2026}
}
```

---

## License
MIT License — see [LICENSE](LICENSE)
