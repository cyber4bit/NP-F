# Quick Start · 快速开始

> 示例药物：血筒（大血藤，*Sargentodoxa cuneata*）
> 示例疾病：类风湿关节炎（Rheumatoid Arthritis, RA）

---

## 0. 推荐阅读顺序

如果你现在需要快速交付项目，建议先看 [`user-guide.md`](./user-guide.md)。

本文档更适合已经准备开始运行流程的用户；
如果你要逐条查看每个手工下载步骤，再去看 [`manual-inputs.md`](./manual-inputs.md)。

---

## 1. 环境安装

### 方式 A：Conda（推荐）

```bash
git clone <repo_url>
cd np-pipeline
conda env create -f env/environment.yml
conda activate np-pipeline
Rscript env/install_r_packages.R
```

### 方式 B：pip + 本地 R

```bash
pip install -r env/requirements.txt
Rscript env/install_r_packages.R
```

### 验证环境

```bash
python scripts/python/00_preflight.py --config config.yaml
```

所有 `[PASS]` = 就绪；`[WARN]` = 可选组件缺失；`[FAIL]` = 必须修复。

---

## 2. 修改配置

打开 `config.yaml`，修改药物和疾病：

```yaml
drug:
  name: "Xuetong"
  chinese_name: "血筒"
  latin_name: "Sargentodoxa cuneata"

disease:
  name: "rheumatoid arthritis"
  abbreviation: "RA"
  mesh_term: "Arthritis, Rheumatoid"
  geo_bulk_accession: "GSE93272"
```

---

## 3. 准备手工数据

**流程不是全自动的。** 以下步骤需要手工从数据库下载文件：

| 步骤 | 需要手工下载 | 放置路径 | 详见 |
|------|-------------|---------|------|
| Step 01 | TCMSP 成分表 | `data/raw/cached/tcmsp_YYYYMMDD.csv` | `docs/manual-inputs.md` |
| Step 02 | STP/PharmMapper/CTD/SEA | `data/raw/cached/*_manual.csv` | `docs/manual-inputs.md` |
| Step 03 | GeneCards/OMIM/DisGeNET | `data/raw/cached/*_manual.csv` | `docs/manual-inputs.md` |
| Step 05 | MCODE 聚类结果 | `data/raw/cached/mcode_clusters.csv` | `docs/manual-inputs.md` |
| Step 07 | eQTLGen + FinnGen GWAS | `data/raw/cached/eqtlgen_*.tsv` | `docs/manual-inputs.md` |
| Step 09 | GEO scRNA 数据 | `data/raw/geo_datasets/GSE*/` | `docs/manual-inputs.md` |

> **首次运行建议**：先只跑 Step 01，检查输出，再逐步补充手工数据。
> 如果你需要一份更像“交付文档”的版本，请优先使用 [`user-guide.md`](./user-guide.md)。

---

## 4. 运行流程

### 逐步运行（推荐新手）

```bash
# Step 01: 成分收集
python scripts/python/01_compounds.py --config config.yaml

# Step 02: 靶点预测（需手工数据）
python scripts/python/02_targets.py --config config.yaml

# Step 03: 疾病靶点（需手工数据）
python scripts/python/03_disease.py --config config.yaml

# Step 04: 韦恩图 + PPI
python scripts/python/04_ppi.py --config config.yaml

# Step 05: Hub 基因筛选
python scripts/python/05_hub_genes.py --config config.yaml

# Step 06: 富集分析（需 R）
python scripts/python/06_enrichment.py --config config.yaml

# Step 07: 孟德尔随机化（需 R + 数据）
Rscript scripts/R/07_mr.R config.yaml

# Step 08: 分子对接
python scripts/python/08_docking.py --config config.yaml

# Step 10: 汇总出图
python scripts/python/10_visualization.py --config config.yaml
```

### Snakemake 全流程

```bash
snakemake --cores 8
```

---

## 5. 输出目录说明

```
data/processed/          中间数据（CSV）
results/figures/         图表（PDF + PNG）
results/tables/          结果表格
results/docking/         对接文件（PDB/PDBQT）
logs/                    运行日志
data/raw/cached/         数据库导出缓存（带日期戳）
```

---

## 6. 如何识别 Placeholder 结果

**重要**：缺少手工数据时，流程会生成 placeholder 继续执行。这些**不是真实结论**。

识别方法：
- CSV 中出现 `placeholder`、`N/A`、`UNKNOWN` 字样
- `method` 列值为 `"placeholder"`
- `query_date` 包含 `YYYY-MM-DD`
- 数据行数异常少（< 5 行成分、< 10 个靶点）

运行 preflight 自动检查：

```bash
python scripts/python/00_preflight.py --config config.yaml
```

---

## 7. 最小成功标志

完成一次**真实**运行的最低要求：

- [ ] `01_compounds_filtered.csv` 有 >= 5 个真实成分（有 SMILES）
- [ ] `02_targets_merged.csv` 有 >= 10 个靶点（来自 >= 2 个数据库）
- [ ] `03_disease_targets.csv` 有 >= 50 个疾病靶点
- [ ] `04_intersection_genes.csv` 有 >= 5 个交集基因
- [ ] `05_hub_genes.csv` 有 >= 3 个 hub 基因
- [ ] `06_enrichment_kegg.csv` 有显著通路（p < 0.05）
- [ ] 所有 `query_date` 已填写真实日期
- [ ] `data/raw/cached/` 中有对应日期戳的原始导出文件
## Optional Step 09b

Place an MSigDB C7 immune GMT file at the path in `config.yaml -> immune.gene_sets_gmt`,
then run:

```bash
python scripts/python/09b_immune.py --config config.yaml
```
