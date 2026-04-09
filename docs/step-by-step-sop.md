# 标准操作规程 · Step-by-Step SOP

> 示例药物：血筒（大血藤，*Sargentodoxa cuneata*）
> 示例疾病：类风湿关节炎（RA）

每个步骤使用统一模板：**输入 → 执行 → 预期输出 → 成功判据 → 常见报错 → 下一步检查项**

---

## Step 00 — Preflight 预检

| 项目 | 内容 |
|------|------|
| **输入** | `config.yaml` |
| **执行** | `python scripts/python/00_preflight.py --config config.yaml` |
| **预期输出** | 终端打印各项 `[PASS]`/`[WARN]`/`[FAIL]` |
| **成功判据** | 所有 `[FAIL]` 项已修复；`[WARN]` 了解即可 |
| **常见报错** | |
| | `[FAIL] Rscript not found` → 安装 R 并加入 PATH |
| | `[FAIL] pandas` → `pip install pandas` |
| **下一步** | 确认 config.yaml 中 drug/disease 设置正确 → Step 01 |

---

## Step 01 — 成分收集与 ADME 筛选

| 项目 | 内容 |
|------|------|
| **输入** | `config.yaml`（drug.name）；可选：TCMSP 手工导出 CSV |
| **执行** | |

```bash
# 方式A：有 TCMSP 导出（推荐）
python scripts/python/01_compounds.py --config config.yaml \
  --tcmsp_export data/raw/cached/tcmsp_20260409.csv

# 方式B：无 TCMSP 导出（PubChem 回退，成分少）
python scripts/python/01_compounds.py --config config.yaml
```

| 项目 | 内容 |
|------|------|
| **预期输出** | `data/processed/01_compounds_filtered.csv` |
| | `data/raw/cached/tcmsp_20260409.csv`（缓存） |
| **成功判据** | CSV 有 >= 5 行；含 `Molecule Name`, `OB`, `DL` 列；OB/DL 有数值 |
| **Placeholder 识别** | 如果用方式 B，OB/DL 全为 NaN → 这是 placeholder，需补 TCMSP |
| **常见报错** | |
| | `No PubChem results` → 药名拼写错误或 PubChem 无此药 |
| | CSV 只有 1-2 行 → PubChem 只返回单体化合物，非中药全方 |
| **下一步检查** | 打开 CSV，确认成分名称合理 → Step 02 |

**血筒示例预期**：TCMSP 导出后经 OB>=30%, DL>=0.18 过滤，约 10-30 个活性成分。

---

## Step 02 — 多数据库靶点预测

| 项目 | 内容 |
|------|------|
| **输入** | `data/processed/01_compounds_filtered.csv` |
| | 手工文件（可选）：`stp_manual.csv`, `pharmmapper_manual.csv`, `ctd_manual.csv`, `sea_manual.csv` |
| **执行** | |

```bash
# 基础运行（自动查 STP + 加载已有手工文件）
python scripts/python/02_targets.py --config config.yaml

# 指定 TCMSP 靶点文件
python scripts/python/02_targets.py --config config.yaml \
  --tcmsp_targets data/raw/cached/tcmsp_targets_manual.csv

# 指定 STP 手工导出
python scripts/python/02_targets.py --config config.yaml \
  --stp_csv data/raw/cached/stp_manual.csv
```

| 项目 | 内容 |
|------|------|
| **预期输出** | `data/processed/02_targets_merged.csv` |
| | `data/raw/cached/02_targets_all_YYYYMMDD.csv`（完整记录） |
| **成功判据** | >= 10 个唯一 `gene_symbol`；`source` 列有 >= 2 个来源 |
| **Placeholder 识别** | 如果只有 STP 自动查询且全部失败 → 0 靶点，需手工补数据 |
| **常见报错** | |
| | `SwissTargetPrediction query failed` → 网站限速，用手工导出 |
| | `No gene_symbol column` → 手工 CSV 列名不对，参见 manual-inputs.md |
| | `mygene not installed` → `pip install mygene` |
| **下一步检查** | 确认 gene_symbol 看起来合理（都是人类基因名） → Step 03 |

**血筒示例预期**：TCMSP + STP 合并后约 100-300 个靶点。

---

## Step 03 — 疾病靶点收集

| 项目 | 内容 |
|------|------|
| **输入** | `config.yaml`（disease.name, mesh_term） |
| | 手工文件：`genecards_manual.csv`, `omim_manual.csv`, `disgenet_manual.csv` |
| **执行** | |

```bash
# 有 DisGeNET API key
export DISGENET_API_KEY=your_key_here
python scripts/python/03_disease.py --config config.yaml

# 无 API key（只用手工文件）
python scripts/python/03_disease.py --config config.yaml
```

| 项目 | 内容 |
|------|------|
| **预期输出** | `data/processed/03_disease_targets.csv` |
| **成功判据** | >= 50 个唯一 gene_symbol；source 列有 >= 2 个来源 |
| **常见报错** | |
| | `No disease targets collected` → 没有手工文件也没有 API key |
| | `No 'gene_symbol' column` → CSV 列名不匹配 |
| **下一步检查** | gene_symbol 数量合理 → Step 04 |

**RA 示例预期**：GeneCards(~500) + OMIM(~80) + DisGeNET(~300) 去重后约 800-1500 个靶点。

---

## Step 04 — 韦恩交集 + PPI 网络

| 项目 | 内容 |
|------|------|
| **输入** | `data/processed/02_targets_merged.csv`, `data/processed/03_disease_targets.csv` |
| **执行** | `python scripts/python/04_ppi.py --config config.yaml` |
| **预期输出** | `data/processed/04_intersection_genes.csv` |
| | `data/processed/04_ppi_edges.csv` |
| | `data/processed/04_ppi_network.sif`（Cytoscape 格式） |
| | `results/figures/04_venn.pdf` |
| **成功判据** | 交集 >= 5 个基因；PPI edges >= 10 条 |
| **常见报错** | |
| | `Empty intersection` → 药物靶点与疾病靶点无重叠，检查 Step 02/03 |
| | `STRING API error` → 网络问题或基因名不对 |
| **下一步检查** | |
| | 1. 查看 venn.pdf 比例是否合理 |
| | 2. 打开 `04_ppi_network.sif` 在 Cytoscape 中可视化 |
| | 3. 如需 MCODE 分析 → 在 Cytoscape 中运行（见 manual-inputs.md） |

**血筒 × RA 示例预期**：交集约 20-80 个基因，PPI edges 约 50-500 条。

---

## Step 05 — Hub 基因筛选（三层筛选）

| 项目 | 内容 |
|------|------|
| **输入** | `data/processed/04_intersection_genes.csv`, `04_ppi_edges.csv` |
| | 可选：`data/raw/cached/mcode_clusters.csv`（MCODE 结果） |
| | 可选：`data/raw/cached/geo_expression_matrix.csv`（GEO 表达矩阵） |
| **执行** | `python scripts/python/05_hub_genes.py --config config.yaml` |
| **预期输出** | `data/processed/05_hub_genes.csv` |
| **成功判据** | >= 3 个 hub 基因；degree, betweenness 列有数值 |
| **三层说明** | |
| | Layer 1（拓扑学）：degree + betweenness 前 20 → 始终执行 |
| | Layer 2（MCODE）：需手工运行 Cytoscape → 可选 |
| | Layer 3（ML）：需 GEO 表达矩阵 → 可选 |
| | 2-of-3 投票 → 最终 hub 列表；层数不足时自动降级 |
| **常见报错** | |
| | `Layer 3 skipped: no GEO expression matrix` → 正常，非必需 |
| **下一步检查** | hub 基因名是否在 RA 文献中有报道 → Step 06 |

**血筒 × RA 示例预期**：hub 基因约 5-15 个，如 TNF, IL6, AKT1, VEGFA, TP53 等。

---

## Step 06 — 富集分析（GO/KEGG/GSEA/GSVA）

| 项目 | 内容 |
|------|------|
| **输入** | `data/processed/05_hub_genes.csv`, `04_intersection_genes.csv` |
| **执行** | `python scripts/python/06_enrichment.py --config config.yaml` |
| **预期输出** | `results/figures/06_enrichment_kegg.pdf` |
| | `results/figures/06_enrichment_go.pdf` |
| | `data/processed/06_enrichment_kegg.csv` |
| | `data/processed/06_enrichment_go.csv` |
| **成功判据** | KEGG CSV 有 >= 5 个显著通路（p.adjust < 0.05） |
| **常见报错** | |
| | `Rscript not found` → R 未安装或不在 PATH |
| | `there is no package called 'clusterProfiler'` → 运行 `Rscript env/install_r_packages.R` |
| | `No significant KEGG pathways` → 基因数太少或无显著富集 |
| **下一步检查** | 查看 KEGG 通路是否与 RA 机制相关 → Step 07 |

**RA 示例预期**：应富集到 TNF signaling, IL-17 signaling, Th17 cell differentiation, NF-kB signaling 等 RA 经典通路。

---

## Step 07 — 孟德尔随机化（MR）

| 项目 | 内容 |
|------|------|
| **输入** | `data/processed/05_hub_genes.csv` |
| | `data/raw/cached/eqtlgen_cis_YYYYMMDD.tsv`（eQTL 工具变量） |
| | `data/raw/cached/finngen_r10_RA_YYYYMMDD.gz`（GWAS 结局） |
| **执行** | `Rscript scripts/R/07_mr.R config.yaml` |
| **预期输出** | `data/processed/07_mr_results.csv` |
| | `results/figures/07_mr_forest.pdf` |
| **成功判据** | CSV 中 `method` 列不是 `"placeholder"`；有 IVW 结果 |
| **Placeholder 识别** | `method = "placeholder"`, `b = NA` → eQTL 或 GWAS 数据缺失 |
| **常见报错** | |
| | `eQTL file not found` → 下载 eQTLGen 数据（见 manual-inputs.md） |
| | `FinnGen GWAS file not found` → 下载 FinnGen 数据 |
| | `Clumping failed` → 无 plink 或远程 API 不可用，可忽略 |
| | `No instruments found for GENE` → 该基因无显著 cis-eQTL |
| **下一步检查** | IVW p < 0.05 的基因 → 有因果推断证据 → Step 08 |

**重要提示**：MR 结果为 placeholder 时，**不可写入论文**。必须下载真实数据重跑。

---

## Step 08 — 分子对接

| 项目 | 内容 |
|------|------|
| **输入** | `data/processed/01_compounds_filtered.csv`, `05_hub_genes.csv` |
| **执行** | `python scripts/python/08_docking.py --config config.yaml` |
| **预期输出** | `data/processed/08_docking_scores.csv` |
| | `results/docking/` 目录下的 PDBQT 文件和 MD guide |
| **成功判据** | `status` 列有 `"docked"`；binding_energy <= -5.0 kcal/mol |
| **常见限制** | |
| | 受体 PDBQT 需手工制备（AutoDockTools4） |
| | 结合位点坐标 (0,0,0) 是占位值，**必须替换** |
| | box size 20x20x20 是通用值，需按蛋白调整 |
| | KNOWN_PDB 字典可能缺少你的 hub 基因 → 需手动添加 PDB ID |
| **常见报错** | |
| | `No PDB ID registered for GENE` → 在脚本中添加映射 |
| | `Receptor PDBQT not found` → 需用 ADT 制备受体 |
| | `AutoDock Vina not found` → 安装 Vina 并加入 PATH |
| **下一步检查** | 结合能 <= -5.0 的配对 → 有对接证据 → Step 10 |

**重要提示**：当前 docking 模块是**半自动框架**。真实对接需要：
1. 手动从 RCSB 选择合适的 PDB 结构（优选共结晶配体的高分辨率结构）
2. 用 AutoDockTools4 去水、加氢、定义结合口袋
3. 替换脚本中 `center = (0,0,0)` 为真实坐标
4. 对接后用 PyMOL 检查 pose 合理性

---

## Step 09 — scRNA-seq / 空间转录组（可选）

| 项目 | 内容 |
|------|------|
| **输入** | `data/processed/05_hub_genes.csv` |
| | `data/raw/geo_datasets/GSE159117/`（10x matrix 文件） |
| **执行** | `Rscript scripts/R/09_scrna.R config.yaml` |
| **前提** | `config.yaml` 中 `scrna.enabled: true` |
| **预期输出** | `results/figures/09_umap_clusters.pdf` |
| | `results/figures/09_hub_target_dotplot.pdf` |
| | `data/processed/09_cell_type_expression.csv` |
| **成功判据** | UMAP 有清晰分群；dotplot 显示 hub 基因在不同细胞类型中的表达 |
| **常见报错** | |
| | `GEO dataset not found` → 下载 10x 数据到指定目录 |
| | `scrna.enabled: false` → config 中打开开关 |
| **下一步检查** | hub 基因在哪些细胞类型高表达 → 生物学解释 → Step 10 |

---

## Step 10 — 汇总可视化

| 项目 | 内容 |
|------|------|
| **输入** | Steps 04-08 的所有输出文件 |
| **执行** | `python scripts/python/10_visualization.py --config config.yaml` |
| **预期输出** | `results/figures/10_summary_figure.pdf`（多页汇总） |
| | 各面板单独 PDF：`10_venn.pdf`, `10_ppi_network.pdf`, `10_hub_heatmap.pdf` 等 |
| **成功判据** | summary_figure.pdf 有 >= 3 页内容 |
| **常见报错** | |
| | `matplotlib-venn not installed` → `pip install matplotlib-venn` |
| | 部分面板为空 → 前置步骤输出缺失（正常，已有面板仍会生成） |

---

## 全流程速查命令

```bash
# 0. 环境检查
python scripts/python/00_preflight.py --config config.yaml

# 1-10. 逐步执行
python scripts/python/01_compounds.py --config config.yaml --tcmsp_export data/raw/cached/tcmsp_20260409.csv
python scripts/python/02_targets.py --config config.yaml
python scripts/python/03_disease.py --config config.yaml
python scripts/python/04_ppi.py --config config.yaml
python scripts/python/05_hub_genes.py --config config.yaml
python scripts/python/06_enrichment.py --config config.yaml
Rscript scripts/R/07_mr.R config.yaml
python scripts/python/08_docking.py --config config.yaml
Rscript scripts/R/09_scrna.R config.yaml
python scripts/python/10_visualization.py --config config.yaml

# 或 Snakemake 一键
snakemake --cores 8
```
## Optional Step 09b - Immune Infiltration

1. Ensure `data/raw/cached/geo_expression_matrix.csv` exists.
2. Place the MSigDB C7 immune GMT file at `config.yaml -> immune.gene_sets_gmt`.
3. Run `python scripts/python/09b_immune.py --config config.yaml`.
4. Review:
   - `data/processed/09b_immune_scores.csv`
   - `results/figures/09b_immune_heatmap.pdf`
