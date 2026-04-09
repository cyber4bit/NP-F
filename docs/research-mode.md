# 研究模式指南 · Research Mode Guide

> 哪些结果能用于论文，哪些不能？

本文档帮助用户区分「模板跑通」和「真实科研结论」，避免将 placeholder 或不充分的结果写入发表论文。

---

## 核心原则

> **跑完流程 ≠ 得到结论。**
> 流程是工具，结论需要判断。

---

## 1. Placeholder 识别规则

### 什么是 Placeholder？

当关键输入数据缺失时，流程会生成占位数据以便后续步骤继续执行。这些占位数据**绝对不可**作为科研结论。

### 识别方法

| 标志 | 位置 | 含义 |
|------|------|------|
| `method = "placeholder"` | `07_mr_results.csv` | MR 未真实执行 |
| `b = NA, se = NA, pval = NA` | `07_mr_results.csv` | 无真实统计量 |
| `status = "receptor_not_prepared"` | `08_docking_scores.csv` | 对接未执行 |
| `binding_energy = None` | `08_docking_scores.csv` | 无对接结果 |
| `OB = NaN, DL = NaN` | `01_compounds_filtered.csv` | PubChem 回退，无 ADME |
| `query_date = "YYYY-MM-DD"` | 任何文件 | 查询日期未记录 |
| `source` 列只有 1 个来源 | `02/03_*.csv` | 靶点来源单一 |
| 行数 < 5 | `01_compounds_filtered.csv` | 成分过少 |
| 行数 < 10 | `02_targets_merged.csv` | 靶点过少 |

### 自动检测

```bash
python scripts/python/00_preflight.py --config config.yaml
```

输出中 `CONTAINS PLACEHOLDERS` 的文件不可用于论文。

---

## 2. 各步骤的「可发表」判据

### Step 01 — 成分收集

| 状态 | 判据 |
|------|------|
| **可发表** | 来自 TCMSP 的真实导出；OB/DL 有数值；>= 5 个成分 |
| **不可发表** | PubChem 回退（OB/DL 全 NaN）；无 TCMSP 缓存文件 |
| **论文写法** | "通过 TCMSP 数据库（访问日期：YYYY-MM-DD）检索大血藤活性成分，以 OB >= 30%, DL >= 0.18 为筛选标准，共获得 X 个活性成分。" |

### Step 02 — 靶点预测

| 状态 | 判据 |
|------|------|
| **可发表** | >= 2 个数据库来源；gene_symbol 已标准化；>= 10 个靶点 |
| **不可发表** | 只有 1 个来源且自动查询失败；靶点 < 10 |
| **论文写法** | "通过 TCMSP、SwissTargetPrediction、PharmMapper 等 X 个数据库预测药物靶点，合并去重后获得 Y 个潜在靶点。" |

### Step 03 — 疾病靶点

| 状态 | 判据 |
|------|------|
| **可发表** | >= 2 个数据库来源；>= 50 个靶点；GeneCards relevance 已过滤 |
| **不可发表** | 只有 1 个来源；靶点 < 50 |
| **论文写法** | "以 'rheumatoid arthritis' 为关键词，从 GeneCards（relevance >= 1.0）、OMIM、DisGeNET 检索疾病相关基因，合并去重后获得 Z 个疾病靶点。" |

### Step 04 — 交集 + PPI

| 状态 | 判据 |
|------|------|
| **可发表** | 交集 >= 5 个基因；STRING confidence 符合 config 设定 |
| **需注意** | STRING confidence 0.4 是 medium，0.9 是 highest — 论文需注明 |

### Step 05 — Hub 基因

| 状态 | 判据 |
|------|------|
| **可发表** | 使用 >= 2 层筛选（拓扑 + MCODE 或 ML）；hub 基因 >= 3 个 |
| **较弱** | 仅拓扑学 1 层（Layer 2/3 跳过）→ 可用但证据强度低 |
| **论文写法** | "通过 degree/betweenness centrality、MCODE 聚类、LASSO + Random Forest 三层筛选，取 2/3 投票交集，确定 X 个核心靶点。" |

### Step 06 — 富集分析

| 状态 | 判据 |
|------|------|
| **可发表** | p.adjust < 0.05 的通路 >= 5 条；R 脚本真实执行 |
| **不可发表** | R 未安装导致脚本未执行（log 中 `Rscript not found`） |

### Step 07 — MR

| 状态 | 判据 |
|------|------|
| **可发表** | 使用真实 eQTL + GWAS 数据；method 不是 "placeholder"；有敏感性分析 |
| **绝对不可发表** | `method = "placeholder"` 或 `b = NA` |
| **论文写法** | "以 eQTLGen cis-eQTLs 为暴露，FinnGen R10 RA GWAS 为结局，P < 5×10⁻⁸, F > 10, r² < 0.001 筛选工具变量，IVW 为主要分析方法，辅以 MR-Egger、Weighted Median 敏感性分析。" |

### Step 08 — 分子对接

| 状态 | 判据 |
|------|------|
| **可发表** | 受体 PDBQT 手工制备；结合位点坐标为真实值（非 0,0,0）；binding energy 有数值 |
| **绝对不可发表** | `status = "receptor_not_prepared"`；center = (0,0,0)；binding_energy = None |
| **论文写法** | "从 RCSB PDB 获取靶蛋白晶体结构（PDB ID: XXXX），使用 AutoDock Vina 1.2.5 进行分子对接，exhaustiveness = 8，结合能 <= -5.0 kcal/mol 判定为良好结合。" |

### Step 09 — scRNA

| 状态 | 判据 |
|------|------|
| **可发表** | 使用真实 GEO 数据；QC 参数合理；细胞分群注释可靠 |
| **不可发表** | `scrna.enabled: false` 或无 GEO 数据 |

---

## 3. 缺失输入时必须中断的步骤

有些步骤缺少数据可以跳过，有些**必须中断**：

| 步骤 | 缺失数据 | 处理方式 |
|------|---------|---------|
| Step 01 | 无 TCMSP 导出 | ⚠ 可继续（PubChem 回退）但成分不全 |
| Step 02 | 无任何靶点数据库 | ❌ **必须中断** — 至少需要 1 个来源 |
| Step 03 | 无任何疾病靶点 | ❌ **必须中断** — 至少需要 GeneCards |
| Step 04 | 交集为空 | ❌ **必须中断** — 检查 Step 02/03 |
| Step 05 | 无 MCODE/GEO | ⚠ 可继续（仅拓扑学）但证据弱 |
| Step 06 | R 未安装 | ❌ **必须中断** — 安装 R |
| Step 07 | 无 eQTL/GWAS | ⚠ 生成 placeholder → **不可发表** |
| Step 08 | 无 Vina/受体 | ⚠ 生成 placeholder → **不可发表** |
| Step 09 | 无 GEO scRNA | ⚠ 跳过（可选模块） |

---

## 4. 数据版本与查询日期记录

**每次真实运行必须记录：**

```yaml
# config.yaml 中填写
target_prediction:
  query_date: "2026-04-09"    # ← 实际查询日期

# 缓存文件带日期戳
data/raw/cached/tcmsp_20260409.csv
data/raw/cached/genecards_manual.csv   # 文件内记录下载日期
```

**论文 Methods 中必须注明：**
- 各数据库的访问日期
- STRING 版本号
- FinnGen release 版本（R10）
- eQTLGen 版本
- 软件版本号（config.yaml → software_versions）

---

## 5. 什么叫「最小可发表运行」

一次可以支撑论文的运行至少需要：

- [ ] Step 01: 真实 TCMSP 成分（有 OB/DL）
- [ ] Step 02: >= 2 个数据库的靶点
- [ ] Step 03: >= 2 个数据库的疾病靶点
- [ ] Step 04: 交集 >= 5 基因 + STRING PPI
- [ ] Step 05: >= 2 层筛选的 hub 基因
- [ ] Step 06: R 真实执行的 GO/KEGG 富集
- [ ] Step 07: 真实 eQTL + GWAS 的 MR 结果（非 placeholder）
- [ ] Step 08: 至少 1 对化合物-靶点的真实对接（binding energy 有值）
- [ ] 所有 query_date 已填写
- [ ] 所有原始数据库导出文件已缓存在 `data/raw/cached/`

---

## 6. 网络药理学方法学的天然局限

即使工程完全跑通，以下结论**仍需谨慎**：

| 结果 | 实际含义 | 不等于 |
|------|---------|--------|
| 靶点预测命中 | 数据库中有记录或算法预测 | 真实作用靶点 |
| PPI 中心性高 | 网络拓扑指标突出 | 关键致病基因 |
| MR 显著 | 血液 eQTL 与疾病 GWAS 有关联 | 药物在靶组织的作用机制 |
| 对接结合能低 | 几何形状互补 | 真实结合（需 SPR/CETSA 验证） |
| 富集到通路 | 基因集与通路基因集重叠 | 药物通过该通路起效 |

**补强建议：**
- MR 补充组织特异性 eQTL（GTEx 滑膜组织）
- 对接后做 100 ns MD 验证稳定性
- 用 scRNA 确认 hub 基因在靶细胞类型中高表达
- 湿实验验证：SPR（结合常数）、CETSA（热稳定性）、Western blot
