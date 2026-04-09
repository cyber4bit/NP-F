# Claude Code 任务提示词
# Network Pharmacology Pipeline — 开源项目搭建

## 背景
这是一个单一药物网络药理学分析的开源 GitHub 项目。
示例药物：血筒（Caulis Sargentodoxae）。
方法论来源：
  - CTS论文 (Xie et al., Phytomedicine 2025): scRNA-seq + SPR/CETSA 验证
  - BBR论文 (Zhao et al., Pharmacogenomics Pers Med 2025): MR + 全公开数据库链路
  - SPI论文 (Zhu & He, Front Pharmacol 2026): bulk + scRNA + 空间转录组 + AI辅助

## 项目文件已存在
- config.yaml         ← 所有参数集中管理
- README.md           ← 中英双语文档
- Snakefile           ← 工作流管理
- scripts/python/01_compounds.py
- scripts/python/04_ppi.py
- scripts/python/08_docking.py
- scripts/R/07_mr.R   ← TwoSampleMR (BBR模板)
- scripts/R/09_scrna.R ← Seurat 4.3 + Harmony + Nebulosa (CTS + SPI模板)
- env/environment.yml

## 需要补全的脚本（请逐一创建）

### 02_targets.py
"""
Step 02: Multi-database target prediction
- PubChem → structure download
- SwissTargetPrediction API (structure upload)
- PharmMapper: log instructions (no public API)
- CTD: log download instructions
- SEA: log instructions
- 合并去重，转换为官方gene symbol (UniProt/MyGene)
- 输出: data/processed/02_targets_merged.csv
"""

### 03_disease.py
"""
Step 03: Disease target collection
- GeneCards API 或手动下载指导
- OMIM download instructions
- DisGeNET API (requires free registration)
- 合并去重
- 输出: data/processed/03_disease_targets.csv
"""

### 05_hub_genes.py
"""
Step 05: Hub gene selection (三层筛选)
Layer 1: 拓扑指标 — degree + betweenness centrality from PPI edges (networkx)
Layer 2: MCODE — 打印Cytoscape使用说明
Layer 3: LASSO + Random Forest (scikit-learn)
         以共现分数或GEO差异基因矩阵作为特征矩阵
         LASSO: alpha search, 10-fold CV
         RF: n_estimators=1000, feature_importance threshold
- 输出: data/processed/05_hub_genes.csv (含 degree, betweenness, lasso_coef, rf_importance)
"""

### 06_enrichment.py  
"""
Step 06: Enrichment analysis
- 调用 R 包 clusterProfiler (通过 rpy2 或 subprocess)
- GO (BP/CC/MF) + KEGG + GSEA
- GSVA (SPI模板: 按基因表达高低30%分组)
- 生成: bubble plot (KEGG), bar plot (GO)
- 输出 figures: results/figures/06_enrichment_kegg.pdf
                results/figures/06_enrichment_go.pdf
"""

### 10_visualization.py
"""
Step 10: Summary figures
1. Venn diagram (已有04_venn, 这里做最终版)
2. PPI network visualization (networkx + matplotlib)
3. Hub gene heatmap
4. MR forest plot (从07_mr_results.csv读取)
5. Docking bar chart (binding energies)
6. CETSA/SPR curve template (根据CTS论文Fig 2G-J格式)
7. Kaplan-Meier template (根据SPI论文Fig 2B格式, survminer)
- 输出: results/figures/10_summary_figure.pdf
"""

## Codex 审查要求
完成脚本后：
- /codex:review scripts/python/05_hub_genes.py  ← 检查LASSO/RF算法实现
- /codex:adversarial-review scripts/R/07_mr.R   ← 检查MR统计假设和工具变量选择
- /codex:result  ← 汇总所有结果

## 代码要求
1. 每个脚本顶部：Input/Output/数据库版本/查询日期占位符注释
2. 所有阈值从 config.yaml 读取，不硬编码
3. 异常处理：数据库不可用时给出清晰的下载指引
4. 数据库原始导出缓存到 data/raw/cached/ 目录
5. logging 替代 print
6. 脚本末尾提示下一步命令
