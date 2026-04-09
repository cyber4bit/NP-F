# data/ 目录说明

## 结构

```
data/
├── raw/
│   ├── cached/              # 数据库导出的原始文件（带日期戳，勿修改）
│   │   ├── tcmsp_20260409.csv
│   │   ├── genecards_manual.csv
│   │   ├── stp_manual.csv
│   │   ├── eqtlgen_cis_20260409.tsv
│   │   └── ...
│   └── geo_datasets/        # GEO 下载的原始数据
│       ├── GSE93272/
│       └── GSE159117/
└── processed/               # 流程生成的中间文件（可通过 snakemake clean 清除）
    ├── 01_compounds_filtered.csv
    ├── 02_targets_merged.csv
    ├── ...
    └── 09_scrna_done.flag
```

## 规则

1. `raw/cached/` 中的文件是**原始证据**，不可修改，带日期戳
2. `processed/` 中的文件由脚本生成，可随时重新生成
3. 手工下载的文件命名规则见 `docs/manual-inputs.md`
4. `.provenance.yaml` 文件记录每个输出的来源信息
