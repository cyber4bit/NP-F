# Manual Input Column Formats

This file documents the expected columns for every manual input placed in `data/raw/cached/`.
Scripts in `scripts/python/` are written against these conventions.

## `tcmsp_YYYY-MM-DD.csv`

| Standard column | Accepted variants |
| --- | --- |
| `Molecule Name` | `mol_name`, `molecule`, `Name`, `化合物`, `成分名称` |
| `OB` | `OB(%)`, `oral bioavailability`, `ob` |
| `DL` | `drug-likeness`, `DL score`, `dl` |
| `MW` | `Molecular Weight`, `mol_weight` |
| `SMILES` | `Canonical SMILES`, `smiles`, `SMILES` |

## `tcmsp_targets_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Gene symbol or target identifier | `gene`, `symbol`, `gene name`, `target name`, `genesymbol` |
| Optional score column | `score`, `probability`, `fit score` |
| Optional UniProt accession | `uniprot`, `uniprot id`, `accession` |

## `stp_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Gene symbol or target label | `Target`, `Gene Name`, `Common name` |
| UniProt accession | `Uniprot ID`, `UniProt`, `accession` |
| Probability | `Probability`, `probability`, `prob` |

## `pharmmapper_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Target or gene column | `target name`, `gene`, `symbol`, `gene name` |
| Optional score | `fit score`, `score` |
| Optional UniProt accession | `uniprot`, `accession` |

## `ctd_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Gene symbol | `Gene Symbol`, `symbol`, `gene_symbol` |
| Optional score or evidence | `score`, `interaction` |

## `sea_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Gene or target column | `Target Name`, `Gene Symbol`, `gene`, `symbol` |
| Optional score | `score`, `probability`, `e-value`, `p-value` |

## `genecards_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Gene symbol | `Gene Symbol`, `Symbol`, `GeneSymbol`, `gene_symbol` |
| Relevance score | `Relevance score`, `Relevance Score`, `Score`, `relevance` |

## `omim_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Gene symbol | `Gene Symbol`, `symbol`, `gene`, `gene_symbol` |
| Optional entry type | `Entry Type`, `Record Type`, `type` |

Rows are kept when the entry type contains `gene` or when a gene symbol is non-empty.
Phenotype MIM numbers are stripped from the gene symbol field automatically.

## `disgenet_manual.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Gene symbol | `geneSymbol`, `gene_symbol`, `Gene Symbol`, `symbol` |
| Score | `score`, `Score` |

## `mcode_clusters.csv`

| Expected content | Accepted column patterns |
| --- | --- |
| Cluster gene column | `gene_symbol`, `Gene`, `Node`, `Name`, `symbol` |

## `geo_expression_matrix.csv`

Final materialized format expected by downstream scripts:

| Requirement | Value |
| --- | --- |
| Rows | samples |
| Columns | HGNC gene symbols |
| Index column | sample identifier |
| Values | numeric expression values |

## `geo_sample_labels.csv`

| Required column | Meaning |
| --- | --- |
| `sample_id` | Must match the expression matrix row index |
| `label` | Case/control or another sample-level endpoint |

## `eqtlgen_cis_YYYY-MM-DD.tsv`

| Required column |
| --- |
| `SNP` |
| `beta` |
| `se` |
| `pval` |
| `effect_allele` |
| `other_allele` |
| `eaf` |
| `gene_symbol` |
| `N` |

## `finngen_r10_*.gz`

| Required column |
| --- |
| `rsid` |
| `beta` |
| `sebeta` |
| `pval` |
| `alt` |
| `ref` |
| `af_alt` |

## `msigdb_c7_immune.gmt`

Used by `scripts/python/09b_immune.py`.

| Requirement | Value |
| --- | --- |
| File format | GMT |
| Content | Immune-related MSigDB C7 gene sets |
| Config key | `immune.gene_sets_gmt` |
