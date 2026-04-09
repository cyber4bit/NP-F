# Codex Implementation Prompt
# Network Pharmacology Pipeline: 血筒 × RA
# Task: Implement gaps to make this a fully runnable project

---

## Your Role

You are implementing a **network pharmacology research pipeline** in Python and R. The framework scaffolding already exists. Your job is to fill in every gap so the pipeline runs end-to-end on real data with no silent failures, no placeholder leakage, and no ambiguity about what is real vs. what is a stub.

Read all existing files before writing anything. Never overwrite working logic—extend it.

---

## Repository State (read these files first)

```
config.yaml                          # global parameters, drug=Xuetong, disease=RA
Snakefile                            # workflow orchestration
scripts/python/00_preflight.py       # env + data checks — COMPLETE
scripts/python/01_compounds.py       # ADME filter — MOSTLY COMPLETE, gaps below
scripts/python/02_targets.py         # target prediction — NEEDS FIXES
scripts/python/03_disease.py         # disease targets — NEEDS FIXES
scripts/python/04_ppi.py             # Venn + STRING PPI — MOSTLY COMPLETE
scripts/python/05_hub_genes.py       # hub gene selection — COMPLETE
scripts/python/06_enrichment.py      # GO/KEGG via R — COMPLETE
scripts/python/08_docking.py         # docking framework — NEEDS RA PDB + fixes
scripts/python/10_visualization.py   # summary figures — COMPLETE
scripts/python/utils_validate.py     # shared validators — COMPLETE
scripts/R/07_mr.R                    # MR analysis — COMPLETE
scripts/R/09_scrna.R                 # scRNA-seq — COMPLETE
docs/quickstart.md
docs/manual-inputs.md
docs/step-by-step-sop.md
docs/research-mode.md
env/environment.yml
env/requirements.txt
env/install_r_packages.R
```

---

## Design Constraints (non-negotiable)

1. **No silent placeholder leakage.** Any script that cannot proceed due to missing data must: (a) log a clear human-readable instruction, (b) write a CSV with a `status` column value of `"placeholder"` or exit with code 1. Never write a file that looks real but contains fake data without that marker.

2. **Provenance on every output.** After writing any `data/processed/*.csv`, call `log_data_provenance()` from `utils_validate.py`. Every output file must have a companion `.provenance.yaml`.

3. **Config-driven, not hardcoded.** All thresholds, file paths, database URLs come from `config.yaml`. Never hardcode a cutoff or path in script logic.

4. **Fail loudly on missing required inputs.** Call `validate_csv()` from `utils_validate.py` at the top of every script's `main()` before doing any work.

5. **Gene symbol normalization is mandatory in Steps 02 and 03.** All gene symbols must be converted to uppercase official HGNC symbols via MyGene.info before writing output.

6. **Column names in output CSVs are fixed contracts.** Downstream scripts depend on exact column names. Do not change them:
   - Step 01 output: `Molecule Name`, `SMILES`, `OB`, `DL`, `MW`, `source`, `query_date`
   - Step 02 output: `gene_symbol`, `source`, `query_date`
   - Step 03 output: `gene_symbol`, `source`, `disease`, `query_date`
   - Step 04 output (intersection): `gene_symbol`
   - Step 04 output (edges): `gene1`, `gene2`, `combined_score`
   - Step 05 output: `gene_symbol`, `degree`, `betweenness`, `lasso_coef`, `rf_importance`
   - Step 07 output: `gene`, `method`, `b`, `se`, `pval`, `or`, `or_lci95`, `or_uci95`
   - Step 08 output: `compound`, `target`, `pdb_id`, `binding_energy_kcal_mol`, `status`

---

## Task 1: Fix `scripts/python/01_compounds.py`

### Gap A — Column normalization for TCMSP exports

TCMSP exports use inconsistent headers across different download methods. Add a function `normalize_tcmsp_columns(df)` that maps any of these variants to the standard names:

| Standard | Accepted variants |
|----------|------------------|
| `Molecule Name` | `mol_name`, `molecule`, `Name`, `化合物`, `成分名称` |
| `OB` | `OB(%)`, `oral bioavailability`, `ob` |
| `DL` | `drug-likeness`, `DL score`, `dl` |
| `MW` | `Molecular Weight`, `mol_weight` |
| `SMILES` | `Canonical SMILES`, `smiles`, `SMILES` |

### Gap B — RDKit ADME fallback

When OB/DL are absent (PubChem fallback path), compute approximate values from SMILES using RDKit:
- MW from `Descriptors.MolWt`
- logP from `Descriptors.MolLogP`
- HBD from `Descriptors.NumHDonors`
- HBA from `Descriptors.NumHAcceptors`
- Set `OB = NaN` and `DL = NaN` (these require TCMSP model — do not fake them)
- Flag the row with `adme_source = "rdkit_partial"` vs `adme_source = "tcmsp"`

Add a guard: if rdkit is not installed, skip silently and log a warning.

### Gap C — SMILES validation

After loading and filtering, for every row with a SMILES value, validate it parses with RDKit (`Chem.MolFromSmiles`). Rows with invalid SMILES get `smiles_valid = False`. Log count of invalid SMILES but do not drop them (downstream scripts may still use them for target lookup by name).

---

## Task 2: Fix `scripts/python/02_targets.py`

### Gap A — STP HTML parsing is fragile

The current STP scraper uses `pd.read_html()` on raw HTML, which breaks whenever the site changes. Replace with:

```python
def query_swiss_target_prediction(smiles, species="Homo sapiens", prob_cutoff=0.5):
    """
    Primary: POST to STP API endpoint, parse JSON if available.
    Fallback: parse HTML table with column detection.
    If both fail: return empty DataFrame and log manual instructions.
    Never raise — always return DataFrame (possibly empty).
    """
```

Column detection must handle: `Target`, `Gene Name`, `Uniprot ID`, `Common name`, `Probability`. Use case-insensitive substring matching.

### Gap B — Standardize all manual export formats

Add a function `standardize_target_df(df, source_name)` that:
1. Finds the gene symbol column (case-insensitive search for: `gene`, `symbol`, `target name`, `gene name`, `genesymbol`)
2. Finds a probability/score column if present
3. Renames both to `gene_symbol` and `score`
4. Adds `source = source_name`
5. Returns cleaned DataFrame

Apply this to every manual export loaded in `load_manual_exports()`.

### Gap C — UniProt ID → gene symbol fallback

When a source provides UniProt IDs but not gene symbols, resolve them via MyGene.info:

```python
def uniprot_to_gene_symbol(uniprot_ids: list[str]) -> dict[str, str]:
    """Query MyGene.info to map UniProt accessions to HGNC symbols."""
```

Apply after merging all sources, before the final normalization step.

### Gap D — Source coverage report in log

After deduplication, log a per-source count AND the total unique gene count. Also log the top 20 gene symbols by source count (genes appearing in most databases = strongest evidence).

---

## Task 3: Fix `scripts/python/03_disease.py`

### Gap A — DisGeNET API update

The current API endpoint `https://www.disgenet.org/api/gda/disease` and auth pattern may be outdated. Update to use the v7 REST API:

```
Base URL: https://www.disgenet.org/api
Auth: Bearer token in header
Endpoint: GET /gda/disease/{disease_id}?source=ALL&format=json
```

Add a helper `search_disgenet_by_name(disease_name, api_key)` that:
1. First searches for the disease ID using `/disease/search?q={disease_name}`
2. Then fetches GDAs for the best-matching disease ID
3. Filters by `score >= config["disease_targets"]["disgenet_score_min"]`

If the API is unavailable or the key is missing, log the manual download instructions and return empty DataFrame (never exit 1 — disease targets can come from GeneCards alone).

### Gap B — GeneCards column detection

GeneCards exports from different browser sessions produce different column headers. Extend `load_manual_exports()` to handle:
- `Gene Symbol`, `Symbol`, `GeneSymbol`, `gene_symbol`
- `Relevance score`, `Relevance Score`, `Score`, `relevance`

Also strip whitespace and trailing/leading quotes from all string columns.

### Gap C — OMIM gene list parsing

OMIM exports contain phenotype entries mixed with gene entries. Add filtering:
- Keep rows where the entry type contains `gene` (case-insensitive) OR where a gene symbol column is non-empty
- Strip phenotype MIM numbers from the gene symbol column if they appear

---

## Task 4: Fix `scripts/python/08_docking.py`

### Gap A — Expand KNOWN_PDB for RA hub genes

Replace the current KNOWN_PDB dict with a comprehensive RA-relevant mapping. Use PDB IDs that have:
- Human protein structure (Homo sapiens)
- Resolution <= 2.5 Å preferred
- Bound co-crystallized ligand preferred (useful for binding site reference)

```python
KNOWN_PDB = {
    # RA core targets
    "TNF":      ["2AZ5", "1TNF"],       # TNF-alpha; 2AZ5 has SPD304 inhibitor bound
    "IL6":      ["1ALU", "4CNI"],       # IL-6; 4CNI with receptor complex
    "IL1B":     ["2MIB", "1HIB"],       # IL-1 beta
    "IL6ST":    ["3L5H"],               # gp130
    "JAK1":     ["3EYH", "4E5W"],       # JAK1 kinase domain; 4E5W with inhibitor
    "JAK2":     ["6E2Q", "2B7A"],       # JAK2; 6E2Q high-res
    "JAK3":     ["1YVJ", "3LXK"],       # JAK3
    "STAT3":    ["6NJS", "3CWG"],       # STAT3 SH2 domain
    "PTPN11":   ["2SHP", "3B7O"],       # SHP2
    "AKT1":     ["4EKL", "3CQW"],       # AKT1 with allosteric inhibitor
    "MAPK14":   ["3HEC", "1A9U"],       # p38-alpha (MAPK14)
    "MAPK1":    ["4GT3"],               # ERK2
    "MAPK3":    ["2ZOQ"],               # ERK1
    "VEGFA":    ["4WPB", "1VPF"],       # VEGF-A
    "TP53":     ["2OCJ", "3ZME"],       # p53 DNA-binding domain
    "IL17A":    ["5HHV"],               # IL-17A
    "TGFB1":    ["3KFD"],              # TGF-beta1
    "PTGS2":    ["5IKR", "1CX2"],       # COX-2 (PTGS2); 5IKR with celecoxib
    "MMP9":     ["4H82", "2OVZ"],       # MMP-9
    "MMP3":     ["1HFS"],               # MMP-3
    "CASP3":    ["2XYG", "1PAU"],       # Caspase-3
    "BCL2":     ["4LVT", "2W3L"],       # BCL-2
    "EGFR":     ["4WKQ", "1IVO"],       # EGFR kinase domain
    "SRC":      ["2SRC", "3EL8"],       # SRC kinase
    "HRAS":     ["4EFL", "3L8Y"],       # H-RAS
    "HSP90AA1": ["5J64", "2YE4"],       # HSP90-alpha
    "MAPK8":    ["3O2M"],               # JNK1
    "RELA":     ["2RAM"],               # NF-kB p65
    "PIK3CA":   ["4JPS"],               # PI3K-alpha
    "PTEN":     ["1D5R"],               # PTEN
    "CDK4":     ["2W9Z"],               # CDK4
    "CDK6":     ["2EUF"],               # CDK6
}
```

### Gap B — Binding site center from PDB

Replace the hardcoded `center = (0, 0, 0)` placeholder with a function that attempts to derive a binding site center from the PDB structure:

```python
def estimate_binding_center_from_pdb(pdb_path: Path) -> tuple[float, float, float] | None:
    """
    If the PDB contains a HETATM ligand record (non-water, non-ion),
    compute the geometric center of that ligand as the docking box center.
    Returns None if no ligand found — caller must warn user.
    Requires: biopython (Bio.PDB) OR manual coordinate specification.
    """
```

If biopython is not installed or no ligand is found:
- Log: `"WARNING: No co-crystallized ligand found in {pdb_id}. Set docking box center manually."`
- Write `status = "needs_manual_center"` in the output CSV for this target
- Do NOT run Vina — skip to next target

This prevents the silent (0,0,0) invalid docking result.

### Gap C — Receptor PDBQT preparation instructions

When `receptor_pdbqt` does not exist, instead of just logging a warning, write a per-target shell script with exact commands:

```python
def write_receptor_prep_script(pdb_path: Path, receptor_pdbqt: Path, out_dir: Path):
    """Write receptor_prep_{gene}.sh with AutoDockTools4 commands."""
    script = out_dir / f"receptor_prep_{pdb_path.stem}.sh"
    script.write_text(f"""\
#!/bin/bash
# Prepare receptor PDBQT for {pdb_path.stem}
# Requires: AutoDockTools4 (MGLTools) in PATH
# Step 1: Remove water and heteroatoms, add hydrogens
prepare_receptor -r {pdb_path} \\
    -o {receptor_pdbqt} \\
    -A hydrogens \\
    -U nphs_lps_waters_deleteAltB
echo "Done: {receptor_pdbqt}"
""")
    log.info(f"Receptor prep script → {script}")
```

---

## Task 5: Add `scripts/python/utils_geo.py` — GEO Data Download Helper

Create a new utility module (not a step script) that other steps can import:

```python
"""
GEO data download utilities.
Used by Step 05 (expression matrix) and Step 09 (scRNA).
"""

def download_geo_matrix(accession: str, out_dir: Path, force: bool = False) -> Path | None:
    """
    Download GEO series matrix file for bulk RNA-seq.
    URL pattern: https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/matrix/
    Returns path to downloaded file, or None if failed.
    """

def parse_geo_series_matrix(matrix_gz: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse GEO series matrix .gz file.
    Returns (expression_df, metadata_df).
    expression_df: rows=probes, cols=samples
    metadata_df: rows=samples, cols=metadata fields
    """

def map_probes_to_genes(expr_df: pd.DataFrame, platform_id: str) -> pd.DataFrame:
    """
    Map probe IDs to gene symbols using GEO platform annotation.
    Falls back to MyGene.info for probe ID mapping.
    """

def extract_sample_labels(metadata_df: pd.DataFrame) -> pd.DataFrame | None:
    """
    Attempt to extract case/control labels from GEO metadata.
    Looks for 'disease state', 'source_name', 'characteristics_ch1' fields.
    Returns DataFrame with columns [sample_id, label] or None if ambiguous.
    """
```

Integrate into Step 05: after `load_geo_matrix()` fails to find a local file, call `download_geo_matrix(cfg["disease"]["geo_bulk_accession"], cache_dir)` and then `parse_geo_series_matrix()`.

---

## Task 6: Integrate `utils_validate.py` into every step script

Every step script's `main()` must:

**At start:**
```python
from utils_validate import validate_csv, log_data_provenance, print_step_summary
# validate inputs
if not validate_csv("data/processed/XX_previous.csv", required_columns=["gene_symbol"], min_rows=1, step_name="Step NN"):
    sys.exit(1)
```

**At end (after writing output):**
```python
log_data_provenance(
    step="Step NN: Description",
    output_path=out_path,
    sources=list_of_sources_used,
    cfg=cfg,
    extra={"n_input": n_input, "n_output": len(result_df), "filters_applied": {...}}
)
print_step_summary("Step NN", out_path, cfg, next_step="python scripts/python/XX_next.py")
```

Do this for: 01, 02, 03, 04, 05, 08.

---

## Task 7: Add `.gitkeep` and directory scaffolding

Create the following empty files so git tracks the directory structure:

```
data/raw/cached/.gitkeep
data/raw/geo_datasets/.gitkeep
data/processed/.gitkeep
results/figures/.gitkeep
results/tables/.gitkeep
results/docking/.gitkeep
logs/.gitkeep
```

Also create `data/raw/cached/COLUMN_FORMATS.md` documenting the expected column format for every manual input file (one table per file, matching the content of `docs/manual-inputs.md` but in a format the scripts themselves can reference).

---

## Task 8: Add `scripts/python/09b_immune.py` — ssGSEA Immune Infiltration

Create a new Python script that wraps the R ssGSEA call (mirrors the pattern in `06_enrichment.py`):

```
Input:  data/raw/cached/geo_expression_matrix.csv (bulk RNA-seq)
        data/processed/05_hub_genes.csv
        config.yaml → immune settings
Output: data/processed/09b_immune_scores.csv
        results/figures/09b_immune_heatmap.pdf
```

The R script it generates should:
1. Load the expression matrix
2. Run ssGSEA using GSVA package with MSigDB C7 immune cell gene sets
3. Compute Spearman correlations between hub gene expression and immune cell infiltration scores
4. Output heatmap: hub genes × immune cell types, color = correlation coefficient
5. Flag hub genes with |r| > 0.3 and p < 0.05 as "immune-correlated"

Add rule `step09b_immune` to Snakefile (optional, not in `rule all`).

---

## Task 9: Fix Snakefile — forest plot output dependency

Current issue: `rule step07_mr` declares `results/figures/07_mr_forest.pdf` as output, but when MR data is missing, `07_mr.R` writes only the placeholder CSV and exits — the PDF is never created, breaking Snakemake.

Fix: change the forest plot from a hard output to a soft output using `ancient()` or make it a separate rule:

```python
rule step07_mr:
    input:  "data/processed/05_hub_genes.csv"
    output:
        results = "data/processed/07_mr_results.csv",
        forest  = touch("results/figures/07_mr_forest.pdf")   # touch ensures file exists
    log: "logs/07_mr.log"
    shell:
        "Rscript scripts/R/07_mr.R config.yaml > {log} 2>&1 || true"
```

And in `07_mr.R`, when data is missing: still create an empty PDF with a single text page saying "MR not executed — missing input data. See logs/07_mr.log."

Apply the same fix to `rule step08_docking` — docking can partially succeed; always write `08_docking_scores.csv` even if all rows have `status = "receptor_not_prepared"`.

---

## Task 10: Add `scripts/python/check_results.py` — Post-Run Quality Report

Create a standalone script (not a pipeline step) that reads all outputs and generates a human-readable quality report:

```
Usage: python scripts/python/check_results.py --config config.yaml

Output: results/tables/quality_report.txt  (and stdout)
```

Report structure:
```
=== Network Pharmacology Pipeline Quality Report ===
Drug:    Xuetong (血筒)
Disease: rheumatoid arthritis (RA)
Run date: 2026-04-09

STEP | FILE | ROWS | SOURCES | PLACEHOLDERS | PUBLISHABLE
01   | 01_compounds_filtered.csv | 18 | TCMSP | None | YES
02   | 02_targets_merged.csv | 247 | TCMSP,STP,PharmMapper | None | YES
03   | 03_disease_targets.csv | 1243 | GeneCards,OMIM,DisGeNET | None | YES
04   | 04_intersection_genes.csv | 52 | — | None | YES
05   | 05_hub_genes.csv | 12 | topo+MCODE+ML | None | YES
06   | 06_enrichment_kegg.csv | 34 | clusterProfiler | None | YES
07   | 07_mr_results.csv | 8 | eQTLGen+FinnGen | DETECTED | NO
08   | 08_docking_scores.csv | 96 | AutoDock Vina | status:receptor_not_prepared | PARTIAL

SUMMARY: 6/8 steps publishable. Fix steps: 07, 08
```

Logic: `publishable = True` iff no placeholder markers detected AND row count >= min threshold from `config.yaml → validation_rules`.

---

## Coding Standards

- All Python files: UTF-8, shebang `#!/usr/bin/env python3`, docstring at top
- Line length: 100 chars
- Logging: use module-level `log = logging.getLogger(__name__)`, not print()
- All file paths: `pathlib.Path`, never string concatenation
- All YAML loading: `yaml.safe_load()`, never `yaml.load()`
- All CSV writing: `df.to_csv(path, index=False)`
- Never use `sys.exit(1)` in library functions — only in `main()`
- R scripts: use `suppressPackageStartupMessages()`, set `options(warn=1)`
- R scripts: always check `file.exists()` before reading, use `tryCatch()` around all I/O

---

## Verification

After implementing, run:

```bash
python scripts/python/00_preflight.py --config config.yaml --skip-r
python scripts/python/check_results.py --config config.yaml
```

The preflight must show no `[FAIL]` on Python dependencies. The quality report must run without crashing even when all data files are missing (it should report each step as "not run" rather than crashing).

---

## What NOT to do

- Do not add new config keys without updating `config.yaml` and `docs/`
- Do not change output column names (contracts listed in Design Constraints)
- Do not add interactive prompts — all scripts must run non-interactively
- Do not hardcode file paths — use `Path("data/processed/")` relative to CWD
- Do not call `plt.show()` anywhere — always save to file
- Do not install packages inside scripts — check availability and log instructions
- Do not remove existing working logic — only add to it
