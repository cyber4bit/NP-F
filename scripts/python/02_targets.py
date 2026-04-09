#!/usr/bin/env python3
"""
Step 02: Multi-database target prediction.

Input:
  - data/processed/01_compounds_filtered.csv
  - Optional manual exports in data/raw/cached/

Output:
  - data/processed/02_targets_merged.csv
  - data/raw/cached/02_targets_all_YYYY-MM-DD.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from datetime import date
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import requests
import yaml

from utils_validate import log_data_provenance, print_step_summary, validate_csv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def load_config(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def strip_string_columns(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = df.copy()
    for column in cleaned.columns:
        if cleaned[column].dtype == object:
            cleaned[column] = (
                cleaned[column]
                .astype(str)
                .str.strip()
                .str.strip('"')
                .str.strip("'")
                .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
            )
    return cleaned


def find_column(columns: list[str], patterns: list[str]) -> str | None:
    lowered = {column: str(column).strip().lower() for column in columns}
    for pattern in patterns:
        for column, lowered_name in lowered.items():
            if pattern in lowered_name:
                return column
    return None


def normalize_with_mygene(values: list[str], scopes: str) -> dict[str, str]:
    """Resolve identifiers to HGNC symbols via mygene, with an HTTP fallback."""
    queries = [value for value in values if value and str(value).strip()]
    if not queries:
        return {}

    mapping: dict[str, str] = {}
    try:
        import mygene

        mg = mygene.MyGeneInfo()
        results = mg.querymany(
            queries,
            scopes=scopes,
            fields="symbol",
            species="human",
            as_dataframe=False,
            returnall=False,
            verbose=False,
        )
        for hit in results:
            query = str(hit.get("query", "")).strip()
            symbol = hit.get("symbol")
            if query and symbol:
                mapping[query] = str(symbol).upper()
        if mapping:
            return mapping
    except ImportError:
        log.warning("mygene package not installed; falling back to MyGene.info HTTP queries.")
    except Exception as exc:
        log.warning("mygene querymany failed; falling back to HTTP queries: %s", exc)

    for query in queries:
        try:
            response = requests.get(
                "https://mygene.info/v3/query",
                params={
                    "q": query,
                    "scopes": scopes,
                    "fields": "symbol",
                    "species": "human",
                    "size": 1,
                },
                timeout=20,
            )
            response.raise_for_status()
            hits = response.json().get("hits", [])
            if hits and hits[0].get("symbol"):
                mapping[query] = str(hits[0]["symbol"]).upper()
        except Exception as exc:
            log.debug("MyGene.info HTTP lookup failed for %s: %s", query, exc)
    return mapping


def normalize_gene_symbols(df: pd.DataFrame, column: str = "gene_symbol") -> pd.DataFrame:
    if column not in df.columns:
        return df
    cleaned = df.copy()
    cleaned[column] = cleaned[column].astype("string").str.strip()
    queries = cleaned[column].dropna().unique().tolist()
    mapping = normalize_with_mygene(
        queries,
        scopes="symbol,alias,retired,reporter,ensembl.gene,entrezgene,uniprot",
    )
    if mapping:
        cleaned[column] = cleaned[column].map(lambda value: mapping.get(str(value), value))
    cleaned[column] = cleaned[column].astype("string").str.upper().str.strip()
    return cleaned


def uniprot_to_gene_symbol(uniprot_ids: list[str]) -> dict[str, str]:
    """Query MyGene.info to map UniProt accessions to HGNC symbols."""
    return normalize_with_mygene(
        uniprot_ids,
        scopes="uniprot,uniprot.Swiss-Prot,uniprot.TrEMBL",
    )


def standardize_target_df(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    """
    Standardize heterogeneous target export formats.

    Returns a DataFrame with the best-effort columns:
      - gene_symbol
      - score
      - uniprot_id
      - source
    """
    cleaned = strip_string_columns(df)
    gene_col = find_column(
        list(cleaned.columns),
        ["gene symbol", "genesymbol", "gene name", "target name", "gene", "symbol"],
    )
    score_col = find_column(list(cleaned.columns), ["probability", "score", "fit score", "pvalue"])
    uniprot_col = find_column(list(cleaned.columns), ["uniprot id", "uniprot", "accession"])

    rename_map: dict[str, str] = {}
    if gene_col:
        rename_map[gene_col] = "gene_symbol"
    if score_col and score_col != gene_col:
        rename_map[score_col] = "score"
    if uniprot_col and uniprot_col not in rename_map:
        rename_map[uniprot_col] = "uniprot_id"

    standardized = cleaned.rename(columns=rename_map).copy()
    keep_columns = [column for column in ["gene_symbol", "score", "uniprot_id"] if column in standardized]
    if not keep_columns:
        log.warning("%s export has no recognizable gene or UniProt columns; skipping.", source_name)
        return pd.DataFrame(columns=["gene_symbol", "score", "uniprot_id", "source"])
    standardized = standardized[keep_columns].copy()
    if "score" in standardized.columns:
        standardized["score"] = pd.to_numeric(standardized["score"], errors="coerce")
    standardized["source"] = source_name
    return standardized


def parse_stp_table(table_df: pd.DataFrame, prob_cutoff: float) -> pd.DataFrame:
    column_patterns = {
        "gene_symbol": ["gene name", "gene", "target", "common name"],
        "uniprot_id": ["uniprot id", "uniprot"],
        "score": ["probability", "prob"],
    }
    rename_map: dict[str, str] = {}
    for target, patterns in column_patterns.items():
        column = find_column(list(table_df.columns), patterns)
        if column:
            rename_map[column] = target

    parsed = table_df.rename(columns=rename_map).copy()
    if "score" in parsed.columns:
        parsed["score"] = pd.to_numeric(parsed["score"], errors="coerce")
        parsed = parsed[parsed["score"].fillna(0) >= prob_cutoff].copy()
    keep = [column for column in ["gene_symbol", "uniprot_id", "score"] if column in parsed.columns]
    if not keep:
        return pd.DataFrame()
    parsed = strip_string_columns(parsed[keep])
    parsed["source"] = "SwissTargetPrediction"
    return parsed


def log_stp_manual_instructions() -> None:
    log.info("Manual SwissTargetPrediction fallback:")
    log.info("  1. Visit https://www.swisstargetprediction.ch/")
    log.info("  2. Paste the compound SMILES and select Homo sapiens")
    log.info("  3. Download the CSV result")
    log.info("  4. Save it to data/raw/cached/stp_manual.csv")


def query_swiss_target_prediction(
    smiles: str,
    species: str = "Homo sapiens",
    prob_cutoff: float = 0.5,
) -> pd.DataFrame:
    """
    Primary: POST to STP endpoint and parse JSON if available.
    Fallback: parse HTML tables with case-insensitive column detection.
    If both fail: return an empty DataFrame and log manual instructions.
    Never raises.
    """
    payload = {"smiles": smiles, "organism": species.replace(" ", "_"), "ioi": "2"}
    endpoints = [
        "https://www.swisstargetprediction.ch/predict.php",
        "http://www.swisstargetprediction.ch/predict.php",
    ]

    for endpoint in endpoints:
        try:
            response = requests.post(endpoint, data=payload, timeout=90)
            response.raise_for_status()
            content_type = response.headers.get("Content-Type", "").lower()
            if "json" in content_type:
                data = response.json()
                if isinstance(data, dict):
                    for key in ["results", "targets", "data"]:
                        if isinstance(data.get(key), list):
                            parsed = parse_stp_table(pd.DataFrame(data[key]), prob_cutoff)
                            if not parsed.empty:
                                return parsed
                elif isinstance(data, list):
                    parsed = parse_stp_table(pd.DataFrame(data), prob_cutoff)
                    if not parsed.empty:
                        return parsed

            tables = pd.read_html(StringIO(response.text))
            for table in tables:
                parsed = parse_stp_table(table, prob_cutoff)
                if not parsed.empty:
                    return parsed
        except ValueError:
            pass
        except Exception as exc:
            log.warning("SwissTargetPrediction request failed via %s: %s", endpoint, exc)

    log.warning("SwissTargetPrediction parsing failed for the submitted SMILES.")
    log_stp_manual_instructions()
    return pd.DataFrame(columns=["gene_symbol", "uniprot_id", "score", "source"])


def load_manual_export(path: Path, source_name: str) -> pd.DataFrame:
    if not validate_csv(path, min_rows=1, step_name=f"Step 02 {source_name}"):
        return pd.DataFrame(columns=["gene_symbol", "score", "uniprot_id", "source"])
    return standardize_target_df(pd.read_csv(path), source_name)


def load_tcmsp_targets(cache_dir: Path, tcmsp_path: str | None) -> pd.DataFrame:
    if tcmsp_path and Path(tcmsp_path).exists():
        df = load_manual_export(Path(tcmsp_path), "TCMSP")
        log.info("Loaded TCMSP target export: %d row(s)", len(df))
        return df
    log.info("TCMSP herb-target instructions:")
    log.info("  1. Visit https://www.tcmsp-e.com/")
    log.info("  2. Search the herb and open the Related Targets tab")
    log.info("  3. Export CSV and save to %s", cache_dir / "tcmsp_targets_manual.csv")
    return pd.DataFrame(columns=["gene_symbol", "score", "uniprot_id", "source"])


def log_pharmmapper_instructions() -> None:
    log.info("PharmMapper manual instructions:")
    log.info("  1. Visit https://www.lilab-ecust.cn/pharmmapper/")
    log.info("  2. Upload a MOL2 structure and restrict to human targets")
    log.info("  3. Save the result to data/raw/cached/pharmmapper_manual.csv")


def log_ctd_instructions(drug_name: str) -> None:
    log.info("CTD manual instructions:")
    log.info("  1. Visit https://ctdbase.org/")
    log.info("  2. Search the compound or herb-related active ingredients for %s", drug_name)
    log.info("  3. Export the human chemical-gene interactions")
    log.info("  4. Save the result to data/raw/cached/ctd_manual.csv")


def log_sea_instructions() -> None:
    log.info("SEA manual instructions:")
    log.info("  1. Visit https://sea.bkslab.org/")
    log.info("  2. Upload SMILES or SDF input")
    log.info("  3. Save the result to data/raw/cached/sea_manual.csv")


def load_manual_exports(cache_dir: Path, stp_override: str | None = None) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    file_map = {
        "stp_manual.csv": "SwissTargetPrediction",
        "pharmmapper_manual.csv": "PharmMapper",
        "ctd_manual.csv": "CTD",
        "sea_manual.csv": "SEA",
    }

    for filename, source_name in file_map.items():
        file_path = cache_dir / filename
        if file_path.exists():
            frame = load_manual_export(file_path, source_name)
            if not frame.empty:
                frames.append(frame)
                log.info("Loaded %s manual export: %d row(s)", source_name, len(frame))

    if stp_override and Path(stp_override).exists():
        override_path = Path(stp_override)
        frame = load_manual_export(override_path, "SwissTargetPrediction")
        if not frame.empty:
            frames.append(frame)
            log.info("Loaded SwissTargetPrediction override export: %d row(s)", len(frame))

    if not frames:
        return pd.DataFrame(columns=["gene_symbol", "score", "uniprot_id", "source"])
    return pd.concat(frames, ignore_index=True)


def aggregate_output(df: pd.DataFrame, today: str) -> pd.DataFrame:
    aggregated = (
        df.groupby("gene_symbol", as_index=False)
        .agg(source=("source", lambda values: "|".join(sorted(set(values)))))
        .sort_values("gene_symbol")
        .reset_index(drop=True)
    )
    aggregated["query_date"] = today
    return aggregated[["gene_symbol", "source", "query_date"]]


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 02: Multi-database target prediction")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--tcmsp_targets", default=None, help="TCMSP target CSV export")
    parser.add_argument("--stp_csv", default=None, help="SwissTargetPrediction CSV export")
    args = parser.parse_args()

    cfg = load_config(args.config)
    today = date.today().isoformat()
    tp_cfg = cfg["target_prediction"]

    out_dir = Path("data/processed")
    cache_dir = Path("data/raw/cached")
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    compounds_path = out_dir / "01_compounds_filtered.csv"
    if not validate_csv(
        compounds_path,
        required_columns=["Molecule Name", "SMILES", "source"],
        min_rows=1,
        step_name="Step 02",
    ):
        sys.exit(1)

    compounds = pd.read_csv(compounds_path)
    compounds = strip_string_columns(compounds)
    log.info("=== Step 02: Target prediction | %d compound(s) ===", len(compounds))

    all_targets: list[pd.DataFrame] = []

    tcmsp_df = load_tcmsp_targets(cache_dir, args.tcmsp_targets)
    if not tcmsp_df.empty:
        all_targets.append(tcmsp_df)

    manual_df = load_manual_exports(cache_dir, args.stp_csv)
    if not manual_df.empty:
        all_targets.append(manual_df)

    if args.stp_csv is None or not Path(args.stp_csv).exists():
        for _, row in compounds.iterrows():
            smiles = row.get("SMILES")
            if pd.isna(smiles) or not str(smiles).strip():
                log.warning("Skipping SwissTargetPrediction for %s: missing SMILES", row.get("Molecule Name"))
                continue
            if "smiles_valid" in row and row.get("smiles_valid") is False:
                log.warning(
                    "Skipping SwissTargetPrediction for %s: invalid SMILES",
                    row.get("Molecule Name"),
                )
                continue
            compound_name = str(row.get("Molecule Name", "compound")).strip()
            log.info("Querying SwissTargetPrediction for %s", compound_name)
            stp_df = query_swiss_target_prediction(
                str(smiles).strip(),
                species=tp_cfg["species"],
                prob_cutoff=float(tp_cfg["min_probability"]),
            )
            if not stp_df.empty:
                stp_df["compound"] = compound_name
                all_targets.append(stp_df)
            time.sleep(5)

    log_pharmmapper_instructions()
    log_ctd_instructions(cfg["drug"]["name"])
    log_sea_instructions()

    if not all_targets:
        log.error("No targets were collected from any source. Follow the manual export instructions.")
        sys.exit(1)

    merged = pd.concat(all_targets, ignore_index=True)
    merged = strip_string_columns(merged)

    if "uniprot_id" in merged.columns:
        needs_gene = merged.get("gene_symbol").isna() if "gene_symbol" in merged.columns else pd.Series(True, index=merged.index)
        uniprot_ids = merged.loc[needs_gene, "uniprot_id"].dropna().astype(str).unique().tolist()
        if uniprot_ids:
            uniprot_mapping = uniprot_to_gene_symbol(uniprot_ids)
            if "gene_symbol" not in merged.columns:
                merged["gene_symbol"] = pd.NA
            merged.loc[needs_gene, "gene_symbol"] = merged.loc[needs_gene, "uniprot_id"].map(
                lambda value: uniprot_mapping.get(str(value), pd.NA)
            )

    merged = normalize_gene_symbols(merged, "gene_symbol")
    merged["gene_symbol"] = merged["gene_symbol"].astype("string").str.upper().str.strip()
    merged = merged.dropna(subset=["gene_symbol", "source"]).copy()
    merged = merged[merged["gene_symbol"].ne("")]

    if merged.empty:
        log.error("Target collection succeeded, but no gene symbols remained after normalization.")
        sys.exit(1)

    source_counts = merged.groupby("source")["gene_symbol"].nunique().sort_values(ascending=False)
    log.info("Per-source unique gene counts:")
    for source_name, count in source_counts.items():
        log.info("  %s: %d", source_name, count)

    support_counts = (
        merged.groupby("gene_symbol")["source"]
        .nunique()
        .sort_values(ascending=False)
        .head(20)
    )
    log.info("Top 20 genes by source support:")
    for gene_symbol, count in support_counts.items():
        log.info("  %s: %d source(s)", gene_symbol, count)

    cache_path = cache_dir / f"02_targets_all_{today}.csv"
    merged.to_csv(cache_path, index=False)
    log.info("Cached merged target evidence -> %s", cache_path)

    result = aggregate_output(merged, today)
    out_path = out_dir / "02_targets_merged.csv"
    result.to_csv(out_path, index=False)
    log.info("Output -> %s (%d unique genes)", out_path, len(result))

    log_data_provenance(
        step="Step 02: Multi-database target prediction",
        output_path=out_path,
        sources=sorted(source_counts.index.tolist()),
        cfg=cfg,
        extra={
            "n_input": int(len(compounds)),
            "n_output": int(len(result)),
            "n_evidence_rows": int(len(merged)),
            "filters_applied": {"min_probability": tp_cfg["min_probability"]},
        },
    )
    print_step_summary(
        "Step 02",
        out_path,
        cfg,
        next_step="python scripts/python/03_disease.py --config config.yaml",
    )


if __name__ == "__main__":
    main()
