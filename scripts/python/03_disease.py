#!/usr/bin/env python3
"""
Step 03: Disease target collection.

Input:
  - config.yaml
  - Optional manual exports in data/raw/cached/

Output:
  - data/processed/03_disease_targets.csv
  - data/raw/cached/disgenet_YYYY-MM-DD.csv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import requests
import yaml

from utils_validate import log_data_provenance, print_step_summary, validate_csv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

DISGENET_SCORE_MIN = 0.1


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


def log_genecards_instructions(disease_name: str, relevance_min: float) -> None:
    log.info("GeneCards manual instructions:")
    log.info("  1. Visit https://www.genecards.org/")
    log.info("  2. Search %s", disease_name)
    log.info("  3. Filter relevance score >= %s", relevance_min)
    log.info("  4. Save the CSV to data/raw/cached/genecards_manual.csv")


def log_omim_instructions(disease_name: str) -> None:
    log.info("OMIM manual instructions:")
    log.info("  1. Visit https://www.omim.org/")
    log.info("  2. Search %s", disease_name)
    log.info("  3. Export gene-focused results to data/raw/cached/omim_manual.csv")


def log_disgenet_instructions() -> None:
    log.info("DisGeNET manual instructions:")
    log.info("  1. Register at https://www.disgenet.org/signup/")
    log.info("  2. Download the GDA result for rheumatoid arthritis")
    log.info("  3. Save it to data/raw/cached/disgenet_manual.csv")


def select_best_disgenet_hit(results: list[dict[str, Any]], disease_name: str) -> str | None:
    disease_name_lower = disease_name.lower()
    best_id: str | None = None
    best_score = -1
    for item in results:
        identifier = item.get("diseaseid") or item.get("diseaseId") or item.get("id")
        label = str(item.get("name") or item.get("disease_name") or item.get("diseaseName") or "")
        if not identifier:
            continue
        score = 0
        if label.lower() == disease_name_lower:
            score += 3
        if disease_name_lower in label.lower():
            score += 2
        if str(identifier).upper().startswith("C"):
            score += 1
        if score > best_score:
            best_score = score
            best_id = str(identifier)
    return best_id


def parse_disgenet_payload(payload: Any) -> pd.DataFrame:
    if isinstance(payload, list):
        return pd.DataFrame(payload)
    if isinstance(payload, dict):
        for key in ["payload", "results", "data"]:
            if isinstance(payload.get(key), list):
                return pd.DataFrame(payload[key])
    return pd.DataFrame()


def search_disgenet_by_name(disease_name: str, api_key: str | None) -> pd.DataFrame:
    """
    Search DisGeNET v7 for a disease name, then fetch filtered gene-disease associations.
    Returns an empty DataFrame on missing credentials or API errors.
    """
    if not api_key:
        log.warning("DisGeNET API key not set.")
        log_disgenet_instructions()
        return pd.DataFrame(columns=["gene_symbol", "score", "source"])

    base_url = "https://www.disgenet.org/api"
    headers = {"Authorization": f"Bearer {api_key}", "Accept": "application/json"}
    try:
        search_response = requests.get(
            f"{base_url}/disease/search",
            headers=headers,
            params={"q": disease_name},
            timeout=30,
        )
        search_response.raise_for_status()
        search_df = parse_disgenet_payload(search_response.json())
        if search_df.empty:
            log.warning("DisGeNET disease search returned no matches for %s", disease_name)
            return pd.DataFrame(columns=["gene_symbol", "score", "source"])

        disease_id = select_best_disgenet_hit(search_df.to_dict("records"), disease_name)
        if not disease_id:
            log.warning("Could not determine a DisGeNET disease identifier for %s", disease_name)
            return pd.DataFrame(columns=["gene_symbol", "score", "source"])

        gda_response = requests.get(
            f"{base_url}/gda/disease/{disease_id}",
            headers=headers,
            params={"source": "ALL", "format": "json"},
            timeout=60,
        )
        gda_response.raise_for_status()
        gda_df = parse_disgenet_payload(gda_response.json())
        if gda_df.empty:
            log.warning("DisGeNET GDA endpoint returned no associations for %s (%s)", disease_name, disease_id)
            return pd.DataFrame(columns=["gene_symbol", "score", "source"])

        gene_col = find_column(list(gda_df.columns), ["gene symbol", "genesymbol", "symbol", "gene"])
        score_col = find_column(list(gda_df.columns), ["score"])
        rename_map: dict[str, str] = {}
        if gene_col:
            rename_map[gene_col] = "gene_symbol"
        if score_col and score_col != gene_col:
            rename_map[score_col] = "score"
        standardized = strip_string_columns(gda_df.rename(columns=rename_map))
        keep_columns = [column for column in ["gene_symbol", "score"] if column in standardized.columns]
        if not keep_columns:
            log.warning("DisGeNET response did not include recognizable gene columns.")
            return pd.DataFrame(columns=["gene_symbol", "score", "source"])
        standardized = standardized[keep_columns].copy()
        if "score" in standardized.columns:
            standardized["score"] = pd.to_numeric(standardized["score"], errors="coerce")
            standardized = standardized[standardized["score"].fillna(0) >= DISGENET_SCORE_MIN].copy()
        standardized["source"] = "DisGeNET"
        return standardized
    except Exception as exc:
        log.warning("DisGeNET API unavailable: %s", exc)
        log_disgenet_instructions()
        return pd.DataFrame(columns=["gene_symbol", "score", "source"])


def standardize_genecards(df: pd.DataFrame, relevance_min: float) -> pd.DataFrame:
    cleaned = strip_string_columns(df)
    gene_col = find_column(list(cleaned.columns), ["gene symbol", "genesymbol", "symbol", "gene_symbol"])
    score_col = find_column(list(cleaned.columns), ["relevance score", "relevance", "score"])
    if not gene_col:
        log.warning("GeneCards export does not contain a recognizable gene symbol column.")
        return pd.DataFrame(columns=["gene_symbol", "source"])

    standardized = cleaned.rename(columns={gene_col: "gene_symbol"})
    if score_col and score_col != gene_col:
        standardized = standardized.rename(columns={score_col: "score"})
        standardized["score"] = pd.to_numeric(standardized["score"], errors="coerce")
        standardized = standardized[standardized["score"].fillna(0) >= relevance_min].copy()
    standardized["source"] = "GeneCards"
    keep = [column for column in ["gene_symbol", "score", "source"] if column in standardized.columns]
    return standardized[keep]


def clean_omim_gene_symbol(value: Any) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().strip('"').strip("'")
    text = re.sub(r"\b\d{5,7}\b", " ", text)
    text = text.replace("///", ";").replace("/", ";").replace(",", ";")
    candidates = [candidate.strip() for candidate in text.split(";") if candidate.strip()]
    for candidate in candidates:
        if re.fullmatch(r"[A-Za-z0-9\-]+", candidate):
            return candidate.upper()
    return candidates[0].upper() if candidates else None


def standardize_omim(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = strip_string_columns(df)
    gene_col = find_column(list(cleaned.columns), ["gene symbol", "genesymbol", "symbol", "gene"])
    entry_type_col = find_column(list(cleaned.columns), ["entry type", "record type", "type"])
    if gene_col:
        cleaned = cleaned.rename(columns={gene_col: "gene_symbol"})
        cleaned["gene_symbol"] = cleaned["gene_symbol"].map(clean_omim_gene_symbol)
    else:
        cleaned["gene_symbol"] = pd.NA

    keep_mask = cleaned["gene_symbol"].notna()
    if entry_type_col:
        keep_mask = keep_mask | cleaned[entry_type_col].astype("string").str.contains(
            "gene",
            case=False,
            na=False,
        )

    standardized = cleaned[keep_mask].copy()
    standardized["source"] = "OMIM"
    return standardized[["gene_symbol", "source"]]


def standardize_disgenet_manual(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = strip_string_columns(df)
    gene_col = find_column(list(cleaned.columns), ["gene symbol", "genesymbol", "gene_symbol", "symbol", "gene"])
    score_col = find_column(list(cleaned.columns), ["score"])
    if not gene_col:
        log.warning("DisGeNET manual export does not contain a recognizable gene symbol column.")
        return pd.DataFrame(columns=["gene_symbol", "source"])

    standardized = cleaned.rename(columns={gene_col: "gene_symbol"})
    if score_col and score_col != gene_col:
        standardized = standardized.rename(columns={score_col: "score"})
        standardized["score"] = pd.to_numeric(standardized["score"], errors="coerce")
        standardized = standardized[standardized["score"].fillna(0) >= DISGENET_SCORE_MIN].copy()
    standardized["source"] = "DisGeNET"
    keep = [column for column in ["gene_symbol", "score", "source"] if column in standardized.columns]
    return standardized[keep]


def load_manual_exports(cache_dir: Path, relevance_min: float) -> list[pd.DataFrame]:
    frames: list[pd.DataFrame] = []

    genecards_path = cache_dir / "genecards_manual.csv"
    if genecards_path.exists() and validate_csv(genecards_path, min_rows=1, step_name="Step 03 GeneCards"):
        frame = standardize_genecards(pd.read_csv(genecards_path), relevance_min)
        if not frame.empty:
            frames.append(frame)
            log.info("Loaded GeneCards export: %d row(s)", len(frame))

    omim_path = cache_dir / "omim_manual.csv"
    if omim_path.exists() and validate_csv(omim_path, min_rows=1, step_name="Step 03 OMIM"):
        frame = standardize_omim(pd.read_csv(omim_path))
        if not frame.empty:
            frames.append(frame)
            log.info("Loaded OMIM export: %d row(s)", len(frame))

    disgenet_path = cache_dir / "disgenet_manual.csv"
    if disgenet_path.exists() and validate_csv(disgenet_path, min_rows=1, step_name="Step 03 DisGeNET"):
        frame = standardize_disgenet_manual(pd.read_csv(disgenet_path))
        if not frame.empty:
            frames.append(frame)
            log.info("Loaded DisGeNET manual export: %d row(s)", len(frame))

    return frames


def aggregate_output(df: pd.DataFrame, disease_name: str, today: str) -> pd.DataFrame:
    aggregated = (
        df.groupby("gene_symbol", as_index=False)
        .agg(source=("source", lambda values: "|".join(sorted(set(values)))))
        .sort_values("gene_symbol")
        .reset_index(drop=True)
    )
    aggregated["disease"] = disease_name
    aggregated["query_date"] = today
    return aggregated[["gene_symbol", "source", "disease", "query_date"]]


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 03: Disease target collection")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--disgenet_key", default=None, help="DisGeNET API key")
    args = parser.parse_args()

    cfg = load_config(args.config)
    disease_name = cfg["disease"]["name"]
    disease_query = cfg["disease"].get("mesh_term") or disease_name
    today = date.today().isoformat()
    relevance_min = float(cfg["disease_targets"]["genecards_relevance_min"])

    global DISGENET_SCORE_MIN
    DISGENET_SCORE_MIN = float(cfg["disease_targets"]["disgenet_score_min"])

    out_dir = Path("data/processed")
    cache_dir = Path("data/raw/cached")
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    log.info("=== Step 03: Disease target collection | %s ===", disease_name)
    log_genecards_instructions(disease_name, relevance_min)
    log_omim_instructions(disease_name)

    api_key = args.disgenet_key or os.environ.get("DISGENET_API_KEY")
    all_targets: list[pd.DataFrame] = []

    disgenet_df = search_disgenet_by_name(disease_query, api_key)
    if not disgenet_df.empty:
        cache_path = cache_dir / f"disgenet_{today}.csv"
        disgenet_df.to_csv(cache_path, index=False)
        log.info("Cached DisGeNET API result -> %s", cache_path)
        all_targets.append(disgenet_df)

    manual_frames = load_manual_exports(cache_dir, relevance_min)
    all_targets.extend(manual_frames)

    if not all_targets:
        log.error("No disease targets were collected. Follow the manual export instructions above.")
        sys.exit(1)

    merged = pd.concat(all_targets, ignore_index=True)
    merged = strip_string_columns(merged)
    if "gene_symbol" not in merged.columns:
        log.error("No gene_symbol column remained after standardizing disease target inputs.")
        sys.exit(1)

    merged = normalize_gene_symbols(merged, "gene_symbol")
    merged["gene_symbol"] = merged["gene_symbol"].astype("string").str.upper().str.strip()
    merged = merged.dropna(subset=["gene_symbol", "source"]).copy()
    merged = merged[merged["gene_symbol"].ne("")]

    if merged.empty:
        log.error("Disease target collection succeeded, but normalization removed all gene symbols.")
        sys.exit(1)

    source_counts = merged.groupby("source")["gene_symbol"].nunique().sort_values(ascending=False)
    log.info("Per-source unique disease gene counts:")
    for source_name, count in source_counts.items():
        log.info("  %s: %d", source_name, count)

    result = aggregate_output(merged, disease_name, today)
    out_path = out_dir / "03_disease_targets.csv"
    result.to_csv(out_path, index=False)
    log.info("Output -> %s (%d unique genes)", out_path, len(result))

    log_data_provenance(
        step="Step 03: Disease target collection",
        output_path=out_path,
        sources=sorted(source_counts.index.tolist()),
        cfg=cfg,
        extra={
            "n_input": int(sum(len(frame) for frame in all_targets)),
            "n_output": int(len(result)),
            "filters_applied": {
                "genecards_relevance_min": relevance_min,
                "disgenet_score_min": DISGENET_SCORE_MIN,
            },
        },
    )
    print_step_summary(
        "Step 03",
        out_path,
        cfg,
        next_step="python scripts/python/04_ppi.py --config config.yaml",
    )


if __name__ == "__main__":
    main()
