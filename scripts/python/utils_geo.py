#!/usr/bin/env python3
"""
GEO data download utilities.

Used by Step 05 (expression matrix) and Step 09 (scRNA).
"""

from __future__ import annotations

import csv
import gzip
import logging
import re
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import requests

log = logging.getLogger(__name__)


def _series_prefix(accession: str) -> str:
    return re.sub(r"\d{3}$", "nnn", accession)


def _parse_geo_line(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))


def _clean_gene_symbol(value: Any) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text in {"---", "NA", "nan"}:
        return None
    text = text.replace("///", ";").replace("//", ";").replace(",", ";")
    candidates = [candidate.strip() for candidate in text.split(";") if candidate.strip()]
    for candidate in candidates:
        if re.fullmatch(r"[A-Za-z0-9\-]+", candidate):
            return candidate.upper()
    return candidates[0].upper() if candidates else None


def _query_mygene_reporters(queries: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not queries:
        return mapping

    try:
        import mygene

        mg = mygene.MyGeneInfo()
        results = mg.querymany(
            queries,
            scopes="reporter,accession",
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
    except Exception as exc:
        log.warning("Reporter-to-gene lookup via mygene failed: %s", exc)

    for query in queries:
        try:
            response = requests.get(
                "https://mygene.info/v3/query",
                params={
                    "q": query,
                    "scopes": "reporter,accession",
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
            log.debug("MyGene.info reporter lookup failed for %s: %s", query, exc)
    return mapping


def download_geo_matrix(accession: str, out_dir: Path, force: bool = False) -> Path | None:
    """
    Download a GEO series matrix file.

    URL pattern:
      https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/matrix/
    """
    if not accession:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{accession}_series_matrix.txt.gz"
    if out_path.exists() and not force:
        log.info("Using cached GEO series matrix: %s", out_path)
        return out_path

    prefix = _series_prefix(accession)
    url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{prefix}/{accession}/matrix/{accession}_series_matrix.txt.gz"
    )
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        out_path.write_bytes(response.content)
        log.info("Downloaded GEO series matrix -> %s", out_path)
        return out_path
    except Exception as exc:
        log.warning("Failed to download GEO series matrix for %s: %s", accession, exc)
        return None


def parse_geo_series_matrix(matrix_gz: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a GEO series matrix file.

    Returns:
      - expression_df: rows = probes, cols = samples
      - metadata_df: rows = samples, cols = metadata fields
    """
    sample_meta: dict[str, list[str]] = {}
    field_counts: dict[str, int] = {}
    series_meta: dict[str, str] = {}
    table_lines: list[str] = []
    in_table = False

    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            lower = line.lower()
            if lower.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if lower.startswith("!series_matrix_table_end"):
                in_table = False
                break
            if in_table:
                table_lines.append(line)
                continue
            if line.startswith("!Sample_"):
                parts = _parse_geo_line(line)
                key = parts[0].replace("!Sample_", "", 1)
                field_counts[key] = field_counts.get(key, 0) + 1
                if field_counts[key] > 1:
                    key = f"{key}_{field_counts[key]}"
                sample_meta[key] = parts[1:]
            elif line.startswith("!Series_"):
                parts = _parse_geo_line(line)
                key = parts[0].replace("!Series_", "", 1)
                value = parts[1] if len(parts) > 1 else ""
                series_meta[key] = value

    if not table_lines:
        raise ValueError(f"No expression table found in {matrix_gz}")

    expr_df = pd.read_csv(StringIO("\n".join(table_lines)), sep="\t")
    first_column = expr_df.columns[0]
    expr_df = expr_df.rename(columns={first_column: "probe_id"}).set_index("probe_id")
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

    sample_ids = expr_df.columns.astype(str).tolist()
    metadata_df = pd.DataFrame(index=sample_ids)
    for key, values in sample_meta.items():
        values = values[: len(sample_ids)] + [pd.NA] * max(0, len(sample_ids) - len(values))
        metadata_df[key] = values[: len(sample_ids)]

    metadata_df.index.name = "sample_id"
    metadata_df.attrs["platform_id"] = str(series_meta.get("platform_id", "")).strip()
    metadata_df.attrs["geo_accession"] = str(series_meta.get("geo_accession", "")).strip()
    return expr_df, metadata_df


def _download_platform_annotation(platform_id: str, out_dir: Path) -> Path | None:
    if not platform_id:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{platform_id}.annot.gz"
    if out_path.exists():
        return out_path

    prefix = _series_prefix(platform_id)
    url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/"
        f"{prefix}/{platform_id}/annot/{platform_id}.annot.gz"
    )
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        out_path.write_bytes(response.content)
        return out_path
    except Exception as exc:
        log.warning("Failed to download GEO platform annotation for %s: %s", platform_id, exc)
        return None


def _parse_platform_annotation(annotation_gz: Path) -> pd.DataFrame:
    table_lines: list[str] = []
    in_table = False
    with gzip.open(annotation_gz, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            lower = line.lower()
            if lower.startswith("!platform_table_begin"):
                in_table = True
                continue
            if lower.startswith("!platform_table_end"):
                break
            if in_table:
                table_lines.append(line)
    if not table_lines:
        raise ValueError(f"No platform table found in {annotation_gz}")
    return pd.read_csv(StringIO("\n".join(table_lines)), sep="\t")


def map_probes_to_genes(expr_df: pd.DataFrame, platform_id: str) -> pd.DataFrame:
    """
    Map probe IDs to gene symbols using the GEO platform annotation.

    Falls back to MyGene.info reporter mapping when annotation is unavailable.
    Returns a DataFrame with rows = genes and cols = samples.
    """
    if expr_df.empty:
        return expr_df

    mapping: dict[str, str] = {}
    annotation_path = _download_platform_annotation(platform_id, Path("data/raw/cached"))
    if annotation_path:
        try:
            annotation_df = _parse_platform_annotation(annotation_path)
            probe_col = None
            for candidate in ["ID", "ID_REF", "Probe Set ID", "ProbeID"]:
                if candidate in annotation_df.columns:
                    probe_col = candidate
                    break
            gene_col = None
            for candidate in ["Gene symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol"]:
                if candidate in annotation_df.columns:
                    gene_col = candidate
                    break
            if probe_col and gene_col:
                for _, row in annotation_df[[probe_col, gene_col]].dropna().iterrows():
                    gene_symbol = _clean_gene_symbol(row[gene_col])
                    if gene_symbol:
                        mapping[str(row[probe_col]).strip()] = gene_symbol
            else:
                log.warning("Platform annotation %s lacked probe/gene columns.", platform_id)
        except Exception as exc:
            log.warning("Failed to parse platform annotation %s: %s", platform_id, exc)

    if not mapping:
        reporter_queries = [str(index).strip() for index in expr_df.index]
        mapping = _query_mygene_reporters(reporter_queries)

    mapped = expr_df.copy()
    mapped["gene_symbol"] = [mapping.get(str(index).strip()) for index in mapped.index]
    mapped = mapped.dropna(subset=["gene_symbol"]).copy()
    if mapped.empty:
        log.warning("Probe-to-gene mapping produced no gene symbols.")
        return pd.DataFrame(columns=expr_df.columns)

    mapped["gene_symbol"] = mapped["gene_symbol"].astype(str).str.upper().str.strip()
    mapped = mapped.drop(columns=[], errors="ignore")
    gene_expr = mapped.groupby("gene_symbol").mean(numeric_only=True)
    gene_expr.index.name = "gene_symbol"
    return gene_expr


def extract_sample_labels(metadata_df: pd.DataFrame) -> pd.DataFrame | None:
    """
    Infer case/control labels from GEO sample metadata.

    Returns a DataFrame with columns [sample_id, label] or None if ambiguous.
    """
    if metadata_df.empty:
        return None

    candidate_columns = [
        column
        for column in metadata_df.columns
        if any(
            token in str(column).lower()
            for token in ["disease state", "source_name", "characteristics_ch1", "title"]
        )
    ]
    if not candidate_columns:
        return None

    labels: list[dict[str, str]] = []
    ambiguous = 0
    for sample_id, row in metadata_df[candidate_columns].fillna("").iterrows():
        text = " | ".join(row.astype(str)).lower()
        has_control = any(token in text for token in ["healthy control", "control", "normal", "hc"])
        has_case = any(
            token in text
            for token in ["rheumatoid arthritis", "arthritis", "ra patient", "patient", "disease"]
        )
        if has_case and not has_control:
            label = "case"
        elif has_control and not has_case:
            label = "control"
        elif "healthy control" in text:
            label = "control"
        elif "rheumatoid arthritis" in text:
            label = "case"
        else:
            ambiguous += 1
            continue
        labels.append({"sample_id": sample_id, "label": label})

    label_df = pd.DataFrame(labels)
    if label_df.empty or label_df["label"].nunique() < 2:
        return None
    if ambiguous > len(metadata_df) * 0.5:
        return None
    return label_df
