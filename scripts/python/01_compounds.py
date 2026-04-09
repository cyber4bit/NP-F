#!/usr/bin/env python3
"""
Step 01: Compound collection and ADME screening.

Input:
  - config.yaml
  - Optional TCMSP CSV export

Output:
  - data/processed/01_compounds_filtered.csv
  - data/raw/cached/tcmsp_YYYY-MM-DD.csv
  - data/raw/cached/pubchem_YYYY-MM-DD.csv
"""

from __future__ import annotations

import argparse
import logging
import math
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

TCMSP_COLUMN_ALIASES = {
    "Molecule Name": {"molecule name", "mol_name", "molecule", "name", "化合物", "成分名称"},
    "OB": {"ob", "ob(%)", "oral bioavailability"},
    "DL": {"dl", "drug-likeness", "dl score"},
    "MW": {"mw", "molecular weight", "mol_weight"},
    "SMILES": {"smiles", "canonical smiles"},
}


def load_config(config_path: str) -> dict[str, Any]:
    with open(config_path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def canonicalize_column_name(column: Any) -> str:
    return str(column).strip().lower()


def normalize_tcmsp_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize common TCMSP export header variants to contract column names."""
    rename_map: dict[str, str] = {}
    used_targets: set[str] = set()
    for column in df.columns:
        canonical = canonicalize_column_name(column)
        for target, aliases in TCMSP_COLUMN_ALIASES.items():
            if target in used_targets:
                continue
            if canonical == target.lower() or canonical in aliases:
                rename_map[column] = target
                used_targets.add(target)
                break
    return df.rename(columns=rename_map)


def import_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        return Chem, Descriptors
    except ImportError:
        return None, None


def pubchem_search_by_name(name: str, max_hits: int = 10) -> list[dict[str, Any]]:
    """Search PubChem by compound name and return normalized rows."""
    cids_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{requests.utils.quote(name)}/cids/JSON"
    )
    try:
        response = requests.get(cids_url, timeout=30)
        response.raise_for_status()
        cids = response.json().get("IdentifierList", {}).get("CID", [])[:max_hits]
        if not cids:
            return []
        properties_url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{','.join(str(cid) for cid in cids)}/property/Title,CanonicalSMILES/JSON"
        )
        properties_response = requests.get(properties_url, timeout=30)
        properties_response.raise_for_status()
        rows = []
        for item in properties_response.json().get("PropertyTable", {}).get("Properties", []):
            rows.append(
                {
                    "pubchem_cid": item.get("CID"),
                    "Molecule Name": item.get("Title") or name,
                    "SMILES": item.get("CanonicalSMILES"),
                }
            )
        return rows
    except Exception as exc:
        log.warning("PubChem name search failed for '%s': %s", name, exc)
        return []


def compute_rdkit_adme(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute partial ADME descriptors from SMILES when TCMSP OB/DL are unavailable.
    OB and DL remain missing by design.
    """
    chem, descriptors = import_rdkit()
    if chem is None or descriptors is None:
        log.warning("RDKit not installed; skipping RDKit ADME fallback.")
        if "adme_source" not in df.columns:
            df["adme_source"] = "unknown"
        return df

    for column in ["MW", "logP", "HBD", "HBA"]:
        if column not in df.columns:
            df[column] = pd.NA
    if "OB" not in df.columns:
        df["OB"] = pd.NA
    if "DL" not in df.columns:
        df["DL"] = pd.NA
    if "adme_source" not in df.columns:
        df["adme_source"] = "unknown"

    needs_partial = (
        df["SMILES"].notna()
        & df["OB"].isna()
        & df["DL"].isna()
        & ~df["adme_source"].isin(["rdkit_partial"])
    )
    for idx in df.index[needs_partial]:
        smiles = str(df.at[idx, "SMILES"]).strip()
        mol = chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        df.at[idx, "MW"] = float(descriptors.MolWt(mol))
        df.at[idx, "logP"] = float(descriptors.MolLogP(mol))
        df.at[idx, "HBD"] = int(descriptors.NumHDonors(mol))
        df.at[idx, "HBA"] = int(descriptors.NumHAcceptors(mol))
        df.at[idx, "OB"] = pd.NA
        df.at[idx, "DL"] = pd.NA
        df.at[idx, "adme_source"] = "rdkit_partial"
    return df


def validate_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Validate SMILES strings with RDKit when available."""
    chem, _ = import_rdkit()
    if "smiles_valid" not in df.columns:
        df["smiles_valid"] = pd.NA
    if chem is None:
        log.warning("RDKit not installed; SMILES validation skipped.")
        return df

    invalid_count = 0
    for idx, smiles in df["SMILES"].items():
        if pd.isna(smiles) or not str(smiles).strip():
            continue
        is_valid = chem.MolFromSmiles(str(smiles).strip()) is not None
        df.at[idx, "smiles_valid"] = bool(is_valid)
        if not is_valid:
            invalid_count += 1
    if invalid_count:
        log.warning("Invalid SMILES detected: %d row(s)", invalid_count)
    else:
        log.info("All non-empty SMILES parsed successfully.")
    return df


def apply_adme_filter(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Apply config-driven filters without fabricating missing ADME values."""
    filters = cfg["compound_filter"]
    filtered = df.copy()

    if "OB" in filtered.columns and filtered["OB"].notna().any():
        ob_series = pd.to_numeric(filtered["OB"], errors="coerce")
        before = len(filtered)
        filtered = filtered[ob_series >= float(filters["OB_threshold"])].copy()
        log.info("OB >= %s retained %d/%d rows", filters["OB_threshold"], len(filtered), before)
    else:
        log.warning("OB values unavailable; OB threshold not applied.")

    if "DL" in filtered.columns and filtered["DL"].notna().any():
        dl_series = pd.to_numeric(filtered["DL"], errors="coerce")
        before = len(filtered)
        filtered = filtered[dl_series >= float(filters["DL_threshold"])].copy()
        log.info("DL >= %s retained %d/%d rows", filters["DL_threshold"], len(filtered), before)
    else:
        log.warning("DL values unavailable; DL threshold not applied.")

    if "MW" in filtered.columns and filtered["MW"].notna().any():
        mw_series = pd.to_numeric(filtered["MW"], errors="coerce")
        before = len(filtered)
        filtered = filtered[mw_series <= float(filters["MW_max"])].copy()
        log.info("MW <= %s retained %d/%d rows", filters["MW_max"], len(filtered), before)

    if "logP" in filtered.columns and filtered["logP"].notna().any():
        lower, upper = filters.get("logP_range", [-math.inf, math.inf])
        logp_series = pd.to_numeric(filtered["logP"], errors="coerce")
        before = len(filtered)
        filtered = filtered[(logp_series >= float(lower)) & (logp_series <= float(upper))].copy()
        log.info("logP in [%s, %s] retained %d/%d rows", lower, upper, len(filtered), before)

    return filtered.reset_index(drop=True)


def finalize_columns(df: pd.DataFrame, source: str, today: str) -> pd.DataFrame:
    """Populate required columns and order the output contract first."""
    working = df.copy()
    required_defaults = {
        "Molecule Name": pd.NA,
        "SMILES": pd.NA,
        "OB": pd.NA,
        "DL": pd.NA,
        "MW": pd.NA,
        "source": source,
        "query_date": today,
    }
    for column, value in required_defaults.items():
        if column not in working.columns:
            working[column] = value
    working["source"] = working["source"].fillna(source)
    working["query_date"] = working["query_date"].fillna(today)

    contract = ["Molecule Name", "SMILES", "OB", "DL", "MW", "source", "query_date"]
    extras = [column for column in working.columns if column not in contract]
    return working[contract + extras]


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 01: Compound collection and ADME screening")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument(
        "--tcmsp_export",
        default=None,
        help="Path to manually downloaded TCMSP CSV export (optional)",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    today = date.today().isoformat()
    drug_name = cfg["drug"]["name"]

    out_dir = Path("data/processed")
    cache_dir = Path("data/raw/cached")
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    log.info("=== Step 01: Compound collection | Drug: %s ===", drug_name)

    if args.tcmsp_export:
        export_path = Path(args.tcmsp_export)
        if not validate_csv(export_path, min_rows=1, step_name="Step 01"):
            sys.exit(1)

    source_name = "TCMSP"
    if args.tcmsp_export and Path(args.tcmsp_export).exists():
        raw_df = pd.read_csv(args.tcmsp_export)
        raw_df = normalize_tcmsp_columns(raw_df)
        raw_df["source"] = "TCMSP"
        raw_df["query_date"] = today
        raw_df["adme_source"] = "tcmsp"
        cache_path = cache_dir / f"tcmsp_{today}.csv"
        raw_df.to_csv(cache_path, index=False)
        log.info("Cached raw TCMSP export -> %s", cache_path)
    else:
        source_name = "PubChem"
        log.warning("No TCMSP export provided. Falling back to PubChem search.")
        log.warning("For a TCM herb study, provide a TCMSP export to obtain real OB/DL values.")
        rows = pubchem_search_by_name(drug_name)
        if not rows:
            log.error("No PubChem results found for '%s'.", drug_name)
            sys.exit(1)
        raw_df = pd.DataFrame(rows)
        raw_df["source"] = "PubChem"
        raw_df["query_date"] = today
        raw_df["OB"] = pd.NA
        raw_df["DL"] = pd.NA
        raw_df["MW"] = pd.NA
        raw_df["adme_source"] = "unknown"
        cache_path = cache_dir / f"pubchem_{today}.csv"
        raw_df.to_csv(cache_path, index=False)
        log.info("Cached raw PubChem results -> %s", cache_path)

    if "SMILES" not in raw_df.columns:
        raw_df["SMILES"] = pd.NA

    raw_df = compute_rdkit_adme(raw_df)
    raw_df = validate_smiles(raw_df)
    filtered_df = apply_adme_filter(raw_df, cfg)
    filtered_df = finalize_columns(filtered_df, source=source_name, today=today)

    out_path = out_dir / "01_compounds_filtered.csv"
    filtered_df.to_csv(out_path, index=False)
    log.info("Output -> %s (%d row(s))", out_path, len(filtered_df))

    log_data_provenance(
        step="Step 01: Compound collection and ADME screening",
        output_path=out_path,
        sources=sorted(filtered_df["source"].dropna().astype(str).unique().tolist()),
        cfg=cfg,
        extra={
            "n_input": int(len(raw_df)),
            "n_output": int(len(filtered_df)),
            "filters_applied": cfg["compound_filter"],
            "invalid_smiles": int((filtered_df.get("smiles_valid") == False).sum()),  # noqa: E712
        },
    )
    print_step_summary(
        "Step 01",
        out_path,
        cfg,
        next_step="python scripts/python/02_targets.py --config config.yaml",
    )

    if filtered_df.empty:
        log.error("No compounds remained after filtering. Check the raw export and thresholds.")
        sys.exit(1)


if __name__ == "__main__":
    main()
