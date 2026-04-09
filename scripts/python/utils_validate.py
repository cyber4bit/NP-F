#!/usr/bin/env python3
"""
Shared validation utilities for NP pipeline.
Used by individual step scripts to check inputs/outputs and detect placeholders.

Usage in scripts:
    from utils_validate import validate_csv, check_placeholders, log_data_provenance
"""

import logging
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger("np_validate")


def load_config(path: str = "config.yaml") -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_validation_rules(cfg: dict) -> dict:
    return cfg.get("validation_rules", {})


# ── Placeholder detection ────────────────────────────────────

def check_placeholders(df: pd.DataFrame, cfg: dict | None = None) -> list[str]:
    """
    Scan DataFrame for placeholder markers.
    Returns list of warning messages (empty = clean).
    """
    markers = ["placeholder", "YYYY-MM-DD", "UNKNOWN", "N/A"]
    if cfg:
        markers = cfg.get("validation_rules", {}).get("placeholder_markers", markers)

    warnings = []
    for col in df.columns:
        if df[col].dtype != object:
            continue
        for marker in markers:
            mask = df[col].str.contains(marker, case=False, na=False)
            count = mask.sum()
            if count > 0:
                warnings.append(
                    f"Column '{col}': {count} rows contain placeholder marker '{marker}'"
                )
    return warnings


def is_placeholder_result(csv_path: str | Path, cfg: dict | None = None) -> bool:
    """Quick check: does this CSV contain any placeholder data?"""
    path = Path(csv_path)
    if not path.exists():
        return True  # missing file is treated as placeholder
    try:
        df = pd.read_csv(path)
        return len(check_placeholders(df, cfg)) > 0
    except Exception:
        return True


# ── CSV validation ───────────────────────────────────────────

def validate_csv(
    path: str | Path,
    required_columns: list[str] | None = None,
    min_rows: int = 1,
    step_name: str = "",
) -> bool:
    """
    Validate a CSV file exists, has required columns, and meets minimum row count.
    Returns True if valid, logs errors and returns False otherwise.
    """
    path = Path(path)
    prefix = f"[{step_name}] " if step_name else ""

    if not path.exists():
        log.error(f"{prefix}File not found: {path}")
        return False

    try:
        df = pd.read_csv(path)
    except Exception as e:
        log.error(f"{prefix}Cannot read {path}: {e}")
        return False

    if len(df) < min_rows:
        log.error(f"{prefix}{path.name}: {len(df)} rows (minimum: {min_rows})")
        return False

    if required_columns:
        # Flexible matching: check if any column contains the keyword
        missing = []
        for req in required_columns:
            found = any(req.lower() in c.lower() for c in df.columns)
            if not found:
                missing.append(req)
        if missing:
            log.error(f"{prefix}{path.name}: missing columns {missing}")
            log.info(f"{prefix}Available columns: {list(df.columns)}")
            return False

    log.info(f"{prefix}{path.name}: OK ({len(df)} rows, {len(df.columns)} cols)")
    return True


# ── Data provenance logging ──────────────────────────────────

def log_data_provenance(
    step: str,
    output_path: str | Path,
    sources: list[str],
    cfg: dict,
    extra: dict | None = None,
):
    """
    Write a provenance sidecar file (.provenance.yaml) alongside an output CSV.
    Records what data went in, when, and with what parameters.
    """
    path = Path(output_path)
    prov_path = path.with_suffix(".provenance.yaml")

    prov = {
        "step": step,
        "output_file": path.name,
        "generated_date": date.today().isoformat(),
        "drug": cfg.get("drug", {}).get("name", "?"),
        "disease": cfg.get("disease", {}).get("name", "?"),
        "sources": sources,
        "config_snapshot": {
            "seed": cfg.get("output", {}).get("seed_global"),
        },
    }
    if extra:
        prov.update(extra)

    with open(prov_path, "w", encoding="utf-8") as f:
        yaml.dump(prov, f, default_flow_style=False, allow_unicode=True)

    log.info(f"Provenance → {prov_path}")


# ── Step output summary ──────────────────────────────────────

def print_step_summary(
    step: str,
    output_path: str | Path,
    cfg: dict,
    next_step: str = "",
):
    """Print standardized step completion summary with placeholder check."""
    path = Path(output_path)
    print("=" * 60)

    if not path.exists():
        print(f"  [{step}] OUTPUT MISSING: {path}")
        print("=" * 60)
        return

    df = pd.read_csv(path)
    placeholders = check_placeholders(df, cfg)

    print(f"  [{step}] Output: {path.name} ({len(df)} rows)")

    if placeholders:
        print(f"  ⚠ PLACEHOLDER DATA DETECTED:")
        for w in placeholders:
            print(f"    - {w}")
        print(f"  ⚠ This output is NOT suitable for publication.")
    else:
        print(f"  ✓ No placeholders detected.")

    if next_step:
        print(f"  NEXT: {next_step}")
    print("=" * 60)
