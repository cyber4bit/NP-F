#!/usr/bin/env python3
"""
Generate a post-run quality report for the network pharmacology pipeline.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from utils_validate import check_placeholders


@dataclass
class StepReport:
    step: str
    file_name: str
    rows: str
    sources: str
    placeholders: str
    publishable: str


def load_config(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def read_provenance(csv_path: Path) -> dict[str, Any] | None:
    prov_path = csv_path.with_suffix(".provenance.yaml")
    if not prov_path.exists():
        return None
    try:
        with open(prov_path, encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except Exception:
        return None


def infer_sources(csv_path: Path, df: pd.DataFrame, provenance: dict[str, Any] | None) -> str:
    if provenance and provenance.get("sources"):
        return ",".join(str(source) for source in provenance["sources"])
    if "source" in df.columns:
        values = sorted({str(value) for value in df["source"].dropna().astype(str) if value})
        return ",".join(values) if values else "-"
    return "-"


def detect_placeholder_note(step: str, df: pd.DataFrame, cfg: dict[str, Any]) -> tuple[str, str]:
    warnings = check_placeholders(df, cfg)
    status_note = "None"
    publishable = "YES"

    if warnings:
        status_note = "DETECTED"
        publishable = "NO"

    if "status" in df.columns:
        statuses = sorted({str(value) for value in df["status"].dropna().astype(str)})
        success_statuses = {"docked", "ok"}
        non_success = [status for status in statuses if status not in success_statuses]
        if non_success:
            status_note = f"status:{','.join(non_success[:3])}"
            publishable = "PARTIAL" if step == "08" else "NO"

    return status_note, publishable


def minimum_rows(cfg: dict[str, Any], step: str) -> int:
    rules = cfg.get("validation_rules", {})
    mapping = {
        "01": int(rules.get("min_compounds", 1)),
        "02": int(rules.get("min_drug_targets", 1)),
        "03": int(rules.get("min_disease_targets", 1)),
        "04": int(rules.get("min_intersection", 1)),
        "05": int(rules.get("min_hub_genes", 1)),
        "06": 1,
        "07": 1,
        "08": 1,
    }
    return mapping.get(step, 1)


def evaluate_step(step: str, csv_path: Path, cfg: dict[str, Any]) -> StepReport:
    if not csv_path.exists():
        return StepReport(step, csv_path.name, "NOT RUN", "-", "not run", "NOT RUN")

    try:
        df = pd.read_csv(csv_path)
    except Exception as exc:
        return StepReport(step, csv_path.name, "ERROR", "-", f"read_error:{exc}", "NO")

    provenance = read_provenance(csv_path)
    sources = infer_sources(csv_path, df, provenance)
    placeholder_note, publishable = detect_placeholder_note(step, df, cfg)

    row_count = len(df)
    if row_count < minimum_rows(cfg, step):
        publishable = "NO" if publishable == "YES" else publishable

    if step == "08" and publishable == "PARTIAL":
        docked_rows = 0
        if "status" in df.columns:
            docked_rows = int((df["status"].astype(str) == "docked").sum())
        if docked_rows == 0:
            publishable = "PARTIAL"

    return StepReport(step, csv_path.name, str(row_count), sources or "-", placeholder_note, publishable)


def build_report(cfg: dict[str, Any]) -> tuple[str, list[StepReport]]:
    processed_dir = Path("data/processed")
    steps = [
        ("01", processed_dir / "01_compounds_filtered.csv"),
        ("02", processed_dir / "02_targets_merged.csv"),
        ("03", processed_dir / "03_disease_targets.csv"),
        ("04", processed_dir / "04_intersection_genes.csv"),
        ("05", processed_dir / "05_hub_genes.csv"),
        ("06", processed_dir / "06_enrichment_kegg.csv"),
        ("07", processed_dir / "07_mr_results.csv"),
        ("08", processed_dir / "08_docking_scores.csv"),
    ]
    reports = [evaluate_step(step, path, cfg) for step, path in steps]

    publishable_steps = [report.step for report in reports if report.publishable == "YES"]
    fix_steps = [report.step for report in reports if report.publishable in {"NO", "PARTIAL"}]

    drug_display = cfg["drug"]["name"]
    if cfg["drug"].get("chinese_name"):
        drug_display = f"{drug_display} ({cfg['drug']['chinese_name']})"
    disease_display = cfg["disease"]["name"]
    if cfg["disease"].get("abbreviation"):
        disease_display = f"{disease_display} ({cfg['disease']['abbreviation']})"

    lines = [
        "=== Network Pharmacology Pipeline Quality Report ===",
        f"Drug:    {drug_display}",
        f"Disease: {disease_display}",
        f"Run date: {date.today().isoformat()}",
        "",
        "STEP | FILE | ROWS | SOURCES | PLACEHOLDERS | PUBLISHABLE",
    ]
    for report in reports:
        lines.append(
            f"{report.step:<4} | {report.file_name} | {report.rows} | "
            f"{report.sources} | {report.placeholders} | {report.publishable}"
        )

    summary = f"SUMMARY: {len(publishable_steps)}/{len(reports)} steps publishable."
    if fix_steps:
        summary += f" Fix steps: {', '.join(fix_steps)}"
    lines.extend(["", summary])
    return "\n".join(lines), reports


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a post-run quality report")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    report_text, _ = build_report(cfg)

    out_dir = Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "quality_report.txt"
    out_path.write_text(report_text + "\n", encoding="utf-8")
    print(report_text)


if __name__ == "__main__":
    main()
