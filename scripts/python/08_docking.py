#!/usr/bin/env python3
"""
Step 08: Molecular docking.
"""

from __future__ import annotations

import argparse
import logging
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import requests
import yaml

from utils_validate import log_data_provenance, print_step_summary, validate_csv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

KNOWN_PDB = {
    "TNF": ["2AZ5", "1TNF"],
    "IL6": ["1ALU", "4CNI"],
    "IL1B": ["2MIB", "1HIB"],
    "IL6ST": ["3L5H"],
    "JAK1": ["3EYH", "4E5W"],
    "JAK2": ["6E2Q", "2B7A"],
    "JAK3": ["1YVJ", "3LXK"],
    "STAT3": ["6NJS", "3CWG"],
    "PTPN11": ["2SHP", "3B7O"],
    "AKT1": ["4EKL", "3CQW"],
    "MAPK14": ["3HEC", "1A9U"],
    "MAPK1": ["4GT3"],
    "MAPK3": ["2ZOQ"],
    "VEGFA": ["4WPB", "1VPF"],
    "TP53": ["2OCJ", "3ZME"],
    "IL17A": ["5HHV"],
    "TGFB1": ["3KFD"],
    "PTGS2": ["5IKR", "1CX2"],
    "MMP9": ["4H82", "2OVZ"],
    "MMP3": ["1HFS"],
    "CASP3": ["2XYG", "1PAU"],
    "BCL2": ["4LVT", "2W3L"],
    "EGFR": ["4WKQ", "1IVO"],
    "SRC": ["2SRC", "3EL8"],
    "HRAS": ["4EFL", "3L8Y"],
    "HSP90AA1": ["5J64", "2YE4"],
    "MAPK8": ["3O2M"],
    "RELA": ["2RAM"],
    "PIK3CA": ["4JPS"],
    "PTEN": ["1D5R"],
    "CDK4": ["2W9Z"],
    "CDK6": ["2EUF"],
}

DEFAULT_BOX_SIZE = (20.0, 20.0, 20.0)
EXCLUDED_HET_RESIDUES = {
    "HOH",
    "WAT",
    "DOD",
    "NA",
    "CL",
    "CA",
    "MG",
    "ZN",
    "MN",
    "K",
    "CU",
    "CO",
    "NI",
    "FE",
    "SO4",
    "PO4",
    "ACT",
    "ACE",
    "EDO",
    "GOL",
    "PEG",
    "DMS",
    "MES",
    "TRS",
    "FMT",
    "BME",
}


def load_config(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def sanitize_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "item"


def fetch_pdb(pdb_id: str, out_dir: Path) -> Path | None:
    pdb_path = out_dir / f"{pdb_id}.pdb"
    if pdb_path.exists():
        return pdb_path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        pdb_path.write_text(response.text, encoding="utf-8")
        log.info("Downloaded PDB %s -> %s", pdb_id, pdb_path)
        return pdb_path
    except Exception as exc:
        log.warning("Failed to download PDB %s: %s", pdb_id, exc)
        return None


def estimate_binding_center_from_pdb(pdb_path: Path) -> tuple[float, float, float] | None:
    """
    Estimate a docking center from the largest non-solvent HETATM residue.

    Returns None when no suitable co-crystallized ligand is present.
    """
    ligands: dict[tuple[str, str, str], list[tuple[float, float, float]]] = {}
    try:
        with open(pdb_path, encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("HETATM"):
                    continue
                res_name = line[17:20].strip().upper()
                if res_name in EXCLUDED_HET_RESIDUES:
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                chain = line[21].strip() or "_"
                res_seq = line[22:26].strip() or "_"
                key = (res_name, chain, res_seq)
                ligands.setdefault(key, []).append((x, y, z))
    except Exception as exc:
        log.warning("Could not parse PDB coordinates from %s: %s", pdb_path, exc)
        return None

    if not ligands:
        return None

    ligand_atoms = max(ligands.values(), key=len)
    xs, ys, zs = zip(*ligand_atoms)
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def write_receptor_prep_script(pdb_path: Path, receptor_pdbqt: Path, out_dir: Path) -> None:
    """Write receptor preparation instructions for AutoDockTools4."""
    script = out_dir / f"receptor_prep_{pdb_path.stem}.sh"
    script.write_text(
        f"""#!/bin/bash
# Prepare receptor PDBQT for {pdb_path.stem}
# Requires: AutoDockTools4 (MGLTools) in PATH
# Step 1: Remove water and heteroatoms, add hydrogens
prepare_receptor -r {pdb_path} \\
    -o {receptor_pdbqt} \\
    -A hydrogens \\
    -U nphs_lps_waters_deleteAltB
echo "Done: {receptor_pdbqt}"
""",
        encoding="utf-8",
    )
    log.info("Receptor prep script -> %s", script)


def prepare_ligand_pdbqt(smiles: str, compound_name: str, out_dir: Path) -> Path | None:
    pdbqt_path = out_dir / f"{sanitize_name(compound_name)}.pdbqt"
    if pdbqt_path.exists():
        return pdbqt_path

    mol2_path = out_dir / f"{sanitize_name(compound_name)}.mol2"
    try:
        subprocess.run(
            ["obabel", f"-:{smiles}", "--gen3D", "-O", str(mol2_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            ["obabel", str(mol2_path), "-O", str(pdbqt_path), "--partialcharge", "gasteiger"],
            check=True,
            capture_output=True,
            text=True,
        )
        return pdbqt_path
    except FileNotFoundError:
        log.warning("Open Babel not found. Install it to prepare ligand PDBQT files.")
        return None
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else str(exc)
        log.warning("Open Babel failed for %s: %s", compound_name, stderr[:300])
        return None


def run_vina(
    ligand_pdbqt: Path,
    receptor_pdbqt: Path,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    out_dir: Path,
    cfg: dict[str, Any],
) -> tuple[float | None, str]:
    dock_cfg = cfg["docking"]
    out_pdbqt = out_dir / f"{ligand_pdbqt.stem}_{receptor_pdbqt.stem}_out.pdbqt"
    command = [
        "vina",
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(box_size[0]),
        "--size_y",
        str(box_size[1]),
        "--size_z",
        str(box_size[2]),
        "--exhaustiveness",
        str(dock_cfg["exhaustiveness"]),
        "--num_modes",
        str(dock_cfg["num_modes"]),
        "--energy_range",
        str(dock_cfg["energy_range"]),
        "--out",
        str(out_pdbqt),
    ]
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        return None, "vina_not_found"
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else str(exc)
        log.warning("Vina docking failed for %s/%s: %s", ligand_pdbqt.stem, receptor_pdbqt.stem, stderr[:300])
        return None, "vina_failed"

    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1]), "docked"
            except ValueError:
                break
    return None, "vina_output_unparsed"


def write_gromacs_guide(out_path: Path, target: str, ligand: str, cfg: dict[str, Any]) -> None:
    md_cfg = cfg["docking"]["md"]
    out_path.write_text(
        f"""# GROMACS MD Setup Guide
# Target: {target}
# Ligand: {ligand}
# Simulation: {md_cfg['simulation_time_ns']} ns at {md_cfg['temperature_K']} K

gmx pdb2gmx -f receptor.pdb -o receptor.gro -water tip3p -ff amber99sb-ildn
python acpype.py -i {sanitize_name(ligand)}.mol2 -c bcc -a gasteiger
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro
gmx mdrun -v -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -r nvt.gro
gmx mdrun -v -deffnm npt
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
""",
        encoding="utf-8",
    )


def resolve_pdb_choice(gene: str, docking_dir: Path) -> tuple[str | None, Path | None, tuple[float, float, float] | None]:
    pdb_ids = KNOWN_PDB.get(gene, [])
    fallback: tuple[str | None, Path | None, tuple[float, float, float] | None] = (None, None, None)
    for pdb_id in pdb_ids:
        pdb_path = fetch_pdb(pdb_id, docking_dir)
        if pdb_path is None:
            continue
        center = estimate_binding_center_from_pdb(pdb_path)
        if fallback[0] is None:
            fallback = (pdb_id, pdb_path, center)
        if center is not None:
            return pdb_id, pdb_path, center
    return fallback


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 08: Molecular docking")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    docking_dir = Path("results/docking")
    processed_dir = Path("data/processed")
    docking_dir.mkdir(parents=True, exist_ok=True)

    compounds_path = processed_dir / "01_compounds_filtered.csv"
    hub_path = processed_dir / "05_hub_genes.csv"
    if not validate_csv(
        compounds_path,
        required_columns=["Molecule Name", "SMILES"],
        min_rows=1,
        step_name="Step 08",
    ):
        sys.exit(1)
    if not validate_csv(
        hub_path,
        required_columns=["gene_symbol"],
        min_rows=1,
        step_name="Step 08",
    ):
        sys.exit(1)

    compounds = pd.read_csv(compounds_path)
    hub_genes = pd.read_csv(hub_path)
    log.info("=== Step 08: Molecular docking ===")
    log.info("Compounds: %d | Targets: %d", len(compounds), len(hub_genes))

    results: list[dict[str, Any]] = []
    box_size = DEFAULT_BOX_SIZE

    for _, gene_row in hub_genes.iterrows():
        gene = str(gene_row["gene_symbol"]).strip().upper()
        pdb_id, pdb_path, center = resolve_pdb_choice(gene, docking_dir)
        if pdb_id is None or pdb_path is None:
            log.warning("No downloadable PDB structure configured for %s.", gene)

        receptor_pdbqt = docking_dir / f"{pdb_id}_receptor.pdbqt" if pdb_id else docking_dir / "UNKNOWN_receptor.pdbqt"
        if pdb_path and not receptor_pdbqt.exists():
            write_receptor_prep_script(pdb_path, receptor_pdbqt, docking_dir)

        if pdb_id and center is None:
            log.warning(
                "WARNING: No co-crystallized ligand found in %s. Set docking box center manually.",
                pdb_id,
            )

        for _, compound_row in compounds.iterrows():
            compound = str(compound_row.get("Molecule Name", "compound")).strip()
            smiles = compound_row.get("SMILES")
            base_result = {
                "compound": compound,
                "target": gene,
                "pdb_id": pdb_id or "UNKNOWN",
                "binding_energy_kcal_mol": pd.NA,
                "status": "pending",
            }

            if pd.isna(smiles) or not str(smiles).strip():
                base_result["status"] = "missing_smiles"
                results.append(base_result)
                continue
            if pdb_id is None or pdb_path is None:
                base_result["status"] = "pdb_not_available"
                results.append(base_result)
                continue
            if not receptor_pdbqt.exists():
                base_result["status"] = "receptor_not_prepared"
                results.append(base_result)
                continue
            if center is None:
                base_result["status"] = "needs_manual_center"
                results.append(base_result)
                continue

            ligand_pdbqt = prepare_ligand_pdbqt(str(smiles).strip(), compound, docking_dir)
            if ligand_pdbqt is None:
                base_result["status"] = "ligand_not_prepared"
                results.append(base_result)
                continue

            energy, status = run_vina(ligand_pdbqt, receptor_pdbqt, center, box_size, docking_dir, cfg)
            base_result["binding_energy_kcal_mol"] = energy if energy is not None else pd.NA
            base_result["status"] = status
            results.append(base_result)

            if energy is not None and energy <= float(cfg["docking"]["binding_energy_cutoff"]):
                guide_path = docking_dir / f"md_guide_{sanitize_name(compound)}_{gene}.sh"
                write_gromacs_guide(guide_path, gene, compound, cfg)

    results_df = pd.DataFrame(
        results,
        columns=["compound", "target", "pdb_id", "binding_energy_kcal_mol", "status"],
    )
    out_path = processed_dir / "08_docking_scores.csv"
    results_df.to_csv(out_path, index=False)
    log.info("Docking scores -> %s (%d row(s))", out_path, len(results_df))

    status_counts = results_df["status"].value_counts(dropna=False).to_dict() if not results_df.empty else {}
    for status, count in status_counts.items():
        log.info("  %s: %d", status, count)

    log_data_provenance(
        step="Step 08: Molecular docking",
        output_path=out_path,
        sources=["RCSB PDB", cfg["docking"]["software"]],
        cfg=cfg,
        extra={
            "n_input": int(len(compounds) * len(hub_genes)),
            "n_output": int(len(results_df)),
            "filters_applied": {
                "binding_energy_cutoff": cfg["docking"]["binding_energy_cutoff"],
                "box_size": box_size,
            },
        },
    )
    print_step_summary(
        "Step 08",
        out_path,
        cfg,
        next_step="python scripts/python/10_visualization.py --config config.yaml",
    )


if __name__ == "__main__":
    main()
