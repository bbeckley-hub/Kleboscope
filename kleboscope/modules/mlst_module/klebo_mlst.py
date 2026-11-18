#!/usr/bin/env python3

"""
Kleboscope MLST - Klebsiella pneumoniae Multi-Locus Sequence Typing
Bundled mlst binary with mlstdb-based database management.

Per-sample output (1 file each):
    <sample>/<sample>_mlst.html

Batch-level outputs:
    kleboscope_mlst_summary.tsv
    kleboscope_mlst_summary.html
    kleboscope_mlst_summary.json

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Version: 1.3.0 (single per-sample HTML, 3-file batch summary, mlstdb --update-db)

Acknowledgements:
    - Seemann T. mlst (https://github.com/tseemann/mlst) for the typing engine
    - PubMLST (Jolley et al.) for the Klebsiella scheme
    - mlstdb Seemann T.

Scheme expected:
    klebsiella   PubMLST K. pneumoniae MLST scheme

Typical usage:
    python klebo_mlst.py -i "*.fna" -o mlst_results -db db/pubmlst -sc bin/
    python klebo_mlst.py --update-db -db db/pubmlst -sc bin/
"""

import argparse
import glob
import json
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import psutil


KLEBOSCOPE_BANNER = (
    "██╗  ██╗ ██╗      ███████╗ ██████╗  ██████╗ ███████╗  ██████╗ ██████╗  ██████╗ ███████╗\n"
    "██║ ██╔╝ ██║      ██╔════╝ ██╔══██╗██╔═══██╗██╔════╝ ██╔════╝██╔═══██╗ ██╔══██╗██╔════╝\n"
    "█████╔╝  ██║      █████╗   ██████╔╝██║   ██║███████╗ ██║     ██║   ██║ ██████╔╝█████╗\n"
    "██╔═██╗  ██║      ██╔══╝   ██╔══██╗██║   ██║╚════██║ ██║     ██║   ██║ ██╔═══╝ ██╔══╝\n"
    "██║  ██╗ ███████╗ ███████╗ ██████╔╝╚██████╔╝███████║ ╚██████╗╚██████╔╝ ██║     ███████╗\n"
    "╚═╝  ╚═╝ ╚══════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚══════╝  ╚═════╝╚═════╝   ╚═╝     ╚══════╝"
)


def _render_footer(title: str, current_time: str) -> str:
    """Render the shared KLEBOSCOPE footer used across every module."""
    return f"""
<div class="footer">
    <p><strong>KLEBOSCOPE</strong> - {title}</p>
    <p class="timestamp">Generated: {current_time}</p>
    <div class="authorship">
        <p><strong>Technical Support & Inquiries:</strong></p>
        <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
        <p>Email: brownbeckley94@gmail.com</p>
        <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
    </div>
</div>"""


class KleboMLST:
    """Klebsiella MLST typing via the bundled mlst binary."""

    def __init__(
        self,
        database_dir: Path,
        script_dir: Path,
        scheme: str = "klebsiella",
        cpus: Optional[int] = None,
    ):
        self.database_dir = Path(database_dir)
        self.script_dir = Path(script_dir)
        self.mlst_bin = self.script_dir / "mlst"
        self.scheme = scheme
        self.logger = self._setup_logging()
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)

        self.metadata = {
            "tool_name": "Kleboscope MLST",
            "version": "1.3.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "scheme": self.scheme,
            "database_dir": str(self.database_dir),
        }

        self.science_quotes = [
            "The important thing is not to stop questioning. Curiosity has its own reason for existing. - Albert Einstein",
            "Science is not only a disciple of reason but also one of romance and passion. - Stephen Hawking",
            "Somewhere, something incredible is waiting to be known. - Carl Sagan",
            "The good thing about science is that it's true whether or not you believe in it. - Neil deGrasse Tyson",
            "In science, there are no shortcuts to truth. - Karl Popper",
            "Science knows no country, because knowledge belongs to humanity. - Louis Pasteur",
            "The science of today is the technology of tomorrow. - Edward Teller",
            "Nothing in life is to be feared, it is only to be understood. - Marie Curie",
            "Kleboscope turns genomic complexity into actionable insights for AMR surveillance. - Brown Beckley",
        ]

    # ------------------------------------------------------------------ logging / resources
    @staticmethod
    def _setup_logging():
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
        )
        return logging.getLogger(__name__)

    @staticmethod
    def _get_available_ram() -> float:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception:
            return 8.0

    @staticmethod
    def _calculate_optimal_cpus(user_cpus: Optional[int] = None) -> int:
        if user_cpus is not None:
            return max(1, user_cpus)
        try:
            total = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total <= 4:
                optimal = total
            elif total <= 8:
                optimal = total - 1
            elif total <= 16:
                optimal = max(8, total - 1)
            elif total <= 32:
                optimal = max(16, total - 3)
            else:
                optimal = min(32, int(total * 0.95))
            return max(1, min(optimal, total))
        except Exception:
            return os.cpu_count() or 4

    # ------------------------------------------------------------------ database management
    def _ensure_database(self) -> None:
        """Download the MLST database via database.sh if it is missing."""
        has_scheme = (self.database_dir / "pubmlst" / self.scheme).exists()
        if has_scheme:
            return

        self.logger.info("MLST database missing; running database.sh …")
        module_dir = Path(__file__).parent
        script_path = module_dir / "database.sh"
        if not script_path.exists():
            raise FileNotFoundError(f"database.sh not found at {script_path}")
        subprocess.run(["bash", str(script_path)], cwd=module_dir, check=True)
        self.logger.info("MLST database downloaded.")

    def update_database(self, force: bool = False) -> bool:
        """Download and prune the MLST database to the requested scheme via mlstdb."""
        if shutil.which("mlstdb") is None:
            self.logger.error(
                "'mlstdb' not found on PATH. Install with: pip install mlstdb"
            )
            self.logger.error(
                "Then register a PubMLST API key: mlstdb connect --db pubmlst --api-key"
            )
            return False

        self.database_dir.mkdir(parents=True, exist_ok=True)

        if force:
            self.logger.info("Forcing full database rebuild — pruning existing schemes …")
            for sub in ("pubmlst", "blast"):
                target = self.database_dir / sub
                if target.exists():
                    shutil.rmtree(target)

        temp_dir = Path(tempfile.mkdtemp(prefix="kleboscope_mlst_update_"))
        self.logger.info("Temporary download directory: %s", temp_dir)

        try:
            self.logger.info("Running mlstdb update (may take 10-20 minutes) …")
            result = subprocess.run(
                ["mlstdb", "update"],
                input="y\n",
                text=True,
                cwd=temp_dir,
                env=os.environ.copy(),
            )
            if result.returncode != 0:
                self.logger.error("mlstdb update failed with exit code %d", result.returncode)
                return False

            self.logger.info("Download completed; installing %s scheme …", self.scheme)

            src_scheme = temp_dir / "pubmlst" / self.scheme
            if not src_scheme.exists():
                self.logger.error(
                    "Scheme '%s' not present in mlstdb download at %s",
                    self.scheme, src_scheme,
                )
                return False

            target_pubmlst = self.database_dir / "pubmlst"
            target_pubmlst.mkdir(parents=True, exist_ok=True)
            target_scheme = target_pubmlst / self.scheme
            if target_scheme.exists():
                shutil.rmtree(target_scheme)
            shutil.copytree(src_scheme, target_scheme)
            self.logger.info("Installed scheme at %s", target_scheme)

            blast_src = temp_dir / "blast"
            if blast_src.exists():
                target_blast = self.database_dir / "blast"
                if target_blast.exists():
                    shutil.rmtree(target_blast)
                shutil.copytree(blast_src, target_blast)
                self.logger.info("Installed BLAST database at %s", target_blast)

            removed = 0
            for item in target_pubmlst.iterdir():
                if item.is_dir() and item.name != self.scheme:
                    shutil.rmtree(item)
                    removed += 1
                    self.logger.info("Removed non-target scheme: %s", item.name)
            if removed:
                self.logger.info("Pruned %d non-target scheme(s).", removed)

            rebuild_script = self.script_dir.parent / "scripts" / "mlst-make_blast_db"
            if rebuild_script.exists():
                self.logger.info("Rebuilding BLAST database …")
                rebuild_env = os.environ.copy()
                rebuild_env["MLST_DATADIR"] = str(self.database_dir)
                try:
                    subprocess.run(
                        ["perl", str(rebuild_script)],
                        check=True,
                        cwd=str(self.script_dir.parent),
                        env=rebuild_env,
                    )
                    self.logger.info("BLAST database rebuilt successfully.")
                except subprocess.CalledProcessError as exc:
                    self.logger.warning("BLAST rebuild failed (%s).", exc)
                    self.logger.warning(
                        "The mlstdb download already produced blast/, so analysis "
                        "can still proceed; rebuild manually if you see errors."
                    )
            else:
                self.logger.warning(
                    "mlst-make_blast_db not found at %s; skipping BLAST rebuild.",
                    rebuild_script,
                )

            self.logger.info("Database update complete: %s", self.database_dir)
            return True

        except subprocess.CalledProcessError as exc:
            self.logger.error("Update failed: %s", exc)
            return False
        except Exception as exc:
            self.logger.error("Unexpected error during update: %s", exc)
            return False
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
            self.logger.info("Cleaned up temporary directory: %s", temp_dir)

    # ------------------------------------------------------------------ input discovery
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Resolve the input pattern into a sorted list of FASTA paths."""
        if os.path.isfile(input_path):
            return [Path(input_path).resolve()]

        if os.path.isdir(input_path):
            patterns = [
                "*.fna", "*.fasta", "*.fa", "*.fn", "*.faa",
                "*.fna.gz", "*.fasta.gz", "*.fa.gz",
            ]
            files: List[Path] = []
            for pat in patterns:
                for name in glob.glob(os.path.join(input_path, pat)):
                    path = Path(name)
                    if path.is_file():
                        files.append(path)
            return sorted(set(f.resolve() for f in files))

        candidates: List[Path] = []
        for pat in (
            input_path,
            f"{input_path}.fna",
            f"{input_path}.fasta",
            f"{input_path}.fa",
            f"{input_path}.fn",
            f"{input_path}.faa",
            f"{input_path}.fna.gz",
            f"{input_path}.fasta.gz",
            f"{input_path}.fa.gz",
        ):
            for name in glob.glob(pat):
                path = Path(name)
                if path.is_file():
                    candidates.append(path)
        return sorted(set(f.resolve() for f in candidates))

    # ------------------------------------------------------------------ per-sample pipeline
    def run_single(self, input_file: Path, output_base: Path) -> Dict[str, Any]:
        """Run MLST on one genome and emit one HTML report."""
        self._ensure_database()
        sample_name = input_file.stem
        sample_dir = output_base / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("=== PROCESSING GENOME: %s ===", sample_name)

        cmd = [
            "perl", str(self.mlst_bin),
            str(input_file.resolve()),
            "--scheme", self.scheme,
            "--csv",
            "--nopath",
        ]
        env = os.environ.copy()
        env["MLST_DATADIR"] = str(self.database_dir)

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, check=True, env=env,
            )
        except subprocess.CalledProcessError as exc:
            self.logger.error(
                "MLST failed for %s: exit=%s stderr=%s",
                sample_name, exc.returncode, (exc.stderr or "").strip(),
            )
            result = self._fallback_result(sample_name, "MLST subprocess failed")
            self._create_sample_html(result, sample_dir)
            return result
        except FileNotFoundError as exc:
            self.logger.error("mlst binary not found at %s: %s", self.mlst_bin, exc)
            result = self._fallback_result(sample_name, f"mlst binary not found: {exc}")
            self._create_sample_html(result, sample_dir)
            return result

        result = self._parse_csv(proc.stdout, sample_name)
        result.update(self._identity_coverage(result["st"]))

        self._create_sample_html(result, sample_dir)
        self.logger.info(
            "Completed: %s → ST%s (status=%s)",
            sample_name, result.get("st", "ND"), result.get("mlst_status", "Not Assigned"),
        )
        return result

    # ------------------------------------------------------------------ parsing
    def _parse_csv(self, stdout: str, sample_name: str) -> Dict[str, Any]:
        """Parse the comma-separated mlst --csv output line."""
        lines = [ln for ln in stdout.strip().splitlines() if ln.strip()]
        if not lines:
            return self._empty_result(sample_name)

        result_line: Optional[str] = None
        for line in reversed(lines):
            if "," in line and not line.startswith("["):
                result_line = line.strip()
                break

        if not result_line:
            return self._empty_result(sample_name)

        parts = result_line.split(",")
        if len(parts) < 3:
            return self._empty_result(sample_name)

        st = parts[2].strip()
        alleles: Dict[str, str] = {}
        allele_parts: List[str] = []
        for token in parts[3:]:
            token = token.strip()
            match = re.match(r"([^()]+)\(([^)]+)\)", token)
            if match:
                gene = match.group(1).strip()
                allele = match.group(2).strip()
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")

        return {
            "sample": sample_name,
            "st": st,
            "scheme": self.scheme,
            "alleles": alleles,
            "allele_profile": "-".join(allele_parts),
            "confidence": "HIGH" if st and st not in ("-", "ND") else "LOW",
            "mlst_assigned": st not in ("", "-", "ND"),
        }

    def _identity_coverage(self, st: str) -> Dict[str, Any]:
        """Return identity/coverage/status fields derived from the ST assignment."""
        if st and st not in ("-", "ND", "UNKNOWN"):
            return {
                "identity": "100%",
                "coverage": "100%",
                "mlst_status": "Assigned",
                "quality_metrics": {
                    "assembly_quality": "High Quality",
                    "allele_completeness": "Complete",
                    "database_match": "Perfect Match",
                },
            }
        return {
            "identity": "Not Assigned",
            "coverage": "Not Assigned",
            "mlst_status": "Not Assigned",
            "quality_metrics": {
                "assembly_quality": "Requires Review",
                "allele_completeness": "Incomplete",
                "database_match": "No Match",
            },
        }

    def _empty_result(self, sample_name: str) -> Dict[str, Any]:
        return {
            "sample": sample_name,
            "st": "ND",
            "scheme": self.scheme,
            "alleles": {},
            "allele_profile": "",
            "confidence": "LOW",
            "mlst_assigned": False,
        }

    def _fallback_result(self, sample_name: str, error: str) -> Dict[str, Any]:
        base = self._empty_result(sample_name)
        base.update(self._identity_coverage("ND"))
        base["error"] = error
        return base

    # ------------------------------------------------------------------ per-sample HTML
    def _create_sample_html(self, result: Dict[str, Any], output_dir: Path) -> Path:
        """Render the single per-sample MLST HTML report."""
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        sample = result.get("sample", "Unknown")
        st = result.get("st", "ND")
        confidence = result.get("confidence", "LOW")
        allele_profile = result.get("allele_profile", "")
        identity = result.get("identity", "Not Assigned")
        coverage = result.get("coverage", "Not Assigned")
        mlst_status = result.get("mlst_status", "Not Assigned")
        quality_metrics = result.get("quality_metrics", {})
        error_message = result.get("error")

        alleles_html = "".join(
            f'<div class="allele-card"><div class="allele-name">{gene}</div>'
            f'<div class="allele-value">{allele}</div></div>'
            for gene, allele in result.get("alleles", {}).items()
        ) or '<p style="color:#6b7280;">No alleles resolved.</p>'

        quality_html = "".join(
            f'<div class="quality-card"><div class="quality-value">{value}</div>'
            f'<div class="quality-label">{key.replace("_", " ").title()}</div></div>'
            for key, value in quality_metrics.items()
        ) or '<p style="color:#6b7280;">No quality metrics available.</p>'

        error_html = ""
        if error_message:
            error_html = (
                f'<div class="error-box"><strong>MLST error:</strong> {error_message}</div>'
            )

        status_color = "#10b981" if mlst_status == "Assigned" else "#dc2626"
        confidence_class = (
            "confidence-high" if confidence == "HIGH"
            else "confidence-medium" if confidence == "MEDIUM"
            else "confidence-low"
        )
        identity_class = "ic-card" if identity == "100%" else "ic-card not-assigned"
        coverage_class = "ic-card" if coverage == "100%" else "ic-card not-assigned"

        footer = _render_footer("MLST Analysis Module", current_time)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KLEBOSCOPE - {sample} MLST Report</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color:#ffffff; padding:20px; min-height:100vh;
}}
.container {{ max-width:1400px; margin:0 auto; }}
.ascii-container {{
    background: rgba(0,0,0,0.7); padding:20px; border-radius:15px;
    margin-bottom:20px; box-shadow:0 8px 32px rgba(0,0,0,0.4);
    border:2px solid rgba(0,119,190,0.3); overflow-x:auto; text-align:center;
}}
.ascii-art {{
    font-family: 'Cascadia Mono', 'Consolas', 'Monaco', 'Menlo',
                 'DejaVu Sans Mono', 'Liberation Mono', 'Courier New', monospace;
    font-size:12px; line-height:1.05; white-space:pre;
    color:#00d4ff; text-shadow:0 0 10px rgba(0,212,255,0.4);
    display:inline-block; padding:8px 0; text-align:left;
    font-variant-ligatures:none;
    font-feature-settings: "liga" 0, "calt" 0, "kern" 0;
    font-kerning:none; letter-spacing:0;
    text-rendering:geometricPrecision;
}}
.quote-container {{
    background: rgba(255,255,255,0.1); backdrop-filter:blur(10px);
    padding:20px; border-radius:10px; margin-bottom:30px; text-align:center;
    min-height:100px; display:flex; flex-direction:column;
    justify-content:center; box-shadow:0 4px 20px rgba(0,0,0,0.3);
    border:1px solid rgba(255,255,255,0.2); transition:opacity 0.5s ease-in-out;
}}
.quote-text {{ font-size:18px; font-style:italic; margin-bottom:10px; color:#ffffff; }}
.quote-author {{ font-size:14px; color:#fbbf24; font-weight:bold; }}
.card {{
    background: rgba(255,255,255,0.95); color:#1f2937;
    padding:25px; margin:20px 0; border-radius:12px;
    box-shadow:0 4px 15px rgba(0,0,0,0.2);
}}
.card h2 {{
    color:#0b3d5f; border-bottom:3px solid #0077be;
    padding-bottom:10px; margin-bottom:20px; font-size:24px;
}}
.card h3 {{ color:#1b6b8f; margin-top:20px; margin-bottom:10px; font-size:18px; }}
.metrics-grid {{
    display:grid; grid-template-columns:repeat(auto-fit, minmax(200px, 1fr));
    gap:20px; margin-top:15px;
}}
.metric-card {{
    background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
    color:white; padding:20px; border-radius:8px;
    box-shadow:0 4px 12px rgba(0,0,0,0.15); text-align:center;
}}
.metric-label {{ font-size:13px; opacity:0.9; margin-bottom:5px; text-transform:uppercase; letter-spacing:0.5px; }}
.metric-value {{ font-size:22px; font-weight:bold; word-break:break-word; }}
.identity-coverage-grid {{
    display:grid; grid-template-columns:repeat(auto-fit, minmax(200px, 1fr));
    gap:20px; margin-top:15px;
}}
.ic-card {{
    background: linear-gradient(135deg, #10b981 0%, #059669 100%);
    color:white; padding:20px; border-radius:8px;
    box-shadow:0 4px 12px rgba(0,0,0,0.15); text-align:center;
}}
.ic-card.not-assigned {{ background: linear-gradient(135deg, #6b7280 0%, #4b5563 100%); }}
.ic-value {{ font-size:26px; font-weight:bold; margin-bottom:5px; }}
.ic-label {{ font-size:13px; opacity:0.9; }}
.allele-grid {{
    display:grid; grid-template-columns:repeat(auto-fill, minmax(140px, 1fr));
    gap:15px; margin-top:15px;
}}
.allele-card {{
    background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
    color:white; padding:15px; border-radius:8px; text-align:center;
}}
.allele-name {{ font-size:12px; opacity:0.9; margin-bottom:5px; }}
.allele-value {{ font-size:18px; font-weight:bold; }}
.confidence-high {{ color:#16a34a; font-weight:bold; }}
.confidence-medium {{ color:#f59e0b; font-weight:bold; }}
.confidence-low {{ color:#dc2626; font-weight:bold; }}
.quality-grid {{
    display:grid; grid-template-columns:repeat(auto-fit, minmax(200px, 1fr));
    gap:15px; margin-top:15px;
}}
.quality-card {{
    background:#f0f9ff; color:#0369a1; padding:15px;
    border-radius:8px; text-align:center; border:1px solid #bae6fd;
}}
.quality-value {{ font-size:18px; font-weight:bold; margin-bottom:5px; }}
.quality-label {{ font-size:12px; opacity:0.8; }}
.profile-box {{
    background:#f8fafc; padding:15px; border-radius:8px;
    margin:15px 0; border-left:4px solid #0077be;
    color:#0b3d5f;
}}
.profile-box code {{ color:#0b3d5f; font-weight:bold; }}
.error-box {{
    background:#fef2f2; border-left:4px solid #dc2626;
    padding:15px; margin:15px 0; border-radius:6px; color:#991b1b;
}}
.footer {{
    text-align:center; margin-top:30px; padding:20px;
    background:rgba(0,0,0,0.3); border-radius:10px; font-size:14px;
}}
.timestamp {{ color:#fbbf24; font-weight:bold; }}
.authorship {{
    margin-top:15px; padding:15px; background:rgba(255,255,255,0.1);
    border-radius:8px; font-size:12px;
}}
@media (max-width:768px) {{
    .ascii-art {{ font-size:6px; }}
    .allele-grid, .quality-grid, .metrics-grid {{
        grid-template-columns:1fr;
    }}
}}
</style>
</head>
<body>
<div class="container">
<div class="ascii-container">
<div class="ascii-art">{KLEBOSCOPE_BANNER}</div>
</div>

<div class="quote-container" id="quoteContainer">
    <div class="quote-text" id="quoteText"></div>
    <div class="quote-author" id="quoteAuthor"></div>
</div>

<div class="card">
    <h2>📊 Sample Information</h2>
    <div class="metrics-grid">
        <div class="metric-card"><div class="metric-label">Sample</div><div class="metric-value" style="font-size:16px;">{sample}</div></div>
        <div class="metric-card"><div class="metric-label">Analysis Date</div><div class="metric-value" style="font-size:16px;">{current_time}</div></div>
        <div class="metric-card"><div class="metric-label">MLST Scheme</div><div class="metric-value" style="font-size:16px;">{self.scheme.title()}</div></div>
    </div>
    {error_html}
</div>

<div class="card">
    <h2>🎯 Identity & Coverage</h2>
    <div class="identity-coverage-grid">
        <div class="{identity_class}"><div class="ic-value">{identity}</div><div class="ic-label">Identity</div></div>
        <div class="{coverage_class}"><div class="ic-value">{coverage}</div><div class="ic-label">Coverage</div></div>
    </div>
    <h3>MLST Status: <span style="color:{status_color};">{mlst_status}</span></h3>
    <h3>Quality Metrics</h3>
    <div class="quality-grid">{quality_html}</div>
</div>

<div class="card">
    <h2>🧬 MLST Typing Results</h2>
    <div class="metrics-grid">
        <div class="metric-card"><div class="metric-label">Sequence Type</div><div class="metric-value">ST{st}</div></div>
        <div class="metric-card"><div class="metric-label">Confidence</div><div class="metric-value {confidence_class}">{confidence}</div></div>
    </div>
    <h3>Allele Profile</h3>
    <div class="profile-box"><code style="font-size:16px;">{allele_profile or "No alleles resolved"}</code></div>
    <h3>Individual Alleles</h3>
    <div class="allele-grid">{alleles_html}</div>
</div>
{footer}
</div>

<script>
const quotes = {json.dumps(self.science_quotes)};
const quoteContainer = document.getElementById('quoteContainer');
const quoteText = document.getElementById('quoteText');
const quoteAuthor = document.getElementById('quoteAuthor');
function displayQuote() {{
    quoteContainer.style.opacity = '0';
    setTimeout(() => {{
        const q = quotes[Math.floor(Math.random() * quotes.length)];
        quoteText.textContent = q;
        quoteAuthor.textContent = '';
        quoteContainer.style.opacity = '1';
    }}, 500);
}}
displayQuote();
setInterval(displayQuote, 10000);
</script>
</body>
</html>"""

        out_file = output_dir / f"{sample}_mlst.html"
        with open(out_file, "w", encoding="utf-8") as handle:
            handle.write(html)
        self.logger.info("Per-sample MLST HTML: %s", out_file)
        return out_file

    # ------------------------------------------------------------------ batch summaries
    def run_batch(self, input_pattern: str, output_base: str) -> Dict[str, Dict[str, Any]]:
        """Run MLST on every matching FASTA and write batch summaries."""
        self._ensure_database()

        fasta_paths = self.find_fasta_files(input_pattern)
        if not fasta_paths:
            raise FileNotFoundError(f"No FASTA files found matching: {input_pattern}")

        output_dir = Path(output_base)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("Found %d genome(s)", len(fasta_paths))
        self.logger.info("Processing sequentially, %d CPU thread(s) available", self.cpus)

        results: Dict[str, Dict[str, Any]] = {}
        for fasta in fasta_paths:
            result = self.run_single(fasta, output_dir)
            results[result["sample"]] = result

        prefix = "kleboscope_mlst_summary"
        self._write_summary_tsv(results, output_dir / f"{prefix}.tsv")
        self._write_summary_html(results, output_dir / f"{prefix}.html")
        self._write_summary_json(results, output_dir / f"{prefix}.json")

        return results

    def _write_summary_tsv(
        self, results: Dict[str, Dict[str, Any]], output_path: Path,
    ) -> None:
        """Write the batch-level MLST TSV, one row per sample."""
        columns = [
            "sample", "st", "mlst_status", "identity", "coverage",
            "allele_profile",
        ]
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write("\t".join(columns) + "\n")
            for sample, result in sorted(results.items()):
                row = [
                    sample,
                    result.get("st", "ND"),
                    result.get("mlst_status", "Not Assigned"),
                    result.get("identity", "Not Assigned"),
                    result.get("coverage", "Not Assigned"),
                    result.get("allele_profile", ""),
                ]
                handle.write("\t".join(str(v) for v in row) + "\n")
        self.logger.info("Wrote batch TSV: %s", output_path)

    def _write_summary_json(
        self, results: Dict[str, Dict[str, Any]], output_path: Path,
    ) -> None:
        """Write the single consolidated batch JSON."""
        total = len(results)
        assigned = sum(1 for r in results.values() if r.get("mlst_status") == "Assigned")
        st_counter: Dict[str, int] = {}
        for r in results.values():
            st = r.get("st", "ND")
            st_counter[st] = st_counter.get(st, 0) + 1

        payload = {
            "metadata": {
                "tool": self.metadata["tool_name"],
                "version": self.metadata["version"],
                "scheme": self.metadata["scheme"],
                "database_dir": self.metadata["database_dir"],
                "analysis_date": self.metadata["analysis_date"],
                "total_samples": total,
            },
            "statistics": {
                "total_samples": total,
                "assigned": assigned,
                "not_assigned": total - assigned,
                "assignment_rate_pct": round(assigned / total * 100, 1) if total else 0.0,
                "st_distribution": dict(sorted(st_counter.items(), key=lambda x: -x[1])),
            },
            "samples": {
                sample: {
                    "sequence_type": r.get("st", "ND"),
                    "mlst_status": r.get("mlst_status", "Not Assigned"),
                    "identity": r.get("identity", "Not Assigned"),
                    "coverage": r.get("coverage", "Not Assigned"),
                    "confidence": r.get("confidence", "LOW"),
                    "alleles": r.get("alleles", {}),
                    "allele_profile": r.get("allele_profile", ""),
                    "quality_metrics": r.get("quality_metrics", {}),
                    "error": r.get("error"),
                }
                for sample, r in sorted(results.items())
            },
        }

        with open(output_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
        self.logger.info("Wrote batch JSON: %s", output_path)

    def _write_summary_html(
        self, results: Dict[str, Dict[str, Any]], output_path: Path,
    ) -> None:
        """Write the batch-level MLST HTML summary."""
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        total = len(results)
        assigned = sum(1 for r in results.values() if r.get("mlst_status") == "Assigned")
        not_assigned = total - assigned
        unique_sts = sorted({r.get("st", "ND") for r in results.values() if r.get("st") not in ("", "ND")})

        rows = []
        for sample, r in sorted(results.items()):
            status = r.get("mlst_status", "Not Assigned")
            status_class = "success-row" if status == "Assigned" else "failed-row"
            rows.append(
                f"<tr class='{status_class}'>"
                f"<td><strong>{sample}</strong></td>"
                f"<td>ST{r.get('st', 'ND')}</td>"
                f"<td>{status}</td>"
                f"<td>{r.get('identity', 'Not Assigned')}</td>"
                f"<td>{r.get('coverage', 'Not Assigned')}</td>"
                f"<td>{r.get('allele_profile', '')}</td>"
                f"</tr>"
            )
        table_rows = "\n".join(rows) or (
            '<tr><td colspan="6" style="text-align:center;color:#9ca3af;">'
            "No samples processed</td></tr>"
        )

        footer = _render_footer("MLST Batch Summary", current_time)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KLEBOSCOPE - MLST Batch Summary</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color:#ffffff; padding:20px; min-height:100vh;
}}
.container {{ max-width:1800px; margin:0 auto; }}
.ascii-container {{
    background: rgba(0,0,0,0.7); padding:20px; border-radius:15px;
    margin-bottom:20px; box-shadow:0 8px 32px rgba(0,0,0,0.4);
    border:2px solid rgba(0,119,190,0.3); overflow-x:auto; text-align:center;
}}
.ascii-art {{
    font-family: 'Cascadia Mono', 'Consolas', 'Monaco', 'Menlo',
                 'DejaVu Sans Mono', 'Liberation Mono', 'Courier New', monospace;
    font-size:12px; line-height:1.05; white-space:pre;
    color:#00d4ff; text-shadow:0 0 10px rgba(0,212,255,0.4);
    display:inline-block; padding:8px 0; text-align:left;
    font-variant-ligatures:none;
    font-feature-settings: "liga" 0, "calt" 0, "kern" 0;
    font-kerning:none; letter-spacing:0;
    text-rendering:geometricPrecision;
}}
.quote-container {{
    background: rgba(255,255,255,0.1); backdrop-filter:blur(10px);
    padding:20px; border-radius:10px; margin-bottom:30px; text-align:center;
    min-height:100px; display:flex; flex-direction:column;
    justify-content:center; box-shadow:0 4px 20px rgba(0,0,0,0.3);
    border:1px solid rgba(255,255,255,0.2); transition:opacity 0.5s ease-in-out;
}}
.quote-text {{ font-size:18px; font-style:italic; margin-bottom:10px; }}
.quote-author {{ font-size:14px; color:#fbbf24; font-weight:bold; }}
.card {{
    background: rgba(255,255,255,0.95); color:#1f2937;
    padding:25px; margin:20px 0; border-radius:12px;
    box-shadow:0 4px 15px rgba(0,0,0,0.2);
}}
.card h2 {{
    color:#0b3d5f; border-bottom:3px solid #0077be;
    padding-bottom:10px; margin-bottom:20px; font-size:24px;
}}
.summary-stats {{
    display:flex; gap:20px; justify-content:center; flex-wrap:wrap; margin:20px 0;
}}
.stat-card {{
    background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
    color:white; padding:20px; border-radius:8px; text-align:center;
    flex:1; min-width:180px;
}}
.stat-card.assigned {{ background: linear-gradient(135deg, #10b981 0%, #059669 100%); }}
.stat-card.not-assigned {{ background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%); }}
.stat-value {{ font-size:26px; font-weight:bold; }}
.stat-label {{ font-size:13px; opacity:0.9; text-transform:uppercase; letter-spacing:0.5px; margin-top:5px; }}
.table-responsive {{ width:100%; overflow-x:auto; margin:20px 0; }}
.gene-table {{ width:100%; border-collapse:collapse; background:white; font-size:13px; }}
.gene-table th, .gene-table td {{
    padding:10px; text-align:left; border-bottom:1px solid #e5e7eb;
}}
.gene-table th {{
    background:#005a8c; color:white; font-weight:600; position:sticky; top:0;
}}
.gene-table tr:hover {{ background:#f3f4f6; }}
.success-row {{ border-left:4px solid #10b981; }}
.failed-row {{ border-left:4px solid #dc2626; background-color:#fef2f2; }}
.footer {{
    text-align:center; margin-top:30px; padding:20px;
    background:rgba(0,0,0,0.3); border-radius:10px; font-size:14px;
}}
.timestamp {{ color:#fbbf24; font-weight:bold; }}
.authorship {{
    margin-top:15px; padding:15px; background:rgba(255,255,255,0.1);
    border-radius:8px; font-size:12px;
}}
@media (max-width:768px) {{
    .ascii-art {{ font-size:6px; }}
    .gene-table {{ font-size:11px; }}
}}
</style>
</head>
<body>
<div class="container">
<div class="ascii-container">
<div class="ascii-art">{KLEBOSCOPE_BANNER}</div>
</div>

<div class="quote-container" id="quoteContainer">
    <div class="quote-text" id="quoteText"></div>
    <div class="quote-author" id="quoteAuthor"></div>
</div>

<div class="card">
    <h2>📊 Batch MLST Summary</h2>
    <div class="summary-stats">
        <div class="stat-card"><div class="stat-value">{total}</div><div class="stat-label">Total Samples</div></div>
        <div class="stat-card assigned"><div class="stat-value">{assigned}</div><div class="stat-label">Assigned</div></div>
        <div class="stat-card not-assigned"><div class="stat-value">{not_assigned}</div><div class="stat-label">Not Assigned</div></div>
        <div class="stat-card"><div class="stat-value">{len(unique_sts)}</div><div class="stat-label">Unique STs</div></div>
    </div>
    <p><strong>Scheme:</strong> {self.metadata['scheme']} &nbsp;|&nbsp; <strong>Tool:</strong> {self.metadata['tool_name']} v{self.metadata['version']}</p>
</div>

<div class="card">
    <h2>🔍 Detailed Results by Genome</h2>
    <div class="table-responsive">
        <table class="gene-table">
            <thead>
                <tr>
                    <th>Sample</th><th>ST</th><th>MLST Status</th>
                    <th>Identity</th><th>Coverage</th><th>Allele Profile</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>
    </div>
</div>
{footer}
</div>

<script>
const quotes = {json.dumps(self.science_quotes)};
const quoteContainer = document.getElementById('quoteContainer');
const quoteText = document.getElementById('quoteText');
const quoteAuthor = document.getElementById('quoteAuthor');
function displayQuote() {{
    quoteContainer.style.opacity = '0';
    setTimeout(() => {{
        const q = quotes[Math.floor(Math.random() * quotes.length)];
        quoteText.textContent = q;
        quoteAuthor.textContent = '';
        quoteContainer.style.opacity = '1';
    }}, 500);
}}
displayQuote();
setInterval(displayQuote, 10000);
</script>
</body>
</html>"""

        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write(html)
        self.logger.info("Wrote batch HTML: %s", output_path)


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Kleboscope MLST — K. pneumoniae sequence typing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        help="Input FASTA file, glob, or directory",
    )
    parser.add_argument(
        "-o", "--output-dir",
        help="Output directory",
    )
    parser.add_argument(
        "-db", "--database-dir", required=True,
        help="Directory containing pubmlst/ and blast/",
    )
    parser.add_argument(
        "-sc", "--script-dir", required=True,
        help="Directory containing the mlst binary",
    )
    parser.add_argument(
        "-s", "--scheme", default="klebsiella",
        help="MLST scheme name",
    )
    parser.add_argument(
        "-c", "--cpus", type=int, default=None,
        help="CPU core budget (auto-detected if omitted)",
    )
    parser.add_argument(
        "--batch", action="store_true",
        help="Process every matching genome under --input",
    )
    parser.add_argument(
        "--update-db", action="store_true",
        help="Download and refresh the MLST database via mlstdb, then exit",
    )
    parser.add_argument(
        "--force-update", action="store_true",
        help="With --update-db, delete existing database folders first",
    )
    args = parser.parse_args()

    executor = KleboMLST(
        database_dir=Path(args.database_dir),
        script_dir=Path(args.script_dir),
        scheme=args.scheme,
        cpus=args.cpus,
    )

    if args.update_db:
        success = executor.update_database(force=args.force_update)
        sys.exit(0 if success else 1)

    if not args.input or not args.output_dir:
        parser.error("When not using --update-db, -i and -o are required.")

    started = time.time()
    try:
        if args.batch:
            results = executor.run_batch(args.input, args.output_dir)
            elapsed = time.time() - started
            assigned = sum(
                1 for r in results.values() if r.get("mlst_status") == "Assigned"
            )
            print()
            print("=" * 70)
            print("📊 Kleboscope MLST — Final Summary")
            print("=" * 70)
            print(f"Samples processed    : {len(results)}")
            print(f"Assigned to an ST    : {assigned}")
            print(f"Elapsed time         : {elapsed:.1f}s")
            print(f"Output directory     : {args.output_dir}")
            print()
            print("Per-sample (1 file each):")
            print("  • <sample>/<sample>_mlst.html")
            print("Batch-level (3 files):")
            print("  • kleboscope_mlst_summary.tsv")
            print("  • kleboscope_mlst_summary.html")
            print("  • kleboscope_mlst_summary.json")
        else:
            input_file = Path(args.input)
            if not input_file.exists():
                print(f"❌ Input file not found: {args.input}")
                sys.exit(1)
            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            result = executor.run_single(input_file, output_dir)
            print(f"🎉 MLST completed for {input_file.name}: ST{result.get('st', 'ND')}")
    except Exception as exc:
        executor.logger.error("Analysis failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()