#!/usr/bin/env python3
"""
Kleboscope Kaptive - K. pneumoniae Capsule (K) and Lipopolysaccharide (O) Locus Typing
Kaptive 3 subcommand interface, Prodigal-free, TSV-native.

Per-sample output (1 file each):
    <sample>/<sample>_kaptive.html

Batch-level outputs:
    <prefix>_summary.tsv
    <prefix>_summary.html
    <prefix>_summary.json

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry

Acknowledgements:
    - Kaptive team (Wyres, Holt et al.) for the K/O locus databases and CLI
    - Kaptive 3 subcommand interface: kaptive assembly <database> <genomes...>

Database keywords expected by Kaptive 3 (as per upstream docs):
    kpsc_k   Klebsiella pneumoniae species complex K locus
    kpsc_o   Klebsiella pneumoniae species complex O locus
"""

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import psutil


KLEBOSCOPE_BANNER = (
    "██╗  ██╗ ██╗      ███████╗ ██████╗  ██████╗ ███████╗  ██████╗ ██████╗  ██████╗ ███████╗\n"
    "██║ ██╔╝ ██║      ██╔════╝ ██╔══██╗██╔═══██╗██╔════╝ ██╔════╝██╔═══██╗ ██╔══██╗██╔════╝\n"
    "█████╔╝  ██║      █████╗   ██████╔╝██║   ██║███████╗ ██║     ██║   ██║ ██████╔╝█████╗\n"
    "██╔═██╗  ██║      ██╔══╝   ██╔══██╗██║   ██║╚════██║ ██║     ██║   ██║ ██╔═══╝ ██╔══╝\n"
    "██║  ██╗ ███████╗ ███████╗ ██████╔╝╚██████╔╝███████║ ╚██████╗╚██████╔╝ ██║     ███████╗\n"
    "╚═╝  ╚═╝ ╚══════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚══════╝  ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝"
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


class KleboKaptive:
    """Klebsiella K/O locus typing via Kaptive 3 subcommand CLI."""

    def __init__(
        self,
        k_db: str = "kpsc_k",
        o_db: str = "kpsc_o",
        cpus: Optional[int] = None,
    ):
        self.logger = self._setup_logging()
        self.k_db = k_db
        self.o_db = o_db
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        self.kaptive_available: Optional[bool] = None
        self.kaptive_version = self._get_kaptive_version()

        self.metadata = {
            "tool_name": "Kleboscope Kaptive",
            "version": "1.3.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "kaptive_version": self.kaptive_version,
            "k_database": self.k_db,
            "o_database": self.o_db,
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
        """Configure the module logger once."""
        import logging
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
        )
        return logging.getLogger(__name__)

    @staticmethod
    def _get_available_ram() -> float:
        """Return available RAM in gigabytes, defaulting to 8 if detection fails."""
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception:
            return 8.0

    @staticmethod
    def _calculate_optimal_cpus(user_cpus: Optional[int] = None) -> int:
        """Return the CPU budget with the same tiered scaling used by the MGE module."""
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

    # ------------------------------------------------------------------ Kaptive discovery
    def _get_kaptive_version(self) -> str:
        """Return the installed Kaptive version string."""
        try:
            result = subprocess.run(
                ["kaptive", "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
            text = result.stdout.strip() or result.stderr.strip()
            match = re.search(r"v?(\d+\.\d+(?:\.\d+)?)", text)
            return match.group(1) if match else text
        except Exception:
            return "Unknown"

    def check_kaptive_installed(self) -> bool:
        """Return True if the Kaptive CLI is available."""
        if self.kaptive_available is not None:
            return self.kaptive_available
        try:
            subprocess.run(
                ["kaptive", "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
            self.kaptive_available = True
            self.logger.info("Kaptive is available: %s", self.kaptive_version)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            self.logger.error("Kaptive check failed: %s", exc)
            self.kaptive_available = False
            return False

    # ------------------------------------------------------------------ Kaptive run + parse
    def _run_kaptive(
        self,
        database: str,
        genome_file: Path,
        out_tsv: Path,
    ) -> Tuple[bool, str]:
        """Run a single Kaptive 3 typing call and return (success, stderr).

        Both the genome path and the output path are resolved to absolute
        paths *before* the subprocess is launched. This prevents Kaptive
        from failing with 'File not found' when the subprocess working
        directory differs from the original launch directory.
        """
        genome_abs = genome_file.resolve()
        out_abs = out_tsv.resolve()
        out_abs.parent.mkdir(parents=True, exist_ok=True)

        if not genome_abs.exists():
            return False, f"Genome file does not exist: {genome_abs}"

        cmd = [
            "kaptive",
            "assembly",
            database,
            str(genome_abs),
            "--out",
            str(out_abs),
        ]
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                cwd=str(out_abs.parent),
            )
            if not out_abs.exists():
                # Defensive: Kaptive 3 normally honours --out exactly, but
                # occasionally writes a sibling .tsv with a different suffix.
                candidates = sorted(
                    out_abs.parent.glob("*.tsv"),
                    key=lambda p: p.stat().st_mtime,
                    reverse=True,
                )
                if candidates:
                    # Rewrite out_abs so the caller reads the file that was
                    # actually written.
                    candidates[0].rename(out_abs)
                else:
                    return False, (
                        f"Kaptive completed but no output file was created "
                        f"at {out_abs}"
                    )
            return True, result.stderr or ""
        except subprocess.CalledProcessError as exc:
            return False, (exc.stderr or str(exc)).strip()
        except FileNotFoundError as exc:
            return False, str(exc)

    def _parse_kaptive_tsv(self, tsv_path: Path, locus_type: str) -> List[Dict[str, Any]]:
        """Parse a Kaptive 3 TSV report into a list of row dictionaries."""
        rows: List[Dict[str, Any]] = []
        if not tsv_path.exists():
            return rows
        try:
            with open(tsv_path, "r", encoding="utf-8") as handle:
                lines = handle.readlines()
        except Exception as exc:
            self.logger.error("Could not read %s: %s", tsv_path, exc)
            return rows

        if len(lines) < 2:
            return rows

        headers = [h.strip() for h in lines[0].rstrip("\n").split("\t")]
        for raw in lines[1:]:
            raw = raw.rstrip("\n")
            if not raw.strip():
                continue
            parts = raw.split("\t")
            if len(parts) < len(headers):
                parts.extend([""] * (len(headers) - len(parts)))
            row = dict(zip(headers, parts))
            row["Locus Type"] = locus_type
            row["Expected genes in locus, details parsed"] = self._parse_gene_details(
                row.get("Expected genes in locus, details", "")
            )
            row["Other genes in locus, details parsed"] = self._parse_gene_details(
                row.get("Other genes in locus, details", "")
            )
            rows.append(row)
        return rows

    @staticmethod
    def _parse_gene_details(details: str) -> List[Dict[str, str]]:
        """Parse Kaptive's 'gene,identity,coverage' semicolon list."""
        genes: List[Dict[str, str]] = []
        if not details:
            return genes
        for chunk in details.split(";"):
            chunk = chunk.strip()
            if not chunk:
                continue
            parts = [p.strip() for p in chunk.split(",")]
            if len(parts) >= 3:
                genes.append({
                    "name": parts[0],
                    "identity": parts[1].replace("%", ""),
                    "coverage": parts[2].replace("%", ""),
                    "raw": chunk,
                })
        return genes

    # ------------------------------------------------------------------ per-sample pipeline
    def run_single_genome(self, genome_file: Path, sample_dir: Path) -> Dict[str, Any]:
        """Run K and O typing for one genome and emit one HTML report.

        Returns a result dict with an explicit `status` of:
            * "success" — both K and O completed without errors
            * "partial" — exactly one of K/O completed
            * "failed"  — neither completed
        A clean run that reports zero loci for either or both loci is still
        "success"; zero-locus genomes are legitimate biological outcomes.
        """
        genome_file = genome_file.resolve()
        sample_name = genome_file.stem
        sample_dir.mkdir(parents=True, exist_ok=True)

        hits: List[Dict[str, Any]] = []
        errors: List[str] = []
        k_ok = False
        o_ok = False

        with tempfile.TemporaryDirectory(prefix=f"kaptive_{sample_name}_") as temp_dir:
            temp_path = Path(temp_dir)

            for locus_type, database in (("K", self.k_db), ("O", self.o_db)):
                out_tsv = temp_path / f"{sample_name}_{locus_type}.tsv"
                ok, stderr = self._run_kaptive(database, genome_file, out_tsv)
                if ok:
                    if locus_type == "K":
                        k_ok = True
                    else:
                        o_ok = True
                    parsed = self._parse_kaptive_tsv(out_tsv, locus_type)
                    hits.extend(parsed)
                else:
                    msg = (
                        f"{locus_type} locus typing failed "
                        f"({database}): {stderr.strip()}"
                    )
                    self.logger.error(msg)
                    errors.append(msg)

        if k_ok and o_ok:
            status = "success"
        elif k_ok or o_ok:
            status = "partial"
        else:
            status = "failed"

        result: Dict[str, Any] = {
            "sample": sample_name,
            "input_file": str(genome_file),
            "hits": hits,
            "hit_count": len(hits),
            "errors": errors,
            "status": status,
        }

        self._create_sample_html(result, sample_dir)
        return result

    def _analyze_sample(self, hits: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Summarise K and O locus calls for a single sample."""
        analysis: Dict[str, Any] = {
            "k_locus": "None",
            "k_type": "None",
            "k_confidence": "",
            "k_identity": "",
            "k_coverage": "",
            "o_locus": "None",
            "o_type": "None",
            "o_confidence": "",
            "o_identity": "",
            "o_coverage": "",
            "total_loci": len(hits),
        }
        for hit in hits:
            locus_type = hit.get("Locus Type", "")
            prefix = "k" if locus_type == "K" else "o" if locus_type == "O" else None
            if prefix is None:
                continue
            analysis[f"{prefix}_locus"] = hit.get("Best match locus", "Unknown") or "Unknown"
            analysis[f"{prefix}_type"] = hit.get("Best match type", "Unknown") or "Unknown"
            analysis[f"{prefix}_confidence"] = hit.get("Match confidence", "")
            analysis[f"{prefix}_identity"] = hit.get("Identity", "")
            analysis[f"{prefix}_coverage"] = hit.get("Coverage", "")
        return analysis

    # ------------------------------------------------------------------ per-sample HTML
    def _create_sample_html(self, result: Dict[str, Any], output_dir: Path) -> Path:
        """Render the single per-sample Kaptive HTML report."""
        analysis = self._analyze_sample(result["hits"])
        random_quote = self._get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        if result["hits"]:
            all_keys = list(result["hits"][0].keys())
            table_keys = [
                k for k in all_keys
                if k not in (
                    "Expected genes in locus, details parsed",
                    "Other genes in locus, details parsed",
                )
            ]
        else:
            table_keys = []

        table_rows = ""
        for hit in result["hits"]:
            row_class = "k-locus-row" if hit.get("Locus Type") == "K" else "o-locus-row"
            cells = []
            for key in table_keys:
                value = hit.get(key, "")
                if isinstance(value, str) and len(value) > 200:
                    value = value[:200] + "..."
                cells.append(f"<td>{value}</td>")
            table_rows += f'<tr class="{row_class}">{"".join(cells)}</tr>\n'

        gene_blocks = ""
        for hit in result["hits"]:
            locus_type = hit.get("Locus Type", "")
            color = "#28a745" if locus_type == "K" else "#17a2b8"
            locus = hit.get("Best match locus", "Unknown")
            expected = hit.get("Expected genes in locus, details parsed", [])
            other = hit.get("Other genes in locus, details parsed", [])
            expected_html = "".join(
                f'<span class="gene-item">{g["name"]}: {g["identity"]}% identity, '
                f'{g["coverage"]}% coverage</span><br>'
                for g in expected
            ) or "<span>None</span>"
            other_html = "".join(
                f'<span class="gene-item">{g["name"]}: {g["identity"]}% identity, '
                f'{g["coverage"]}% coverage</span><br>'
                for g in other
            ) or "<span>None</span>"
            gene_blocks += f"""
            <div class="gene-details-container" style="margin-top:20px;padding:15px;background:#f8f9fa;border-radius:8px;border-left:4px solid {color};">
                <h4 style="margin-top:0;color:{color};">{locus} ({locus_type} locus) - Gene Details</h4>
                <div style="margin-top:10px;">
                    <strong>Expected genes in locus:</strong>
                    <div style="margin:5px 0 10px 20px;">{expected_html}</div>
                    <strong>Other genes in locus:</strong>
                    <div style="margin:5px 0 0 20px;">{other_html}</div>
                </div>
            </div>
            """

        error_html = ""
        if result["errors"]:
            items = "".join(f"<li>{e}</li>" for e in result["errors"])
            error_html = f'<div class="error-box"><strong>Kaptive errors:</strong><ul>{items}</ul></div>'

        footer = _render_footer("Kaptive K/O Analysis Module", current_time)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KLEBOSCOPE - {result['sample']} Kaptive Report</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color:#ffffff; padding:20px; min-height:100vh;
}}
.container {{ max-width:1600px; margin:0 auto; }}
.ascii-container {{
    background: rgba(0,0,0,0.7); padding:20px; border-radius:15px;
    margin-bottom:20px; box-shadow:0 8px 32px rgba(0,0,0,0.4);
    border:2px solid rgba(0,119,190,0.3); overflow-x:auto; text-align:center;
}}
.ascii-art {{
    font-family: 'Cascadia Mono', 'Consolas', 'Monaco', 'Menlo',
                 'DejaVu Sans Mono', 'Liberation Mono', 'Courier New', monospace;
    font-size:12px; line-height:1.05; white-space:pre;
    color:#00d4ff; text-shadow:0 0 8px rgba(0,212,255,0.4);
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
.summary-stats {{
    display:flex; justify-content:space-around; flex-wrap:wrap; margin:20px 0;
}}
.stat-card {{
    background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
    color:white; padding:20px; border-radius:12px; text-align:center;
    margin:10px; flex:1; min-width:180px;
    box-shadow:0 4px 15px rgba(0,0,0,0.2);
}}
.k-stat-card {{ background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%); }}
.o-stat-card {{ background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%); }}
.stat-label {{ font-size:13px; opacity:0.9; text-transform:uppercase; letter-spacing:0.5px; }}
.stat-value {{ font-size:24px; font-weight:bold; }}
.table-responsive {{ width:100%; overflow-x:auto; margin:20px 0; }}
.gene-table {{
    width:100%; border-collapse:collapse; background:white; font-size:13px;
}}
.gene-table th, .gene-table td {{
    padding:10px; text-align:left; border-bottom:1px solid #e5e7eb;
    white-space:nowrap;
}}
.gene-table th {{
    background:#005a8c; color:white; font-weight:600; position:sticky; top:0;
}}
.gene-table tr:hover {{ background:#f3f4f6; }}
.k-locus-row {{ background-color:#d4edda; }}
.o-locus-row {{ background-color:#d1ecf1; }}
.gene-item {{
    display:inline-block; background:#e9ecef; padding:3px 8px;
    margin:2px; border-radius:12px; font-size:0.85em;
}}
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
    <h2>📊 Sample Information</h2>
    <div class="summary-stats">
        <div class="stat-card"><div class="stat-label">Sample</div><div class="stat-value" style="font-size:16px;">{result['sample']}</div></div>
        <div class="stat-card"><div class="stat-label">Total loci</div><div class="stat-value">{analysis['total_loci']}</div></div>
        <div class="stat-card k-stat-card"><div class="stat-label">K locus</div><div class="stat-value" style="font-size:18px;">{analysis['k_locus']}</div></div>
        <div class="stat-card o-stat-card"><div class="stat-label">O locus</div><div class="stat-value" style="font-size:18px;">{analysis['o_locus']}</div></div>
    </div>
    <p><strong>K type:</strong> {analysis['k_type']} &nbsp;|&nbsp; <strong>K confidence:</strong> {analysis['k_confidence'] or '—'}</p>
    <p><strong>O type:</strong> {analysis['o_type']} &nbsp;|&nbsp; <strong>O confidence:</strong> {analysis['o_confidence'] or '—'}</p>
    <p><strong>Kaptive version:</strong> {self.metadata['kaptive_version']} &nbsp;|&nbsp; <strong>K database:</strong> {self.metadata['k_database']} &nbsp;|&nbsp; <strong>O database:</strong> {self.metadata['o_database']}</p>
    {error_html}
</div>

<div class="card">
    <h2>🔬 Kaptive Typing Results</h2>
    <div class="table-responsive">
        <table class="gene-table">
            <thead>
                <tr>{"".join(f"<th>{key}</th>" for key in table_keys)}</tr>
            </thead>
            <tbody>
                {table_rows or '<tr><td colspan="1" style="text-align:center;color:#9ca3af;">No Kaptive rows returned</td></tr>'}
            </tbody>
        </table>
    </div>
</div>

<div class="card">
    <h2>🧬 Detailed Gene Information</h2>
    {gene_blocks or '<p style="color:#6b7280;">No gene-level details returned.</p>'}
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

        out_file = output_dir / f"{result['sample']}_kaptive.html"
        with open(out_file, "w", encoding="utf-8") as handle:
            handle.write(html)
        self.logger.info("Per-sample Kaptive HTML: %s", out_file)
        return out_file

    def _get_random_quote(self) -> str:
        """Return a random science quote for initial render."""
        import random
        return random.choice(self.science_quotes)

    # ------------------------------------------------------------------ batch orchestration
    def process_multiple_genomes(
        self,
        genome_pattern: str,
        output_base: str = "kaptive_results",
        workers: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Run Kaptive on every matching assembly and write batch summaries."""
        if not self.check_kaptive_installed():
            raise RuntimeError("Kaptive is not available on PATH")

        fasta_paths = self._resolve_inputs(genome_pattern)
        if not fasta_paths:
            raise FileNotFoundError(f"No FASTA files found matching: {genome_pattern}")

        output_dir = Path(output_base)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("Found %d genome(s)", len(fasta_paths))
        self.logger.info(
            "Processing sequentially, %d CPU thread(s) available to Kaptive",
            self.cpus,
        )

        results: Dict[str, Dict[str, Any]] = {}
        for fasta in fasta_paths:
            sample_dir = output_dir / fasta.stem
            self.logger.info("Processing %s ...", fasta.name)
            result = self.run_single_genome(fasta, sample_dir)
            results[result["sample"]] = result
            self.logger.info(
                "  %s: %d loci | status=%s",
                result["sample"], result["hit_count"], result["status"],
            )

        self._write_summary_tsv(results, output_dir / "kleboscope_kaptive_summary.tsv")
        self._write_summary_json(results, output_dir / "kleboscope_kaptive_summary.json")
        self._write_summary_html(
            results,
            output_dir / "kleboscope_kaptive_summary.html",
        )

        return results

    @staticmethod
    def _resolve_inputs(pattern: str) -> List[Path]:
        """Resolve a glob (or explicit file) into a sorted list of FASTA paths."""
        matches = glob.glob(pattern)
        if not matches:
            for ext in (".fna", ".fasta", ".fa", ".faa"):
                matches.extend(glob.glob(f"{pattern}{ext}"))
        return sorted({Path(m).resolve() for m in matches if Path(m).is_file()})

    # ------------------------------------------------------------------ batch TSV
    def _write_summary_tsv(
        self,
        results: Dict[str, Dict[str, Any]],
        output_path: Path,
    ) -> None:
        """Write the batch-level Kaptive TSV, one row per genome."""
        columns = [
            "sample", "status", "k_locus", "k_type", "k_confidence",
            "k_identity", "k_coverage", "o_locus", "o_type", "o_confidence",
            "o_identity", "o_coverage",
        ]
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write("\t".join(columns) + "\n")
            for sample, result in sorted(results.items()):
                analysis = self._analyze_sample(result["hits"])
                row = [sample, result["status"]] + [
                    analysis[key] for key in columns[2:]
                ]
                handle.write("\t".join(str(v) for v in row) + "\n")
        self.logger.info("Wrote batch TSV: %s", output_path)

    # ------------------------------------------------------------------ batch JSON
    def _write_summary_json(
        self,
        results: Dict[str, Dict[str, Any]],
        output_path: Path,
    ) -> None:
        """Write the single consolidated batch JSON."""
        payload: Dict[str, Any] = {
            "metadata": {
                "tool": self.metadata["tool_name"],
                "version": self.metadata["version"],
                "kaptive_version": self.metadata["kaptive_version"],
                "k_database": self.metadata["k_database"],
                "o_database": self.metadata["o_database"],
                "analysis_date": self.metadata["analysis_date"],
                "total_genomes": len(results),
            },
            "statistics": self._calculate_statistics(results),
            "samples": {},
        }
        for sample, result in sorted(results.items()):
            analysis = self._analyze_sample(result["hits"])
            payload["samples"][sample] = {
                "status": result["status"],
                "total_loci": result["hit_count"],
                "errors": result["errors"],
                "summary": analysis,
                "hits": result["hits"],
            }
        with open(output_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
        self.logger.info("Wrote batch JSON: %s", output_path)

    def _calculate_statistics(self, results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Compute cross-genome K/O locus statistics for the batch JSON."""
        stats: Dict[str, Any] = {
            "total_genomes": len(results),
            "successful": sum(1 for r in results.values() if r["status"] == "success"),
            "partial": sum(1 for r in results.values() if r["status"] == "partial"),
            "failed": sum(1 for r in results.values() if r["status"] == "failed"),
            "genomes_with_k": 0,
            "genomes_with_o": 0,
            "k_type_counts": {},
            "o_type_counts": {},
            "k_locus_counts": {},
            "o_locus_counts": {},
        }
        for result in results.values():
            analysis = self._analyze_sample(result["hits"])
            if analysis["k_locus"] != "None":
                stats["genomes_with_k"] += 1
                stats["k_locus_counts"][analysis["k_locus"]] = (
                    stats["k_locus_counts"].get(analysis["k_locus"], 0) + 1
                )
                stats["k_type_counts"][analysis["k_type"]] = (
                    stats["k_type_counts"].get(analysis["k_type"], 0) + 1
                )
            if analysis["o_locus"] != "None":
                stats["genomes_with_o"] += 1
                stats["o_locus_counts"][analysis["o_locus"]] = (
                    stats["o_locus_counts"].get(analysis["o_locus"], 0) + 1
                )
                stats["o_type_counts"][analysis["o_type"]] = (
                    stats["o_type_counts"].get(analysis["o_type"], 0) + 1
                )
        return stats

    # ------------------------------------------------------------------ batch HTML
    def _write_summary_html(
        self,
        results: Dict[str, Dict[str, Any]],
        output_path: Path,
    ) -> None:
        """Write the batch-level Kaptive HTML summary."""
        stats = self._calculate_statistics(results)
        random_quote = self._get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        rows = []
        for sample, result in sorted(results.items()):
            analysis = self._analyze_sample(result["hits"])
            if result["status"] == "success":
                row_class = "success-row"
            elif result["status"] == "partial":
                row_class = "partial-row"
            else:
                row_class = "failed-row"
            rows.append(
                f"<tr class='{row_class}'>"
                f"<td><strong>{sample}</strong></td>"
                f"<td>{result['status']}</td>"
                f"<td>{analysis['k_locus']}</td>"
                f"<td>{analysis['k_type']}</td>"
                f"<td>{analysis['k_confidence'] or '—'}</td>"
                f"<td>{analysis['o_locus']}</td>"
                f"<td>{analysis['o_type']}</td>"
                f"<td>{analysis['o_confidence'] or '—'}</td>"
                f"</tr>"
            )
        table_rows = "\n".join(rows)

        footer = _render_footer("Kaptive K/O Batch Summary", current_time)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KLEBOSCOPE - Kaptive Batch Summary</title>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color:#ffffff; padding:20px; min-height:100vh;
}}
.container {{ max-width:1600px; margin:0 auto; }}
.ascii-container {{
    background: rgba(0,0,0,0.7); padding:20px; border-radius:15px;
    margin-bottom:20px; box-shadow:0 8px 32px rgba(0,0,0,0.4);
    border:2px solid rgba(0,119,190,0.3); overflow-x:auto; text-align:center;
}}
.ascii-art {{
    font-family: 'Cascadia Mono', 'Consolas', 'Monaco', 'Menlo',
                 'DejaVu Sans Mono', 'Liberation Mono', 'Courier New', monospace;
    font-size:12px; line-height:1.05; white-space:pre;
    color:#00d4ff; text-shadow:0 0 8px rgba(0,212,255,0.4);
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
.k-stat-card {{ background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%); }}
.o-stat-card {{ background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%); }}
.table-responsive {{ width:100%; overflow-x:auto; margin:20px 0; }}
.gene-table {{
    width:100%; border-collapse:collapse; background:white; font-size:13px;
}}
.gene-table th, .gene-table td {{
    padding:10px; text-align:left; border-bottom:1px solid #e5e7eb;
}}
.gene-table th {{
    background:#005a8c; color:white; font-weight:600; position:sticky; top:0;
}}
.gene-table tr:hover {{ background:#f3f4f6; }}
.success-row {{ border-left:4px solid #28a745; }}
.partial-row {{ border-left:4px solid #f59e0b; }}
.failed-row {{ border-left:4px solid #dc3545; background-color:#f8d7da; }}
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
    <h2>📊 Batch Kaptive Summary</h2>
    <div class="summary-stats">
        <div class="stat-card"><div class="stat-label">Total genomes</div><div class="stat-value" style="font-size:2em;">{stats['total_genomes']}</div></div>
        <div class="stat-card"><div class="stat-label">Successful</div><div class="stat-value" style="font-size:2em;">{stats['successful']}</div></div>
        <div class="stat-card k-stat-card"><div class="stat-label">Genomes with K locus</div><div class="stat-value" style="font-size:2em;">{stats['genomes_with_k']}</div></div>
        <div class="stat-card o-stat-card"><div class="stat-label">Genomes with O locus</div><div class="stat-value" style="font-size:2em;">{stats['genomes_with_o']}</div></div>
    </div>
    <p><strong>Kaptive version:</strong> {self.metadata['kaptive_version']} &nbsp;|&nbsp; <strong>K database:</strong> {self.metadata['k_database']} &nbsp;|&nbsp; <strong>O database:</strong> {self.metadata['o_database']}</p>
</div>

<div class="card">
    <h2>🔍 Detailed Results by Genome</h2>
    <div class="table-responsive">
        <table class="gene-table">
            <thead>
                <tr>
                    <th>Genome</th><th>Status</th>
                    <th>K locus</th><th>K type</th><th>K confidence</th>
                    <th>O locus</th><th>O type</th><th>O confidence</th>
                </tr>
            </thead>
            <tbody>
                {table_rows or '<tr><td colspan="8" style="text-align:center;color:#9ca3af;">No samples processed</td></tr>'}
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
        description="Kleboscope Kaptive - K. pneumoniae K/O locus typing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FASTA file, glob, or prefix without extension",
    )
    parser.add_argument(
        "-o", "--output", default="kaptive_results",
        help="Output directory",
    )
    parser.add_argument(
        "--k-db", default="kpsc_k",
        help="Kaptive 3 K-locus database keyword or path",
    )
    parser.add_argument(
        "--o-db", default="kpsc_o",
        help="Kaptive 3 O-locus database keyword or path",
    )
    parser.add_argument(
        "-c", "--cpus", type=int, default=None,
        help="CPU core budget (auto-detected if omitted)",
    )
    args = parser.parse_args()

    started = time.time()
    executor = KleboKaptive(k_db=args.k_db, o_db=args.o_db, cpus=args.cpus)

    print("=" * 70)
    print("🧬 Kleboscope Kaptive K/O Locus Analysis")
    print("=" * 70)
    print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print(f"Affiliation: University of Ghana Medical School")
    print(f"Kaptive: {executor.kaptive_version} | K db: {args.k_db} | O db: {args.o_db}")
    print("=" * 70)

    try:
        results = executor.process_multiple_genomes(args.input, args.output)
        elapsed = time.time() - started
        total_loci = sum(r["hit_count"] for r in results.values())
        successful = sum(1 for r in results.values() if r["status"] == "success")
        partial = sum(1 for r in results.values() if r["status"] == "partial")

        print()
        print("=" * 70)
        print("📊 Kaptive Analysis Final Summary")
        print("=" * 70)
        print(f"Genomes processed    : {len(results)}")
        print(f"Successful analyses  : {successful}")
        print(f"Partial analyses     : {partial}")
        print(f"Total loci detected  : {total_loci}")
        print(f"Elapsed time         : {elapsed:.1f}s")
        print(f"Output directory     : {args.output}")
        print()
        print("Per-sample (1 file each):")
        print("  • <sample>/<sample>_kaptive.html")
        print("Batch-level (3 files):")
        print("  • kleboscope_kaptive_summary.tsv")
        print("  • kleboscope_kaptive_summary.html")
        print("  • kleboscope_kaptive_summary.json")
    except Exception as exc:
        executor.logger.error("Analysis failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()