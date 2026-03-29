#!/usr/bin/env python3
"""
Kleboscope Kaptive K/O Locus Analysis - K. pneumoniae Capsule and Lipopolysaccharide Typing
Comprehensive Kaptive analysis for K. pneumoniae with beautiful HTML reporting
Author: Brown Beckley
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-03-20
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from typing import List, Dict, Any
import argparse
import re
from datetime import datetime
import json
import random
from collections import defaultdict

class KleboscopeKaptive:
    def __init__(self, k_db="kp_k", o_db="kp_o"):
        self.logger = self._setup_logging()
        self.k_db = k_db
        self.o_db = o_db
        self.kaptive_available = None  # will be set after check

        # ASCII Art for Kleboscope (blue theme)
        self.ascii_art = r"""
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
"""

        self.metadata = {
            "tool_name": "Kleboscope Kaptive K/O Analysis",
            "version": "1.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "database_version": "2026-01-21 (kp_k/kp_o)"
        }

        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "Kleboscope turns genomic complexity into actionable insights for AMR surveillance.", "author": "Brown Beckley"}
        ]

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)

    def _get_kaptive_version(self) -> str:
        try:
            result = subprocess.run(["kaptive", "--version"], capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            match = re.search(r'v?(\d+\.\d+\.\d+)', version_line)
            return match.group(1) if match else version_line
        except Exception as e:
            self.logger.warning(f"Could not determine Kaptive version: {e}")
            return "Unknown"

    def check_kaptive_installed(self) -> bool:
        """Check Kaptive once and store result."""
        if self.kaptive_available is not None:
            return self.kaptive_available
        try:
            subprocess.run(["kaptive", "--version"], capture_output=True, text=True, check=True)
            self.logger.info("Kaptive is available")
            self.kaptive_available = True
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Kaptive check failed: {e}")
            self.logger.info("Please install Kaptive: conda install -c bioconda kaptive")
            self.kaptive_available = False
            return False

    def run_kaptive_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        """Run Kaptive for a single genome. output_dir should already be the genome‑specific directory."""
        genome_name = Path(genome_file).stem
        # output_dir is the genome‑specific directory (e.g., kaptive_results/GCF_xxx)
        os.makedirs(output_dir, exist_ok=True)

        k_tsv = os.path.join(output_dir, f"{genome_name}_K.tsv")
        o_tsv = os.path.join(output_dir, f"{genome_name}_O.tsv")
        combined_tsv = os.path.join(output_dir, f"{genome_name}_combined.tsv")

        self.logger.info(f"Running Kaptive K locus analysis on {genome_name}")
        try:
            cmd_k = ["kaptive", "assembly", self.k_db, genome_file, "-o", k_tsv, "--verbose"]
            subprocess.run(cmd_k, capture_output=True, text=True, check=True)

            self.logger.info(f"Running Kaptive O locus analysis on {genome_name}")
            cmd_o = ["kaptive", "assembly", self.o_db, genome_file, "-o", o_tsv, "--verbose"]
            subprocess.run(cmd_o, capture_output=True, text=True, check=True)

            # Combine results
            self._combine_k_o_results(k_tsv, o_tsv, combined_tsv)

            # Parse combined TSV fully
            hits = self._parse_full_kaptive_output(combined_tsv)

            # Create reports
            self._create_kaptive_html_report(genome_name, hits, output_dir)
            self._create_kaptive_json_report(genome_name, hits, output_dir)

            return {
                'genome': genome_name,
                'output_file': combined_tsv,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Kaptive failed for {genome_name}: {e}")
            return {
                'genome': genome_name,
                'output_file': combined_tsv,
                'hits': [],
                'hit_count': 0,
                'status': 'failed',
                'error': str(e)
            }

    def _combine_k_o_results(self, k_tsv: str, o_tsv: str, combined_tsv: str):
        try:
            with open(combined_tsv, 'w') as out:
                # Write K results if exists
                if os.path.exists(k_tsv):
                    with open(k_tsv, 'r') as f:
                        lines = f.readlines()
                        if lines:
                            # Header
                            out.write(lines[0])
                            for line in lines[1:]:
                                out.write(line)
                # Write O results if exists
                if os.path.exists(o_tsv):
                    with open(o_tsv, 'r') as f:
                        lines = f.readlines()
                        if lines:
                            # If we already wrote K header, skip O header
                            if not os.path.exists(k_tsv):
                                out.write(lines[0])
                            for line in lines[1:]:
                                out.write(line)
            self.logger.info(f"Combined K and O results to {combined_tsv}")
        except Exception as e:
            self.logger.error(f"Error combining K/O results: {e}")

    def _parse_full_kaptive_output(self, combined_file: str) -> List[Dict]:
        """Parse the combined TSV into a list of dicts, preserving all fields."""
        hits = []
        try:
            with open(combined_file, 'r') as f:
                lines = f.readlines()
            if not lines or len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            for line in lines[1:]:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < len(headers):
                    parts.extend([''] * (len(headers) - len(parts)))
                hit = dict(zip(headers, parts))
                # Determine locus type from Best match locus (KL vs OL)
                locus = hit.get('Best match locus', '')
                if locus.startswith('KL'):
                    hit['Locus Type'] = 'K'
                elif locus.startswith('OL'):
                    hit['Locus Type'] = 'O'
                else:
                    hit['Locus Type'] = 'Unknown'
                # Parse detailed gene list from Expected genes in locus, details
                gene_details = hit.get('Expected genes in locus, details', '')
                hit['Gene Details Parsed'] = self._parse_gene_details(gene_details)
                # Also parse other genes in locus if needed
                other_genes = hit.get('Other genes in locus, details', '')
                hit['Other Genes Parsed'] = self._parse_gene_details(other_genes)
                hits.append(hit)
        except Exception as e:
            self.logger.error(f"Error parsing {combined_file}: {e}")
        self.logger.info(f"Parsed {len(hits)} Kaptive hits from {combined_file}")
        return hits

    def _parse_gene_details(self, gene_details_str: str) -> List[Dict]:
        genes = []
        if not gene_details_str:
            return genes
        try:
            for part in gene_details_str.split(';'):
                if part:
                    sub = part.split(',')
                    if len(sub) >= 3:
                        genes.append({
                            'name': sub[0],
                            'identity': sub[1].replace('%', ''),
                            'coverage': sub[2].replace('%', ''),
                            'raw': part
                        })
        except Exception as e:
            self.logger.warning(f"Could not parse gene details: {e}")
        return genes

    def _create_kaptive_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """AcinetoScope-style per-genome HTML report: full table with sorting, plus gene details below."""
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Prepare all columns from the first hit (excluding parsed lists)
        if hits:
            all_keys = list(hits[0].keys())
            # Remove parsed list fields from table display
            table_keys = [k for k in all_keys if k not in ['Gene Details Parsed', 'Other Genes Parsed']]
        else:
            table_keys = []

        # Build table rows
        table_rows = ""
        for hit in hits:
            row_class = "k-locus-row" if hit.get('Locus Type') == 'K' else "o-locus-row" if hit.get('Locus Type') == 'O' else ""
            table_rows += f'<tr class="{row_class}">'
            for key in table_keys:
                val = hit.get(key, '')
                # Truncate very long values for readability
                if isinstance(val, str) and len(val) > 200:
                    val = val[:200] + '...'
                table_rows += f"<td>{val}</td>"
            table_rows += "</tr>\n"

        # Build gene details sections for each hit
        gene_details_html = ""
        for idx, hit in enumerate(hits):
            locus = hit.get('Best match locus', 'Unknown')
            locus_type = hit.get('Locus Type', '')
            color = "#28a745" if locus_type == 'K' else "#17a2b8" if locus_type == 'O' else "#6c757d"
            expected_genes = hit.get('Gene Details Parsed', [])
            other_genes = hit.get('Other Genes Parsed', [])

            gene_details_html += f"""
            <div class="gene-details-container" style="margin-top: 20px; padding: 15px; background: #f8f9fa; border-radius: 8px; border-left: 4px solid {color};">
                <h4 style="margin-top: 0; color: {color};">{locus} ({locus_type} Locus) - Gene Details:</h4>
                <div style="margin-top: 10px;">
                    <strong>Expected genes in locus:</strong>
                    <div style="margin: 5px 0 10px 20px;">
            """
            if expected_genes:
                for gene in expected_genes:
                    gene_details_html += f'<span class="gene-item">{gene["name"]}: {gene["identity"]}% identity, {gene["coverage"]}% coverage</span><br>'
            else:
                gene_details_html += "<span>None</span>"
            gene_details_html += """
                    </div>
                    <strong>Other genes in locus:</strong>
                    <div style="margin: 5px 0 0 20px;">
            """
            if other_genes:
                for gene in other_genes:
                    gene_details_html += f'<span class="gene-item">{gene["name"]}: {gene["identity"]}% identity, {gene["coverage"]}% coverage</span><br>'
            else:
                gene_details_html += "<span>None</span>"
            gene_details_html += """
                    </div>
                </div>
            </div>
            """

        # HTML content as a single f-string
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KLEBOSCOPE - Kaptive K/O Analysis Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1800px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid rgba(0,119,190,0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.2;
            white-space: pre;
            color: #0077be;
            text-shadow: 0 0 10px rgba(0,119,190,0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255,255,255,0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .card {{
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.2);
        }}
        .card h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .k-stat-card {{ background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%); }}
        .o-stat-card {{ background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%); }}
        .controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        .search-box input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #3b82f6;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
        }}
        .export-buttons button.print {{ background: #10b981; }}
        .gene-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .gene-table th, .gene-table td {{
            padding: 12px 8px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
            white-space: nowrap;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: relative;
        }}
        .gene-table th:hover {{
            background: linear-gradient(135deg, #2563eb 0%, #1e3a8a 100%);
        }}
        .gene-table th.sortable::after {{
            content: "↕";
            position: absolute;
            right: 8px;
            opacity: 0.6;
        }}
        .gene-table tr:hover {{ background-color: #f8f9fa; }}
        .k-locus-row {{ background-color: #d4edda; }}
        .o-locus-row {{ background-color: #d1ecf1; }}
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        .gene-item {{
            display: inline-block;
            background: #e9ecef;
            padding: 3px 8px;
            margin: 2px;
            border-radius: 12px;
            font-size: 0.85em;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.3);
            border-radius: 10px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255,255,255,0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 1200px) {{
            .gene-table {{ font-size: 12px; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>

        <div class="card">
            <h2>📊 K. pneumoniae K/O Locus Summary</h2>
            <div class="summary-stats">
                <div class="stat-card"><div class="metric-label">Total Loci</div><div class="metric-value">{len(hits)}</div></div>
                <div class="stat-card k-stat-card"><div class="metric-label">K Locus (Capsule)</div><div class="metric-value">{sum(1 for h in hits if h.get('Locus Type') == 'K')}</div></div>
                <div class="stat-card o-stat-card"><div class="metric-label">O Locus (LPS)</div><div class="metric-value">{sum(1 for h in hits if h.get('Locus Type') == 'O')}</div></div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Kaptive Version:</strong> {self.metadata.get('kaptive_version', 'Unknown')}</p>
            <p><strong>Database:</strong> {self.metadata['database_version']}</p>
        </div>

        <div class="card">
            <h2>🔬 Complete Kaptive Analysis Results</h2>
            <div class="controls">
                <div class="search-box"><input type="text" id="search" placeholder="Search..." onkeyup="searchTable()"></div>
                <div class="export-buttons">
                    <button onclick="exportToCSV()">📥 Export CSV</button>
                    <button onclick="printReport()" class="print">🖨️ Print</button>
                </div>
            </div>
            <div class="table-container">
                <table class="gene-table" id="kaptive-table">
                    <thead>
                        <tr>
"""
        # Generate column headers with sortable class
        for i, key in enumerate(table_keys):
            html += f'<th class="sortable" onclick="sortTable({i})">{key}</th>\n'
        html += """
                        </thead>
                    <tbody>
                        {table_rows}
                    </tbody>
                 </table>
            </div>

            <h3 style="margin-top: 30px;">🧬 Detailed Gene Information</h3>
            {gene_details_html}
        </div>

        <div class="card">
            <h2>📁 Generated Files</h2>
            <ul>
                <li><strong>Individual TSV files</strong>: <code>{genome_name}_K.tsv</code>, <code>{genome_name}_O.tsv</code>, <code>{genome_name}_combined.tsv</code></li>
                <li><strong>HTML report</strong>: <code>{genome_name}_kaptive_report.html</code></li>
                <li><strong>JSON report</strong>: <code>{genome_name}_kaptive_report.json</code></li>
            </ul>
        </div>

        <div class="footer">
            <p><strong>KLEBOSCOPE</strong> - Kaptive K/O Analysis Module</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>

    <script>
        function sortTable(n) {{
            const table = document.getElementById("kaptive-table");
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));
            rows.sort((a, b) => {{
                let x = a.children[n].textContent.trim();
                let y = b.children[n].textContent.trim();
                // Try numeric conversion
                let xn = parseFloat(x);
                let yn = parseFloat(y);
                if (!isNaN(xn) && !isNaN(yn)) {{
                    return xn - yn;
                }}
                return x.localeCompare(y);
            }});
            rows.forEach(row => tbody.appendChild(row));
        }}
        function searchTable() {{
            const input = document.getElementById("search");
            const filter = input.value.toLowerCase();
            const table = document.getElementById("kaptive-table");
            const rows = table.getElementsByTagName("tr");
            let visible = 0;
            for (let i = 1; i < rows.length; i++) {{
                const row = rows[i];
                const text = row.textContent.toLowerCase();
                if (text.includes(filter)) {{
                    row.style.display = "";
                    visible++;
                }} else {{
                    row.style.display = "none";
                }}
            }}
        }}
        function exportToCSV() {{
            const table = document.getElementById("kaptive-table");
            const rows = table.getElementsByTagName("tr");
            let csv = [];
            for (let i = 0; i < rows.length; i++) {{
                const row = rows[i];
                if (row.style.display !== "none") {{
                    const cells = row.getElementsByTagName("th").length ? row.getElementsByTagName("th") : row.getElementsByTagName("td");
                    const rowData = [];
                    for (let c = 0; c < cells.length; c++) {{
                        rowData.push(cells[c].textContent.trim().replace(/,/g, ';'));
                    }}
                    csv.push(rowData.join(','));
                }}
            }}
            const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'kaptive_report.csv';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }}
        function printReport() {{
            window.print();
        }}
        const quotes = {json.dumps(self.science_quotes)};
        let currentQuote = 0;
        function rotateQuote() {{
            document.getElementById('quoteText').innerHTML = '"' + quotes[currentQuote].text + '"';
            document.getElementById('quoteAuthor').innerHTML = '— ' + quotes[currentQuote].author;
            currentQuote = (currentQuote + 1) % quotes.length;
        }}
        setInterval(rotateQuote, 10000);
    </script>
</body>
</html>
"""
        # Replace the two large placeholders (the table_rows and gene_details_html) that are inside the f-string
        # but because they contain nested f‑string braces, we need to insert them after the main f-string.
        html = html.replace("{table_rows}", table_rows)
        html = html.replace("{gene_details_html}", gene_details_html)

        out_file = os.path.join(output_dir, f"{genome_name}_kaptive_report.html")
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write(html)
        self.logger.info(f"K. pneumoniae Kaptive HTML report generated: {out_file}")

    def _create_kaptive_json_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        json_data = {
            'metadata': {
                'genome': genome_name,
                'analysis_date': self.metadata['analysis_date'],
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'kaptive_version': self.metadata.get('kaptive_version', 'Unknown'),
                'database_version': self.metadata['database_version']
            },
            'summary': {
                'total_loci': len(hits),
                'k_loci': sum(1 for h in hits if h.get('Locus Type') == 'K'),
                'o_loci': sum(1 for h in hits if h.get('Locus Type') == 'O'),
                'k_types': list(set(h.get('Best match type', '') for h in hits if h.get('Locus Type') == 'K')),
                'o_types': list(set(h.get('Best match type', '') for h in hits if h.get('Locus Type') == 'O'))
            },
            'hits': hits
        }
        out_file = os.path.join(output_dir, f"{genome_name}_kaptive_report.json")
        with open(out_file, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2)
        self.logger.info(f"K. pneumoniae Kaptive JSON report generated: {out_file}")

    def get_random_quote(self):
        return random.choice(self.science_quotes)

    def process_single_genome(self, genome_file: str, output_base: str = "kaptive_results") -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        # Create a genome-specific directory directly under output_base
        genome_dir = os.path.join(output_base, genome_name)
        os.makedirs(genome_dir, exist_ok=True)
        self.logger.info(f"=== PROCESSING GENOME: {genome_name} ===")
        if not self.check_kaptive_installed():
            self.logger.error("Kaptive not available!")
            return {'genome': genome_name, 'hits': [], 'hit_count': 0, 'status': 'failed', 'error': 'Kaptive not available'}
        # Pass genome_dir (the sample directory) to run_kaptive_single_genome
        result = self.run_kaptive_single_genome(genome_file, genome_dir)
        status = "✓" if result['status'] == 'success' else "✗"
        self.logger.info(f"{status} {genome_name}: {result['hit_count']} K/O loci")
        return result

    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "kaptive_results") -> Dict[str, Any]:
        fasta_exts = ['.fna', '.fasta', '.fa', '.faa']
        genome_files = []
        genome_files.extend(glob.glob(genome_pattern))
        if not genome_files:
            for ext in fasta_exts:
                genome_files.extend(glob.glob(f"{genome_pattern}{ext}"))
        genome_files = list(set([f for f in genome_files if os.path.exists(f)]))
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        self.logger.info(f"Found {len(genome_files)} K. pneumoniae genomes")
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        for genome_file in genome_files:
            try:
                result = self.process_single_genome(genome_file, output_base)
                all_results[result['genome']] = result
            except Exception as e:
                self.logger.error(f"Failed: {genome_file} - {e}")
                all_results[Path(genome_file).stem] = {'genome': Path(genome_file).stem, 'hits': [], 'hit_count': 0, 'status': 'failed', 'error': str(e)}
        self.create_kaptive_summary(all_results, output_base)
        self.logger.info("=== K. PNEUMONIAE K/O ANALYSIS COMPLETE ===")
        return all_results

    

    def create_kaptive_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating K. pneumoniae Kaptive summary files...")
        # TSV summary with full details
        summary_file = os.path.join(output_base, "klebo_kaptive_summary.tsv")
        all_hits = []
        for gn, res in all_results.items():
            for hit in res.get('hits', []):
                hit_copy = hit.copy()
                hit_copy['Genome'] = gn
                all_hits.append(hit_copy)
        if all_hits:
            # Collect all headers from the first hit
            headers = list(all_hits[0].keys())
            # Remove parsed list fields
            output_headers = [h for h in headers if h not in ['Gene Details Parsed', 'Other Genes Parsed']]
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write('\t'.join(output_headers) + '\n')
                for hit in all_hits:
                    row = []
                    for h in output_headers:
                        val = hit.get(h, '')
                        if isinstance(val, list):
                            # For gene lists, join into a string
                            val = ';'.join([g.get('raw', '') for g in val]) if val else ''
                        row.append(str(val).replace('\t', ' '))
                    f.write('\t'.join(row) + '\n')
        self.logger.info(f"✓ Created summary TSV: {summary_file}")

        # JSON summary
        json_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'kaptive_version': self.metadata.get('kaptive_version', 'Unknown'),
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results)
            },
            'genome_summaries': {},
            'statistics': self._calculate_statistics(all_results)
        }
        for gn, res in all_results.items():
            json_summary['genome_summaries'][gn] = {
                'status': res.get('status', 'unknown'),
                'hit_count': res.get('hit_count', 0),
                'k_loci': sum(1 for h in res.get('hits', []) if h.get('Locus Type') == 'K'),
                'o_loci': sum(1 for h in res.get('hits', []) if h.get('Locus Type') == 'O'),
                'k_types': list(set(h.get('Best match type', '') for h in res.get('hits', []) if h.get('Locus Type') == 'K')),
                'o_types': list(set(h.get('Best match type', '') for h in res.get('hits', []) if h.get('Locus Type') == 'O')),
                'error': res.get('error', '') if res.get('status') == 'failed' else ''
            }
        json_file = os.path.join(output_base, "klebo_kaptive_summary.json")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(json_summary, f, indent=2)
        self.logger.info(f"✓ Created summary JSON: {json_file}")

        # HTML summary with clean table (no gene lists)
        self._create_summary_html_report(all_results, output_base)

    def _calculate_statistics(self, all_results: Dict[str, Any]) -> Dict:
        stats = {
            'total_genomes': len(all_results),
            'successful_genomes': sum(1 for r in all_results.values() if r.get('status') == 'success'),
            'total_k_loci': 0,
            'total_o_loci': 0,
            'genomes_with_k': 0,
            'genomes_with_o': 0,
            'genomes_with_both': 0,
            'k_type_distribution': defaultdict(int),
            'o_type_distribution': defaultdict(int),
            'k_locus_distribution': defaultdict(int),
            'o_locus_distribution': defaultdict(int)
        }
        for r in all_results.values():
            if r.get('status') == 'success':
                k_hits = [h for h in r.get('hits', []) if h.get('Locus Type') == 'K']
                o_hits = [h for h in r.get('hits', []) if h.get('Locus Type') == 'O']
                stats['total_k_loci'] += len(k_hits)
                stats['total_o_loci'] += len(o_hits)
                if k_hits:
                    stats['genomes_with_k'] += 1
                    for h in k_hits:
                        stats['k_type_distribution'][h.get('Best match type', 'Unknown')] += 1
                        stats['k_locus_distribution'][h.get('Best match locus', 'Unknown')] += 1
                if o_hits:
                    stats['genomes_with_o'] += 1
                    for h in o_hits:
                        stats['o_type_distribution'][h.get('Best match type', 'Unknown')] += 1
                        stats['o_locus_distribution'][h.get('Best match locus', 'Unknown')] += 1
                if k_hits and o_hits:
                    stats['genomes_with_both'] += 1
        stats['k_type_distribution'] = dict(stats['k_type_distribution'])
        stats['o_type_distribution'] = dict(stats['o_type_distribution'])
        stats['k_locus_distribution'] = dict(stats['k_locus_distribution'])
        stats['o_locus_distribution'] = dict(stats['o_locus_distribution'])
        return stats

    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        stats = self._calculate_statistics(all_results)
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Build table rows: Genome, Status, K Locus, K Type, K Identity, K Coverage, K Confidence,
        # O Locus, O Type, O Identity, O Coverage, O Confidence
        table_rows = ""
        for gn, res in all_results.items():
            if res.get('status') == 'success':
                k_locus = "None"
                k_type = "None"
                k_identity = ""
                k_coverage = ""
                k_confidence = ""
                o_locus = "None"
                o_type = "None"
                o_identity = ""
                o_coverage = ""
                o_confidence = ""
                for h in res.get('hits', []):
                    if h.get('Locus Type') == 'K':
                        k_locus = h.get('Best match locus', 'Unknown')
                        k_type = h.get('Best match type', 'Unknown')
                        k_identity = h.get('Identity', '')
                        k_coverage = h.get('Coverage', '')
                        k_confidence = h.get('Match confidence', '')
                    elif h.get('Locus Type') == 'O':
                        o_locus = h.get('Best match locus', 'Unknown')
                        o_type = h.get('Best match type', 'Unknown')
                        o_identity = h.get('Identity', '')
                        o_coverage = h.get('Coverage', '')
                        o_confidence = h.get('Match confidence', '')
                status = "Success"
                row_class = "success-row"
            else:
                k_locus = "N/A"
                k_type = "N/A"
                k_identity = ""
                k_coverage = ""
                k_confidence = ""
                o_locus = "N/A"
                o_type = "N/A"
                o_identity = ""
                o_coverage = ""
                o_confidence = ""
                status = "Failed"
                row_class = "failed-row"

            table_rows += f"""
                <tr class="{row_class}">
                    <td><strong>{gn}</strong></td>
                    <td>{status}</td>
                    <td class="type-cell">{k_locus}</td>
                    <td class="type-cell">{k_type}</td>
                    <td class="num-cell">{k_identity}</td>
                    <td class="num-cell">{k_coverage}</td>
                    <td>{k_confidence}</td>
                    <td class="type-cell">{o_locus}</td>
                    <td class="type-cell">{o_type}</td>
                    <td class="num-cell">{o_identity}</td>
                    <td class="num-cell">{o_coverage}</td>
                    <td>{o_confidence}</td>
                  </tr>
            """

        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>KLEBOSCOPE - Kaptive Summary Report</title>
    <meta charset="UTF-8">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1600px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid rgba(0,119,190,0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            white-space: pre;
            color: #0077be;
            text-shadow: 0 0 10px rgba(0,119,190,0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255,255,255,0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .card {{
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.2);
        }}
        .card h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            background: white;
            font-size: 13px;
        }}
        .summary-table th, .summary-table td {{
            padding: 10px 8px;
            text-align: left;
            border-bottom: 1px solid #e5e7eb;
            vertical-align: top;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            font-weight: 600;
            cursor: pointer;
        }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #e0f2fe; }}
        .success-row {{ border-left: 4px solid #28a745; }}
        .failed-row {{ border-left: 4px solid #dc3545; background-color: #f8d7da; }}
        .type-cell {{ font-weight: bold; color: #1e40af; }}
        .num-cell {{ font-family: monospace; }}
        .summary-stats {{
            display: flex;
            gap: 20px;
            justify-content: center;
            flex-wrap: wrap;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            flex: 1;
            min-width: 150px;
        }}
        .k-stat-card {{ background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%); }}
        .o-stat-card {{ background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%); }}
        .controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        .search-box input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #3b82f6;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }}
        .export-buttons button.print {{ background: #10b981; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,0,0,0.3);
            border-radius: 10px;
        }}
        .timestamp {{ color: #fbbf24; }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255,255,255,0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 1200px) {{
            .summary-table {{ font-size: 11px; }}
        }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="quote-container" id="quoteContainer">
            <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
            <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
        </div>
    </div>

    <div class="card">
        <h2>📊 Overall Kaptive Summary</h2>
        <div class="summary-stats">
            <div class="stat-card"><div class="stat-label">Total Genomes</div><div class="stat-value" style="font-size:2em;">{stats['total_genomes']}</div></div>
            <div class="stat-card"><div class="stat-label">Successful Analyses</div><div class="stat-value" style="font-size:2em;">{stats['successful_genomes']}</div></div>
            <div class="stat-card k-stat-card"><div class="stat-label">Genomes with K Locus</div><div class="stat-value" style="font-size:2em;">{stats['genomes_with_k']}</div></div>
            <div class="stat-card o-stat-card"><div class="stat-label">Genomes with O Locus</div><div class="stat-value" style="font-size:2em;">{stats['genomes_with_o']}</div></div>
        </div>
        <p><strong>Kaptive Version:</strong> {self.metadata.get('kaptive_version', 'Unknown')} | <strong>Database:</strong> {self.metadata['database_version']}</p>
    </div>

    <div class="card">
        <h2>🔍 Detailed Results by Genome</h2>
        <div class="controls">
            <div class="search-box"><input type="text" id="search" placeholder="Search genomes, locus, types..." onkeyup="searchTable()"></div>
            <div class="export-buttons"><button onclick="exportToCSV()">📥 Export CSV</button><button onclick="printReport()" class="print">🖨️ Print</button></div>
        </div>
        <div class="table-container">
            <table class="summary-table" id="summary-table">
                <thead>
                    <tr>
                        <th onclick="sortTable(0)">Genome</th>
                        <th onclick="sortTable(1)">Status</th>
                        <th onclick="sortTable(2)">K Locus</th>
                        <th onclick="sortTable(3)">K Type</th>
                        <th onclick="sortTable(4)">K Identity</th>
                        <th onclick="sortTable(5)">K Coverage</th>
                        <th onclick="sortTable(6)">K Confidence</th>
                        <th onclick="sortTable(7)">O Locus</th>
                        <th onclick="sortTable(8)">O Type</th>
                        <th onclick="sortTable(9)">O Identity</th>
                        <th onclick="sortTable(10)">O Coverage</th>
                        <th onclick="sortTable(11)">O Confidence</th>
                    </thead>
                <tbody>
                    {table_rows}
                </tbody>
             </table>
        </div>
    </div>

    <div class="card">
        <h2>📁 Generated Files</h2>
        <ul>
            <li><strong>klebo_kaptive_summary.tsv</strong> – Complete Kaptive data for all genomes</li>
            <li><strong>klebo_kaptive_summary.json</strong> – JSON summary</li>
            <li><strong>Individual genome HTML reports</strong> – Detailed per genome</li>
            <li><strong>Individual genome JSON reports</strong> – JSON per genome</li>
            <li><strong>This summary report</strong> – Cross‑genome overview</li>
        </ul>
    </div>

    <div class="footer">
        <p><strong>KLEBOSCOPE</strong> – Kaptive K/O Analysis Summary</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</div>

<script>
    function sortTable(n) {{
        const table = document.getElementById("summary-table");
        const tbody = table.querySelector('tbody');
        const rows = Array.from(tbody.querySelectorAll('tr'));
        const isNumber = [4,5,9,10].includes(n); // identity and coverage columns
        rows.sort((a, b) => {{
            let x = a.children[n].textContent.trim();
            let y = b.children[n].textContent.trim();
            if (isNumber) {{
                x = parseFloat(x) || 0;
                y = parseFloat(y) || 0;
                return x - y;
            }}
            return x.localeCompare(y);
        }});
        rows.forEach(row => tbody.appendChild(row));
    }}
    function searchTable() {{
        const input = document.getElementById("search");
        const filter = input.value.toLowerCase();
        const table = document.getElementById("summary-table");
        const rows = table.getElementsByTagName("tr");
        let visible = 0;
        for (let i = 1; i < rows.length; i++) {{
            const row = rows[i];
            const text = row.textContent.toLowerCase();
            if (text.includes(filter)) {{
                row.style.display = "";
                visible++;
            }} else {{
                row.style.display = "none";
            }}
        }}
    }}
    function exportToCSV() {{
        const table = document.getElementById("summary-table");
        const rows = table.getElementsByTagName("tr");
        let csv = [];
        for (let i = 0; i < rows.length; i++) {{
            const row = rows[i];
            if (row.style.display !== "none") {{
                const cells = row.getElementsByTagName("th").length ? row.getElementsByTagName("th") : row.getElementsByTagName("td");
                const rowData = [];
                for (let c = 0; c < cells.length; c++) {{
                    rowData.push(cells[c].textContent.trim().replace(/,/g, ';'));
                }}
                csv.push(rowData.join(','));
            }}
        }}
        const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'kaptive_summary.csv';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }}
    function printReport() {{ window.print(); }}
    const quotes = {json.dumps(self.science_quotes)};
    let currentQuote = 0;
    function rotateQuote() {{
        document.getElementById('quoteText').innerHTML = '"' + quotes[currentQuote].text + '"';
        document.getElementById('quoteAuthor').innerHTML = '— ' + quotes[currentQuote].author;
        currentQuote = (currentQuote + 1) % quotes.length;
    }}
    setInterval(rotateQuote, 10000);
</script>
</body>
</html>"""
        # Write the HTML
        out_file = os.path.join(output_base, "klebo_kaptive_summary.html")
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write(html)
        self.logger.info(f"✓ Created summary HTML report: {out_file}")

def main():
    parser = argparse.ArgumentParser(description='Kleboscope Kaptive K/O Analysis - K. pneumoniae Capsule and Lipopolysaccharide Typing')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file or pattern (e.g., "*.fna")')
    parser.add_argument('-o', '--output', default='kaptive_results', help='Output directory')
    parser.add_argument('--k-db', default='kp_k', help='Kaptive K locus database name (default: kp_k)')
    parser.add_argument('--o-db', default='kp_o', help='Kaptive O locus database name (default: kp_o)')
    args = parser.parse_args()

    print("\n" + "="*80)
    print("🧬 Kleboscope Kaptive K/O Locus Analysis")
    print("="*80)
    print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
    print(f"Affiliation: University of Ghana Medical School")
    print("="*80)

    executor = KleboscopeKaptive(k_db=args.k_db, o_db=args.o_db)
    # Get Kaptive version once
    executor.metadata['kaptive_version'] = executor._get_kaptive_version()

    try:
        results = executor.process_multiple_genomes(args.input, args.output)
        executor.logger.info("\n" + "="*50)
        executor.logger.info("📊 K. pneumoniae Kaptive Analysis FINAL SUMMARY")
        executor.logger.info("="*50)
        total_k = sum(r.get('hit_count', 0) for r in results.values())
        successful = sum(1 for r in results.values() if r.get('status') == 'success')
        executor.logger.info(f"Genomes processed: {len(results)}")
        executor.logger.info(f"Successful analyses: {successful}")
        executor.logger.info(f"Total K/O loci detected: {total_k}")
        executor.logger.info(f"Results saved to: {args.output}")
        random_quote = executor.get_random_quote()
        executor.logger.info(f"\n💡 {random_quote['text']} - {random_quote['author']}")
    except Exception as e:
        executor.logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()