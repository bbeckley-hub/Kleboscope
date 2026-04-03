#!/usr/bin/env python3
"""
Kleboscope Kaptive K/O Locus Analysis - K. pneumoniae Capsule and Lipopolysaccharide Typing
Comprehensive Kaptive analysis for K. pneumoniae with beautiful HTML reporting
Author: Brown Beckley
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Version: 2.0.0
"""

import subprocess
import sys
import os
import glob
import logging
import json
import random
import argparse
import re
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime
from collections import defaultdict

class KleboscopeKaptive:
    def __init__(self, k_db="kp_k", o_db="kp_o"):
        self.logger = self._setup_logging()
        self.k_db = k_db
        self.o_db = o_db
        self.kaptive_available = None

        self.metadata = {
            "tool_name": "Kleboscope Kaptive K/O Analysis",
            "version": "2.0.0",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "database_version": f"{k_db} / {o_db}"
        }

        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "Kleboscope turns genomic complexity into actionable insights.", "author": "Brown Beckley"}
        ]

        self.ascii_art = r"""
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
"""

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)

    # =========================================================================
    # CORE LOGIC
    # =========================================================================
    def _get_kaptive_version(self) -> str:
        try:
            result = subprocess.run(["kaptive", "--version"], capture_output=True, text=True, check=True)
            match = re.search(r'v?(\d+\.\d+\.\d+)', result.stdout.strip())
            return match.group(1) if match else result.stdout.strip()
        except Exception:
            return "Unknown"

    def check_kaptive_installed(self) -> bool:
        if self.kaptive_available is not None: return self.kaptive_available
        try:
            subprocess.run(["kaptive", "--version"], capture_output=True, text=True, check=True)
            self.kaptive_available = True
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("Kaptive check failed. Please install Kaptive.")
            self.kaptive_available = False
            return False

    def run_kaptive_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        os.makedirs(output_dir, exist_ok=True)

        k_tsv = os.path.join(output_dir, f"{genome_name}_K.tsv")
        o_tsv = os.path.join(output_dir, f"{genome_name}_O.tsv")
        combined_tsv = os.path.join(output_dir, f"{genome_name}_combined.tsv")

        try:
            subprocess.run(["kaptive", "assembly", self.k_db, genome_file, "-o", k_tsv], capture_output=True, text=True, check=True)
            subprocess.run(["kaptive", "assembly", self.o_db, genome_file, "-o", o_tsv], capture_output=True, text=True, check=True)

            self._combine_k_o_results(k_tsv, o_tsv, combined_tsv)
            hits = self._parse_full_kaptive_output(combined_tsv)

            self._create_kaptive_html_report(genome_name, hits, output_dir)
            self._create_kaptive_json_report(genome_name, hits, output_dir)

            return {'genome': genome_name, 'output_file': combined_tsv, 'hits': hits, 'hit_count': len(hits), 'status': 'success'}
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Kaptive failed for {genome_name}: {e}")
            return {'genome': genome_name, 'output_file': combined_tsv, 'hits': [], 'hit_count': 0, 'status': 'failed', 'error': str(e)}

    def _combine_k_o_results(self, k_tsv: str, o_tsv: str, combined_tsv: str):
        try:
            with open(combined_tsv, 'w') as out:
                if os.path.exists(k_tsv):
                    with open(k_tsv, 'r') as f: out.writelines(f.readlines())
                if os.path.exists(o_tsv):
                    with open(o_tsv, 'r') as f:
                        lines = f.readlines()
                        if lines:
                            if not os.path.exists(k_tsv): out.write(lines[0])
                            out.writelines(lines[1:])
        except Exception as e:
            self.logger.error(f"Error combining K/O results: {e}")

    def _parse_full_kaptive_output(self, combined_file: str) -> List[Dict]:
        hits = []
        try:
            with open(combined_file, 'r') as f: lines = f.readlines()
            if len(lines) < 2: return hits
            
            headers = lines[0].strip().split('\t')
            for line in lines[1:]:
                parts = line.strip().split('\t')
                parts.extend([''] * (len(headers) - len(parts)))
                hit = dict(zip(headers, parts))
                
                locus = hit.get('Best match locus', '')
                hit['Locus Type'] = 'K' if locus.startswith('KL') else 'O' if locus.startswith('O') else 'Unknown'
                hit['Gene Details Parsed'] = self._parse_gene_details(hit.get('Expected genes in locus, details', ''))
                hit['Other Genes Parsed'] = self._parse_gene_details(hit.get('Other genes in locus, details', ''))
                hits.append(hit)
        except Exception as e:
            self.logger.error(f"Error parsing {combined_file}: {e}")
        return hits

    def _parse_gene_details(self, details_str: str) -> List[Dict]:
        genes = []
        if not details_str: return genes
        for part in details_str.split(';'):
            if part:
                sub = part.split(',')
                if len(sub) >= 3:
                    genes.append({'name': sub[0], 'identity': sub[1].replace('%', ''), 'coverage': sub[2].replace('%', ''), 'raw': part})
        return genes

    def _calculate_statistics(self, all_results: Dict[str, Any]) -> Dict:
        stats = {
            'total_genomes': len(all_results), 'successful_genomes': sum(1 for r in all_results.values() if r.get('status') == 'success'),
            'genomes_with_k': 0, 'genomes_with_o': 0
        }
        for r in all_results.values():
            if r.get('status') == 'success':
                if any(h.get('Locus Type') == 'K' for h in r.get('hits', [])): stats['genomes_with_k'] += 1
                if any(h.get('Locus Type') == 'O' for h in r.get('hits', [])): stats['genomes_with_o'] += 1
        return stats

    # =========================================================================
    # HTML TEMPLATING
    # =========================================================================
    def _get_base_html(self, title: str, body_content: str) -> str:
        quote = random.choice(self.science_quotes)
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%); font-family: 'Segoe UI', Tahoma, Geneva, sans-serif; color: #fff; padding: 20px; }}
        .container {{ max-width: 1600px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px; border: 2px solid rgba(0,119,190,0.3); }}
        .ascii-art {{ font-family: monospace; font-size: 10px; color: #0077be; white-space: pre; overflow-x: auto; }}
        .quote-container {{ background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px; text-align: center; border: 1px solid rgba(255,255,255,0.2); margin-bottom: 20px; }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .card {{ background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px; margin-bottom: 20px; }}
        .card h2 {{ color: #1e3a8a; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; }}
        .summary-stats {{ display: flex; gap: 20px; justify-content: space-around; flex-wrap: wrap; }}
        .stat-card {{ background: linear-gradient(135deg, #0077be 0%, #005a8c 100%); color: white; padding: 20px; border-radius: 8px; text-align: center; flex: 1; }}
        .k-stat-card {{ background: linear-gradient(135deg, #28a745 0%, #1e7e34 100%); }}
        .o-stat-card {{ background: linear-gradient(135deg, #17a2b8 0%, #117a8b 100%); }}
        .gene-table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 13px; background: white; }}
        .gene-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; }}
        .gene-table td {{ padding: 12px; border-bottom: 1px solid #e5e7eb; }}
        .gene-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .k-locus-row {{ background-color: #d4edda; }}
        .o-locus-row {{ background-color: #d1ecf1; }}
        .gene-item {{ display: inline-block; background: #e9ecef; padding: 3px 8px; margin: 2px; border-radius: 12px; font-size: 0.85em; }}
        .footer {{ text-align: center; padding: 20px; background: rgba(0,0,0,0.8); border-radius: 10px; margin-top: 30px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
            <div class="quote-container"><div class="quote-text">"{quote['text']}"</div><div class="quote-author">— {quote['author']}</div></div>
        </div>
        {body_content}
        <div class="footer"><p>KLEBOSCOPE Kaptive | Generated: {self.metadata['analysis_date']}</p></div>
    </div>
</body>
</html>"""

    def _create_kaptive_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        table_keys = [k for k in hits[0].keys() if k not in ['Gene Details Parsed', 'Other Genes Parsed']] if hits else []
        
        table_rows = ""
        gene_details_html = ""
        
        for hit in hits:
            locus_type = hit.get('Locus Type', '')
            row_class = "k-locus-row" if locus_type == 'K' else "o-locus-row" if locus_type == 'O' else ""
            table_rows += f'<tr class="{row_class}">' + "".join([f"<td>{str(hit.get(k, ''))[:200]}</td>" for k in table_keys]) + "</tr>"
            
            color = "#28a745" if locus_type == 'K' else "#17a2b8"
            genes_html = lambda genes: "".join([f'<span class="gene-item">{g["name"]}: {g["identity"]}% ID, {g["coverage"]}% COV</span>' for g in genes]) if genes else "<span>None</span>"
            
            gene_details_html += f"""
            <div style="margin-top: 15px; padding: 15px; background: #f8f9fa; border-left: 4px solid {color};">
                <h4 style="color: {color};">{hit.get('Best match locus', 'Unknown')} ({locus_type} Locus)</h4>
                <p><strong>Expected genes:</strong> {genes_html(hit.get('Gene Details Parsed', []))}</p>
                <p><strong>Other genes:</strong> {genes_html(hit.get('Other Genes Parsed', []))}</p>
            </div>
            """

        headers_html = "".join([f"<th>{k}</th>" for k in table_keys])
        content = f"""
        <div class="card">
            <h2>📊 K. pneumoniae K/O Locus Summary ({genome_name})</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total Loci</h3><h2>{len(hits)}</h2></div>
                <div class="stat-card k-stat-card"><h3>K Locus</h3><h2>{sum(1 for h in hits if h.get('Locus Type') == 'K')}</h2></div>
                <div class="stat-card o-stat-card"><h3>O Locus</h3><h2>{sum(1 for h in hits if h.get('Locus Type') == 'O')}</h2></div>
            </div>
        </div>
        <div class="card">
            <h2>🔬 Complete Kaptive Results</h2>
            <div style="overflow-x:auto;"><table class="gene-table"><thead><tr>{headers_html}</tr></thead><tbody>{table_rows}</tbody></table></div>
            <h3 style="margin-top: 30px;">🧬 Detailed Gene Information</h3>
            {gene_details_html}
        </div>
        """
        with open(os.path.join(output_dir, f"{genome_name}_kaptive_report.html"), 'w', encoding='utf-8') as f:
            f.write(self._get_base_html(f"Kaptive Report - {genome_name}", content))

    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        stats = self._calculate_statistics(all_results)
        
        table_rows = ""
        for gn, res in all_results.items():
            if res.get('status') == 'success':
                k, o = {}, {}
                for h in res.get('hits', []):
                    if h.get('Locus Type') == 'K': k = h
                    elif h.get('Locus Type') == 'O': o = h
                
                table_rows += f"""<tr>
                    <td><strong>{gn}</strong></td><td>Success</td>
                    <td>{k.get('Best match locus', 'None')}</td><td>{k.get('Best match type', 'None')}</td><td>{k.get('Identity', '')}</td><td>{k.get('Match confidence', '')}</td>
                    <td>{o.get('Best match locus', 'None')}</td><td>{o.get('Best match type', 'None')}</td><td>{o.get('Identity', '')}</td><td>{o.get('Match confidence', '')}</td>
                </tr>"""
            else:
                table_rows += f'<tr style="background:#f8d7da; border-left:4px solid #dc3545"><td><strong>{gn}</strong></td><td>Failed</td><td colspan="10">Error processing</td></tr>'

        content = f"""
        <div class="card">
            <h2>📊 Batch Kaptive Summary</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total Genomes</h3><h2>{stats['total_genomes']}</h2></div>
                <div class="stat-card k-stat-card"><h3>Genomes with K Locus</h3><h2>{stats['genomes_with_k']}</h2></div>
                <div class="stat-card o-stat-card"><h3>Genomes with O Locus</h3><h2>{stats['genomes_with_o']}</h2></div>
            </div>
        </div>
        <div class="card">
            <h2>🔍 Detailed Results</h2>
            <div style="overflow-x:auto;">
                <table class="gene-table">
                    <thead><tr><th>Genome</th><th>Status</th><th>K Locus</th><th>K Type</th><th>K Identity</th><th>K Confidence</th><th>O Locus</th><th>O Type</th><th>O Identity</th><th>O Confidence</th></tr></thead>
                    <tbody>{table_rows}</tbody>
                </table>
            </div>
        </div>
        """
        with open(os.path.join(output_base, "klebo_kaptive_summary.html"), 'w', encoding='utf-8') as f:
            f.write(self._get_base_html("Kaptive Batch Summary", content))

    # =========================================================================
    # BATCH PROCESSORS
    # =========================================================================
    def _create_kaptive_json_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        json_data = {'metadata': {'genome': genome_name}, 'hits': hits}
        with open(os.path.join(output_dir, f"{genome_name}_kaptive_report.json"), 'w') as f: json.dump(json_data, f, indent=2)

    def create_kaptive_summary(self, all_results: Dict[str, Any], output_base: str):
        self._create_summary_html_report(all_results, output_base)
        with open(os.path.join(output_base, "klebo_kaptive_summary.json"), 'w') as f: json.dump(all_results, f, indent=2)

    def process_single_genome(self, genome_file: str, output_base: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        # Create a genome-specific directory directly under output_base
        genome_dir = os.path.join(output_base, genome_name)
        os.makedirs(genome_dir, exist_ok=True)
        
        self.logger.info(f"=== PROCESSING GENOME: {genome_name} ===")
        result = self.run_kaptive_single_genome(genome_file, genome_dir)
        
        status = "✓" if result['status'] == 'success' else "✗"
        self.logger.info(f"{status} {genome_name}: {result['hit_count']} K/O loci")
        return result

    def process_multiple_genomes(self, genome_pattern: str, output_base: str) -> Dict[str, Any]:
        genome_files = list(set(glob.glob(genome_pattern) + glob.glob(f"{genome_pattern}.fasta") + glob.glob(f"{genome_pattern}.fna")))
        if not genome_files: raise FileNotFoundError(f"No files found for {genome_pattern}")
        
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        for g in genome_files:
            all_results[Path(g).stem] = self.process_single_genome(g, output_base)
            
        self.create_kaptive_summary(all_results, output_base)
        return all_results

def main():
    parser = argparse.ArgumentParser(description='Kleboscope Kaptive - Orchestrator Compliant')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA pattern')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--k-db', default='kp_k', help='Kaptive K database')
    parser.add_argument('--o-db', default='kp_o', help='Kaptive O database')
    args = parser.parse_args()

    executor = KleboscopeKaptive(k_db=args.k_db, o_db=args.o_db)
    executor.metadata['kaptive_version'] = executor._get_kaptive_version()

    try:
        executor.process_multiple_genomes(args.input, args.output)
    except Exception as e:
        executor.logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()