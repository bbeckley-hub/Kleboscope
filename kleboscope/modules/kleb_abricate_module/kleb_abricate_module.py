#!/usr/bin/env python3
"""
Kleboscope ABRicate Standalone Module
Comprehensive ABRicate analysis for Klebsiella pneumoniae with HTML, TSV, and JSON reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
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
#import re
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any

class AbricateExecutor:
    """ABRicate executor for Klebsiella with comprehensive reporting"""
    
    def __init__(self, cpus: int = 1):
        self.logger = self._setup_logging()
        self.cpus = cpus
        self.logger.info(f"Initialized with {self.cpus} threads allocated by orchestrator.")
        
        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf', 'bacmet2'
        ]
        
        # Load externalized gene sets
        self._load_gene_dictionaries()
        
        self.metadata = {
            "tool_name": "Kleboscope ABRicate",
            "version": "2.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Kleboscope turns genomic complexity into actionable insights for AMR surveillance.", "author": "Brown Beckley"}
        ]
        
        self.ascii_art = """
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
"""

    def _load_gene_dictionaries(self):
        db_path = Path(__file__).parent / "kleb_genes.json"
        try:
            with open(db_path, 'r') as f:
                gene_dbs = json.load(f)
                self.critical_resistance_genes = set(gene_dbs.get("critical_resistance_genes", []))
                self.high_risk_virulence_genes = set(gene_dbs.get("high_risk_virulence_genes", []))
                self.beta_lactamase_genes = set(gene_dbs.get("beta_lactamase_genes", []))
            self.logger.info("Successfully loaded gene dictionaries from kleb_genes.json")
        except FileNotFoundError:
            self.logger.error("kleb_genes.json not found. Creating empty sets.")
            self.critical_resistance_genes = set()
            self.high_risk_virulence_genes = set()
            self.beta_lactamase_genes = set()

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)

    # =========================================================================
    # CORE LOGIC
    # =========================================================================
    def check_abricate_installed(self) -> bool:
        try:
            result = subprocess.run(['abricate', '--version'], capture_output=True, text=True, check=True)
            self.logger.info("✓ ABRicate installed: %s", result.stdout.strip())
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found.")
            return False
            
    def setup_abricate_databases(self):
        available_dbs = []
        try:
            check = subprocess.run(['abricate', '--list'], capture_output=True, text=True, check=True)
            for db in self.required_databases:
                if db in check.stdout: available_dbs.append(db)
            self.required_databases = available_dbs
        except Exception as e:
            self.logger.error("Error setting up databases: %s", e)

    def run_abricate_single_db(self, genome_file: str, database: str, output_dir: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"abricate_{database}.txt")
        cmd = ['abricate', genome_file, '--db', database, '--minid', '80', '--mincov', '80', '--threads', '1']
        try:
            with open(output_file, 'w') as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            hits = self._parse_abricate_output(output_file)
            self._create_database_html_report(genome_name, database, hits, output_dir)
            return {'database': database, 'genome': genome_name, 'output_file': output_file, 'hits': hits, 'hit_count': len(hits), 'status': 'success'}
        except subprocess.CalledProcessError:
            return {'database': database, 'genome': genome_name, 'output_file': output_file, 'hits': [], 'hit_count': 0, 'status': 'failed'}

    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        hits = []
        try:
            with open(abricate_file, 'r') as f: lines = [l.strip() for l in f.readlines() if l.strip()]
            if not lines: return hits
            headers = lines[0].replace('#', '').split('\t')
            for line in lines[1:]:
                parts = line.split('\t')
                parts.extend([''] * (len(headers) - len(parts))) # pad if short
                hit = dict(zip(headers, parts))
                hits.append({
                    'genome': hit.get('FILE', ''), 'gene': hit.get('GENE', ''),
                    'coverage_percent': hit.get('%COVERAGE', ''), 'identity_percent': hit.get('%IDENTITY', ''),
                    'database': hit.get('DATABASE', ''), 'accession': hit.get('ACCESSION', ''),
                    'product': hit.get('PRODUCT', '')
                })
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
        return hits

    def analyze_klebsiella_genes(self, all_hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'critical_resistance_genes': [], 'high_risk_virulence_genes': [],
            'beta_lactamase_genes': [], 'other_genes': [], 'resistance_classes': {},
            'total_hits': len(all_hits)
        }
        for hit in all_hits:
            gene, product = hit['gene'], hit['product']
            if any(c in gene for c in self.critical_resistance_genes):
                analysis['critical_resistance_genes'].append(hit)
            elif any(v in gene for v in self.high_risk_virulence_genes):
                analysis['high_risk_virulence_genes'].append(hit)
            elif any(b in gene for b in self.beta_lactamase_genes):
                analysis['beta_lactamase_genes'].append(hit)
            else:
                analysis['other_genes'].append(hit)
                
            res_class = self._classify_resistance(product)
            if res_class:
                if res_class not in analysis['resistance_classes']: analysis['resistance_classes'][res_class] = []
                analysis['resistance_classes'][res_class].append(hit)
                
        analysis['total_critical_resistance'] = len(analysis['critical_resistance_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        analysis['total_beta_lactamase'] = len(analysis['beta_lactamase_genes'])
        return analysis

    def _classify_resistance(self, product: str) -> str:
        p = product.lower()
        if any(t in p for t in ['carbapenem', 'kpc', 'ndm', 'imp', 'vim', 'oxa']): return 'Carbapenem resistance'
        if any(t in p for t in ['esbl', 'ctx-m', 'shv', 'tem']): return 'ESBL/AmpC resistance'
        if any(t in p for t in ['colistin', 'mcr']): return 'Colistin resistance'
        return 'Other resistance'

    # =========================================================================
    # HTML TEMPLATING & REPORT GENERATION
    # =========================================================================
    def _get_base_html(self, title: str, body_content: str) -> str:
        """Centralized HTML wrapper to eliminate repetitive CSS/JS boilerplate."""
        quote = random.choice(self.science_quotes)
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%); font-family: 'Segoe UI', Tahoma, Geneva, sans-serif; color: #fff; padding: 20px; min-height: 100vh; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px; border: 2px solid rgba(0,119,190,0.3); }}
        .ascii-art {{ font-family: monospace; font-size: 10px; color: #0077be; white-space: pre; overflow-x: auto; }}
        .quote-container {{ background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px; text-align: center; border: 1px solid rgba(255,255,255,0.2); }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{ background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px; margin-bottom: 20px; }}
        .report-section h2 {{ color: #1e3a8a; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; }}
        .metrics-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; }}
        .metric-card {{ background: linear-gradient(135deg, #0077be 0%, #005a8c 100%); color: white; padding: 20px; border-radius: 8px; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .summary-table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px; }}
        .summary-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; }}
        .summary-table td {{ padding: 12px; border-bottom: 1px solid #e5e7eb; }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .high-risk {{ background-color: #fef3c7; }}
        .present {{ background-color: #d1fae5; }}
        .risk-badge {{ background: #dc2626; color: white; padding: 5px 10px; border-radius: 15px; margin: 2px; display: inline-block; }}
        .warning-badge {{ background: #f59e0b; color: black; padding: 5px 10px; border-radius: 15px; margin: 2px; display: inline-block; }}
        .footer {{ text-align: center; padding: 20px; background: rgba(0,0,0,0.3); border-radius: 10px; margin-top: 30px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
            <div class="quote-container">
                <div class="quote-text">"{quote['text']}"</div>
                <div class="quote-author">— {quote['author']}</div>
            </div>
        </div>
        {body_content}
        <div class="footer">
            <p><strong>KLEBOSCOPE</strong> - ABRicate Analysis Module</p>
            <p>Generated: {self.metadata['analysis_date']}</p>
        </div>
    </div>
</body>
</html>"""

    def _build_html_table(self, hits: List[Dict], row_class_key: str = "present") -> str:
        if not hits: return "<p>No genes detected.</p>"
        html = '<div style="overflow-x: auto;"><table class="summary-table"><thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead><tbody>'
        for hit in hits:
            html += f'<tr class="{row_class_key}"><td><strong>{hit["gene"]}</strong></td><td>{hit["product"]}</td><td>{hit["database"]}</td><td>{hit["coverage_percent"]}%</td><td>{hit["identity_percent"]}%</td></tr>'
        return html + "</tbody></table></div>"

    def _create_database_html_report(self, genome_name: str, database: str, hits: List[Dict], output_dir: str):
        content = f"""
        <div class="report-section">
            <h2>📊 Database Information: {database.upper()}</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{len(hits)}</div><div>Total Hits</div></div>
                <div class="metric-card"><div class="metric-value">{genome_name}</div><div>Genome</div></div>
            </div>
        </div>
        <div class="report-section"><h2>🔍 Genes Detected</h2>{self._build_html_table(hits)}</div>
        """
        with open(os.path.join(output_dir, f"abricate_{database}_report.html"), 'w') as f:
            f.write(self._get_base_html(f"KLEBOSCOPE - {database.upper()}", content))

    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        all_hits = [h for r in results.values() for h in r['hits']]
        analysis = self.analyze_klebsiella_genes(all_hits)
        
        content = f"""
        <div class="report-section">
            <h2>📊 Klebsiella AMR/Virulence Summary ({genome_name})</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{analysis['total_hits']}</div><div>Total Genes</div></div>
                <div class="metric-card"><div class="metric-value">{analysis['total_critical_resistance']}</div><div>Critical Res</div></div>
                <div class="metric-card"><div class="metric-value">{analysis['total_high_risk_virulence']}</div><div>High-Risk Virulence</div></div>
            </div>
        </div>
        """
        if analysis['critical_resistance_genes']:
            content += f'<div class="report-section"><h2>🔴 Critical Resistance Genes</h2>{self._build_html_table(analysis["critical_resistance_genes"], "critical")}</div>'
        if analysis['high_risk_virulence_genes']:
            content += f'<div class="report-section"><h2>🟡 High-Risk Virulence Genes</h2>{self._build_html_table(analysis["high_risk_virulence_genes"], "high-risk")}</div>'
        if analysis['beta_lactamase_genes']:
            content += f'<div class="report-section"><h2>🔵 Beta-Lactamase Genes</h2>{self._build_html_table(analysis["beta_lactamase_genes"][:1000], "present")}</div>'
        
        with open(os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html"), 'w') as f:
            f.write(self._get_base_html(f"KLEBOSCOPE - Comprehensive ({genome_name})", content))

    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        unique_genomes = list(set(h['genome'] for h in hits))
        content = f"""
        <div class="report-section">
            <h2>📊 Batch Summary: {database.upper()}</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{len(hits)}</div><div>Total Hits</div></div>
                <div class="metric-card"><div class="metric-value">{len(unique_genomes)}</div><div>Genomes</div></div>
            </div>
        </div>
        """
        with open(os.path.join(output_base, f"klebo_{database}_summary_report.html"), 'w') as f:
            f.write(self._get_base_html(f"KLEBOSCOPE - Batch {database.upper()}", content))

    # =========================================================================
    # PROCESS RUNNERS
    # =========================================================================
    def process_single_genome(self, genome_file: str, output_base: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        results = {db: self.run_abricate_single_db(genome_file, db, output_base) for db in self.required_databases}
        self.create_comprehensive_html_report(genome_name, results, output_base)
        return {'genome': genome_name, 'results': results, 'total_hits': sum(r['hit_count'] for r in results.values())}

    def process_multiple_genomes(self, genome_pattern: str, output_base: str) -> Dict[str, Any]:
        if not self.check_abricate_installed(): raise RuntimeError("ABRicate not installed")
        self.setup_abricate_databases()
        
        genome_files = list(set(glob.glob(genome_pattern) + glob.glob(f"{genome_pattern}.*")))
        if not genome_files: raise FileNotFoundError("No FASTA files found.")
        
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        
        if len(genome_files) > 1 and self.cpus > 1:
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                futures = {executor.submit(self.process_single_genome, g, output_base): g for g in genome_files}
                for future in as_completed(futures):
                    res = future.result()
                    all_results[res['genome']] = res
        else:
            for g in genome_files:
                res = self.process_single_genome(g, output_base)
                all_results[res['genome']] = res
                
        # Generate Batch Summaries
        for db in self.required_databases:
            hits = [h for r in all_results.values() for h in r['results'][db]['hits']]
            if hits: self._create_database_summary_html(db, hits, output_base)
            
        return all_results

def main():
    parser = argparse.ArgumentParser(description='Kleboscope ABRicate Analysis - Orchestrator Compliant')
    parser.add_argument('-i', '--input', required=True, help='File pattern or path for genomes')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of CPU cores')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    args = parser.parse_args()
    
    executor = AbricateExecutor(cpus=args.threads)
    try:
        results = executor.process_multiple_genomes(args.input, args.output)
        executor.logger.info("Analysis Complete. Results saved to: %s", args.output)
    except Exception as e:
        executor.logger.error("Analysis failed: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()