#!/usr/bin/env python3
"""
Kleboscope AMRfinderPlus - K. pneumoniae Antimicrobial Resistance
Comprehensive AMR analysis with interactive per‑genome reports and clean batch summary
Author: Brown Beckley
Affiliation: University of Ghana Medical School
Version: 2.0.0 
Uses BUNDLED AMRFinderPlus 4.2.7 with database 2026-01-21.1
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import argparse
#import re
from datetime import datetime
import json
import random
from collections import defaultdict

class KleboAMRfinderPlus:
    def __init__(self, cpus: int = 1):
        self.logger = self._setup_logging()
        self.cpus = cpus
        self.logger.info(f"Initialized with {self.cpus} threads allocated by orchestrator.")
        
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.bundled_amrfinder = "amrfinder"
        db_base_dir = Path(self.module_dir) / "data" / "amrfinder_db"
        
        if (db_base_dir / "latest").exists():
            self.bundled_database = str(db_base_dir / "latest")
        elif db_base_dir.exists():
            db_folders = sorted([d for d in db_base_dir.iterdir() if d.is_dir()])
            self.bundled_database = str(db_folders[-1]) if db_folders else None
        else:
            self.bundled_database = None
        
        self._load_gene_dictionaries()
        
        self.high_risk_genes = set.union(
            self.critical_carbapenemases, self.critical_esbls,
            self.critical_colistin, self.critical_aminoglycoside, self.high_risk_resistance
        )
        self.critical_risk_genes = set.union(self.critical_carbapenemases, self.critical_colistin)

        self.metadata = {
            "tool_name": "Kleboscope AMRfinderPlus", "version": "2.0.0",
            "authors": ["Brown Beckley"], "email": "brownbeckley94@gmail.com",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7", "database_version": "2026-01-21.1"
        }
        
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning.", "author": "Albert Einstein"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Kleboscope turns genomic complexity into actionable insights.", "author": "Brown Beckley"}
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
        db_path = Path(__file__).parent / "klebo_amr_genes.json"
        try:
            with open(db_path, 'r') as f:
                dbs = json.load(f)
                self.critical_carbapenemases = set(dbs.get("critical_carbapenemases", []))
                self.critical_esbls = set(dbs.get("critical_esbls", []))
                self.critical_colistin = set(dbs.get("critical_colistin", []))
                self.critical_aminoglycoside = set(dbs.get("critical_aminoglycoside", []))
                self.high_risk_resistance = set(dbs.get("high_risk_resistance", []))
            self.logger.info("Loaded gene dictionaries from klebo_amr_genes.json")
        except FileNotFoundError:
            self.logger.error("klebo_amr_genes.json not found. Creating empty sets.")
            self.critical_carbapenemases, self.critical_esbls, self.critical_colistin, self.critical_aminoglycoside, self.high_risk_resistance = set(), set(), set(), set(), set()

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)
    
    def check_amrfinder_installed(self) -> bool:
        try:
            res = subprocess.run([self.bundled_amrfinder, '--version'], capture_output=True, text=True, check=True)
            self.logger.info(f"AMRfinderPlus version: {res.stdout.strip()}")
            return True
        except Exception: return False

    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str, threads_per_job: int) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        cmd = [self.bundled_amrfinder, '-n', genome_file, '--output', output_file, '--threads', str(threads_per_job), '--plus', '--organism', 'Klebsiella_pneumoniae']
        if os.path.exists(self.bundled_database): cmd.extend(['--database', self.bundled_database])
        
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = self._parse_amrfinder_output(output_file)
            return {'genome': genome_name, 'output_file': output_file, 'hits': hits, 'hit_count': len(hits), 'status': 'success'}
        except subprocess.CalledProcessError as e:
            self.logger.error(f"AMRfinderPlus failed for {genome_name}: {e.stderr}")
            return {'genome': genome_name, 'output_file': output_file, 'hits': [], 'hit_count': 0, 'status': 'failed'}

    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        hits = []
        try:
            with open(amrfinder_file, 'r') as f: 
                lines = [l.strip('\n\r') for l in f.readlines() if l.strip()]
            if len(lines) < 2: return hits
            headers = lines[0].split('\t')
            for line in lines[1:]:
                parts = line.split('\t')
                parts.extend([''] * (len(headers) - len(parts)))
                if len(parts) >= len(headers):
                    h = dict(zip(headers, parts))
                    hits.append({
                        'gene_symbol': h.get('Element symbol', ''), 'sequence_name': h.get('Element name', ''),
                        'class': h.get('Class', ''), 'subclass': h.get('Subclass', ''),
                        'coverage': h.get('% Coverage of reference', '').replace('%', ''),
                        'identity': h.get('% Identity to reference', '').replace('%', ''),
                        'scope': h.get('Scope', ''), 'contig_id': h.get('Contig id', ''),
                        'start': h.get('Start', ''), 'stop': h.get('Stop', ''), 'element_type': h.get('Type', ''),
                        'accession': h.get('Closest reference accession', '')
                    })
        except Exception as e:
            self.logger.error(f"Error parsing {amrfinder_file}: {e}")
        return hits

    def _analyze_kpneumoniae_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'total_genes': len(hits), 'resistance_classes': {}, 'high_risk_genes': 0, 'critical_risk_genes': 0,
            'high_risk_list': [], 'critical_risk_list': [], 'resistance_mechanisms': defaultdict(list)
        }
        for hit in hits:
            gene = hit.get('gene_symbol', '')
            if not gene: continue
            
            # Mechanism
            if any(k in gene for k in ('KPC', 'NDM', 'IMP', 'VIM', 'OXA-48', 'GES-2', 'SME')): analysis['resistance_mechanisms']['carbapenemase'].append(gene)
            elif 'mcr' in gene.lower(): analysis['resistance_mechanisms']['colistin_resistance'].append(gene)
            elif any(k in gene for k in ('CTX-M', 'SHV', 'TEM')): analysis['resistance_mechanisms']['esbl'].append(gene)
            else: analysis['resistance_mechanisms']['other_amr'].append(gene)

            # Risks
            if gene in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene not in analysis['critical_risk_list']: analysis['critical_risk_list'].append(gene)
            elif gene in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene not in analysis['high_risk_list']: analysis['high_risk_list'].append(gene)
                
            # Classes
            cls = hit.get('class', '')
            if cls:
                if cls not in analysis['resistance_classes']: analysis['resistance_classes'][cls] = []
                if gene not in analysis['resistance_classes'][cls]: analysis['resistance_classes'][cls].append(gene)
                
        return analysis

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
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px; border: 2px solid rgba(0,119,190,0.3); }}
        .ascii-art {{ font-family: monospace; font-size: 10px; color: #0077be; white-space: pre; overflow-x: auto; }}
        .quote-container {{ background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px; text-align: center; border: 1px solid rgba(255,255,255,0.2); }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .card {{ background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px; margin-bottom: 20px; }}
        .card h2 {{ color: #1e3a8a; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; }}
        .summary-stats {{ display: flex; gap: 20px; justify-content: space-around; flex-wrap: wrap; }}
        .stat-card {{ background: linear-gradient(135deg, #0077be 0%, #005a8c 100%); color: white; padding: 20px; border-radius: 8px; text-align: center; flex: 1; }}
        .critical-stat {{ background: linear-gradient(135deg, #dc3545 0%, #c82333 100%); }}
        .gene-table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px; background: white; }}
        .gene-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; }}
        .gene-table td {{ padding: 12px; border-bottom: 1px solid #e5e7eb; }}
        .gene-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .critical-row {{ background-color: #f8d7da; font-weight: bold; }}
        .high-risk-row {{ background-color: #fff3cd; }}
        .risk-badge {{ display: inline-block; background: #dc2626; color: white; padding: 4px 8px; border-radius: 12px; font-size: 0.8em; margin: 2px; }}
        .warning-badge {{ display: inline-block; background: #f59e0b; color: black; padding: 4px 8px; border-radius: 12px; font-size: 0.8em; margin: 2px; }}
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
        <div class="footer"><p>KLEBOSCOPE AMRfinderPlus | Generated: {self.metadata['analysis_date']}</p></div>
    </div>
</body>
</html>"""

    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        analysis = self._analyze_kpneumoniae_amr_results(hits)
        content = f"""
        <div class="card">
            <h2>📊 K. pneumoniae AMR Summary ({genome_name})</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total AMR Genes</h3><p style="font-size: 2em;">{analysis['total_genes']}</p></div>
                <div class="stat-card"><h3>High Risk Genes</h3><p style="font-size: 2em;">{analysis['high_risk_genes']}</p></div>
                <div class="stat-card critical-stat"><h3>Critical Risk</h3><p style="font-size: 2em;">{analysis['critical_risk_genes']}</p></div>
            </div>
        </div>
        """
        if hits:
            table = '<table class="gene-table"><thead><tr><th>Gene</th><th>Sequence Name</th><th>Class</th><th>Subclass</th><th>Coverage</th><th>Identity</th></tr></thead><tbody>'
            for h in hits:
                g = h.get('gene_symbol', '')
                cls = "critical-row" if g in self.critical_risk_genes else "high-risk-row" if g in self.high_risk_genes else ""
                table += f'<tr class="{cls}"><td><strong>{g}</strong></td><td>{h.get("sequence_name","")}</td><td>{h.get("class","")}</td><td>{h.get("subclass","")}</td><td>{h.get("coverage","")}%</td><td>{h.get("identity","")}%</td></tr>'
            content += f'<div class="card"><h2>🔬 Detailed AMR Genes</h2><div style="overflow-x:auto;">{table}</div></div>'
            
        with open(os.path.join(output_dir, f"{genome_name}_amrfinder_report.html"), 'w') as f:
            f.write(self._get_base_html(f"AMR Report - {genome_name}", content))

    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        total_genomes = len(all_results)
        total_hits = sum(r['hit_count'] for r in all_results.values())
        
        genes_per_genome = {}
        gene_frequency = defaultdict(set)
        for gn, res in all_results.items():
            genes = [h.get('gene_symbol', '') for h in res['hits'] if h.get('gene_symbol')]
            genes_per_genome[gn] = set(genes)
            for g in genes: gene_frequency[g].add(gn)
            
        t1_rows = "".join([f"<tr><td>{gn}</td><td>{len(genes)}</td><td>{', '.join(genes)}</td></tr>" for gn, genes in genes_per_genome.items()])
        
        t2_rows = ""
        for g, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            t2_rows += f"<tr><td>{g}</td><td>{len(genomes)}</td><td>{', '.join(genomes)}</td><td>Standard</td></tr>"
            
        content = f"""
        <div class="card">
            <h2>📊 Batch AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Genomes</h3><p style="font-size: 2em;">{total_genomes}</p></div>
                <div class="stat-card"><h3>Total Hits</h3><p style="font-size: 2em;">{total_hits}</p></div>
            </div>
        </div>
        <div class="card">
            <h2>Genes by Genome</h2>
            <table class="summary-table"><thead><tr><th>Genome</th><th>Count</th><th>Genes Detected</th></tr></thead><tbody>{t1_rows}</tbody></table>
        </div>
        <div class="card">
            <h2>Gene Frequency</h2>
            <table class="summary-table"><thead><tr><th>Gene</th><th>Frequency</th><th>Genomes</th><th>Risk Level</th></tr></thead><tbody>{t2_rows}</tbody></table>
        </div>
        """
        with open(os.path.join(output_base, "klebo_amrfinder_summary_report.html"), 'w', encoding='utf-8') as f:
            f.write(self._get_base_html("AMR Batch Summary", content))
    # =========================================================================
    # PROCESS RUNNERS
    # =========================================================================
    def process_single_genome(self, genome_file: str, output_base: str, threads: int) -> Dict[str, Any]:
        result = self.run_amrfinder_single_genome(genome_file, output_base, threads)
        if result['status'] == 'success':
            self._create_amrfinder_html_report(result['genome'], result['hits'], output_base)
        return result

    def process_multiple_genomes(self, genome_pattern: str, output_base: str) -> Dict[str, Any]:
        if not self.check_amrfinder_installed(): raise RuntimeError("AMRfinderPlus not installed")
        
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa", f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        genome_files = []
        for p in fasta_patterns:
            genome_files.extend(glob.glob(p))
        genome_files = list(set(genome_files))    
        if not genome_files: raise FileNotFoundError(f"No files found for {genome_pattern}")
        
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        
        # Parallel Execution Logic
        max_concurrent = min(self.cpus, len(genome_files))
        threads_per_job = max(1, self.cpus // max_concurrent)
        
        self.logger.info(f"Processing {len(genome_files)} genomes. Concurrent jobs: {max_concurrent}, Threads per job: {threads_per_job}")
        
        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            futures = {executor.submit(self.process_single_genome, g, output_base, threads_per_job): g for g in genome_files}
            for future in as_completed(futures):
                res = future.result()
                all_results[res['genome']] = res
                self.logger.info(f"✓ Completed {res['genome']}: {res['hit_count']} hits")
                
        self._create_summary_html_report(all_results, output_base)
        return all_results

def main():
    parser = argparse.ArgumentParser(description='Kleboscope AMRfinderPlus - Orchestrator Compliant')
    parser.add_argument('pattern', nargs='?', help='File pattern (Legacy Support)')
    parser.add_argument('-i', '--input', help='File pattern for genomes')
    parser.add_argument('--cpus', '-t', type=int, default=1, help='CPU cores allocated by orchestrator')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    args = parser.parse_args()
    
    pattern = args.input if args.input else args.pattern
    if not pattern:
        print("Error: Must provide an input pattern using -i")
        sys.exit(1)

    executor = KleboAMRfinderPlus(cpus=args.cpus)
    try:
        executor.process_multiple_genomes(pattern, args.output)
        executor.logger.info(f"Analysis saved to {args.output}")
    except Exception as e:
        executor.logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()