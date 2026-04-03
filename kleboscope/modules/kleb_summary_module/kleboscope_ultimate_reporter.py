#!/usr/bin/env python3
"""
Kleboscope Ultimate Reporter – Comprehensive Gene‑Centric Analysis for K. pneumoniae
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.0.0
"""

#import os
#import sys
import json
import re
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Any
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup

# =============================================================================
# HTML PARSER
# =============================================================================
class KleboHTMLParser:
    """Parse all Kleboscope HTML reports."""

    def __init__(self):
        self.abricate_databases = [
            'card', 'resfinder', 'argannot', 'vfdb', 'plasmidfinder',
            'megares', 'ncbi', 'ecoh', 'ecoli_vf', 'bacmet2'
        ]
        self.db_name_mapping = {f"klebo_{db}": db for db in self.abricate_databases}

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id).strip()
        for ext in ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']:
            if sample.endswith(ext): sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample: sample = Path(sample).name
        return sample

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables): return pd.DataFrame()
            
            table = tables[table_index]
            rows = table.find_all('tr')
            if not rows: return pd.DataFrame()
            
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                row_data = [col.get_text().strip() for col in cols]
                row_data.extend([''] * (len(headers) - len(row_data)))
                data.append(row_data[:len(headers)])
                
            if not data: return pd.DataFrame()
            df = pd.DataFrame(data, columns=headers)
            df.columns = [c.replace('\n', ' ').strip() for c in df.columns]
            return df
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()

    def parse_mlst_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing MLST: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f: html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty: return {}

            sample_col = next((c for c in df.columns if 'sample' in c.lower()), df.columns[0] if len(df.columns) > 0 else None)
            st_col = next((c for c in df.columns if 'st' in c.lower() and 'sample' not in c.lower()), None)

            results = {}
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(row[sample_col]) if sample_col and row[sample_col] else ''
                if not sample: continue
                st = 'ND'
                if st_col and pd.notna(row.get(st_col)):
                    st_val = str(row[st_col]).strip()
                    if st_val.lower() not in ['', 'nan', 'none', 'nd', 'unknown']:
                        st = st_val[2:] if st_val.startswith('ST') else st_val
                results[sample] = {'ST': st}
            return results
        except Exception: return {}

    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f: html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty: return {}

            sample_col = next((c for c in df.columns if 'filename' in c.lower() or 'sample' in c.lower()), df.columns[0])
            results = {}
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(row[sample_col]) if row[sample_col] else ''
                if not sample: continue
                
                qc_data = {}
                for col in df.columns:
                    if col == sample_col: continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND': qc_data[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try: qc_data[col] = float(cleaned)
                        except: qc_data[col] = str(val)
                results[sample] = qc_data
            return results
        except Exception: return {}

    def parse_kaptive_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing Kaptive: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f: html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty: return {}

            genome_col = next((c for c in df.columns if 'genome' in c.lower()), df.columns[0])
            results = {}
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(row[genome_col]) if row[genome_col] else ''
                if not sample: continue

                k_locus = row.get('K Locus', 'ND')
                o_locus = row.get('O Locus', 'ND')
                
                if 'unknown' in str(k_locus).lower():
                    m = re.search(r'KL(\d+)', str(k_locus), re.I)
                    if m: k_locus = f"KL{m.group(1)}"
                if 'unknown' in str(o_locus).lower():
                    m = re.search(r'OCL(\d+)', str(o_locus), re.I)
                    if m: o_locus = f"OC{m.group(1)}"

                results[sample] = {
                    'K_Locus': k_locus, 'O_Locus': o_locus,
                    'K_Identity': row.get('K Identity', 'ND'), 'K_Coverage': row.get('K Coverage', 'ND'),
                    'O_Identity': row.get('O Identity', 'ND'), 'O_Coverage': row.get('O Coverage', 'ND')
                }
            return results
        except Exception: return {}

    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f: html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2: return {}, {}

            df_genomes = self.parse_html_table(str(tables[0]), 0)
            genes_by_genome = {}
            if not df_genomes.empty:
                sample_col = next((c for c in df_genomes.columns if 'genome' in c.lower()), df_genomes.columns[0])
                gene_col = next((c for c in df_genomes.columns if 'genes' in c.lower() or 'detected' in c.lower()), None)
                if sample_col and gene_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col]) if row[sample_col] else ''
                        if sample: genes_by_genome[sample] = [g.strip() for g in str(row[gene_col]).split(',') if g.strip()]

            df_freq = self.parse_html_table(str(tables[1]), 0)
            gene_freq = {}
            if not df_freq.empty:
                for _, row in df_freq.iterrows():
                    gene = row.get('Gene', '')
                    if not gene: continue
                    count_match = re.search(r'(\d+)', str(row.get('Frequency', '0')))
                    count = int(count_match.group(1)) if count_match else 0
                    pct = (count / total_samples * 100) if total_samples else 0
                    genomes = [self.normalize_sample_id(g.strip()) for g in str(row.get('Genomes', '')).split(',') if g.strip()]
                    gene_freq[gene] = {
                        'count': count, 'percentage': round(pct, 2),
                        'frequency_display': f"{count} ({pct:.1f}%)",
                        'genomes': genomes, 'risk_level': row.get('Risk Level', 'Standard'), 'database': 'amrfinder'
                    }
            return genes_by_genome, gene_freq
        except Exception: return {}, {}

    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f: html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2: return {}, {}

            db_name = next((v for k, v in self.db_name_mapping.items() if k in str(file_path.name).lower()), 'unknown')

            df_genomes = self.parse_html_table(str(tables[0]), 0)
            genes_by_genome = {}
            if not df_genomes.empty:
                sample_col = next((c for c in df_genomes.columns if 'genome' in c.lower() or 'sample' in c.lower()), df_genomes.columns[0])
                gene_col = next((c for c in df_genomes.columns if 'genes' in c.lower() or 'detected' in c.lower()), None)
                if sample_col and gene_col:
                    for _, row in df_genomes.iterrows():
                        sample = self.normalize_sample_id(row[sample_col]) if row[sample_col] else ''
                        if sample: genes_by_genome[sample] = [g.strip() for g in str(row[gene_col]).split(',') if g.strip()]

            df_freq = self.parse_html_table(str(tables[1]), 0)
            gene_freq = {}
            if not df_freq.empty:
                for _, row in df_freq.iterrows():
                    gene_full = str(row.get('Gene', '')).strip()
                    if not gene_full: continue
                    gene = re.sub(r'^\([^)]+\)', '', gene_full).strip() or gene_full
                    
                    count_match = re.search(r'(\d+)', str(row.get('Frequency', '0')))
                    count = int(count_match.group(1)) if count_match else 0
                    pct = (count / total_samples * 100) if total_samples else 0
                    genomes = [self.normalize_sample_id(g.strip()) for g in str(row.get('Genomes', '')).split(',') if g.strip()]
                    
                    gene_freq[gene] = {
                        'count': count, 'percentage': round(pct, 2),
                        'frequency_display': f"{count} ({pct:.1f}%)",
                        'genomes': genomes, 'database': db_name, 'full_name': gene_full
                    }
            return genes_by_genome, gene_freq
        except Exception: return {}, {}

# =============================================================================
# DATA ANALYZER
# =============================================================================
class KleboDataAnalyzer:
    def __init__(self):
        self.critical_resistance_genes = {
            'blaKPC', 'blaNDM', 'blaIMP', 'blaVIM', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaGES', 'blaIMI', 'blaSME', 'blaDHA', 'blaCMY', 'blaACT',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9',
            'tet(X)', 'tet(X1)', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)',
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'npmA',
            'fosA3', 'fosA4', 'fosA5', 'fosA6', 'fosA7', 'fosC2'
        }
        self.high_risk_virulence_genes = {
            'ybtA', 'ybtE', 'ybtP', 'ybtQ', 'ybtS', 'ybtT', 'ybtU', 'ybtX', 'irp1', 'irp2',
            'clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH', 'clbI', 'clbJ',
            'clbK', 'clbL', 'clbM', 'clbN', 'clbO', 'clbP', 'clbQ', 'clbR', 'clbS',
            'iroB', 'iroC', 'iroD', 'iroE', 'iroN', 'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
            'rmpA', 'rmpA2', 'rmpC', 'rmpD', 'peg-344', 'allS', 'kfuA', 'kfuB', 'kfuC',
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'mrkA', 'mrkB', 'mrkC', 'mrkD', 'mrkE', 'mrkF', 'mrkH', 'mrkI',
            'ecpA', 'ecpB', 'ecpC', 'ecpD', 'ecpE', 'ecpR',
            'tssA', 'tssB', 'tssC', 'tssD', 'tssE', 'tssF', 'tssG', 'tssH', 'tssI', 'tssJ',
            'tssK', 'tssL', 'tssM', 'tssN', 'tssO', 'tssP', 'tssQ', 'tssR', 'tssS', 'tssT',
            'hcp', 'vgrG', 'icmF', 'impA', 'impB', 'impC', 'impD', 'impE', 'impF', 'impG', 'impH',
            'siderophore', 'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepG', 'fecA', 'fecB', 'fecC', 'fecD', 'fecE',
            'adhesin', 'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH',
            'toxin', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'cnf1', 'cnf2', 'sat', 'pic', 'vat',
            'cdtA', 'cdtB', 'cdtC', 'cfa', 'cfaA', 'cfaB', 'cfaC', 'cfaD', 'cfaE'
        }
        self.beta_lactamase_genes = {
            'blaTEM', 'blaSHV', 'blaCTX-M', 'blaOXA-1', 'blaOXA-2', 'blaOXA-9', 'blaOXA-10',
            'blaCMY-2', 'blaDHA-1', 'blaACC', 'blaMIR', 'blaACT', 'blaFOX'
        }
        self.icekp_markers = {'ybt', 'clb', 'iro', 'rmp'}
        self.virulence_plasmid_markers = {'iro', 'iuc', 'rmp', 'rmpA2'}

    def categorize_gene(self, gene: str) -> str:
        g = gene.lower()
        if any(crit.lower() in g for crit in self.critical_resistance_genes): return 'Critical Resistance'
        elif any(vir.lower() in g for vir in self.high_risk_virulence_genes): return 'High-Risk Virulence'
        elif any(bla.lower() in g for bla in self.beta_lactamase_genes): return 'Beta-Lactamase'
        else: return 'Other'

    def create_gene_centric_tables(self, gene_freqs: Dict[str, Dict], total_samples: int) -> Dict:
        gene_centric = {'all_genes': [], 'by_category': defaultdict(list), 'database_stats': {}}
        all_genes = {}
        for db, genes in gene_freqs.items():
            for gene, data in genes.items():
                if gene not in all_genes:
                    all_genes[gene] = {'count': 0, 'percentage': 0, 'genomes': set(), 'databases': set(), 'risk_level': 'Standard'}
                all_genes[gene]['count'] += data['count']
                all_genes[gene]['percentage'] += data.get('percentage', 0)
                all_genes[gene]['genomes'].update(data['genomes'])
                all_genes[gene]['databases'].add(data.get('database', db))
                if data.get('risk_level') and data['risk_level'] != 'Standard': all_genes[gene]['risk_level'] = data['risk_level']

        for gene, data in all_genes.items():
            data['genomes'] = list(data['genomes'])
            data['databases'] = list(data['databases'])
            data['category'] = self.categorize_gene(gene)
            data['gene'] = gene
            data['frequency_display'] = f"{data['count']} ({data['count']/total_samples*100:.1f}%)" if total_samples else "0 (0%)"
            gene_centric['by_category'][data['category']].append(data)
            gene_centric['all_genes'].append(data)

        for cat in gene_centric['by_category']: gene_centric['by_category'][cat].sort(key=lambda x: x['count'], reverse=True)
        gene_centric['all_genes'].sort(key=lambda x: x['count'], reverse=True)

        for db, genes in gene_freqs.items():
            gene_centric['database_stats'][db] = {
                'total_genes': len(genes), 'total_occurrences': sum(g['count'] for g in genes.values()),
                'critical_genes': sum(1 for g in genes if self.categorize_gene(g) == 'Critical Resistance')
            }
        return gene_centric

    def create_cross_genome_patterns(self, samples_data: Dict, gene_centric: Dict) -> Dict:
        patterns = {
            'st_distribution': Counter(), 'k_locus_distribution': Counter(), 'o_locus_distribution': Counter(),
            'st_k_combinations': defaultdict(list), 'st_o_combinations': defaultdict(list), 'ko_combinations': defaultdict(list),
            'st_ko_combinations': defaultdict(list), 'high_risk_combinations': [],
            'icekp_marker_presence': defaultdict(list), 'virulence_plasmid_marker_presence': defaultdict(list)
        }
        sample_genes = defaultdict(set)
        for g_data in gene_centric.get('all_genes', []):
            for gen in g_data['genomes']: sample_genes[gen].add(g_data['gene'])

        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            k = data.get('kaptive', {}).get('K_Locus', 'ND')
            o = data.get('kaptive', {}).get('O_Locus', 'ND')
            if st != 'ND': patterns['st_distribution'][st] += 1
            if k != 'ND': patterns['k_locus_distribution'][k] += 1
            if o != 'ND': patterns['o_locus_distribution'][o] += 1
            if st != 'ND' and k != 'ND': patterns['st_k_combinations'][f"{st}-{k}"].append(sample)
            if st != 'ND' and o != 'ND': patterns['st_o_combinations'][f"{st}-{o}"].append(sample)
            if k != 'ND' and o != 'ND': patterns['ko_combinations'][f"{k}:{o}"].append(sample)
            if st != 'ND' and k != 'ND' and o != 'ND': patterns['st_ko_combinations'][f"{st}-{k}:{o}"].append(sample)

            genes = sample_genes.get(sample, set())
            crit_res = [g for g in genes if self.categorize_gene(g) == 'Critical Resistance']
            hv = [g for g in genes if self.categorize_gene(g) == 'High-Risk Virulence']
            if crit_res and hv: patterns['high_risk_combinations'].append({'sample': sample, 'st': st, 'k_locus': k, 'critical_resistance': crit_res, 'high_risk_virulence': hv})
            
            icekp = [g for g in genes if any(m in g.lower() for m in self.icekp_markers)]
            if icekp: patterns['icekp_marker_presence'][sample] = icekp
            
            vp = [g for g in genes if any(m in g.lower() for m in self.virulence_plasmid_markers)]
            if vp: patterns['virulence_plasmid_marker_presence'][sample] = vp
            
        return patterns

# =============================================================================
# HTML GENERATOR (Restored Interactive Tabs)
# =============================================================================
class KleboHTMLGenerator:
    """Generate ultimate HTML report for K. pneumoniae."""

    def __init__(self, analyzer: KleboDataAnalyzer):
        self.analyzer = analyzer

    def generate_main_report(self, integrated_data: Dict, output_dir: Path) -> str:
        print("\n🎨 Generating Kleboscope Ultimate HTML report...")
        html = self._create_ultimate_html(integrated_data)
        output_file = output_dir / "kleboscope_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    def _create_ultimate_html(self, data: Dict) -> str:
        metadata = data.get('metadata', {})
        samples = data.get('samples', {})
        patterns = data.get('patterns', {})
        gene_centric = data.get('gene_centric', {})
        total_samples = len(samples)

        css = """
        <style>
        :root {
            --summary-color: #4CAF50; --samples-color: #2196F3; --mlst-color: #FF9800;
            --qc-color: #9C27B0; --kaptive-color: #E91E63; --amr-color: #F44336;
            --virulence-color: #E91E63; --plasmids-color: #673AB7; --patterns-color: #FF5722;
            --databases-color: #607D8B; --aiguide-color: #3F51B5; --export-color: #3F51B5;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: 'Segoe UI', Tahoma, Geneva, sans-serif; line-height: 1.6; color: #333; background: #f5f5f5; }
        .container { max-width: 100%; margin: 0 auto; padding: 20px; }
        .main-header { background: linear-gradient(135deg, #1e5a3a 0%, #2c7a4d 100%); color: white; padding: 30px; border-radius: 15px; box-shadow: 0 10px 30px rgba(0,0,0,0.2); margin-bottom: 30px; text-align: center; }
        .main-header h1 { font-size: 2.8em; margin-bottom: 10px; }
        .metadata-bar { background: rgba(255,255,255,0.1); padding: 15px; border-radius: 10px; margin: 20px 0; display: flex; justify-content: space-around; flex-wrap: wrap; gap: 15px; }
        .metadata-item { display: flex; align-items: center; gap: 8px; font-size: 0.95em; }
        .dashboard-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin-bottom: 30px; }
        .dashboard-card { background: white; padding: 25px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); text-align: center; cursor: pointer; border-left: 5px solid; transition: transform 0.3s ease; }
        .dashboard-card:hover { transform: translateY(-10px); }
        .card-summary { border-left-color: var(--summary-color); }
        .card-samples { border-left-color: var(--samples-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-qc { border-left-color: var(--qc-color); }
        .card-kaptive { border-left-color: var(--kaptive-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-number { font-size: 3em; font-weight: bold; margin: 15px 0; background: linear-gradient(90deg, #1e5a3a, #2c7a4d); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
        .tab-navigation { display: flex; gap: 5px; margin-bottom: 20px; flex-wrap: wrap; background: white; padding: 15px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); position: sticky; top: 10px; z-index: 100; }
        .tab-button { padding: 12px 25px; background: #f5f5f5; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; color: #666; }
        .tab-button.active { color: white; }
        .tab-button.summary.active { background: var(--summary-color); }
        .tab-button.samples.active { background: var(--samples-color); }
        .tab-button.mlst.active { background: var(--mlst-color); }
        .tab-button.qc.active { background: var(--qc-color); }
        .tab-button.kaptive.active { background: var(--kaptive-color); }
        .tab-button.amr.active { background: var(--amr-color); }
        .tab-button.virulence.active { background: var(--virulence-color); }
        .tab-button.plasmids.active { background: var(--plasmids-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.databases.active { background: var(--databases-color); }
        .tab-button.aiguide.active { background: var(--aiguide-color); }
        .tab-button.export.active { background: var(--export-color); }
        .tab-content { display: none; background: white; padding: 30px; border-radius: 15px; box-shadow: 0 10px 30px rgba(0,0,0,0.1); margin-bottom: 30px; animation: fadeIn 0.5s ease; }
        .tab-content.active { display: block; }
        @keyframes fadeIn { from { opacity: 0; transform: translateY(20px); } to { opacity: 1; transform: translateY(0); } }
        .section-header { margin-bottom: 25px; padding-bottom: 15px; border-bottom: 3px solid; font-size: 1.8em; }
        .summary-header { border-color: var(--summary-color); }
        .samples-header { border-color: var(--samples-color); }
        .mlst-header { border-color: var(--mlst-color); }
        .qc-header { border-color: var(--qc-color); }
        .kaptive-header { border-color: var(--kaptive-color); }
        .amr-header { border-color: var(--amr-color); }
        .virulence-header { border-color: var(--virulence-color); }
        .plasmids-header { border-color: var(--plasmids-color); }
        .patterns-header { border-color: var(--patterns-color); }
        .databases-header { border-color: var(--databases-color); }
        .export-header { border-color: var(--export-color); }
        .data-table { width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.95em; }
        .data-table th { background: #2c3e50; color: white; padding: 15px; text-align: left; position: sticky; top: 0; }
        .data-table td { padding: 12px 15px; border-bottom: 1px solid #e0e0e0; }
        .data-table tr:hover { background: #f8f9fa; }
        .master-scrollable-container { width: 100%; overflow-x: auto; border: 1px solid #e0e0e0; border-radius: 8px; margin: 20px 0; }
        .search-box { width: 100%; padding: 12px; margin-bottom: 20px; border: 2px solid #e0e0e0; border-radius: 8px; font-size: 1em; }
        .badge { display: inline-block; padding: 5px 15px; border-radius: 20px; font-size: 0.85em; font-weight: 600; margin: 2px; white-space: nowrap; }
        .badge-low { background: #4CAF50; color: white; }
        .badge-medium { background: #FF9800; color: black; }
        .badge-high { background: #F44336; color: white; }
        .badge-critical { background: #9C27B0; color: white; }
        .alert-box { padding: 20px; border-radius: 10px; margin: 20px 0; border-left: 5px solid; background: #d1ecf1; color: #0c5460; border-color: #17a2b8; }
        .action-buttons { display: flex; gap: 10px; margin: 20px 0; flex-wrap: wrap; }
        .action-btn { padding: 10px 20px; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; color: white; }
        .btn-primary { background: #2c7a4d; } .btn-success { background: #28a745; } .btn-danger { background: #dc3545; } .btn-warning { background: #ffc107; color: black; }
        .genome-list { max-height: 150px; overflow-y: auto; padding: 5px; background: #f8f9fa; border-radius: 5px; border: 1px solid #e0e0e0; }
        .genome-tag { display: inline-block; background: #e0f2f1; color: #2c7a4d; padding: 3px 10px; border-radius: 12px; font-size: 0.85em; border: 1px solid #b2dfdb; margin: 2px; }
        .category-chip { display: inline-block; padding: 4px 12px; border-radius: 15px; font-size: 0.8em; font-weight: 600; margin: 2px; }
        .chip-critical-resistance { background: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; }
        .chip-high-risk-virulence { background: #fff3cd; color: #856404; border: 1px solid #ffeaa7; }
        .chip-beta-lactamase { background: #d1ecf1; color: #0c5460; border: 1px solid #bee5eb; }
        .chip-other { background: #f5f5f5; color: #212121; border: 1px solid #e0e0e0; }
        .footer { text-align: center; padding: 30px; color: white; margin-top: 40px; border-radius: 15px; background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%); }
        </style>
        """

        js = """
        <script>
        function switchTab(tabName) {
            document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(button => button.classList.remove('active'));
            document.getElementById(tabName + '-tab').classList.add('active');
            document.querySelector('.tab-button.' + tabName).classList.add('active');
            window.location.hash = tabName;
        }
        function searchTable(tableId, searchId) {
            const filter = document.getElementById(searchId).value.toUpperCase();
            const rows = document.getElementById(tableId).getElementsByTagName('tr');
            for (let i = 1; i < rows.length; i++) {
                rows[i].style.display = rows[i].innerText.toUpperCase().indexOf(filter) > -1 ? '' : 'none';
            }
        }
        function exportTableToCSV(tableId, filename) {
            const rows = document.getElementById(tableId).querySelectorAll('tr');
            const csv = [];
            for (let i = 0; i < rows.length; i++) {
                const row = [], cols = rows[i].querySelectorAll('td, th');
                for (let j = 0; j < cols.length; j++) row.push('"' + cols[j].innerText.replace(/"/g, '""') + '"');
                csv.push(row.join(','));
            }
            const blob = new Blob([csv.join('\\n')], {type: 'text/csv'});
            const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = filename;
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
        }
        document.addEventListener('DOMContentLoaded', () => {
            const hash = window.location.hash.substring(1);
            if (hash) document.querySelector(`.tab-button.${hash}`)?.click();
            else document.querySelector('.tab-button').click();
        });
        </script>
        """

        html_parts = [f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><title>Kleboscope Ultimate Report</title>{css}{js}</head>
<body>
    <div class="container">
        <div class="main-header">
            <h1>🦠 Kleboscope Ultimate Report</h1>
            <p>Gene-Centric Analysis for Klebsiella pneumoniae</p>
            <div class="metadata-bar">
                <div class="metadata-item"><span>Generated: {metadata.get('analysis_date', 'Unknown')}</span></div>
                <div class="metadata-item"><span>Samples: {total_samples}</span></div>
            </div>
        </div>
        
        <div class="dashboard-grid">
            <div class="dashboard-card card-summary" onclick="switchTab('summary')"><div class="card-number">{total_samples}</div><div>Total Samples</div></div>
            <div class="dashboard-card card-mlst" onclick="switchTab('mlst')"><div class="card-number">{len(patterns.get('st_distribution', {}))}</div><div>Unique STs</div></div>
            <div class="dashboard-card card-kaptive" onclick="switchTab('kaptive')"><div class="card-number">{len(patterns.get('k_locus_distribution', {}))}</div><div>Capsule Types</div></div>
            <div class="dashboard-card card-amr" onclick="switchTab('amr')"><div class="card-number">{len(gene_centric.get('all_genes', []))}</div><div>Unique Genes</div></div>
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')"><div class="card-number">{len(gene_centric.get('by_category', {}).get('High-Risk Virulence', []))}</div><div>Virulence Genes</div></div>
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')"><div class="card-number">{len(patterns.get('high_risk_combinations', []))}</div><div>High-Risk Combos</div></div>
        </div>
        
        <div class="tab-navigation">
            <button class="tab-button summary active" onclick="switchTab('summary')">Summary</button>
            <button class="tab-button samples" onclick="switchTab('samples')">Sample Overview</button>
            <button class="tab-button mlst" onclick="switchTab('mlst')">MLST</button>
            <button class="tab-button qc" onclick="switchTab('qc')">FASTA QC</button>
            <button class="tab-button kaptive" onclick="switchTab('kaptive')">Kaptive</button>
            <button class="tab-button amr" onclick="switchTab('amr')">AMR Genes</button>
            <button class="tab-button virulence" onclick="switchTab('virulence')">Virulence Genes</button>
            <button class="tab-button plasmids" onclick="switchTab('plasmids')">Plasmids</button>
            <button class="tab-button patterns" onclick="switchTab('patterns')">Patterns</button>
            <button class="tab-button databases" onclick="switchTab('databases')">Databases</button>
            <button class="tab-button aiguide" onclick="switchTab('aiguide')">AI Guide</button>
            <button class="tab-button export" onclick="switchTab('export')">Export</button>
        </div>
        
        <div id="summary-tab" class="tab-content active">{self._generate_summary_section(data)}</div>
        <div id="samples-tab" class="tab-content">{self._generate_sample_overview_section(data)}</div>
        <div id="mlst-tab" class="tab-content">{self._generate_mlst_section(data)}</div>
        <div id="qc-tab" class="tab-content">{self._generate_qc_section(data)}</div>
        <div id="kaptive-tab" class="tab-content">{self._generate_kaptive_section(data)}</div>
        <div id="amr-tab" class="tab-content">{self._generate_amr_section(data)}</div>
        <div id="virulence-tab" class="tab-content">{self._generate_virulence_section(data)}</div>
        <div id="plasmids-tab" class="tab-content">{self._generate_plasmids_section(data)}</div>
        <div id="patterns-tab" class="tab-content">{self._generate_patterns_section(data)}</div>
        <div id="databases-tab" class="tab-content">{self._generate_databases_section(data)}</div>
        <div id="aiguide-tab" class="tab-content">{self._generate_aiguide_section(data)}</div>
        <div id="export-tab" class="tab-content">{self._generate_export_section()}</div>
        
        <div class="footer"><p>Kleboscope Ultimate Reporter v4.0.0 | Generated on {metadata.get('analysis_date', 'Unknown')}</p></div>
    </div>
</body></html>"""]
        return ''.join(html_parts)

    # def _generate_summary_section(self, data: Dict) -> str:
    #     # 1. Extract data dictionaries
    #     samples = data.get('samples', {})
    #     patterns = data.get('patterns', {})
    #     gene_centric = data.get('gene_centric', {})
        
    #     # 2. Calculate metrics
    #     total_samples = len(samples)
    #     unique_sts = len(patterns.get('st_distribution', {}))
    #     unique_capsules = len(patterns.get('k_locus_distribution', {}))
    #     total_unique_genes = len(gene_centric.get('all_genes', []))
    #     high_risk_virulence = len(gene_centric.get('by_category', {}).get('High-Risk Virulence', []))
    #     high_risk_combos = len(patterns.get('high_risk_combinations', []))
        
    #     # Count the number of samples that contain these specific markers
    #     icekp_markers = len(patterns.get('icekp_marker_presence', {}))
    #     virulence_plasmids = len(patterns.get('virulence_plasmid_marker_presence', {}))

        # 3. Generate HTML with the Dashboard Grid
        # return f"""
        # <h2 class="section-header summary-header">Executive Summary</h2>
        # <p style="margin-bottom: 20px; font-size: 1.1em;">
        #     <strong>{total_samples}</strong> genomes analyzed. Gene-centric view shows each gene with all genomes that carry it.
        # </p>
        
        # <div class="dashboard-grid">
        #     <div class="dashboard-card card-summary">
        #         <div class="card-number">{total_samples}</div>
        #         <div>Total Samples</div>
        #     </div>
        #     <div class="dashboard-card card-mlst">
        #         <div class="card-number">{unique_sts}</div>
        #         <div>Unique STs</div>
        #     </div>
        #     <div class="dashboard-card card-kaptive">
        #         <div class="card-number">{unique_capsules}</div>
        #         <div>Capsule Types</div>
        #     </div>
        #     <div class="dashboard-card card-amr">
        #         <div class="card-number">{total_unique_genes}</div>
        #         <div>Unique Genes</div>
        #     </div>
        #     <div class="dashboard-card card-virulence">
        #         <div class="card-number">{high_risk_virulence}</div>
        #         <div>High-Risk Virulence</div>
        #     </div>
        #     <div class="dashboard-card card-patterns">
        #         <div class="card-number">{high_risk_combos}</div>
        #         <div>High-Risk Combos</div>
        #     </div>
        #     <div class="dashboard-card card-qc">
        #         <div class="card-number">{icekp_markers}</div>
        #         <div>ICEKp Markers</div>
        #     </div>
        #     <div class="dashboard-card card-plasmids">
        #         <div class="card-number">{virulence_plasmids}</div>
        #         <div>Virulence Plasmid Markers</div>
        #     </div>
        # </div>
        # """
    def _generate_summary_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        patterns = data.get('patterns', {})
        gene_centric = data.get('gene_centric', {})
        total = len(samples)
        st_count = len(patterns.get('st_distribution', {}))
        k_count = len(patterns.get('k_locus_distribution', {}))
        amr_genes = len(gene_centric.get('all_genes', []))
        hv_count = len(gene_centric.get('by_category', {}).get('High-Risk Virulence', []))
        high_risk = len(patterns.get('high_risk_combinations', []))
        icekp_samples = len(patterns.get('icekp_marker_presence', {}))
        vp_samples = len(patterns.get('virulence_plasmid_marker_presence', {}))
        return f"""
        <h2 class="section-header summary-header"><i class="fas fa-chart-pie"></i> Executive Summary</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>K. pneumoniae Analysis Overview</h3><p><strong>{total}</strong> genomes analyzed. Gene‑centric view shows each gene with all genomes that carry it.</p></div></div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:15px;margin:20px 0;">
            <div class="metric-card" style="background:#2c7a4d;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{total}</div><div>Total Samples</div></div>
            <div class="metric-card" style="background:#FF9800;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{st_count}</div><div>Unique STs</div></div>
            <div class="metric-card" style="background:#E91E63;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{k_count}</div><div>Capsule Types</div></div>
            <div class="metric-card" style="background:#F44336;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{amr_genes}</div><div>Unique Genes</div></div>
            <div class="metric-card" style="background:#E91E63;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{hv_count}</div><div>High‑Risk Virulence</div></div>
            <div class="metric-card" style="background:#FF5722;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{high_risk}</div><div>High‑Risk Combos</div></div>
            <div class="metric-card" style="background:#673AB7;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{icekp_samples}</div><div>ICEKp Markers</div></div>
            <div class="metric-card" style="background:#3F51B5;color:white;padding:15px;border-radius:8px;"><div style="font-size:24px;font-weight:bold;">{vp_samples}</div><div>Virulence Plasmid Markers</div></div>
        </div>
        """
    
    def _generate_sample_overview_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        rows = []
        for sample, d in sorted(samples.items()):
            st = d.get('mlst', {}).get('ST', 'ND')
            k = d.get('kaptive', {}).get('K_Locus', 'ND')
            o = d.get('kaptive', {}).get('O_Locus', 'ND')
            vc = sum(len(genes) for db, genes in d.get('abricate_genes', {}).items() if db in ['vfdb', 'ecoli_vf'])
            rows.append(f"<tr><td><strong>{sample}</strong></td><td>{st}</td><td>{k}</td><td>{o}</td><td>{vc}</td></tr>")
        return f"""<h2 class="section-header samples-header">Sample Overview</h2>
        <input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Search...">
        <div class="master-scrollable-container"><table id="samples-table" class="data-table"><thead><tr><th>Sample</th><th>ST</th><th>K Locus</th><th>O Locus</th><th>Virulence Genes</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_mlst_section(self, data: Dict) -> str:
        st_dist = data.get('patterns', {}).get('st_distribution', {})
        rows = [f"<tr><td><strong>ST{st}</strong></td><td>{cnt}</td></tr>" for st, cnt in sorted(st_dist.items(), key=lambda x: x[1], reverse=True)]
        return f"""<h2 class="section-header mlst-header">MLST Distribution</h2>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>ST</th><th>Count</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_qc_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        metric_list = sorted(list(set(m for d in samples.values() for m in d.get('qc', {}).keys())))
        rows = []
        for sample, d in sorted(samples.items()):
            qc = d.get('qc', {})
            r = [f"<td><strong>{sample}</strong></td>"] + [f"<td>{qc.get(m, 'ND')}</td>" for m in metric_list]
            rows.append(f"<tr>{''.join(r)}</tr>")
        return f"""<h2 class="section-header qc-header">FASTA QC</h2>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Sample</th>{''.join([f'<th>{m}</th>' for m in metric_list])}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_kaptive_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        rows = [f"<tr><td><strong>{s}</strong></td><td>{d.get('kaptive', {}).get('K_Locus', 'ND')}</td><td>{d.get('kaptive', {}).get('O_Locus', 'ND')}</td></tr>" for s, d in sorted(samples.items()) if d.get('kaptive')]
        return f"""<h2 class="section-header kaptive-header">Kaptive</h2>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Sample</th><th>K Locus</th><th>O Locus</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_amr_section(self, data: Dict) -> str:
        gene_freqs = data.get('gene_frequencies', {})
        amr_rows = []
        for db_name, db_genes in gene_freqs.items():
            if db_name in ['vfdb', 'ecoli_vf', 'plasmidfinder']: continue
            for gene, gdata in db_genes.items():
                amr_rows.append((gene, self.analyzer.categorize_gene(gene), db_name.upper(), gdata.get('frequency_display', ''), gdata.get('genomes', [])))
        
        rows = []
        for r in sorted(amr_rows, key=lambda x: len(x[4]), reverse=True):
            genome_tags = ''.join([f"<span class='genome-tag'>{g}</span>" for g in r[4]])
            cat_class = r[1].lower().replace(' ', '-')
            rows.append(f"<tr><td><strong>{r[0]}</strong></td><td><span class='category-chip chip-{cat_class}'>{r[1]}</span></td><td>{r[2]}</td><td>{r[3]}</td><td><div class='genome-list'>{genome_tags}</div></td></tr>")
        
        return f"""<h2 class="section-header amr-header">AMR Genes</h2>
        <input type="text" class="search-box" id="search-amr" onkeyup="searchTable('amr-table','search-amr')" placeholder="🔍 Search AMR...">
        <div class="master-scrollable-container"><table id="amr-table" class="data-table"><thead><tr><th>Gene</th><th>Category</th><th>Database</th><th>Frequency</th><th>Genomes</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_virulence_section(self, data: Dict) -> str:
        gene_freqs = data.get('gene_frequencies', {})
        vir_rows = []
        for db_name in ['vfdb', 'ecoli_vf']:
            if db_name in gene_freqs:
                for gene, gdata in gene_freqs[db_name].items():
                    vir_rows.append((gene, self.analyzer.categorize_gene(gene), db_name.upper(), gdata.get('frequency_display', ''), gdata.get('genomes', [])))
        
        rows = []
        for r in sorted(vir_rows, key=lambda x: len(x[4]), reverse=True):
            genome_tags = ''.join([f"<span class='genome-tag'>{g}</span>" for g in r[4]])
            cat_class = r[1].lower().replace(' ', '-')
            rows.append(f"<tr><td><strong>{r[0]}</strong></td><td><span class='category-chip chip-{cat_class}'>{r[1]}</span></td><td>{r[2]}</td><td>{r[3]}</td><td><div class='genome-list'>{genome_tags}</div></td></tr>")
        
        return f"""<h2 class="section-header virulence-header">Virulence Genes</h2>
        <input type="text" class="search-box" id="search-vir" onkeyup="searchTable('vir-table','search-vir')" placeholder="🔍 Search Virulence...">
        <div class="master-scrollable-container"><table id="vir-table" class="data-table"><thead><tr><th>Gene</th><th>Category</th><th>Database</th><th>Frequency</th><th>Genomes</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_plasmids_section(self, data: Dict) -> str:
        plasmid_genes = [g for g in data.get('gene_centric', {}).get('all_genes', []) if 'plasmidfinder' in g.get('databases', [])]
        
        rows = []
        for g in plasmid_genes:
            genome_tags = ''.join([f"<span class='genome-tag'>{gen}</span>" for gen in g['genomes']])
            rows.append(f"<tr><td><strong>{g['gene']}</strong></td><td>{g['frequency_display']}</td><td><div class='genome-list'>{genome_tags}</div></td></tr>")
            
        return f"""<h2 class="section-header plasmids-header">Plasmid Replicons</h2>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Replicon</th><th>Frequency</th><th>Genomes</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_patterns_section(self, data: Dict) -> str:
        high_risk = data.get('patterns', {}).get('high_risk_combinations', [])
        rows = [f"<tr><td><strong>{c['sample']}</strong></td><td>{c['st']}</td><td>{c['k_locus']}</td><td>{', '.join(c['critical_resistance'])}</td><td>{', '.join(c['high_risk_virulence'])}</td></tr>" for c in high_risk]
        return f"""<h2 class="section-header patterns-header">High-Risk Patterns</h2>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Sample</th><th>ST</th><th>K Locus</th><th>Critical Res</th><th>High-Risk Virulence</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_databases_section(self, data: Dict) -> str:
        db_stats = data.get('gene_centric', {}).get('database_stats', {})
        rows = [f"<tr><td><strong>{db.upper()}</strong></td><td>{s.get('total_genes', 0)}</td><td>{s.get('total_occurrences', 0)}</td></tr>" for db, s in db_stats.items()]
        return f"""<h2 class="section-header databases-header">Database Coverage</h2>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Database</th><th>Unique Genes</th><th>Occurrences</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>"""

    def _generate_aiguide_section(self, data: Dict) -> str:
        return """<h2 class="section-header aiguide-header">AI Assistant Guide</h2>
        <div class="alert-box">Upload the generated <code>kleboscope_ultimate_report.json</code> file to an LLM (like ChatGPT or Claude) to ask detailed questions about your dataset.</div>"""

    def _generate_export_section(self) -> str:
        return """<h2 class="section-header export-header">Export Data</h2>
        <div style="display:flex;gap:20px;flex-wrap:wrap;">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','samples.csv')">Export Samples</button>
            <button class="action-btn btn-danger" onclick="exportTableToCSV('amr-table','amr.csv')">Export AMR</button>
            <button class="action-btn btn-warning" onclick="exportTableToCSV('vir-table','virulence.csv')">Export Virulence</button>
        </div>"""


# =============================================================================
# CSV / JSON EXPORTER
# =============================================================================
class KleboscopeUltimateReporter:
    def __init__(self, input_dir: Path, output_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = KleboHTMLParser()
        self.analyzer = KleboDataAnalyzer()
        self.generator = KleboHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "Kleboscope Ultimate Reporter",
            "version": "4.0.0",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }

    def run(self):
        print("=" * 80)
        print("🧬 Kleboscope Ultimate Reporter")
        print("=" * 80)
        
        html_files = { 'mlst': [], 'qc': [], 'kaptive': [], 'amrfinder': [], 'abricate': defaultdict(list) }
        for f in self.input_dir.glob("**/*.html"):
            name = f.name.lower()
            if 'mlst' in name: html_files['mlst'].append(f)
            elif 'qc' in name: html_files['qc'].append(f)
            elif 'kaptive' in name: html_files['kaptive'].append(f)
            elif 'amrfinder' in name and 'summary' in name: 
                html_files['amrfinder'].append(f)
            else:
                for db in self.parser.abricate_databases:
                    if db in name and 'summary' in name:
                        html_files['abricate'][db].append(f)
                        break
                        
        if not any(html_files.values()):
            print("❌ No HTML files found in the input directory!")
            return False
            
        integrated = {'metadata': self.metadata, 'samples': {}}
        all_samples = set()
        
        if html_files['mlst']: 
            for s, d in self.parser.parse_mlst_report(html_files['mlst'][0]).items():
                all_samples.add(s); integrated['samples'].setdefault(s, {})['mlst'] = d
        if html_files['qc']:
            for s, d in self.parser.parse_qc_report(html_files['qc'][0]).items():
                all_samples.add(s); integrated['samples'].setdefault(s, {})['qc'] = d
        if html_files['kaptive']:
            for s, d in self.parser.parse_kaptive_report(html_files['kaptive'][0]).items():
                all_samples.add(s); integrated['samples'].setdefault(s, {})['kaptive'] = d

        gene_freqs = {}
        if html_files['amrfinder']:
            _, f = self.parser.parse_amrfinder_report(html_files['amrfinder'][0], len(all_samples))
            if f: gene_freqs['amrfinder'] = f
            
        for db, files in html_files['abricate'].items():
            if files:
                _, f = self.parser.parse_abricate_database_report(files[0], len(all_samples))
                if f: gene_freqs[db] = f

        # --- THE CRITICAL FIX IS RIGHT HERE --- #
        integrated['gene_frequencies'] = gene_freqs

        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(gene_freqs, len(all_samples))
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(integrated['samples'], integrated['gene_centric'])

        with open(self.output_dir / "kleboscope_ultimate_report.json", 'w') as f: json.dump(integrated, f, indent=2, default=str)
        self.generate_csv_reports(integrated)
        self.generator.generate_main_report(integrated, self.output_dir)
        print("✅ Ultimate Reports generated successfully!")
        return True

    def generate_csv_reports(self, integrated: Dict):
        rows = []
        for sample, d in integrated['samples'].items():
            rows.append({
                'Sample': sample, 'ST': d.get('mlst', {}).get('ST', 'ND'),
                'K_Locus': d.get('kaptive', {}).get('K_Locus', 'ND'), 'O_Locus': d.get('kaptive', {}).get('O_Locus', 'ND'),
                'Virulence_Count': sum(len(genes) for db, genes in d.get('abricate_genes', {}).items() if db in ['vfdb', 'ecoli_vf'])
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)

        amr_rows = []
        for g in integrated['gene_centric'].get('all_genes', []):
            if 'plasmidfinder' in g.get('databases', []) or any(db in ['vfdb', 'ecoli_vf'] for db in g.get('databases', [])): continue
            amr_rows.append({'Gene': g['gene'], 'Category': g['category'], 'Database': ', '.join(g['databases']), 'Count': g['count'], 'Genomes': ';'.join(g['genomes'])})
        if amr_rows: pd.DataFrame(amr_rows).to_csv(self.output_dir / "amr_genes.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description='Kleboscope Ultimate Reporter - Orchestrator Compliant')
    parser.add_argument('-i', '--input', required=True, dest='input_dir', help='Input directory containing HTML reports')
    parser.add_argument('-o', '--output', required=True, dest='output_dir', help='Output directory for final reports')
    args = parser.parse_args()
    
    reporter = KleboscopeUltimateReporter(Path(args.input_dir), Path(args.output_dir))
    reporter.run()

if __name__ == "__main__":
    main()