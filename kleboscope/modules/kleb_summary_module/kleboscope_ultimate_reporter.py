#!/usr/bin/env python3
"""
Kleboscope Ultimate Reporter – Comprehensive Gene‑Centric Analysis for K. pneumoniae
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0 – Per‑Database AMR & Virulence Tables, Enhanced Filter Buttons, Full Database Coverage
Date: 2026-03-25
"""

import os
import sys
import json
import re
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
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
        # All expected ABRicate databases
        self.abricate_databases = [
            'card', 'resfinder', 'argannot', 'vfdb', 'plasmidfinder',
            'megares', 'ncbi', 'ecoh', 'ecoli_vf', 'bacmet2'
        ]
        # Database name mapping (from filenames)
        self.db_name_mapping = {
            'klebo_card': 'card',
            'klebo_resfinder': 'resfinder',
            'klebo_argannot': 'argannot',
            'klebo_vfdb': 'vfdb',
            'klebo_plasmidfinder': 'plasmidfinder',
            'klebo_megares': 'megares',
            'klebo_ncbi': 'ncbi',
            'klebo_ecoh': 'ecoh',
            'klebo_ecoli_vf': 'ecoli_vf',
            'klebo_bacmet2': 'bacmet2',
            'acineto_card': 'card',  # fallback
            'acineto_resfinder': 'resfinder',
            'acineto_argannot': 'argannot',
            'acineto_vfdb': 'vfdb',
            'acineto_plasmidfinder': 'plasmidfinder',
            'acineto_megares': 'megares',
            'acineto_ncbi': 'ncbi',
            'acineto_ecoh': 'ecoh',
            'acineto_ecoli_vf': 'ecoli_vf',
            'acineto_bacmet2': 'bacmet2'
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        """Remove .fna, .fasta, etc., for consistent matching."""
        sample = str(sample_id).strip()
        for ext in ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Parse an HTML table into a pandas DataFrame."""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            table = tables[table_index]
            rows = table.find_all('tr')
            if not rows:
                return pd.DataFrame()
            # headers
            headers = []
            ths = rows[0].find_all(['th', 'td'])
            for th in ths:
                headers.append(th.get_text().strip())
            # data rows
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                row_data = [col.get_text().strip() for col in cols]
                # pad if necessary
                while len(row_data) < len(headers):
                    row_data.append('')
                if len(row_data) > len(headers):
                    row_data = row_data[:len(headers)]
                data.append(row_data)
            if not data:
                return pd.DataFrame()
            df = pd.DataFrame(data, columns=headers)
            df.columns = [c.replace('\n', ' ').strip() for c in df.columns]
            return df
        except Exception as e:
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()

    # --------------------------------------------------------------------------
    # MLST
    # --------------------------------------------------------------------------
    def parse_mlst_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse MLST summary HTML."""
        print(f"  🧬 Parsing MLST: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}

            # Find sample column
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower():
                    sample_col = col
                    break
            if not sample_col and len(df.columns) > 0:
                sample_col = df.columns[0]

            # Find ST column
            st_col = None
            for col in df.columns:
                if col.lower() == 'st' or 'st' in col.lower() and 'sample' not in col.lower():
                    st_col = col
                    break

            results = {}
            for _, row in df.iterrows():
                sample_raw = row[sample_col] if sample_col else ''
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                st = 'ND'
                if st_col and pd.notna(row.get(st_col)):
                    st_val = str(row[st_col]).strip()
                    if st_val and st_val.lower() not in ['', 'nan', 'none', 'nd', 'unknown']:
                        if st_val.startswith('ST'):
                            st = st_val[2:]
                        else:
                            st = st_val
                results[sample] = {'ST': st}
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing MLST: {e}")
            return {}

    # --------------------------------------------------------------------------
    # FASTA QC – parse ALL columns
    # --------------------------------------------------------------------------
    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse FASTA QC summary HTML – extract all available metrics."""
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}

            # Find sample column
            sample_col = None
            for col in df.columns:
                if 'filename' in col.lower() or 'sample' in col.lower() or col == df.columns[0]:
                    sample_col = col
                    break
            if not sample_col:
                return {}

            results = {}
            for _, row in df.iterrows():
                sample_raw = row[sample_col]
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                # Extract all available columns
                qc_data = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc_data[col] = 'ND'
                    else:
                        # try to convert to numeric if it looks like a number
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc_data[col] = float(cleaned)
                        except:
                            qc_data[col] = str(val)
                results[sample] = qc_data
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing QC: {e}")
            return {}

    # --------------------------------------------------------------------------
    # Kaptive
    # --------------------------------------------------------------------------
    def parse_kaptive_report(self, file_path: Path) -> Dict[str, Dict]:
        """Parse Kaptive K/O summary HTML."""
        print(f"  🧬 Parsing Kaptive: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}

            genome_col = None
            for col in df.columns:
                if 'genome' in col.lower():
                    genome_col = col
                    break
            if not genome_col and len(df.columns) > 0:
                genome_col = df.columns[0]

            results = {}
            for _, row in df.iterrows():
                sample_raw = row[genome_col]
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)

                k_locus = row.get('K Locus', 'ND') if 'K Locus' in df.columns else 'ND'
                o_locus = row.get('O Locus', 'ND') if 'O Locus' in df.columns else 'ND'
                # Clean up K locus: if "unknown (KLxx)" -> "KLxx"
                if 'unknown' in str(k_locus).lower():
                    k_match = re.search(r'KL(\d+)', str(k_locus), re.I)
                    if k_match:
                        k_locus = f"KL{k_match.group(1)}"
                if 'unknown' in str(o_locus).lower():
                    o_match = re.search(r'OCL(\d+)', str(o_locus), re.I)
                    if o_match:
                        o_locus = f"OC{o_match.group(1)}"

                results[sample] = {
                    'K_Locus': k_locus,
                    'O_Locus': o_locus,
                    'K_Identity': row.get('K Identity', 'ND') if 'K Identity' in df.columns else 'ND',
                    'K_Coverage': row.get('K Coverage', 'ND') if 'K Coverage' in df.columns else 'ND',
                    'O_Identity': row.get('O Identity', 'ND') if 'O Identity' in df.columns else 'ND',
                    'O_Coverage': row.get('O Coverage', 'ND') if 'O Coverage' in df.columns else 'ND',
                }
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing Kaptive: {e}")
            return {}

    # --------------------------------------------------------------------------
    # AMRfinder
    # --------------------------------------------------------------------------
    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse AMRfinderPlus batch summary HTML."""
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}

            # Table 0: Genes by Genome
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            genes_by_genome = {}
            if not df_genomes.empty:
                sample_col = None
                for col in df_genomes.columns:
                    if 'genome' in col.lower():
                        sample_col = col
                        break
                if not sample_col and len(df_genomes.columns) > 0:
                    sample_col = df_genomes.columns[0]
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample_raw = row[sample_col]
                        if not sample_raw:
                            continue
                        sample = self.normalize_sample_id(sample_raw)
                        # Find genes column (often "Genes Detected")
                        genes = []
                        for col in df_genomes.columns:
                            if 'genes' in col.lower() or 'detected' in col.lower():
                                gene_str = str(row[col])
                                genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                                break
                        genes_by_genome[sample] = genes

            # Table 1: Gene Frequency
            df_freq = self.parse_html_table(str(tables[1]), 0)
            gene_freq = {}
            if not df_freq.empty:
                for _, row in df_freq.iterrows():
                    gene = row.get('Gene', '')
                    if not gene:
                        continue
                    freq_str = row.get('Frequency', '0')
                    count_match = re.search(r'(\d+)', freq_str)
                    count = int(count_match.group(1)) if count_match else 0
                    percentage = (count / total_samples * 100) if total_samples else 0
                    genomes_str = row.get('Genomes', '')
                    genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()] if genomes_str else []
                    risk = row.get('Risk Level', 'Standard')
                    gene_freq[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",
                        'genomes': genomes,
                        'risk_level': risk,
                        'database': 'amrfinder'
                    }
            print(f"    ✓ Parsed {len(genes_by_genome)} samples, {len(gene_freq)} genes")
            return genes_by_genome, gene_freq
        except Exception as e:
            print(f"    ❌ Error parsing AMRfinder: {e}")
            return {}, {}

    # --------------------------------------------------------------------------
    # ABRicate database reports (all)
    # --------------------------------------------------------------------------
    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        """Parse any ABRicate database HTML report."""
        print(f"  🧬 Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}

            # Determine database name from filename
            db_name = 'unknown'
            filename = str(file_path.name).lower()
            for key, value in self.db_name_mapping.items():
                if key in filename:
                    db_name = value
                    break

            # Table 0: Genes by Genome
            df_genomes = self.parse_html_table(str(tables[0]), 0)
            genes_by_genome = {}
            if not df_genomes.empty:
                sample_col = None
                for col in df_genomes.columns:
                    if 'genome' in col.lower() or 'sample' in col.lower() or 'id' in col.lower():
                        sample_col = col
                        break
                if not sample_col and len(df_genomes.columns) > 0:
                    sample_col = df_genomes.columns[0]
                if sample_col:
                    for _, row in df_genomes.iterrows():
                        sample_raw = row[sample_col]
                        if not sample_raw:
                            continue
                        sample = self.normalize_sample_id(sample_raw)
                        genes = []
                        # Find genes column
                        for col in df_genomes.columns:
                            if 'genes' in col.lower() or 'detected' in col.lower():
                                gene_str = str(row[col])
                                genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                                break
                        genes_by_genome[sample] = genes

            # Table 1: Gene Frequency
            df_freq = self.parse_html_table(str(tables[1]), 0)
            gene_freq = {}
            if not df_freq.empty:
                for _, row in df_freq.iterrows():
                    gene_full = str(row.get('Gene', '')).strip()
                    if not gene_full:
                        continue
                    # Clean gene name (remove category prefix like "(AGly)")
                    gene = re.sub(r'^\([^)]+\)', '', gene_full).strip()
                    if not gene:
                        gene = gene_full
                    freq_str = row.get('Frequency', '0')
                    count_match = re.search(r'(\d+)', freq_str)
                    count = int(count_match.group(1)) if count_match else 0
                    percentage = (count / total_samples * 100) if total_samples else 0
                    genomes_str = row.get('Genomes', '')
                    genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()] if genomes_str else []
                    gene_freq[gene] = {
                        'count': count,
                        'percentage': round(percentage, 2),
                        'frequency_display': f"{count} ({percentage:.1f}%)",
                        'genomes': genomes,
                        'database': db_name,
                        'full_name': gene_full
                    }
            print(f"    ✓ {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_freq)} genes")
            return genes_by_genome, gene_freq
        except Exception as e:
            print(f"    ❌ Error parsing ABRicate report: {e}")
            import traceback
            traceback.print_exc()
            return {}, {}


# =============================================================================
# DATA ANALYZER
# =============================================================================
class KleboDataAnalyzer:
    """Analyze K. pneumoniae data, categorize genes, find patterns."""

    def __init__(self):
        # Critical resistance genes (carbapenemases, mcr, tet(X), 16S rRNA methyltransferases, etc.)
        self.critical_resistance_genes = {
            'blaKPC', 'blaNDM', 'blaIMP', 'blaVIM', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaGES', 'blaIMI', 'blaSME', 'blaDHA', 'blaCMY', 'blaACT',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9',
            'tet(X)', 'tet(X1)', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)',
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'npmA',
            'fosA3', 'fosA4', 'fosA5', 'fosA6', 'fosA7', 'fosC2'
        }

        # High-risk virulence genes – focus on ICEKp and virulence plasmid associated loci
        self.high_risk_virulence_genes = {
            # Yersiniabactin (ybt) – ICEKp associated
            'ybtA', 'ybtE', 'ybtP', 'ybtQ', 'ybtS', 'ybtT', 'ybtU', 'ybtX', 'irp1', 'irp2',
            # Colibactin (clb) – ICEKp associated
            'clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH', 'clbI', 'clbJ',
            'clbK', 'clbL', 'clbM', 'clbN', 'clbO', 'clbP', 'clbQ', 'clbR', 'clbS',
            # Salmochelin (iro) – ICEKp and virulence plasmid
            'iroB', 'iroC', 'iroD', 'iroE', 'iroN',
            # Aerobactin (iuc) – virulence plasmid associated
            'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
            # Hypermucoidy (rmp, rmpA2) – virulence plasmid associated
            'rmpA', 'rmpA2', 'rmpC', 'rmpD',
            # Additional hypervirulence markers
            'peg-344', 'allS', 'kfuA', 'kfuB', 'kfuC',
            # Fimbrial adhesins
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'mrkA', 'mrkB', 'mrkC', 'mrkD', 'mrkE', 'mrkF', 'mrkH', 'mrkI',
            'ecpA', 'ecpB', 'ecpC', 'ecpD', 'ecpE', 'ecpR',
            # Type VI secretion system
            'tssA', 'tssB', 'tssC', 'tssD', 'tssE', 'tssF', 'tssG', 'tssH', 'tssI', 'tssJ',
            'tssK', 'tssL', 'tssM', 'tssN', 'tssO', 'tssP', 'tssQ', 'tssR', 'tssS', 'tssT',
            # Additional virulence categories
            'hcp', 'vgrG', 'icmF', 'impA', 'impB', 'impC', 'impD', 'impE', 'impF', 'impG', 'impH',
            'siderophore', 'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepG', 'fecA', 'fecB', 'fecC', 'fecD', 'fecE',
            'adhesin', 'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH',
            'toxin', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'cnf1', 'cnf2', 'sat', 'pic', 'vat',
            'cdtA', 'cdtB', 'cdtC', 'cfa', 'cfaA', 'cfaB', 'cfaC', 'cfaD', 'cfaE'
        }

        # Beta-lactamases (non-critical)
        self.beta_lactamase_genes = {
            'blaTEM', 'blaSHV', 'blaCTX-M', 'blaOXA-1', 'blaOXA-2', 'blaOXA-9', 'blaOXA-10',
            'blaCMY-2', 'blaDHA-1', 'blaACC', 'blaMIR', 'blaACT', 'blaFOX'
        }

        # ICEKp marker keywords (for pattern discovery)
        self.icekp_markers = {'ybt', 'clb', 'iro', 'rmp'}
        # Virulence plasmid marker keywords
        self.virulence_plasmid_markers = {'iro', 'iuc', 'rmp', 'rmpA2'}

    def categorize_gene(self, gene: str) -> str:
        g = gene.lower()
        if any(crit.lower() in g for crit in self.critical_resistance_genes):
            return 'Critical Resistance'
        elif any(vir.lower() in g for vir in self.high_risk_virulence_genes):
            return 'High-Risk Virulence'
        elif any(bla.lower() in g for bla in self.beta_lactamase_genes):
            return 'Beta-Lactamase'
        else:
            return 'Other'

    def create_gene_centric_tables(self, gene_freqs: Dict[str, Dict], total_samples: int) -> Dict:
        """Build gene-centric data structure (merged by gene)."""
        gene_centric = {
            'all_genes': [],
            'by_category': defaultdict(list),
            'database_stats': {}
        }
        all_genes = {}
        for db, genes in gene_freqs.items():
            for gene, data in genes.items():
                if gene not in all_genes:
                    all_genes[gene] = {
                        'count': 0,
                        'percentage': 0,
                        'frequency_display': '',
                        'genomes': set(),
                        'databases': set(),
                        'risk_level': 'Standard'
                    }
                all_genes[gene]['count'] += data['count']
                all_genes[gene]['percentage'] += data.get('percentage', 0)
                all_genes[gene]['frequency_display'] = data.get('frequency_display', f"{data['count']} ({data['count']/total_samples*100:.1f}%)")
                all_genes[gene]['genomes'].update(data['genomes'])
                all_genes[gene]['databases'].add(data.get('database', db))
                if data.get('risk_level') and data['risk_level'] != 'Standard':
                    all_genes[gene]['risk_level'] = data['risk_level']

        for gene, data in all_genes.items():
            data['genomes'] = list(data['genomes'])
            data['databases'] = list(data['databases'])
            data['category'] = self.categorize_gene(gene)
            data['gene'] = gene
            gene_centric['by_category'][data['category']].append(data)
            gene_centric['all_genes'].append(data)

        for cat in gene_centric['by_category']:
            gene_centric['by_category'][cat].sort(key=lambda x: x['count'], reverse=True)
        gene_centric['all_genes'].sort(key=lambda x: x['count'], reverse=True)

        for db, genes in gene_freqs.items():
            gene_centric['database_stats'][db] = {
                'total_genes': len(genes),
                'total_occurrences': sum(g['count'] for g in genes.values()),
                'critical_genes': sum(1 for g in genes if self.categorize_gene(g) == 'Critical Resistance')
            }
        return gene_centric

    def create_cross_genome_patterns(self, samples_data: Dict, gene_centric: Dict) -> Dict:
        patterns = {
            'st_distribution': Counter(),
            'k_locus_distribution': Counter(),
            'o_locus_distribution': Counter(),
            'st_k_combinations': defaultdict(list),
            'st_o_combinations': defaultdict(list),
            'ko_combinations': defaultdict(list),
            'st_ko_combinations': defaultdict(list),
            'high_risk_combinations': [],
            'icekp_marker_presence': defaultdict(list),
            'virulence_plasmid_marker_presence': defaultdict(list)
        }

        # Build gene set per sample
        sample_genes = defaultdict(set)
        for gene_data in gene_centric.get('all_genes', []):
            for genome in gene_data['genomes']:
                sample_genes[genome].add(gene_data['gene'])

        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            k = data.get('kaptive', {}).get('K_Locus', 'ND')
            o = data.get('kaptive', {}).get('O_Locus', 'ND')
            if st != 'ND':
                patterns['st_distribution'][st] += 1
            if k != 'ND':
                patterns['k_locus_distribution'][k] += 1
            if o != 'ND':
                patterns['o_locus_distribution'][o] += 1
            if st != 'ND' and k != 'ND':
                patterns['st_k_combinations'][f"{st}-{k}"].append(sample)
            if st != 'ND' and o != 'ND':
                patterns['st_o_combinations'][f"{st}-{o}"].append(sample)
            if k != 'ND' and o != 'ND':
                patterns['ko_combinations'][f"{k}:{o}"].append(sample)
            if st != 'ND' and k != 'ND' and o != 'ND':
                patterns['st_ko_combinations'][f"{st}-{k}:{o}"].append(sample)

            genes = sample_genes.get(sample, set())
            crit_res = [g for g in genes if self.categorize_gene(g) == 'Critical Resistance']
            hv = [g for g in genes if self.categorize_gene(g) == 'High-Risk Virulence']
            if crit_res and hv:
                patterns['high_risk_combinations'].append({
                    'sample': sample,
                    'st': st,
                    'k_locus': k,
                    'critical_resistance': crit_res,
                    'high_risk_virulence': hv
                })

            # ICEKp associated markers
            icekp_found = []
            for gene in genes:
                if any(marker in gene.lower() for marker in self.icekp_markers):
                    icekp_found.append(gene)
            if icekp_found:
                patterns['icekp_marker_presence'][sample] = icekp_found

            # Virulence plasmid associated markers
            vp_found = []
            for gene in genes:
                if any(marker in gene.lower() for marker in self.virulence_plasmid_markers):
                    vp_found.append(gene)
            if vp_found:
                patterns['virulence_plasmid_marker_presence'][sample] = vp_found

        return patterns


# =============================================================================
# HTML GENERATOR (adapted from Acinetoscope & StaphScope ultimate reporter)
# =============================================================================
class KleboHTMLGenerator:
    """Generate ultimate HTML report for K. pneumoniae."""

    def __init__(self, analyzer: KleboDataAnalyzer):
        self.analyzer = analyzer
        self.tab_colors = {
            'summary': "#4CAF50",
            'samples': '#2196F3',
            'mlst': '#FF9800',
            'qc': '#9C27B0',
            'kaptive': '#E91E63',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'plasmids': '#673AB7',
            'patterns': '#FF5722',
            'databases': '#607D8B',
            'aiguide': '#3F51B5',
            'export': '#3F51B5'
        }

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

        # CSS (identical to Acinetoscope but with green/teal theme)
        css = """
        <style>
        :root {
            --summary-color: #4CAF50;
            --samples-color: #2196F3;
            --mlst-color: #FF9800;
            --qc-color: #9C27B0;
            --kaptive-color: #E91E63;
            --amr-color: #F44336;
            --virulence-color: #E91E63;
            --plasmids-color: #673AB7;
            --patterns-color: #FF5722;
            --databases-color: #607D8B;
            --aiguide-color: #3F51B5;
            --export-color: #3F51B5;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
            overflow-x: auto;
        }
        
        .container {
            width: 100%;
            max-width: 100%;
            margin: 0 auto;
            padding: 20px;
            overflow-x: hidden;
        }
        
        .main-header {
            background: linear-gradient(135deg, #1e5a3a 0%, #2c7a4d 100%);
            color: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-bottom: 30px;
            text-align: center;
        }
        
        .main-header h1 {
            font-size: 2.8em;
            margin-bottom: 10px;
            color: white;
        }
        
        .metadata-bar {
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            justify-content: space-around;
            flex-wrap: wrap;
            gap: 15px;
            backdrop-filter: blur(10px);
        }
        
        .metadata-item {
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 0.95em;
        }
        
        .dashboard-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        
        .dashboard-card {
            background: white;
            padding: 25px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            text-align: center;
            transition: all 0.3s ease;
            cursor: pointer;
            border-left: 5px solid;
            position: relative;
            overflow: hidden;
        }
        
        .dashboard-card:hover {
            transform: translateY(-10px);
            box-shadow: 0 15px 30px rgba(0,0,0,0.2);
        }
        
        .card-summary { border-left-color: var(--summary-color); }
        .card-samples { border-left-color: var(--samples-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-qc { border-left-color: var(--qc-color); }
        .card-kaptive { border-left-color: var(--kaptive-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-plasmids { border-left-color: var(--plasmids-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-databases { border-left-color: var(--databases-color); }
        .card-aiguide { border-left-color: var(--aiguide-color); }
        .card-export { border-left-color: var(--export-color); }
        
        .card-number {
            font-size: 3em;
            font-weight: bold;
            margin: 15px 0;
            background: linear-gradient(90deg, #1e5a3a, #2c7a4d);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .tab-navigation {
            display: flex;
            gap: 5px;
            margin-bottom: 20px;
            flex-wrap: wrap;
            background: white;
            padding: 15px;
            border-radius: 12px;
            box-shadow: 0 5px 20px rgba(0,0,0,0.1);
            position: sticky;
            top: 10px;
            z-index: 100;
        }
        
        .tab-button {
            padding: 12px 25px;
            background: #f5f5f5;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            color: #666;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            gap: 8px;
            position: relative;
            overflow: hidden;
        }
        
        .tab-button::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 50%;
            right: 50%;
            height: 3px;
            background: currentColor;
            transition: all 0.3s ease;
        }
        
        .tab-button:hover::after {
            left: 10%;
            right: 10%;
        }
        
        .tab-button.active {
            color: white;
        }
        
        .tab-button.active::after {
            left: 10%;
            right: 10%;
        }
        
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
        
        .tab-content {
            display: none;
            background: white;
            padding: 30px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            margin-bottom: 30px;
            animation: fadeIn 0.5s ease;
            width: 100%;
            overflow: hidden;
        }
        
        .tab-content.active {
            display: block;
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .section-header {
            color: #2c3e50;
            margin-bottom: 25px;
            padding-bottom: 15px;
            border-bottom: 3px solid;
            font-size: 1.8em;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
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
        .aiguide-header { border-color: var(--aiguide-color); }
        .export-header { border-color: var(--export-color); }
        
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }
        
        .data-table th {
            background: #2c3e50;
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
            position: sticky;
            top: 0;
            white-space: nowrap;
            z-index: 10;
        }
        
        .data-table td {
            padding: 12px 15px;
            border-bottom: 1px solid #e0e0e0;
            vertical-align: top;
            word-wrap: break-word;
            word-break: break-word;
            white-space: normal;
        }
        
        .data-table tr:hover {
            background: #f8f9fa;
        }
        
        .master-scrollable-container {
            width: 100%;
            max-width: 100%;
            overflow-x: auto;
            overflow-y: visible;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            margin: 20px 0;
            position: relative;
        }
        
        .master-scrollable-container table {
            min-width: 100%;
            width: auto;
        }
        
        .master-scrollable-container table.data-table {
            margin: 0;
        }
        
        .search-box {
            width: 100%;
            max-width: 100%;
            padding: 12px;
            margin-bottom: 20px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 1em;
            transition: all 0.3s ease;
        }
        
        .search-box:focus {
            outline: none;
            border-color: #2c7a4d;
            box-shadow: 0 0 0 3px rgba(44,122,77,0.1);
        }
        
        .badge {
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            margin: 2px;
            white-space: nowrap;
        }
        
        .badge-low { background: #4CAF50; color: white; }
        .badge-medium { background: #FF9800; color: black; }
        .badge-high { background: #F44336; color: white; }
        .badge-critical { background: #9C27B0; color: white; }
        
        .alert-box {
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            align-items: center;
            gap: 20px;
            border-left: 5px solid;
        }
        
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        
        .action-buttons {
            display: flex;
            gap: 10px;
            margin: 20px 0;
            flex-wrap: wrap;
        }
        
        .action-btn {
            padding: 10px 20px;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            display: flex;
            align-items: center;
            gap: 8px;
            transition: all 0.3s ease;
            white-space: nowrap;
        }
        
        .action-btn:hover {
            transform: translateY(-2px);
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }
        
        .btn-primary { background: #2c7a4d; color: white; }
        .btn-success { background: #28a745; color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        
        .database-section {
            margin: 30px 0;
            padding: 25px;
            border-radius: 12px;
            background: #f8f9fa;
            box-shadow: 0 3px 15px rgba(0,0,0,0.08);
            overflow: hidden;
        }
        
        .genome-list {
            display: block;
            max-height: 200px;
            overflow-y: auto;
            padding: 5px;
            background: #f8f9fa;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
        }
        
        .genome-tag {
            display: inline-block;
            background: #e0f2f1;
            color: #2c7a4d;
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.85em;
            border: 1px solid #b2dfdb;
            margin: 2px;
            word-break: break-all;
            white-space: normal;
        }
        
        .footer {
            text-align: center;
            padding: 30px;
            color: white;
            margin-top: 40px;
            border-radius: 15px;
            background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
        }
        
        .category-chip {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 15px;
            font-size: 0.8em;
            font-weight: 600;
            margin: 2px;
            white-space: nowrap;
        }
        
        .chip-critical-resistance { background: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; }
        .chip-high-risk-virulence { background: #fff3cd; color: #856404; border: 1px solid #ffeaa7; }
        .chip-beta-lactamase { background: #d1ecf1; color: #0c5460; border: 1px solid #bee5eb; }
        .chip-other { background: #f5f5f5; color: #212121; border: 1px solid #e0e0e0; }
        
        .frequency-display {
            font-weight: 600;
            color: #2c3e50;
        }
        
        @media (max-width: 768px) {
            .container {
                padding: 10px;
            }
            .main-header {
                padding: 20px;
            }
            .main-header h1 {
                font-size: 2em;
            }
            .tab-button {
                padding: 10px 15px;
                font-size: 0.9em;
            }
            .dashboard-grid {
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            }
            .data-table {
                font-size: 0.85em;
            }
            .metadata-bar {
                flex-direction: column;
                gap: 10px;
            }
        }
        
        .col-gene { min-width: 200px; }
        .col-category { min-width: 180px; }
        .col-database { min-width: 120px; }
        .col-frequency { min-width: 120px; }
        .col-sample { min-width: 250px; }
        .col-genomes { min-width: 400px; }
        .col-risk { min-width: 100px; }
        .col-st { min-width: 100px; }
        .col-k-locus { min-width: 100px; }
        .col-o-locus { min-width: 100px; }
        </style>
        """

        # JavaScript (identical to Acinetoscope)
        js = """
        <script>
        // Tab switching
        function switchTab(tabName) {
            document.querySelectorAll('.tab-content').forEach(tab => {
                tab.classList.remove('active');
            });
            document.querySelectorAll('.tab-button').forEach(button => {
                button.classList.remove('active');
            });
            document.getElementById(tabName + '-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }
        
        // Search functionality
        function searchTable(tableId, searchId) {
            const input = document.getElementById(searchId);
            const filter = input.value.toUpperCase();
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            
            for (let i = 1; i < rows.length; i++) {
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {
                    const cell = cells[j];
                    if (cell) {
                        const txtValue = cell.textContent || cell.innerText;
                        if (txtValue.toUpperCase().indexOf(filter) > -1) {
                            found = true;
                            break;
                        }
                    }
                }
                rows[i].style.display = found ? '' : 'none';
            }
        }
        
        // Print current section
        function printSection(sectionId) {
            const content = document.getElementById(sectionId);
            const printWindow = window.open('', '_blank');
            printWindow.document.write('<html><head><title>Print Section</title>');
            printWindow.document.write('<style>' + document.querySelector('style').textContent + '</style>');
            printWindow.document.write('</head><body>');
            printWindow.document.write(content.innerHTML);
            printWindow.document.write('</body></html>');
            printWindow.document.close();
            printWindow.print();
        }
        
        // Export table to CSV
        function exportTableToCSV(tableId, filename) {
            const table = document.getElementById(tableId);
            const rows = table.querySelectorAll('tr');
            const csv = [];
            for (let i = 0; i < rows.length; i++) {
                const row = [], cols = rows[i].querySelectorAll('td, th');
                for (let j = 0; j < cols.length; j++) {
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }
                csv.push(row.join(','));
            }
            const csvFile = new Blob([csv.join('\\n')], {type: 'text/csv'});
            const downloadLink = document.createElement('a');
            downloadLink.download = filename;
            downloadLink.href = window.URL.createObjectURL(csvFile);
            downloadLink.style.display = 'none';
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }
        
        // Initialize from URL hash
        document.addEventListener('DOMContentLoaded', function() {
            const hash = window.location.hash.substring(1);
            if (hash) {
                const tabButton = document.querySelector(`.tab-button.${hash}`);
                if (tabButton) {
                    tabButton.click();
                }
            } else {
                document.querySelector('.tab-button').click();
            }
        });
        </script>
        """

        # Build HTML content
        html_parts = []
        html_parts.append(f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Kleboscope Ultimate Report – K. pneumoniae</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
    <div class="container">
        <div class="main-header">
            <h1><i class="fas fa-bacterium"></i> Kleboscope Ultimate Analysis Report</h1>
            <p>Gene‑Centric, Cross‑Genome Analysis for Klebsiella pneumoniae</p>
            <div class="metadata-bar">
                <div class="metadata-item"><i class="fas fa-calendar"></i><span>Generated: {metadata.get('analysis_date', 'Unknown')}</span></div>
                <div class="metadata-item"><i class="fas fa-database"></i><span>Samples: {total_samples}</span></div>
                <div class="metadata-item"><i class="fas fa-user-md"></i><span>Tool: Kleboscope Ultimate v1.0.0</span></div>
                <div class="metadata-item"><i class="fas fa-university"></i><span>University of Ghana Medical School</span></div>
            </div>
        </div>
        
        <div class="dashboard-grid">
            <div class="dashboard-card card-summary" onclick="switchTab('summary')">
                <div class="card-number">{total_samples}</div><div class="card-label">Total Samples</div>
                <i class="fas fa-vial fa-2x" style="color: var(--summary-color); margin-top: 10px;"></i>
            </div>
            <div class="dashboard-card card-mlst" onclick="switchTab('mlst')">
                <div class="card-number">{len(patterns.get('st_distribution', {}))}</div><div class="card-label">Unique STs</div>
                <i class="fas fa-code-branch fa-2x" style="color: var(--mlst-color); margin-top: 10px;"></i>
            </div>
            <div class="dashboard-card card-kaptive" onclick="switchTab('kaptive')">
                <div class="card-number">{len(patterns.get('k_locus_distribution', {}))}</div><div class="card-label">Capsule Types</div>
                <i class="fas fa-shield-alt fa-2x" style="color: var(--kaptive-color); margin-top: 10px;"></i>
            </div>
            <div class="dashboard-card card-amr" onclick="switchTab('amr')">
                <div class="card-number">{len(gene_centric.get('all_genes', []))}</div><div class="card-label">Unique Genes</div>
                <i class="fas fa-biohazard fa-2x" style="color: var(--amr-color); margin-top: 10px;"></i>
            </div>
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')">
                <div class="card-number">{len(gene_centric.get('by_category', {}).get('High-Risk Virulence', []))}</div><div class="card-label">Virulence Genes</div>
                <i class="fas fa-virus fa-2x" style="color: var(--virulence-color); margin-top: 10px;"></i>
            </div>
            <div class="dashboard-card card-patterns" onclick="switchTab('patterns')">
                <div class="card-number">{len(patterns.get('high_risk_combinations', []))}</div><div class="card-label">High‑Risk Combos</div>
                <i class="fas fa-project-diagram fa-2x" style="color: var(--patterns-color); margin-top: 10px;"></i>
            </div>
        </div>
        
        <div class="tab-navigation">
            <button class="tab-button summary active" onclick="switchTab('summary')"><i class="fas fa-chart-pie"></i> Summary</button>
            <button class="tab-button samples" onclick="switchTab('samples')"><i class="fas fa-list-alt"></i> Sample Overview</button>
            <button class="tab-button mlst" onclick="switchTab('mlst')"><i class="fas fa-code-branch"></i> MLST</button>
            <button class="tab-button qc" onclick="switchTab('qc')"><i class="fas fa-chart-line"></i> FASTA QC</button>
            <button class="tab-button kaptive" onclick="switchTab('kaptive')"><i class="fas fa-shield-alt"></i> Kaptive</button>
            <button class="tab-button amr" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR Genes</button>
            <button class="tab-button virulence" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence Genes</button>
            <button class="tab-button plasmids" onclick="switchTab('plasmids')"><i class="fas fa-dna"></i> Plasmids</button>
            <button class="tab-button patterns" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
            <button class="tab-button databases" onclick="switchTab('databases')"><i class="fas fa-database"></i> Databases</button>
            <button class="tab-button aiguide" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
            <button class="tab-button export" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
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
        <div id="export-tab" class="tab-content">{self._generate_export_section(data)}</div>
        
        <div class="footer">
            <h3>Kleboscope Ultimate Reporter v1.0.0</h3>
            <p>University of Ghana Medical School | Brown Beckley <brownbeckley94@gmail.com></p>
            <p>Generated on {metadata.get('analysis_date', 'Unknown')}</p>
            <p><strong>Critical Genes Tracked:</strong> Carbapenemases (KPC, NDM, OXA-48) • Colistin (mcr) • Tigecycline (tetX) • ICEKp Markers (ybt, clb, iro, rmp) • Virulence Plasmid Markers (iro, iuc, rmp, rmpA2) • Biocides & Heavy Metals (qac, sil, mer, ars, pco) • Adhesins (fim, mrk, ecp) • Secretion Systems (tss) • Siderophores • Toxins</p>
        </div>
    </div>
</body>
</html>""")
        return ''.join(html_parts)

    # --------------------------------------------------------------------------
    # Section generators
    # --------------------------------------------------------------------------
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
        patterns = data.get('patterns', {})

        # Main sample table
        rows = []
        for sample, d in sorted(samples.items()):
            st = d.get('mlst', {}).get('ST', 'ND')
            k = d.get('kaptive', {}).get('K_Locus', 'ND')
            o = d.get('kaptive', {}).get('O_Locus', 'ND')
            # Count virulence genes from all virulence databases (vfdb, ecoli_vf)
            ab_genes = d.get('abricate_genes', {})
            virulence_count = 0
            for db, genes in ab_genes.items():
                if db in ['vfdb', 'ecoli_vf']:
                    virulence_count += len(genes)
            rows.append(f"<tr><td class='col-sample'><strong>{sample}</strong></td><td class='col-st'>{st}</td><td class='col-k-locus'>{k}</td><td class='col-o-locus'>{o}</td><td class='col-frequency'>{virulence_count}</td></tr>")
        main_table = f"""
        <h2 class="section-header samples-header"><i class="fas fa-list-alt"></i> Sample Overview</h2>
        <input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Search by sample name, ST, K locus...">
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','sample_overview.csv')"><i class="fas fa-download"></i> Export CSV</button><button class="action-btn btn-success" onclick="document.getElementById('search-samples').value=''; searchTable('samples-table','search-samples')"><i class="fas fa-sync"></i> Clear Search</button></div>
        <div class="master-scrollable-container"><table id="samples-table" class="data-table"><thead><tr><th>Sample</th><th>ST</th><th>K Locus</th><th>O Locus</th><th>Virulence Genes Count</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

        # Combination tables
        ko_combos = patterns.get('ko_combinations', {})
        st_k_combos = patterns.get('st_k_combinations', {})
        st_o_combos = patterns.get('st_o_combinations', {})
        st_ko_combos = patterns.get('st_ko_combinations', {})

        def make_combo_table(combo_dict, title, col1, col2, search_id):
            if not combo_dict:
                return f"<div class='alert-box alert-warning'>No {title} data available.</div>"
            rows = []
            for combo, samples_list in sorted(combo_dict.items(), key=lambda x: len(x[1]), reverse=True):
                sample_list = ', '.join(samples_list[:1000]) + (' ...' if len(samples_list) > 1000 else '')
                rows.append(f"<tr><td><strong>{combo}</strong></td><td>{len(samples_list)}</td><td>{sample_list}</td></tr>")
            table = f"""
            <h3>{title}</h3>
            <input type="text" class="search-box" id="{search_id}" onkeyup="searchTable('{search_id}-table', '{search_id}')" placeholder="🔍 Search by {col1} or {col2}...">
            <div class="master-scrollable-container"><table id="{search_id}-table" class="data-table"><thead><tr><th>{col1}</th><th>Frequency</th><th>Samples</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
            """
            return table

        combos_html = f"""
        <h3>Combination Analysis</h3>
        {make_combo_table(ko_combos, "K:O Capsule Types", "K Locus:O Locus", "K:O", "search-ko")}
        {make_combo_table(st_k_combos, "ST-K Locus Combinations", "ST", "K Locus", "search-stk")}
        {make_combo_table(st_o_combos, "ST-O Locus Combinations", "ST", "O Locus", "search-sto")}
        {make_combo_table(st_ko_combos, "ST-K:O Combinations", "ST", "K:O", "search-stko")}
        """
        return main_table + combos_html

    def _generate_mlst_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        st_dist = patterns.get('st_distribution', {})
        if not st_dist:
            return """
            <h2 class="section-header mlst-header"><i class="fas fa-code-branch"></i> MLST Distribution</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No MLST Data</h3><p>The MLST report was missing or contained no data.</p></div></div>
            """
        rows = []
        total = sum(st_dist.values())
        for st, cnt in sorted(st_dist.items(), key=lambda x: x[1], reverse=True):
            pct = cnt / total * 100 if total else 0
            rows.append(f"<tr><td class='col-st'><strong>ST{st}</strong></td><td class='col-frequency'>{cnt}</td><td class='col-frequency'>{pct:.1f}%</td></tr>")
        return f"""
        <h2 class="section-header mlst-header"><i class="fas fa-code-branch"></i> MLST Distribution</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>MLST Analysis</h3><p>Each ST is shown with its frequency. Click on filter buttons to highlight specific clones.</p></div></div>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>ST</th><th>Count</th><th>Percentage</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_qc_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        # Check if any sample has QC data
        if not any(d.get('qc') for d in samples.values()):
            return """
            <h2 class="section-header qc-header"><i class="fas fa-chart-line"></i> FASTA Quality Control Metrics</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No QC Data</h3><p>The FASTA QC report was missing or contained no data.</p></div></div>
            """
        # Collect all possible QC metrics from the first sample
        all_metrics = set()
        for d in samples.values():
            qc = d.get('qc', {})
            all_metrics.update(qc.keys())
        # Sort metrics for consistent columns
        metric_list = sorted(all_metrics)
        rows = []
        for sample, d in sorted(samples.items()):
            qc = d.get('qc', {})
            row = [f"<td class='col-sample'><strong>{sample}</strong>"]
            for metric in metric_list:
                val = qc.get(metric, 'ND')
                if isinstance(val, float):
                    val = f"{val:,.0f}" if val > 1000 else f"{val:.1f}"
                row.append(f"<td class='col-frequency'>{val}")
            rows.append("<tr>" + "".join(row) + "</tr>")
        # Build header with data-sort attributes
        header = "<thead><tr><th data-sort='string'>Sample</th>" + \
                 "".join([f"<th data-sort='number'>{m}</th>" for m in metric_list]) + \
                 "</tr></thead>"
        return f"""
        <h2 class="section-header qc-header"><i class="fas fa-chart-line"></i> FASTA Quality Control Metrics</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Assembly Quality</h3><p>Click on column headers to sort. Metrics like N50, total length, and GC% indicate assembly completeness.</p></div></div>
        <input type="text" class="search-box" id="search-qc" onkeyup="searchTable('qc-table','search-qc')" placeholder="🔍 Search by sample name...">
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="exportTableToCSV('qc-table','fasta_qc.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <div class="master-scrollable-container"><table id="qc-table" class="data-table">{header}<tbody>{"".join(rows)}</tbody></table></div>
        <script>
        function sortTable(table, col, type) {{
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.rows);
            const isAscending = table.getAttribute('data-sort-dir') !== 'asc';
            rows.sort((a, b) => {{
                let aVal = a.cells[col].innerText.trim();
                let bVal = b.cells[col].innerText.trim();
                if (type === 'number') {{
                    aVal = parseFloat(aVal.replace(/,/g, '')) || 0;
                    bVal = parseFloat(bVal.replace(/,/g, '')) || 0;
                    return isAscending ? aVal - bVal : bVal - aVal;
                }} else {{
                    return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                }}
            }});
            tbody.append(...rows);
            table.setAttribute('data-sort-dir', isAscending ? 'asc' : 'desc');
        }}
        document.addEventListener('DOMContentLoaded', function() {{
            const table = document.getElementById('qc-table');
            const headers = table.querySelectorAll('th');
            headers.forEach((header, idx) => {{
                header.style.cursor = 'pointer';
                header.addEventListener('click', () => {{
                    const type = header.getAttribute('data-sort') || 'string';
                    sortTable(table, idx, type);
                }});
            }});
        }});
        </script>
        """

    def _generate_kaptive_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        if not any(d.get('kaptive') for d in samples.values()):
            return """
            <h2 class="section-header kaptive-header"><i class="fas fa-shield-alt"></i> Kaptive Capsule Typing</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No Kaptive Data</h3><p>The Kaptive report was missing or contained no data.</p></div></div>
            """
        rows = []
        for sample, d in sorted(samples.items()):
            k = d.get('kaptive', {})
            if k:
                rows.append(f"<tr><td class='col-sample'><strong>{sample}</strong></td><td class='col-k-locus'>{k.get('K_Locus', 'ND')}</td><td class='col-o-locus'>{k.get('O_Locus', 'ND')}</td><td class='col-frequency'>{k.get('K_Identity', 'ND')}</td><td class='col-frequency'>{k.get('K_Coverage', 'ND')}</td></tr>")
        if not rows:
            return """
            <h2 class="section-header kaptive-header"><i class="fas fa-shield-alt"></i> Kaptive Capsule Typing</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No Kaptive Data</h3><p>No capsule typing data was found in the report.</p></div></div>
            """
        return f"""
        <h2 class="section-header kaptive-header"><i class="fas fa-shield-alt"></i> Kaptive Capsule Typing</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Capsule & O‑antigen Typing</h3><p>Kaptive identifies K (capsule) and O (lipopolysaccharide) loci. High identity/coverage indicates reliable typing.</p></div></div>
        <input type="text" class="search-box" id="search-kaptive" onkeyup="searchTable('kaptive-table','search-kaptive')" placeholder="🔍 Search by sample, K locus...">
        <div class="master-scrollable-container"><table id="kaptive-table" class="data-table"><thead><tr><th>Sample</th><th>K Locus</th><th>O Locus</th><th>K Identity</th><th>K Coverage</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    # --------------------------------------------------------------------------
    # UPDATED AMR SECTION – CORRECT ITERATION OVER ALL DATABASES
    # --------------------------------------------------------------------------
    def _generate_amr_section(self, data: Dict) -> str:
        """Generate AMR genes section – per database, from all AMR databases."""
        gene_freqs = data.get('gene_frequencies', {})
        total_samples = len(data.get('samples', {}))

        # Collect all AMR genes from all databases, excluding virulence and plasmid
        all_amr_rows = []

        for db_name, db_genes in gene_freqs.items():
            # Skip virulence and plasmid databases
            if db_name in ['vfdb', 'ecoli_vf', 'plasmidfinder']:
                continue

            for gene, gdata in db_genes.items():
                # Determine risk level (AMRfinder has it, others default)
                risk = gdata.get('risk_level', 'Standard')
                # Frequency display
                freq_disp = gdata.get('frequency_display', f"{gdata['count']} ({gdata.get('percentage', 0):.1f}%)")
                all_amr_rows.append({
                    'gene': gene,
                    'category': self.analyzer.categorize_gene(gene),
                    'database': db_name.upper(),
                    'frequency_display': freq_disp,
                    'risk_level': risk,
                    'genomes': gdata.get('genomes', [])
                })

        if not all_amr_rows:
            return """
            <h2 class="section-header amr-header"><i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No AMR Data</h3><p>No AMR genes were detected or the required reports were missing.</p></div></div>
            """

        # Sort by count descending
        all_amr_rows.sort(key=lambda x: len(x['genomes']), reverse=True)

        rows = []
        for row in all_amr_rows:
            genome_tags = ''.join([f'<span class="genome-tag">{gen}</span>' for gen in row['genomes']])
            risk_badge = "badge-critical" if row['risk_level'] == "CRITICAL" else "badge-high" if row['risk_level'] == "HIGH" else "badge-medium"
            rows.append(f"""
            <tr>
                <td class="col-gene"><strong>{row['gene']}</strong></td>
                <td class="col-category"><span class="category-chip chip-{row['category'].lower().replace(' ', '-')}">{row['category']}</span></td>
                <td class="col-database">{row['database']}</td>
                <td class="col-risk"><span class="badge {risk_badge}">{row['risk_level']}</span></td>
                <td class="col-frequency"><span class="frequency-display">{row['frequency_display']}</span></td>
                <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
            </tr>
            """)

        return f"""
        <h2 class="section-header amr-header"><i class="fas fa-biohazard"></i> Antimicrobial Resistance Genes</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>AMR Gene Analysis</h3><p>Each AMR gene is shown per database with all genomes that contain it. Frequency displayed as <strong>count (percentage%)</strong>. Critical resistance genes are highlighted.</p></div></div>
        <input type="text" class="search-box" id="search-amr" onkeyup="searchTable('amr-table','search-amr')" placeholder="🔍 Search AMR genes by name...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('amr-table','amr_genes.csv')"><i class="fas fa-download"></i> Export All AMR Genes</button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-amr').value='Carbapenem'; searchTable('amr-table','search-amr')"><i class="fas fa-skull-crossbones"></i> Show Carbapenemases</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-amr').value='OXA'; searchTable('amr-table','search-amr')"><i class="fas fa-dna"></i> Show OXA genes</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-amr').value='KPC'; searchTable('amr-table','search-amr')"><i class="fas fa-biohazard"></i> Show KPC genes</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-amr').value='NDM'; searchTable('amr-table','search-amr')"><i class="fas fa-vial"></i> Show NDM genes</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='qac'; searchTable('amr-table','search-amr')"><i class="fas fa-flask"></i> Show Biocides (qac)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='sil'; searchTable('amr-table','search-amr')"><i class="fas fa-coins"></i> Show Silver (sil)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='mer'; searchTable('amr-table','search-amr')"><i class="fas fa-skull"></i> Show Mercury (mer)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='ars'; searchTable('amr-table','search-amr')"><i class="fas fa-skull"></i> Show Arsenic (ars)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-amr').value='pco'; searchTable('amr-table','search-amr')"><i class="fas fa-coins"></i> Show Copper (pco)</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-amr').value=''; searchTable('amr-table','search-amr')"><i class="fas fa-sync"></i> Clear Search</button>
        </div>
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #F44336;">
            <strong><i class="fas fa-info-circle"></i> Role of each gene family:</strong><br>
            • <strong>Carbapenemases</strong> (KPC, NDM, OXA) – confer resistance to carbapenems, last‑line antibiotics.<br>
            • <strong>ESBLs</strong> (CTX‑M, SHV, TEM) – hydrolyze extended‑spectrum cephalosporins.<br>
            • <strong>AmpC</strong> – broad‑spectrum β‑lactamases, often plasmid‑borne.<br>
            • <strong>Colistin Resistance</strong> (mcr, pmr, lpx) – resistance to polymyxins, last‑resort for MDR.<br>
            • <strong>Tigecycline Resistance</strong> (tetX, efflux pumps) – resistance to tigecycline, an important option.<br>
            • <strong>Biofilm Formation</strong> (ompA, csu, bfm) – promotes persistence and device‑related infections.<br>
            • <strong>Efflux Pumps</strong> (ade, acr, mex) – multidrug efflux, contributes to MDR.<br>
            • <strong>Other</strong> – includes aminoglycosides, fluoroquinolones, sulfonamides, etc.
        </div>
        <div class="master-scrollable-container"><table id="amr-table" class="data-table"><thead>
            <tr><th class="col-gene">Gene</th><th class="col-category">Category</th><th class="col-database">Database</th><th class="col-risk">Risk Level</th><th class="col-frequency">Frequency</th><th class="col-genomes">Genomes</th></tr>
        </thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    # --------------------------------------------------------------------------
    # UPDATED VIRULENCE SECTION – COLLECT FROM VFDB AND ECOLI_VF ONLY
    # --------------------------------------------------------------------------
    def _generate_virulence_section(self, data: Dict) -> str:
        """Generate virulence genes section – per database, from VFDB/ecoli_vf only."""
        gene_freqs = data.get('gene_frequencies', {})
        total_samples = len(data.get('samples', {}))

        # Collect all virulence genes from VFDB and ecoli_vf (and any other that may contain virulence)
        vir_rows = []
        for db_name in ['vfdb', 'ecoli_vf']:
            if db_name in gene_freqs:
                for gene, gdata in gene_freqs[db_name].items():
                    vir_rows.append({
                        'gene': gene,
                        'category': self.analyzer.categorize_gene(gene),
                        'database': db_name.upper(),
                        'frequency_display': gdata.get('frequency_display', f"{gdata['count']} ({gdata.get('percentage', 0):.1f}%)"),
                        'genomes': gdata.get('genomes', [])
                    })

        if not vir_rows:
            return """
            <h2 class="section-header virulence-header"><i class="fas fa-virus"></i> Virulence Genes</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No Virulence Data</h3><p>No virulence genes were detected or the VFDB report was missing.</p></div></div>
            """

        # Sort by count descending
        vir_rows.sort(key=lambda x: len(x['genomes']), reverse=True)

        rows = []
        for row in vir_rows:
            genome_tags = ''.join([f'<span class="genome-tag">{gen}</span>' for gen in row['genomes']])
            rows.append(f"""
            <tr>
                <td class="col-gene"><strong>{row['gene']}</strong></td>
                <td class="col-category"><span class="category-chip chip-{row['category'].lower().replace(' ', '-')}">{row['category']}</span></td>
                <td class="col-database">{row['database']}</td>
                <td class="col-frequency"><span class="frequency-display">{row['frequency_display']}</span></td>
                <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
            </tr>
            """)

        return f"""
        <h2 class="section-header virulence-header"><i class="fas fa-virus"></i> Virulence Genes</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Virulence Gene Analysis</h3><p>Each virulence gene is shown per database with all genomes that contain it. Use the filter buttons below to focus on specific virulence families.</p></div></div>
        <input type="text" class="search-box" id="search-virulence" onkeyup="searchTable('virulence-table','search-virulence')" placeholder="🔍 Search virulence genes...">
        <div class="action-buttons">
            <button class="action-btn btn-primary" onclick="exportTableToCSV('virulence-table','virulence_genes.csv')"><i class="fas fa-download"></i> Export All Virulence Genes</button>
            <button class="action-btn btn-success" onclick="document.getElementById('search-virulence').value='ybt'; searchTable('virulence-table','search-virulence')"><i class="fas fa-dna"></i> Show Yersiniabactin (ybt)</button>
            <button class="action-btn btn-danger" onclick="document.getElementById('search-virulence').value='clb'; searchTable('virulence-table','search-virulence')"><i class="fas fa-dna"></i> Show Colibactin (clb)</button>
            <button class="action-btn btn-warning" onclick="document.getElementById('search-virulence').value='iro'; searchTable('virulence-table','search-virulence')"><i class="fas fa-dna"></i> Show Salmochelin (iro)</button>
            <button class="action-btn btn-info" onclick="document.getElementById('search-virulence').value='iuc'; searchTable('virulence-table','search-virulence')"><i class="fas fa-dna"></i> Show Aerobactin (iuc)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='rmp'; searchTable('virulence-table','search-virulence')"><i class="fas fa-dna"></i> Show Hypermucoidy (rmp)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='fim'; searchTable('virulence-table','search-virulence')"><i class="fas fa-brush"></i> Show Type 1 Fimbriae (fim)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='mrk'; searchTable('virulence-table','search-virulence')"><i class="fas fa-brush"></i> Show Type 3 Fimbriae (mrk)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='ecp'; searchTable('virulence-table','search-virulence')"><i class="fas fa-brush"></i> Show ECP Pili (ecp)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='tss'; searchTable('virulence-table','search-virulence')"><i class="fas fa-syringe"></i> Show T6SS (tss)</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='siderophore'; searchTable('virulence-table','search-virulence')"><i class="fas fa-droplet"></i> Show Siderophores</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='adhesin'; searchTable('virulence-table','search-virulence')"><i class="fas fa-hand-holding-heart"></i> Show Adhesins</button>
            <button class="action-btn btn-secondary" onclick="document.getElementById('search-virulence').value='toxin'; searchTable('virulence-table','search-virulence')"><i class="fas fa-skull"></i> Show Toxins</button>
            <button class="action-btn btn-light" onclick="document.getElementById('search-virulence').value=''; searchTable('virulence-table','search-virulence')"><i class="fas fa-sync"></i> Clear Search</button>
        </div>
        <div style="margin: 10px 0 20px 0; background: #f8f9fa; padding: 15px; border-radius: 8px; font-size: 0.9em; border-left: 4px solid #E91E63;">
            <strong><i class="fas fa-info-circle"></i> Role of each virulence family:</strong><br>
            • <strong>ybt</strong> – Yersiniabactin, siderophore for iron acquisition, linked to ICEKp.<br>
            • <strong>clb</strong> – Colibactin, genotoxin involved in DNA damage.<br>
            • <strong>iro</strong> – Salmochelin, siderophore for iron scavenging.<br>
            • <strong>iuc</strong> – Aerobactin, siderophore associated with hypervirulence.<br>
            • <strong>rmp</strong> – Regulator of mucoid phenotype, hypermucoviscosity.<br>
            • <strong>fim</strong> – Type 1 fimbriae, adherence to host surfaces.<br>
            • <strong>mrk</strong> – Type 3 fimbriae, biofilm formation and adhesion.<br>
            • <strong>ecp</strong> – E. coli common pilus, adherence.<br>
            • <strong>tss</strong> – Type VI secretion system, interbacterial competition and host interaction.<br>
            • <strong>siderophore</strong> – General iron‑scavenging systems.<br>
            • <strong>adhesin</strong> – Surface factors promoting attachment.<br>
            • <strong>toxin</strong> – Cytotoxins and hemolysins.
        </div>
        <div class="master-scrollable-container"><table id="virulence-table" class="data-table"><thead>
            <tr><th class="col-gene">Gene</th><th class="col-category">Category</th><th class="col-database">Database</th><th class="col-frequency">Frequency</th><th class="col-genomes">Genomes</th></tr>
        </thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    # --------------------------------------------------------------------------
    # Remaining unchanged sections (plasmids, patterns, databases, export, etc.)
    # --------------------------------------------------------------------------
    def _generate_plasmids_section(self, data: Dict) -> str:
        # unchanged from original
        gene_centric = data.get('gene_centric', {})
        all_genes = gene_centric.get('all_genes', [])
        # Filter: genes that have 'plasmidfinder' in their databases list
        plasmid_genes = [g for g in all_genes if 'plasmidfinder' in g.get('databases', [])]
        if not plasmid_genes:
            return """
            <h2 class="section-header plasmids-header"><i class="fas fa-dna"></i> Plasmid Replicons</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No Plasmid Data</h3><p>No plasmid replicons detected or the PlasmidFinder report was missing.</p></div></div>
            """
        rows = []
        for g in plasmid_genes:
            genome_tags = ''.join([f'<span class="genome-tag">{gen}</span>' for gen in g['genomes']])
            db_str = ', '.join([db.upper() for db in g['databases']])
            rows.append(f"""
            <tr>
                <td class='col-gene'><strong>{g['gene']}</strong></td>
                <td class='col-database'>{db_str}</td>
                <td class='col-frequency'><span class='frequency-display'>{g['frequency_display']}</span></td>
                <td class='col-genomes'><div class='genome-list'>{genome_tags}</div></td>
            </tr>
            """)
        return f"""
        <h2 class="section-header plasmids-header"><i class="fas fa-dna"></i> Plasmid Replicons</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Plasmid Analysis</h3><p>Plasmid replicons indicate the presence of mobile genetic elements that can carry resistance and virulence genes. Each replicon is shown with all genomes that contain it.</p></div></div>
        <input type="text" class="search-box" id="search-plasmids" onkeyup="searchTable('plasmids-table','search-plasmids')" placeholder="🔍 Search by replicon name...">
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="exportTableToCSV('plasmids-table','plasmid_replicons.csv')"><i class="fas fa-download"></i> Export Plasmid Replicons</button></div>
        <div class="master-scrollable-container"><table id="plasmids-table" class="data-table"><thead><tr><th>Plasmid Replicon</th><th>Database</th><th>Frequency</th><th>Genomes</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_patterns_section(self, data: Dict) -> str:
        # unchanged from original
        patterns = data.get('patterns', {})
        high_risk = patterns.get('high_risk_combinations', [])
        icekp = patterns.get('icekp_marker_presence', {})
        vp = patterns.get('virulence_plasmid_marker_presence', {})
        html = "<h2 class='section-header patterns-header'><i class='fas fa-project-diagram'></i> Pattern Discovery</h2>"
        if high_risk:
            rows = []
            for combo in high_risk:
                rows.append(f"<tr><td class='col-sample'><strong>{combo['sample']}</strong></td><td class='col-st'>{combo['st']}</td><td class='col-k-locus'>{combo['k_locus']}</td><td class='col-gene'>{', '.join(combo['critical_resistance'])}</td><td class='col-gene'>{', '.join(combo['high_risk_virulence'])}</td></tr>")
            html += f"""
            <div class='alert-box alert-danger'><i class='fas fa-radiation fa-2x'></i><div><h3>⚠️ High‑Risk Combinations</h3><p>{len(high_risk)} samples carry both critical resistance and high‑risk virulence genes.</p></div></div>
            <div class='master-scrollable-container'><table class='data-table'><thead><tr><th>Sample</th><th>ST</th><th>K Locus</th><th>Critical Resistance</th><th>High‑Risk Virulence</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
            """
        else:
            html += "<div class='alert-box alert-success'><i class='fas fa-check-circle'></i><div>No high‑risk combinations detected.</div></div>"
        if icekp:
            rows = []
            for sample, markers in icekp.items():
                rows.append(f"<tr><td class='col-sample'><strong>{sample}</strong></td><td class='col-gene'>{', '.join(markers)}</td></tr>")
            html += f"""
            <h3>ICEKp Associated Markers</h3>
            <div class='master-scrollable-container'><table class='data-table'><thead><tr><th>Sample</th><th>Markers (ybt, clb, iro, rmp)</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
            """
        if vp:
            rows = []
            for sample, markers in vp.items():
                rows.append(f"<tr><td class='col-sample'><strong>{sample}</strong></td><td class='col-gene'>{', '.join(markers)}</td></tr>")
            html += f"""
            <h3>Virulence Plasmid Associated Markers</h3>
            <div class='master-scrollable-container'><table class='data-table'><thead><tr><th>Sample</th><th>Markers (iro, iuc, rmp, rmpA2)</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
            """
        return html

    def _generate_databases_section(self, data: Dict) -> str:
        # unchanged from original
        gene_centric = data.get('gene_centric', {})
        patterns = data.get('patterns', {})
        database_coverage = patterns.get('database_coverage', {})
        database_stats = gene_centric.get('database_stats', {})

        if not database_stats and not database_coverage:
            return """
            <h2 class="section-header databases-header"><i class="fas fa-database"></i> Database Coverage</h2>
            <div class="alert-box alert-warning"><i class="fas fa-exclamation-circle fa-2x"></i><div><h3>No Database Statistics</h3><p>No database coverage data available.</p></div></div>
            """

        # Build database coverage table
        all_dbs = set(database_stats.keys()) | set(database_coverage.keys())
        if not all_dbs:
            all_dbs = list(database_stats.keys()) + list(database_coverage.keys())
        all_dbs = sorted(all_dbs)

        rows = []
        for db in all_dbs:
            stats = database_stats.get(db, {})
            coverage = database_coverage.get(db, {})
            total_samples = coverage.get('total_samples', 0)
            samples_with_hits = coverage.get('samples_with_hits', stats.get('total_occurrences', 0))
            coverage_pct = coverage.get('coverage_percentage', (samples_with_hits / total_samples * 100) if total_samples else 0)
            rows.append(f"""
            <tr>
                <td><strong>{db.upper()}</strong></td>
                <td>{stats.get('total_genes', 0)}</td>
                <td>{stats.get('total_occurrences', 0)}</td>
                <td>{samples_with_hits} / {total_samples}</td>
                <td>{coverage_pct:.1f}%</td>
                <td>{stats.get('critical_genes', 0)}</td>
            </tr>
            """)

        return f"""
        <h2 class="section-header databases-header"><i class="fas fa-database"></i> Database Coverage and Performance</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Database Statistics</h3><p>This table shows how many unique genes and total occurrences were detected per database, along with coverage across samples and critical gene counts.</p></div></div>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Database</th>
                        <th>Unique Genes</th>
                        <th>Total Occurrences</th>
                        <th>Samples with Hits</th>
                        <th>Coverage (%)</th>
                        <th>Critical Genes</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(rows)}
                </tbody>
            </table>
        </div>
        """

    def _generate_aiguide_section(self, data: Dict) -> str:
        # unchanged from original
        return """
        <h2 class="section-header aiguide-header"><i class="fas fa-robot"></i> AI Assistant Guide</h2>
        <div class="alert-box alert-info"><i class="fas fa-lightbulb fa-2x"></i><div><h3>How to Use AI with This Report</h3><p>You can upload this HTML report (or its JSON data) to AI tools like ChatGPT, Claude, or Gemini to ask detailed questions about your K. pneumoniae dataset. Below are example prompts to get you started.</p></div></div>
        
        <div style="margin: 20px 0;">
            <div class="database-section">
                <h4><i class="fas fa-chart-line"></i> General Questions</h4>
                <ul style="margin-left: 20px;">
                    <li>What are the most common STs in this dataset?</li>
                    <li>Which K and O capsule types are dominant?</li>
                    <li>How many samples carry carbapenemase genes? Which ones?</li>
                    <li>Are there any samples with both colistin resistance (mcr) and carbapenemases?</li>
                    <li>What is the overall quality of the assemblies (N50, total length, GC%)?</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-virus"></i> Virulence Analysis</h4>
                <ul style="margin-left: 20px;">
                    <li>Which virulence genes are most prevalent in hypervirulent clones (ST23, ST86)?</li>
                    <li>List all samples carrying the ICEKp markers (ybt, clb, iro, rmp).</li>
                    <li>Which samples contain the aerobactin (iuc) system? Are they associated with specific STs?</li>
                    <li>Show me the co-occurrence of rmpA and iuc in the same genome.</li>
                    <li>Which samples carry fim, mrk, or ecp adhesins?</li>
                    <li>Identify isolates with Type VI secretion system (T6SS) genes.</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-biohazard"></i> Resistance Patterns</h4>
                <ul style="margin-left: 20px;">
                    <li>What is the prevalence of blaKPC, blaNDM, and blaOXA-48-like genes?</li>
                    <li>Which STs are most associated with carbapenem resistance?</li>
                    <li>Are there samples with high-level aminoglycoside resistance (armA, rmtB)?</li>
                    <li>Which resistance genes co-occur most frequently? Provide examples.</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-project-diagram"></i> Pattern Discovery</h4>
                <ul style="margin-left: 20px;">
                    <li>Identify samples with high‑risk combinations (critical resistance + hypervirulence).</li>
                    <li>Which ST‑K:O combinations are associated with MDR or hypervirulence?</li>
                    <li>Show me the distribution of ICEKp markers across different STs.</li>
                </ul>
            </div>
            
            <div class="database-section">
                <h4><i class="fas fa-upload"></i> How to Upload Data</h4>
                <p>If you are using a text-based AI (like this one), you can copy the relevant tables or the JSON data and paste them into the conversation. For example, you can say: "I'm sending you the sample overview table from the Kleboscope report. Can you analyze the ST distribution?"</p>
                <p>Alternatively, you can save the <strong>kleboscope_ultimate_report.json</strong> file and ask the AI to parse it. Many AI tools accept file uploads.</p>
            </div>
        </div>
        """

    def _generate_export_section(self, data: Dict) -> str:
        # unchanged from original
        return """
        <h2 class="section-header export-header"><i class="fas fa-download"></i> Export Data</h2>
        <div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i><div><h3>Export Options</h3><p>Download data as CSV files for further analysis.</p></div></div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:20px;margin:20px 0;">
            <div class="dashboard-card" onclick="exportTableToCSV('samples-table','sample_overview.csv')"><i class="fas fa-file-csv fa-3x"></i><div>Sample Overview CSV</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('amr-table','amr_genes.csv')"><i class="fas fa-biohazard fa-3x"></i><div>AMR Genes CSV</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('virulence-table','virulence_genes.csv')"><i class="fas fa-virus fa-3x"></i><div>Virulence Genes CSV</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('plasmids-table','plasmid_replicons.csv')"><i class="fas fa-dna fa-3x"></i><div>Plasmid Replicons CSV</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('qc-table','fasta_qc.csv')"><i class="fas fa-chart-line fa-3x"></i><div>FASTA QC CSV</div></div>
            <div class="dashboard-card" onclick="window.location.href='kleboscope_ultimate_report.json'"><i class="fas fa-file-code fa-3x"></i><div>Complete JSON</div></div>
        </div>
        """


# =============================================================================
# MAIN REPORTER
# =============================================================================
class KleboscopeUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "KLEBOSCOPE_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = KleboHTMLParser()
        self.analyzer = KleboDataAnalyzer()
        self.generator = KleboHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "Kleboscope Ultimate Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def find_html_files(self) -> Dict[str, List[Path]]:
        """Find all HTML report files."""
        html_files = {
            'mlst': [],
            'qc': [],
            'kaptive': [],
            'amrfinder': [],
            'abricate': defaultdict(list)
        }
        all_html = list(self.input_dir.glob("**/*.html"))
        if not all_html:
            print("  ⚠️ No HTML files found in the directory!")
            return html_files

        print(f"  📁 Found {len(all_html)} HTML files")
        for f in all_html:
            name = f.name.lower()
            if 'mlst' in name:
                html_files['mlst'].append(f)
            elif 'fasta_qc' in name or 'qc' in name:
                html_files['qc'].append(f)
            elif 'kaptive' in name:
                html_files['kaptive'].append(f)
            elif 'amrfinder' in name:
                html_files['amrfinder'].append(f)
            else:
                # ABRicate database files – match by database name
                for db in self.parser.abricate_databases:
                    if db in name:
                        html_files['abricate'][db].append(f)
                        break
        return html_files

    def integrate_all_data(self, html_files: Dict) -> Dict:
        integrated = {'metadata': self.metadata, 'samples': {}}
        all_samples = set()

        # Parse MLST
        if html_files['mlst']:
            mlst_data = self.parser.parse_mlst_report(html_files['mlst'][0])
            for s, d in mlst_data.items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['mlst'] = d

        # Parse FASTA QC
        if html_files['qc']:
            qc_data = self.parser.parse_qc_report(html_files['qc'][0])
            for s, d in qc_data.items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['qc'] = d

        # Parse Kaptive
        if html_files['kaptive']:
            kaptive_data = self.parser.parse_kaptive_report(html_files['kaptive'][0])
            for s, d in kaptive_data.items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['kaptive'] = d

        total_samples = len(all_samples)

        # Parse AMRfinder
        amr_gene_freq = {}
        if html_files['amrfinder']:
            amr_genes, amr_gene_freq = self.parser.parse_amrfinder_report(html_files['amrfinder'][0], total_samples)
            for s, genes in amr_genes.items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['amr_genes'] = genes

        # Parse all ABRicate databases
        abricate_gene_freq = {}
        if html_files['abricate']:
            for db, files in html_files['abricate'].items():
                if files:
                    genes_by_sample, gene_freq = self.parser.parse_abricate_database_report(files[0], total_samples)
                    # Store per-sample genes
                    for s, genes in genes_by_sample.items():
                        all_samples.add(s)
                        integrated['samples'].setdefault(s, {}).setdefault('abricate_genes', {})[db] = genes
                    # Store gene frequencies
                    for gene, data in gene_freq.items():
                        abricate_gene_freq.setdefault(db, {})[gene] = data

        # Ensure all samples have entries
        for s in all_samples:
            integrated['samples'].setdefault(s, {})

        # Combine gene frequencies
        gene_freqs = {}
        if amr_gene_freq:
            gene_freqs['amrfinder'] = amr_gene_freq
        if abricate_gene_freq:
            for db, freq in abricate_gene_freq.items():
                gene_freqs[db] = freq

        # Store in integrated for HTML generator access
        integrated['gene_frequencies'] = gene_freqs

        # Gene‑centric tables (for other parts)
        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(gene_freqs, total_samples)

        # Patterns
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(
            integrated['samples'], integrated['gene_centric']
        )

        # Also compute database coverage for the databases tab
        if integrated['samples'] and 'gene_centric' in integrated:
            gene_centric = integrated['gene_centric']
            all_samples_set = set(integrated['samples'].keys())
            coverage = {}
            for db, genes in gene_freqs.items():
                samples_with_hits = set()
                for gdata in genes.values():
                    samples_with_hits.update(gdata['genomes'])
                coverage[db] = {
                    'samples_with_hits': len(samples_with_hits),
                    'total_samples': len(all_samples_set),
                    'coverage_percentage': (len(samples_with_hits) / len(all_samples_set) * 100) if all_samples_set else 0,
                }
            integrated['patterns']['database_coverage'] = coverage

        return integrated

    def generate_json_report(self, integrated: Dict) -> Path:
        out = self.output_dir / "kleboscope_ultimate_report.json"
        with open(out, 'w', encoding='utf-8') as f:
            json.dump(integrated, f, indent=2, default=str)
        return out

    def generate_csv_reports(self, integrated: Dict):
        # Sample overview CSV
        rows = []
        for sample, d in integrated['samples'].items():
            rows.append({
                'Sample': sample,
                'ST': d.get('mlst', {}).get('ST', 'ND'),
                'K_Locus': d.get('kaptive', {}).get('K_Locus', 'ND'),
                'O_Locus': d.get('kaptive', {}).get('O_Locus', 'ND'),
                'N50': d.get('qc', {}).get('N50', 'ND'),
                'GC%': d.get('qc', {}).get('GC%', 'ND'),
                'Virulence_Genes_Count': sum(len(genes) for db, genes in d.get('abricate_genes', {}).items() if db in ['vfdb', 'ecoli_vf'])
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)

        # AMR genes CSV (all non‑virulence, non‑plasmid genes)
        amr_rows = []
        for g in integrated['gene_centric'].get('all_genes', []):
            if 'plasmidfinder' in g.get('databases', []):
                continue
            if any(db in ['vfdb', 'ecoli_vf'] for db in g.get('databases', [])):
                continue
            amr_rows.append({
                'Gene': g['gene'],
                'Category': g['category'],
                'Database': ', '.join(g['databases']),
                'Risk': g['risk_level'],
                'Count': g['count'],
                'Frequency': g['frequency_display'],
                'Genomes': ';'.join(g['genomes'])
            })
        if amr_rows:
            pd.DataFrame(amr_rows).to_csv(self.output_dir / "amr_genes.csv", index=False)

        # Virulence genes CSV
        vir_rows = []
        for g in integrated['gene_centric'].get('all_genes', []):
            if g['category'] == 'High-Risk Virulence' or any(db in ['vfdb', 'ecoli_vf'] for db in g.get('databases', [])):
                vir_rows.append({
                    'Gene': g['gene'],
                    'Database': ', '.join(g['databases']),
                    'Count': g['count'],
                    'Frequency': g['frequency_display'],
                    'Genomes': ';'.join(g['genomes'])
                })
        if vir_rows:
            pd.DataFrame(vir_rows).to_csv(self.output_dir / "virulence_genes.csv", index=False)

        # Plasmid replicons CSV
        plasmid_rows = []
        for g in integrated['gene_centric'].get('all_genes', []):
            if 'plasmidfinder' in g.get('databases', []):
                plasmid_rows.append({
                    'Plasmid_Replicon': g['gene'],
                    'Database': ', '.join(g['databases']),
                    'Count': g['count'],
                    'Frequency': g['frequency_display'],
                    'Genomes': ';'.join(g['genomes'])
                })
        if plasmid_rows:
            pd.DataFrame(plasmid_rows).to_csv(self.output_dir / "plasmid_replicons.csv", index=False)

        # FASTA QC CSV
        qc_rows = []
        for sample, d in integrated['samples'].items():
            qc = d.get('qc', {})
            row = {'Sample': sample}
            row.update(qc)
            qc_rows.append(row)
        if qc_rows:
            pd.DataFrame(qc_rows).to_csv(self.output_dir / "fasta_qc.csv", index=False)

    def run(self):
        print("=" * 80)
        print("🧬 Kleboscope Ultimate Reporter – K. pneumoniae Gene‑Centric Analysis v1.0.0")
        print("=" * 80)
        html_files = self.find_html_files()
        if not any(html_files.values()):
            print("❌ No HTML files found! Please ensure the input directory contains Kleboscope reports.")
            return False
        integrated = self.integrate_all_data(html_files)
        if not integrated['samples']:
            print("❌ No data integrated.")
            return False
        self.generate_json_report(integrated)
        self.generate_csv_reports(integrated)
        self.generator.generate_main_report(integrated, self.output_dir)
        print("✅ All reports generated successfully.")
        return True


def main():
    parser = argparse.ArgumentParser(description='Kleboscope Ultimate Reporter – Gene‑Centric K. pneumoniae Analysis')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing Kleboscope HTML reports')
    args = parser.parse_args()
    reporter = KleboscopeUltimateReporter(Path(args.input_dir))
    reporter.run()

if __name__ == "__main__":
    main()