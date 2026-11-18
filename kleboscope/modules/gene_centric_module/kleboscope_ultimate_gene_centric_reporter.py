#!/usr/bin/env python3
"""
Kleboscope Ultimate Reporter - Gene-centric analysis for K. pneumoniae
Full  architecture with Kleboscope blue background + per-tab colors.
Parses MLST / Kaptive / FASTA QC from TSV files.

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 2.3.0
"""

import os
import sys
import json
import re
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Any
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup


# =============================================================================
# PARSER - TSV for MLST/Kaptive/QC, HTML for the rest
# =============================================================================
class KleboHTMLParser:
    """Parse every Kleboscope report into structured dictionaries."""

    def __init__(self):
        self.abricate_databases = [
            'card', 'resfinder', 'argannot', 'vfdb', 'plasmidfinder',
            'megares', 'ncbi', 'ecoh', 'ecoli_vf', 'bacmet2'
        ]
        self.db_name_mapping = {
            'klebo_card': 'card', 'klebo_resfinder': 'resfinder',
            'klebo_argannot': 'argannot', 'klebo_vfdb': 'vfdb',
            'klebo_plasmidfinder': 'plasmidfinder', 'klebo_megares': 'megares',
            'klebo_ncbi': 'ncbi', 'klebo_ecoh': 'ecoh',
            'klebo_ecoli_vf': 'ecoli_vf', 'klebo_bacmet2': 'bacmet2',
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        """Strip extension and directory prefix from a sample identifier."""
        sample = str(sample_id).strip()
        for ext in ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample

    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
        """Return the Nth HTML table from a string as a DataFrame."""
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            tables = soup.find_all('table')
            if not tables or table_index >= len(tables):
                return pd.DataFrame()
            table = tables[table_index]
            rows = table.find_all('tr')
            if not rows:
                return pd.DataFrame()
            headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
            data = []
            for row in rows[1:]:
                cols = row.find_all(['td', 'th'])
                row_data = [col.get_text().strip() for col in cols]
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
            print(f"    Table parsing error: {e}")
            return pd.DataFrame()

    def parse_mlst_tsv(self, file_path: Path) -> Dict[str, Dict]:
        """Parse the MLST summary TSV into a per-sample dict."""
        print(f"  Parsing MLST TSV: {file_path.name}")
        try:
            df = pd.read_csv(file_path, sep='\t', dtype=str).fillna('ND')
            results = {}
            for _, row in df.iterrows():
                sample_raw = row.get('sample', '')
                if not sample_raw or sample_raw == 'ND':
                    continue
                sample = self.normalize_sample_id(sample_raw)
                results[sample] = {
                    'ST':              str(row.get('st', 'ND')).strip(),
                    'MLST_Status':     str(row.get('mlst_status', 'ND')).strip(),
                    'Identity':        str(row.get('identity', 'ND')).strip(),
                    'Coverage':        str(row.get('coverage', 'ND')).strip(),
                    'Allele_Profile':  str(row.get('allele_profile', 'ND')).strip(),
                    'Clonal_Complex':  str(row.get('clonal_complex', 'ND')).strip(),
                    'Classification':  str(row.get('classification', 'ND')).strip(),
                    'Outbreak_Potential': str(row.get('outbreak_potential', 'ND')).strip(),
                }
            print(f"    Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    Error parsing MLST TSV: {e}")
            return {}

    def parse_kaptive_tsv(self, file_path: Path) -> Dict[str, Dict]:
        """Parse the Kaptive summary TSV into per-sample K/O typing."""
        print(f"  Parsing Kaptive TSV: {file_path.name}")
        try:
            df = pd.read_csv(file_path, sep='\t', dtype=str).fillna('ND')
            results = {}
            for _, row in df.iterrows():
                sample_raw = row.get('sample', '')
                if not sample_raw or sample_raw == 'ND':
                    continue
                sample = self.normalize_sample_id(sample_raw)
                results[sample] = {
                    'Status':       str(row.get('status', 'ND')).strip(),
                    'K_Locus':      str(row.get('k_locus', 'ND')).strip(),
                    'K_Type':       str(row.get('k_type', 'ND')).strip(),
                    'K_Confidence': str(row.get('k_confidence', 'ND')).strip(),
                    'K_Identity':   str(row.get('k_identity', 'ND')).strip(),
                    'K_Coverage':   str(row.get('k_coverage', 'ND')).strip(),
                    'O_Locus':      str(row.get('o_locus', 'ND')).strip(),
                    'O_Type':       str(row.get('o_type', 'ND')).strip(),
                    'O_Confidence': str(row.get('o_confidence', 'ND')).strip(),
                    'O_Identity':   str(row.get('o_identity', 'ND')).strip(),
                    'O_Coverage':   str(row.get('o_coverage', 'ND')).strip(),
                }
            print(f"    Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    Error parsing Kaptive TSV: {e}")
            return {}

    def parse_qc_tsv(self, file_path: Path) -> Dict[str, Dict]:
        """Parse the FASTA QC summary TSV into per-sample QC metrics."""
        print(f"  Parsing FASTA QC TSV: {file_path.name}")
        try:
            df = pd.read_csv(file_path, sep='\t', dtype=str).fillna('ND')
            results = {}
            for _, row in df.iterrows():
                sample_raw = row.get('Filename', '')
                if not sample_raw or sample_raw == 'ND':
                    continue
                sample = self.normalize_sample_id(sample_raw)
                qc_data = {}
                for col in df.columns:
                    if col == 'Filename':
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc_data[col] = 'ND'
                    else:
                        cleaned = str(val).replace('%', '').replace(',', '').strip()
                        try:
                            qc_data[col] = float(cleaned)
                        except Exception:
                            qc_data[col] = str(val)
                results[sample] = qc_data
            print(f"    Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    Error parsing QC TSV: {e}")
            return {}

    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict, Dict]:
        """Parse AMRfinder summary HTML into (genes_by_sample, gene_frequencies)."""
        print(f"  Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}

            df_genomes = self.parse_html_table(str(tables[0]), 0)
            genes_by_genome = {}
            if not df_genomes.empty:
                sample_col = next(
                    (c for c in df_genomes.columns if 'genome' in c.lower()),
                    df_genomes.columns[0],
                )
                for _, row in df_genomes.iterrows():
                    sample_raw = row[sample_col]
                    if not sample_raw:
                        continue
                    sample = self.normalize_sample_id(sample_raw)
                    genes = []
                    for col in df_genomes.columns:
                        if 'genes' in col.lower() or 'detected' in col.lower():
                            gene_str = str(row[col])
                            genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                            break
                    genes_by_genome[sample] = genes

            df_freq = self.parse_html_table(str(tables[1]), 0)
            gene_freq = {}
            if not df_freq.empty:
                for _, row in df_freq.iterrows():
                    gene = row.get('Gene', '')
                    if not gene:
                        continue
                    freq_str = row.get('Frequency', '0')
                    m = re.search(r'(\d+)', freq_str)
                    count = int(m.group(1)) if m else 0
                    pct = (count / total_samples * 100) if total_samples else 0
                    genomes_str = row.get('Genomes', '')
                    genomes = [self.normalize_sample_id(g.strip())
                               for g in genomes_str.split(',') if g.strip()] if genomes_str else []
                    gene_freq[gene] = {
                        'count': count,
                        'percentage': round(pct, 2),
                        'frequency_display': f"{count} ({pct:.1f}%)",
                        'genomes': genomes,
                        'risk_level': row.get('Risk Level', 'Standard'),
                        'database': 'amrfinder',
                    }
            print(f"    Parsed {len(genes_by_genome)} samples, {len(gene_freq)} genes")
            return genes_by_genome, gene_freq
        except Exception as e:
            print(f"    Error parsing AMRfinder: {e}")
            return {}, {}

    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict, Dict]:
        """Parse one ABRicate summary HTML file."""
        print(f"  Parsing ABRicate: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}, {}

            db_name = 'unknown'
            filename = str(file_path.name).lower()
            for key, value in self.db_name_mapping.items():
                if key in filename:
                    db_name = value
                    break

            df_genomes = self.parse_html_table(str(tables[0]), 0)
            genes_by_genome = {}
            if not df_genomes.empty:
                sample_col = next(
                    (c for c in df_genomes.columns
                     if 'genome' in c.lower() or 'sample' in c.lower()),
                    df_genomes.columns[0],
                )
                for _, row in df_genomes.iterrows():
                    sample_raw = row[sample_col]
                    if not sample_raw:
                        continue
                    sample = self.normalize_sample_id(sample_raw)
                    genes = []
                    for col in df_genomes.columns:
                        if 'genes' in col.lower() or 'detected' in col.lower():
                            gene_str = str(row[col])
                            genes = [g.strip() for g in gene_str.split(',') if g.strip()]
                            break
                    genes_by_genome[sample] = genes

            df_freq = self.parse_html_table(str(tables[1]), 0)
            gene_freq = {}
            if not df_freq.empty:
                for _, row in df_freq.iterrows():
                    gene_full = str(row.get('Gene', '')).strip()
                    if not gene_full:
                        continue
                    gene = re.sub(r'^\([^)]+\)', '', gene_full).strip() or gene_full
                    freq_str = row.get('Frequency', '0')
                    m = re.search(r'(\d+)', freq_str)
                    count = int(m.group(1)) if m else 0
                    pct = (count / total_samples * 100) if total_samples else 0
                    genomes_str = row.get('Genomes', '')
                    genomes = [self.normalize_sample_id(g.strip())
                               for g in genomes_str.split(',') if g.strip()] if genomes_str else []
                    gene_freq[gene] = {
                        'count': count,
                        'percentage': round(pct, 2),
                        'frequency_display': f"{count} ({pct:.1f}%)",
                        'genomes': genomes,
                        'database': db_name,
                        'full_name': gene_full,
                    }
            print(f"    {db_name.upper()}: {len(genes_by_genome)} samples, {len(gene_freq)} genes")
            return genes_by_genome, gene_freq
        except Exception as e:
            print(f"    Error parsing ABRicate report: {e}")
            return {}, {}

    def parse_mutation_summary_html(self, file_path: Path) -> Dict[str, Any]:
        """Parse the mutation summary HTML."""
        print(f"  Parsing mutations: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            mutation_table = None
            for table in soup.find_all('table'):
                if (table.find(string=re.compile(r'Gene', re.I))
                        and table.find(string=re.compile(r'Mutation', re.I))):
                    mutation_table = table
                    break
            if not mutation_table:
                return {}
            header_row = None
            thead = mutation_table.find('thead')
            if thead:
                header_row = thead.find('tr')
            if not header_row:
                header_row = mutation_table.find('tr')
            if not header_row:
                return {}
            headers = [cell.get_text().strip() for cell in header_row.find_all(['th', 'td'])]
            col_idx = {}
            for idx, h in enumerate(headers):
                hl = h.lower()
                if 'gene' in hl: col_idx['gene'] = idx
                elif 'mutation' in hl: col_idx['mutation'] = idx
                elif 'count' in hl: col_idx['count'] = idx
                elif 'genome' in hl: col_idx['genomes'] = idx
                elif 'class' in hl: col_idx['class'] = idx
                elif 'subclass' in hl: col_idx['subclass'] = idx
            for req in ['gene', 'mutation', 'count', 'genomes']:
                if req not in col_idx:
                    return {}

            tbody = mutation_table.find('tbody')
            rows = tbody.find_all('tr') if tbody else mutation_table.find_all('tr')[1:]
            mutations_list = []
            genome_counts = defaultdict(int)
            for row in rows:
                cells = row.find_all('td')
                if len(cells) <= max(col_idx.values()):
                    continue
                gene = cells[col_idx['gene']].get_text().strip()
                mutation = cells[col_idx['mutation']].get_text().strip()
                count_str = cells[col_idx['count']].get_text().strip()
                cm = re.search(r'(\d+)', count_str)
                count = int(cm.group(1)) if cm else 0
                genomes_str = cells[col_idx['genomes']].get_text().strip()
                genomes = [g.strip() for g in genomes_str.split(',') if g.strip()]
                if not genomes:
                    continue
                for g in genomes:
                    genome_counts[g] += 1
                mutations_list.append({
                    'gene': gene,
                    'mutation': mutation,
                    'class': cells[col_idx['class']].get_text().strip() if 'class' in col_idx else '',
                    'subclass': cells[col_idx['subclass']].get_text().strip() if 'subclass' in col_idx else '',
                    'count': count,
                    'genomes': genomes,
                })
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            print(f"    Parsed {len(mutations_list)} mutations across {len(genome_counts)} genomes")
            return {'mutations': mutations_list, 'genome_mutation_counts': dict(genome_counts)}
        except Exception as e:
            print(f"    Error parsing mutations: {e}")
            return {}

    def parse_mge_summary_html(self, file_path: Path) -> Dict[str, Any]:
        """Parse the MGE summary HTML."""
        print(f"  Parsing MGE summary: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            per_sample = {}
            mge_table = None
            for table in soup.find_all('table'):
                header_text = ' '.join(th.get_text().strip() for th in table.find_all('th'))
                if 'Sample' in header_text and 'mobileOG' in header_text:
                    mge_table = table
                    break
            if mge_table:
                df = pd.read_html(str(mge_table))[0]
                df.columns = [str(c).strip() for c in df.columns]
                for _, row in df.iterrows():
                    sample = self.normalize_sample_id(row.get('Sample', ''))
                    if not sample:
                        continue
                    entry = {}
                    for col in df.columns:
                        if col == 'Sample':
                            continue
                        v = row.get(col, 0)
                        try:
                            entry[col] = int(str(v).replace(',', '').strip())
                        except Exception:
                            entry[col] = v
                    per_sample[sample] = entry
            note = ''
            for tag in soup.find_all(string=re.compile(r'Interpretation note', re.I)):
                parent = tag.find_parent()
                if parent:
                    raw = parent.get_text(" ", strip=True)
                    raw = re.sub(r'^Interpretation note\s*:?\s*', '', raw, flags=re.I)
                    if len(raw) > 40:
                        note = raw
                        break
            print(f"    Parsed {len(per_sample)} samples")
            return {'per_sample': per_sample, 'note': note}
        except Exception as e:
            print(f"    Error parsing MGE summary: {e}")
            return {'per_sample': {}, 'note': ''}


# =============================================================================
# ANALYZER
# =============================================================================
class KleboDataAnalyzer:
    """Aggregate gene frequencies, typing patterns, MGE stats, and compare payloads."""

    def __init__(self):
        self.critical_resistance_genes = {
            'blaKPC', 'blaNDM', 'blaIMP', 'blaVIM', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaGES', 'blaIMI', 'blaSME', 'blaDHA', 'blaCMY', 'blaACT',
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9',
            'tet(X)', 'tet(X1)', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)',
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'npmA',
            'fosA3', 'fosA4', 'fosA5', 'fosA6', 'fosA7', 'fosC2',
        }
        self.high_risk_virulence_genes = {
            'ybtA', 'ybtE', 'ybtP', 'ybtQ', 'ybtS', 'ybtT', 'ybtU', 'ybtX', 'irp1', 'irp2',
            'clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH',
            'clbI', 'clbJ', 'clbK', 'clbL', 'clbM', 'clbN', 'clbO', 'clbP', 'clbQ', 'clbR', 'clbS',
            'iroB', 'iroC', 'iroD', 'iroE', 'iroN',
            'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
            'rmpA', 'rmpA2', 'rmpC', 'rmpD',
            'peg-344', 'allS', 'kfuA', 'kfuB', 'kfuC',
            'fimA', 'fimB', 'fimC', 'fimD', 'fimE', 'fimF', 'fimG', 'fimH',
            'mrkA', 'mrkB', 'mrkC', 'mrkD', 'mrkE', 'mrkF', 'mrkH', 'mrkI',
            'ecpA', 'ecpB', 'ecpC', 'ecpD', 'ecpE', 'ecpR',
            'tssA', 'tssB', 'tssC', 'tssD', 'tssE', 'tssF', 'tssG', 'tssH', 'tssI', 'tssJ',
            'tssK', 'tssL', 'tssM', 'tssN', 'tssO', 'tssP', 'tssQ', 'tssR', 'tssS', 'tssT',
            'hcp', 'vgrG', 'icmF', 'impA', 'impB', 'impC', 'impD', 'impE', 'impF', 'impG', 'impH',
            'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS',
            'fepA', 'fepB', 'fepC', 'fepD', 'fepG', 'fecA', 'fecB', 'fecC', 'fecD', 'fecE',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH',
            'hlyA', 'hlyB', 'hlyC', 'hlyD', 'cnf1', 'cnf2', 'sat', 'pic', 'vat',
            'cdtA', 'cdtB', 'cdtC', 'cfa', 'cfaA', 'cfaB', 'cfaC', 'cfaD', 'cfaE',
        }
        self.beta_lactamase_genes = {
            'blaTEM', 'blaSHV', 'blaCTX-M', 'blaOXA-1', 'blaOXA-2', 'blaOXA-9', 'blaOXA-10',
            'blaCMY-2', 'blaDHA-1', 'blaACC', 'blaMIR', 'blaACT', 'blaFOX',
        }
        self.icekp_markers = {'ybt', 'clb', 'iro', 'rmp'}
        self.virulence_plasmid_markers = {'iro', 'iuc', 'rmp', 'rmpA2'}

    def categorize_gene(self, gene: str) -> str:
        g = gene.lower()
        if any(crit.lower() in g for crit in self.critical_resistance_genes):
            return 'Critical Resistance'
        if any(vir.lower() in g for vir in self.high_risk_virulence_genes):
            return 'High-Risk Virulence'
        if any(bla.lower() in g for bla in self.beta_lactamase_genes):
            return 'Beta-Lactamase'
        return 'Other'

    def create_gene_centric_tables(self, gene_freqs: Dict, total_samples: int) -> Dict:
        gene_centric = {
            'all_genes': [],
            'by_category': defaultdict(list),
            'database_stats': {},
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'bacmet_databases': {},
        }
        all_genes = {}
        for db, genes in gene_freqs.items():
            if db == 'bacmet2':
                target = gene_centric['bacmet_databases']
            elif db in ('vfdb', 'ecoli_vf'):
                target = gene_centric['virulence_databases']
            elif db == 'plasmidfinder':
                target = gene_centric['plasmid_databases']
            else:
                target = gene_centric['amr_databases']

            gene_list = []
            for gene, data in genes.items():
                if gene not in all_genes:
                    all_genes[gene] = {
                        'count': 0, 'percentage': 0, 'frequency_display': '',
                        'genomes': set(), 'databases': set(), 'risk_level': 'Standard',
                    }
                all_genes[gene]['count'] += data['count']
                all_genes[gene]['percentage'] += data.get('percentage', 0)
                all_genes[gene]['frequency_display'] = data.get(
                    'frequency_display',
                    f"{data['count']} ({data['count']/total_samples*100:.1f}%)",
                )
                all_genes[gene]['genomes'].update(data['genomes'])
                all_genes[gene]['databases'].add(data.get('database', db))
                if data.get('risk_level') and data['risk_level'] != 'Standard':
                    all_genes[gene]['risk_level'] = data['risk_level']

                gene_list.append({
                    'gene': gene,
                    'database': db.upper(),
                    'frequency_display': data.get(
                        'frequency_display',
                        f"{data['count']} ({data.get('percentage', 0):.1f}%)",
                    ),
                    'count': data['count'],
                    'percentage': data.get('percentage', 0),
                    'genomes': list(data['genomes']),
                    'risk_level': data.get('risk_level', 'Standard'),
                })
            if gene_list:
                gene_list.sort(key=lambda x: x['count'], reverse=True)
                target[db] = gene_list

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
            'virulence_plasmid_marker_presence': defaultdict(list),
            'gene_cooccurrence': defaultdict(Counter),
        }
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
                patterns['st_k_combinations'][f"ST{st}-{k}"].append(sample)
            if st != 'ND' and o != 'ND':
                patterns['st_o_combinations'][f"ST{st}-{o}"].append(sample)
            if k != 'ND' and o != 'ND':
                patterns['ko_combinations'][f"{k}:{o}"].append(sample)
            if st != 'ND' and k != 'ND' and o != 'ND':
                patterns['st_ko_combinations'][f"ST{st}-{k}:{o}"].append(sample)

            genes = sample_genes.get(sample, set())
            crit = [g for g in genes if self.categorize_gene(g) == 'Critical Resistance']
            hv = [g for g in genes if self.categorize_gene(g) == 'High-Risk Virulence']
            if crit and hv:
                patterns['high_risk_combinations'].append({
                    'sample': sample, 'st': st, 'k_locus': k,
                    'critical_resistance': crit, 'high_risk_virulence': hv,
                })
            icekp = [g for g in genes if any(m in g.lower() for m in self.icekp_markers)]
            if icekp:
                patterns['icekp_marker_presence'][sample] = icekp
            vp = [g for g in genes if any(m in g.lower() for m in self.virulence_plasmid_markers)]
            if vp:
                patterns['virulence_plasmid_marker_presence'][sample] = vp

            genes_list = list(genes)
            for i, g1 in enumerate(genes_list):
                for g2 in genes_list[i + 1:]:
                    patterns['gene_cooccurrence'][g1][g2] += 1
        return patterns

    def compute_mge_stats(self, mge_data: Dict) -> Dict:
        per_sample = mge_data.get('per_sample', {})
        if not per_sample:
            return {'category_totals': [], 'grand_total': 0, 'n_samples': 0}
        all_cols = set()
        for s in per_sample.values():
            all_cols.update(s.keys())
        numeric_cols = [c for c in all_cols
                        if all(isinstance(s.get(c, 0), (int, float))
                               for s in per_sample.values())]
        category_totals = []
        n = len(per_sample)
        grand_total = 0
        for col in numeric_cols:
            vals = {s: per_sample[s].get(col, 0) for s in per_sample}
            total = sum(vals.values())
            grand_total += total
            max_sample = max(vals, key=vals.get) if vals else 'ND'
            category_totals.append({
                'category': col,
                'total': total,
                'mean': round(total / n, 2) if n else 0,
                'max': vals[max_sample] if vals else 0,
                'max_sample': max_sample,
            })
        category_totals.sort(key=lambda x: x['total'], reverse=True)
        for c in category_totals:
            c['pct_of_all'] = round((c['total'] / grand_total) * 100, 2) if grand_total else 0
        return {'category_totals': category_totals, 'grand_total': grand_total, 'n_samples': n}

    def build_compare_payload(self, samples_data: Dict, mutation_data: Dict) -> Tuple[Dict, Dict]:
        sample_mutations = defaultdict(list)
        for m in mutation_data.get('mutations', []):
            label = f"{m.get('gene', '')} {m.get('mutation', '')}".strip()
            for genome in m.get('genomes', []):
                sample_mutations[genome].append(label)

        compare_data = {}
        for sample, sd in samples_data.items():
            t = sd.get('mlst', {})
            k = sd.get('kaptive', {})
            ab = sd.get('abricate_genes', {})
            amr_genes = sorted(set(
                sd.get('amr_genes', [])
                + ab.get('card', [])
                + ab.get('resfinder', [])
                + ab.get('argannot', [])
                + ab.get('megares', [])
                + ab.get('ncbi', [])
            ))
            compare_data[sample] = {
                'typing': {
                    'ST': t.get('ST', 'ND'),
                    'K': k.get('K_Locus', 'ND'),
                    'O': k.get('O_Locus', 'ND'),
                },
                'genes': {
                    'AMR': amr_genes,
                    'Virulence': sorted(set(ab.get('vfdb', []) + ab.get('ecoli_vf', []))),
                    'BACMET': sorted(set(ab.get('bacmet2', []))),
                    'Plasmids': sorted(set(ab.get('plasmidfinder', []))),
                    'Mutations': sorted(set(sample_mutations.get(sample, []))),
                },
            }

        samples = sorted(compare_data.keys())
        n = len(samples)
        if n < 2:
            return compare_data, {'samples': samples, 'flat': [], 'n': n}

        cats = ('AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations')
        gene_sets = {s: {c: set(compare_data[s]['genes'][c]) for c in cats} for s in samples}
        typing_data = {s: compare_data[s]['typing'] for s in samples}
        typing_keys = list(typing_data[samples[0]].keys()) if samples else []
        nk = len(typing_keys)

        flat = []
        for i in range(n):
            for j in range(i + 1, n):
                s1, s2 = samples[i], samples[j]
                matches = sum(1 for key in typing_keys if typing_data[s1].get(key) == typing_data[s2].get(key))
                typing_pct = (matches / nk * 100) if nk else 0
                cat_scores = []
                for cat in cats:
                    sa, sb = gene_sets[s1][cat], gene_sets[s2][cat]
                    union = sa | sb
                    if not union:
                        continue
                    cat_scores.append(len(sa & sb) / len(union) * 100)
                gene_pct = (sum(cat_scores) / len(cat_scores)) if cat_scores else 0
                flat.append(round(typing_pct * 0.5 + gene_pct * 0.5, 1))
        return compare_data, {'samples': samples, 'flat': flat, 'n': n}


# =============================================================================
# HTML GENERATOR -  architecture, Kleboscope blue background, per-tab colors
# =============================================================================
class KleboHTMLGenerator:
    """Render the full interactive report with per-tab colors on blue gradient."""

    TAB_COLORS = {
        'summary':         '#4CAF50',
        'sample_overview': '#2196F3',
        'mlst':            '#FF9800',
        'qc':              '#607D8B',
        'kaptive':         '#9C27B0',
        'combinations':    '#009688',
        'amr':             '#F44336',
        'virulence':       '#E91E63',
        'bacmet':          '#FF5722',
        'plasmids':        '#673AB7',
        'mutations':       '#00BCD4',
        'patterns':        '#3F51B5',
        'highrisk':        '#DC143C',
        'mge':             '#16A085',
        'compare':         '#0d9488',
        'databases':       '#795548',
        'aiguide':         '#00BCD4',
        'citation':        '#FFC107',
        'funding':         '#FFB300',
        'export':          '#9E9E9E',
    }

    def __init__(self, analyzer: KleboDataAnalyzer):
        self.analyzer = analyzer
        self._samples_data = {}

    def _credit_bar(self, icon: str, title: str, body: str) -> str:
        return f'''
        <div class="scientific-note">
            <div style="font-size: 1.4em;">{icon}</div>
            <div>
                <strong style="color: var(--primary-dark); font-size: 1.05em;">{title}</strong><br>
                <span style="font-size: 0.92em; color: #4b5563;">{body}</span>
            </div>
        </div>'''

    def _alert(self, kind: str, icon: str, body: str) -> str:
        return f'''
        <div class="alert-box alert-{kind}">
            <i class="fas {icon} fa-2x"></i>
            <div>{body}</div>
        </div>'''

    def _section_header(self, title: str, icon: str, tab_id: str) -> str:
        return (
            f'<div class="section-header {tab_id}-header">'
            f'<h2><i class="fas {icon}"></i> {title}</h2>'
            f'<button class="print-section-btn" onclick="printSection(\'{tab_id}-tab\')">'
            f'<i class="fas fa-print"></i> Print</button></div>'
        )

    def _search_pair(self, table_id: str, placeholder: str) -> str:
        return (
            f'<input type="text" class="search-box" id="search-{table_id}" '
            f'onkeyup="searchTable(\'{table_id}\',\'search-{table_id}\')" '
            f'placeholder="🔍 {placeholder}">'
            f'<input type="text" class="search-box" id="highlight-{table_id}" '
            f'onkeyup="highlightGenome(\'{table_id}\',\'highlight-{table_id}\')" '
            f'placeholder="🔍 Highlight genomes containing text">'
        )

    def _grouping_controls(self, table_id: str) -> str:
        groups = [('ST', 'ST'), ('K', 'K-locus'), ('O', 'O-locus'),
                  ('ST-K', 'ST-K'), ('ST-O', 'ST-O'), ('ST-K:O', 'ST-K:O')]
        buttons = ''.join(
            f'<button class="group-btn" data-group="{g}" '
            f'onclick="groupGenomesByTyping(\'{table_id}\',\'{g}\')">{label}</button>'
            for g, label in groups
        )
        return (
            f'<div class="grouping-controls">'
            f'<strong><i class="fas fa-layer-group"></i> Group genomes by:</strong>'
            f'{buttons}'
            f'<button class="group-btn" onclick="resetGenomeList(\'{table_id}\')">Reset</button>'
            f'</div>'
        )

    def _action_buttons(self, table_id: str, extra: List[Tuple[str, str, str, str]] = None) -> str:
        html = (
            f'<div class="action-buttons">'
            f'<button class="action-btn btn-primary" onclick="exportTableToCSV(\'{table_id}\', \'{table_id}.csv\')">'
            f'<i class="fas fa-download"></i> Export</button>'
        )
        if extra:
            for label, search_val, css, icon in extra:
                html += (
                    f'<button class="action-btn {css}" '
                    f'onclick="document.getElementById(\'search-{table_id}\').value=\'{search_val}\';'
                    f'searchTable(\'{table_id}\',\'search-{table_id}\')">'
                    f'<i class="fas {icon}"></i> {label}</button>'
                )
        html += (
            f'<button class="action-btn btn-light" '
            f'onclick="document.getElementById(\'search-{table_id}\').value=\'\';'
            f'searchTable(\'{table_id}\',\'search-{table_id}\')">'
            f'<i class="fas fa-sync"></i> Clear</button></div>'
        )
        return html

    def _gene_table(self, table_id: str, title: str, genes: List[Dict], total_samples: int) -> str:
        rows = []
        for g in genes:
            tags = ''.join(f'<span class="genome-tag">{gid}</span>' for gid in g['genomes'])
            pct = (g['count'] / total_samples * 100) if total_samples else 0
            freq = f"{g['count']} ({pct:.1f}%)"
            rows.append(
                f'<tr><td><strong>{g["gene"]}</strong></td>'
                f'<td>{g["database"]}</td>'
                f'<td>{freq}</td>'
                f'<td><div class="genome-list">{tags}</div></td></tr>'
            )
        return f'''
        <h3 class="sub-heading">{title}</h3>
        {self._grouping_controls(table_id)}
        {self._search_pair(table_id, f'Search {title.lower()}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">Gene</th>
                    <th data-sort="string">Database</th>
                    <th data-sort="string">Frequency</th>
                    <th data-sort="string">Genomes</th>
                </tr></thead>
                <tbody>{''.join(rows) or '<tr><td colspan="4" style="text-align:center;color:#9ca3af;">No genes</td></tr>'}</tbody>
            </table>
        </div>'''

    def _combo_table(self, table_id: str, title: str, combo_dict: Dict) -> str:
        if not combo_dict:
            return ''
        rows = []
        for combo, samples in sorted(combo_dict.items(), key=lambda x: len(x[1]), reverse=True):
            tags = ''.join(f'<span class="genome-tag">{s}</span>' for s in samples)
            rows.append(
                f'<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td>'
                f'<td><div class="genome-list">{tags}</div></td></tr>'
            )
        return f'''
        <h3 class="sub-heading">{title}</h3>
        {self._search_pair(table_id, f'Search {title.lower()}...')}
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">Combination</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="string">Samples</th>
                </tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>'''

    def _distribution_table(self, table_id: str, dist: Counter, label: str) -> str:
        total = sum(dist.values())
        rows = []
        for val, cnt in dist.most_common():
            if val in ('Not Assigned', 'ND', ''):
                continue
            pct = (cnt / total * 100) if total else 0
            rows.append(f'<tr><td><strong>{val}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td></tr>')
        return f'''
        <div class="master-scrollable-container">
            <table id="{table_id}" class="data-table">
                <thead><tr>
                    <th data-sort="string">{label}</th>
                    <th data-sort="number">Count</th>
                    <th data-sort="number">Percentage</th>
                </tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>'''

    # ------------------------------------------------------------------ CSS
    def _get_css(self) -> str:
        """Kleboscope blue background + per-tab colors."""
        vars_css = '\n'.join(f'    --{tab}-color: {color};' for tab, color in self.TAB_COLORS.items())
        active_css = '\n'.join(
            f'.tab-button.active.{tab} {{ background: var(--{tab}-color); color: white; border-color: var(--{tab}-color); }}'
            for tab in self.TAB_COLORS
        )
        header_css = '\n'.join(
            f'.{tab}-header {{ border-color: var(--{tab}-color); }}'
            for tab in self.TAB_COLORS
        )
        tab_hover_css = '\n'.join(
            f'.tab-button.{tab}:hover:not(.active) {{ color: var(--{tab}-color); }}'
            for tab in self.TAB_COLORS
        )
        dash_card_css = '\n'.join(
            f'.dashboard-card.card-{tab} {{ border-left-color: var(--{tab}-color); }}'
            for tab in self.TAB_COLORS
        )
        return f"""
<style>
:root {{
    --primary: #0077be;
    --primary-dark: #005a8c;
    --primary-light: #3b82f6;
    --success: #10b981;
    --warning: #f59e0b;
    --danger: #dc2626;
    --bg-light: #f8f9fa;
    --text: #1f2937;
{vars_css}
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    color: white; padding: 20px; min-height: 100vh;
}}
.container {{ max-width: 1800px; margin: 0 auto; }}

.main-header {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 100%);
    color: white; padding: 30px; border-radius: 15px;
    box-shadow: 0 10px 30px rgba(0,0,0,0.3); margin-bottom: 24px;
    text-align: center; border: 2px solid rgba(0,119,190,0.3);
}}
.main-header h1 {{ font-size: 2.4em; margin-bottom: 8px; }}
.main-header p {{ opacity: 0.9; margin-bottom: 12px; }}
.metadata-bar {{
    background: rgba(255,255,255,0.1); padding: 12px; border-radius: 10px;
    margin-top: 16px; display: flex; justify-content: center;
    flex-wrap: wrap; gap: 22px; backdrop-filter: blur(10px);
}}
.metadata-item {{ display: flex; align-items: center; gap: 8px; font-size: 0.92em; }}

.dashboard-grid {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 15px; margin-bottom: 24px;
}}
.dashboard-card {{
    background: white; color: #1f2937; padding: 20px;
    border-radius: 12px; text-align: center; cursor: pointer;
    box-shadow: 0 4px 15px rgba(0,0,0,0.15);
    transition: transform 0.2s, box-shadow 0.2s;
    border-left: 5px solid #0077be;
}}
.dashboard-card:hover {{
    transform: translateY(-4px);
    box-shadow: 0 8px 25px rgba(0,0,0,0.25);
}}
.dashboard-card .card-icon {{ font-size: 1.6em; opacity: 0.8; margin-bottom: 6px; }}
.dashboard-card .card-number {{
    font-size: 2.2em; font-weight: bold; margin: 4px 0; color: #1f2937;
}}
.dashboard-card .card-label {{
    font-size: 0.82em; text-transform: uppercase;
    letter-spacing: 0.5px; color: #6b7280;
}}
{dash_card_css}

.tab-navigation {{
    display: flex; gap: 5px; flex-wrap: wrap;
    background: rgba(255,255,255,0.95); padding: 12px; border-radius: 12px;
    box-shadow: 0 5px 20px rgba(0,0,0,0.2); margin-bottom: 20px;
    position: sticky; top: 10px; z-index: 100;
}}
.tab-button {{
    padding: 10px 16px; background: #e5e7eb; border: 2px solid #e5e7eb;
    border-radius: 8px; cursor: pointer; font-weight: 600;
    color: #4b5563; font-size: 0.85em;
    display: flex; align-items: center; gap: 6px;
    transition: background 0.2s, color 0.2s, border-color 0.2s;
}}
{tab_hover_css}
{active_css}

.tab-content {{
    display: none; background: rgba(255,255,255,0.97);
    color: #1f2937; padding: 28px; border-radius: 15px;
    box-shadow: 0 10px 30px rgba(0,0,0,0.2); margin-bottom: 24px;
}}
.tab-content.active {{ display: block; animation: fadeIn 0.3s ease; }}
@keyframes fadeIn {{ from {{ opacity: 0; transform: translateY(10px); }} to {{ opacity: 1; transform: translateY(0); }} }}

.section-header {{
    display: flex; justify-content: space-between; align-items: center;
    border-bottom: 3px solid #0077be;
    padding-bottom: 12px; margin-bottom: 22px;
}}
.section-header h2 {{
    color: #1f2937; font-size: 1.6em;
    display: flex; align-items: center; gap: 10px;
}}
{header_css}
.sub-heading {{
    color: #005a8c; font-size: 1.2em;
    margin: 24px 0 12px 0; padding-bottom: 6px;
    border-bottom: 1px solid #e5e7eb;
}}

.data-table {{
    width: 100%; border-collapse: collapse; font-size: 0.9em;
    background: white; border-radius: 8px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.05); margin: 12px 0;
}}
.data-table th {{
    background: #005a8c; color: white;
    padding: 12px 10px; text-align: left; font-weight: 600;
    position: sticky; top: 0; cursor: pointer; white-space: nowrap;
}}
.data-table th:hover {{ background: #003f63; }}
.data-table td {{
    padding: 10px; border-bottom: 1px solid #e5e7eb;
    vertical-align: top; word-wrap: break-word;
}}
.data-table tr:hover {{ background: #f3f4f6; }}
.master-scrollable-container {{
    width: 100%; overflow-x: auto;
    border: 1px solid #e5e7eb; border-radius: 8px; margin: 14px 0;
}}

.genome-list {{
    display: flex; flex-wrap: wrap; gap: 4px;
    max-height: 220px; overflow-y: auto;
    padding: 6px; background: #f8f9fa; border-radius: 5px;
}}
.genome-tag {{
    display: inline-block; background: #e6f4fb; color: #005a8c;
    padding: 3px 9px; border-radius: 12px; font-size: 0.78em;
    border: 1px solid #b3e0f4; white-space: nowrap;
}}
.genome-tag.highlight {{
    background: #fef08a !important; color: #1f2937 !important;
    border-color: #facc15;
}}
.genome-group {{ width: 100%; margin-bottom: 8px; }}
.genome-group-header {{
    font-weight: bold; background: #d1d5db; padding: 4px 10px;
    border-radius: 4px; margin: 6px 0 4px; font-size: 0.82em;
    display: inline-block;
}}
.genome-group-tags {{ display: flex; flex-wrap: wrap; gap: 4px; margin-left: 8px; }}

.search-box {{
    width: 100%; padding: 10px 14px; margin-bottom: 10px;
    border: 2px solid #e5e7eb; border-radius: 8px; font-size: 0.95em;
}}
.search-box:focus {{
    outline: none; border-color: #0077be;
    box-shadow: 0 0 0 3px rgba(0,119,190,0.15);
}}
.action-buttons {{ display: flex; gap: 8px; margin: 12px 0; flex-wrap: wrap; }}
.action-btn {{
    padding: 8px 16px; border: none; border-radius: 8px;
    cursor: pointer; font-weight: 600; font-size: 0.85em;
    display: flex; align-items: center; gap: 6px;
}}
.action-btn:hover {{ transform: translateY(-1px); }}
.btn-primary {{ background: #0077be; color: white; }}
.btn-success {{ background: #10b981; color: white; }}
.btn-warning {{ background: #f59e0b; color: white; }}
.btn-danger  {{ background: #dc2626; color: white; }}
.btn-info    {{ background: #3b82f6; color: white; }}
.btn-secondary {{ background: #6b7280; color: white; }}
.btn-light {{ background: #f8f9fa; color: #1f2937; border: 1px solid #d1d5db; }}

.grouping-controls {{
    background: #f8f9fa; padding: 12px; border-radius: 8px;
    margin: 12px 0; display: flex; flex-wrap: wrap;
    gap: 8px; align-items: center;
    border-left: 4px solid #0077be;
}}
.group-btn {{
    background: white; border: 1px solid #0077be;
    color: #0077be; padding: 6px 12px;
    border-radius: 20px; cursor: pointer; font-size: 0.82em;
}}
.group-btn:hover, .group-btn.active {{ background: #0077be; color: white; }}

.metrics-grid {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 15px; margin: 16px 0;
}}
.metric-card {{
    background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
    color: white; padding: 18px; border-radius: 10px;
    text-align: center; box-shadow: 0 4px 12px rgba(0,0,0,0.12);
}}
.metric-value {{ font-size: 1.9em; font-weight: bold; }}
.metric-label {{
    font-size: 0.82em; opacity: 0.92;
    text-transform: uppercase; letter-spacing: 0.5px; margin-top: 4px;
}}

.alert-box {{
    padding: 18px 22px; border-radius: 10px; margin: 16px 0;
    display: flex; align-items: flex-start; gap: 16px;
    border-left: 5px solid; font-size: 0.94em;
}}
.alert-info    {{ background: #e0f2fe; color: #0c4a6e; border-color: #0077be; }}
.alert-success {{ background: #d1fae5; color: #065f46; border-color: #10b981; }}
.alert-warning {{ background: #fef3c7; color: #92400e; border-color: #f59e0b; }}
.alert-danger  {{ background: #fee2e2; color: #7f1d1d; border-color: #dc2626; }}

.scientific-note {{
    background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%);
    border-left: 5px solid #0077be;
    padding: 15px 20px; border-radius: 8px; margin-bottom: 18px;
    display: flex; align-items: flex-start; gap: 14px;
    color: #1f2937; font-size: 0.92em;
}}

.print-section-btn {{
    background: #0077be; color: white; border: none;
    border-radius: 6px; padding: 7px 14px; cursor: pointer;
    font-size: 0.85em; display: flex; align-items: center; gap: 5px;
}}
.print-section-btn:hover {{ background: #005a8c; }}

.badge {{
    display: inline-block; padding: 3px 10px; border-radius: 12px;
    font-size: 0.8em; font-weight: 600; margin: 2px;
}}
.badge-danger  {{ background: #dc2626; color: white; }}
.badge-warning {{ background: #f59e0b; color: white; }}
.badge-info    {{ background: #0077be; color: white; }}

.vir-details {{ margin-top: 4px; }}
.vir-details summary {{
    cursor: pointer; color: #0077be; font-weight: 600;
    font-size: 0.85em; list-style: none;
}}
.vir-details summary::-webkit-details-marker {{ display: none; }}
.vir-details summary::before {{ content: '▸ '; }}
.vir-details[open] summary::before {{ content: '▾ '; }}
.vir-details .vir-gene-tags {{
    margin-top: 4px; display: flex; flex-wrap: wrap; gap: 3px;
}}

/* -style white/black citations */
.citation-item {{
    background: #fafbfc; border-radius: 6px; margin-bottom: 10px;
    padding: 12px 14px; list-style: none;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    border-left: 4px solid #0077be;
    transition: box-shadow 0.2s, transform 0.2s;
}}
.citation-item:hover {{
    box-shadow: 0 4px 12px rgba(0,0,0,0.10);
    transform: translateX(2px);
}}
.citation-body {{
    display: flex; flex-direction: column; gap: 8px;
    font-size: 0.92em; line-height: 1.55;
}}
.citation-name {{ font-size: 1em; font-weight: 700; }}
.citation-actions {{
    display: flex; gap: 8px; flex-wrap: wrap; align-items: center;
}}
.citation-link {{
    display: inline-flex; align-items: center; gap: 4px;
    padding: 4px 14px; border-radius: 16px; font-size: 0.82em;
    font-weight: 600; color: white; text-decoration: none;
    background: #0077be;
}}
.citation-link:hover {{ opacity: 0.88; }}
.citation-actions .copy-btn {{
    background: #6b7280; color: white; border: none;
    padding: 4px 14px; border-radius: 16px; cursor: pointer;
    font-size: 0.82em; font-weight: 600;
}}
.citation-actions .copy-btn:hover {{ background: #4b5563; }}

.accordion {{ margin: 20px 0; }}
.accordion-item {{
    background: #f8f9fa; border: 1px solid #dee2e6;
    margin-bottom: 10px; border-radius: 8px; overflow: hidden;
}}
.accordion-header {{
    background: #e9ecef; padding: 12px 20px; cursor: pointer;
    font-weight: bold; color: #1e3a8a;
    display: flex; justify-content: space-between; align-items: center;
}}
.accordion-header:hover {{ background: #dee2e6; }}
.accordion-content {{
    padding: 15px 20px; border-top: 1px solid #dee2e6; background: white;
}}

/* Footer - original Kleboscope */
.footer {{
    text-align: center; padding: 30px; color: white;
    margin-top: 40px; border-radius: 15px;
    background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
}}
.footer h3 {{ margin-bottom: 8px; }}
.footer p {{ margin: 6px 0; font-size: 0.9em; }}
.footer a {{ color: #fbbf24; text-decoration: none; }}
.footer a:hover {{ text-decoration: underline; }}

/* Compare tab */
.mode-toggle {{ display: flex; gap: 8px; margin-bottom: 18px; }}
.mode-btn {{
    padding: 10px 20px; border: 2px solid #cbd5e1; background: white;
    color: #475569; border-radius: 8px; font-weight: 600;
    cursor: pointer; font-size: 0.9em;
    display: flex; align-items: center; gap: 8px;
}}
.mode-btn:hover {{ border-color: #0077be; color: #0077be; }}
.mode-btn.active {{ background: #0077be; border-color: #0077be; color: white; }}
.compare-mode {{ display: none; }}
.compare-mode.active {{ display: block; }}
.compare-panel {{
    background: #f8f9fa; border-radius: 12px;
    padding: 20px; margin-bottom: 20px;
}}
.compare-select-row {{
    display: flex; gap: 12px; align-items: flex-end;
    flex-wrap: wrap; margin-bottom: 16px;
}}
.compare-select-row label {{ display: flex; flex-direction: column; gap: 4px; }}
.compare-select-row label span {{
    font-size: 0.82em; font-weight: 600; color: #475569;
}}
.compare-select-row select {{
    padding: 10px 14px; border: 2px solid #cbd5e1; border-radius: 8px;
    background: white; font-size: 0.92em; min-width: 220px;
}}
.compare-visuals-grid {{
    display: grid; grid-template-columns: 200px 1fr; gap: 18px;
    margin: 16px 0; align-items: center;
}}
.visual-card {{
    background: white; border-radius: 12px; padding: 18px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06); border: 1px solid #e2e8f0;
}}
.similarity-gauge {{ width: 180px; height: 180px; }}
.metrics-card {{ display: flex; flex-direction: column; gap: 16px; }}
.metric-block {{ display: flex; flex-direction: column; gap: 6px; }}
.metric-track {{
    background: #f1f5f9; height: 22px; border-radius: 11px;
    overflow: hidden; position: relative;
}}
.metric-fill {{ height: 100%; border-radius: 11px; transition: width 0.8s ease; min-width: 4px; }}
.compare-section {{
    background: white; border-radius: 12px; padding: 20px;
    margin-bottom: 18px; box-shadow: 0 2px 8px rgba(0,0,0,0.06);
}}
.compare-row {{
    display: grid; grid-template-columns: 160px 1fr 1fr;
    gap: 12px; padding: 10px 14px; border-radius: 6px;
    align-items: center; font-size: 0.9em;
}}
.compare-row.header {{
    background: #f1f5f9; font-weight: 700;
    border-bottom: 2px solid #cbd5e1;
}}
.compare-row.diff {{ background: #fef2f2; border-left: 3px solid #dc2626; }}
.compare-row.match {{ background: #f0fdf4; border-left: 3px solid #10b981; }}
.compare-row .row-label {{
    font-weight: 600; color: #475569; font-size: 0.82em;
    text-transform: uppercase;
}}
.typing-pill {{
    display: inline-block; padding: 4px 12px; border-radius: 12px;
    font-size: 0.82em; font-weight: 600; border: 1.5px solid;
    font-family: monospace;
}}
.gene-proportion-bar {{
    display: flex; height: 12px; border-radius: 6px;
    overflow: hidden; margin-bottom: 12px; background: #f1f5f9;
}}
.gene-prop-segment {{ height: 100%; min-width: 2px; }}
.gene-prop-segment.shared {{ background: #10b981; }}
.gene-prop-segment.only-a {{ background: #3b82f6; }}
.gene-prop-segment.only-b {{ background: #ec4899; }}
.cluster-card {{
    background: white; border-radius: 10px; padding: 14px 18px;
    margin-bottom: 10px; border-left: 5px solid #0077be;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
}}
.cluster-header {{
    display: flex; align-items: center; gap: 12px;
    margin-bottom: 10px; flex-wrap: wrap;
}}
.cluster-badge {{
    color: white; padding: 3px 12px; border-radius: 12px;
    font-size: 0.82em; font-weight: 700; background: #0077be;
}}
.cluster-stats {{ font-size: 0.85em; color: #64748b; }}
.cluster-samples {{ display: flex; flex-wrap: wrap; gap: 5px; }}
.cluster-sample-pill {{
    background: #e6f4fb; color: #005a8c;
    padding: 3px 10px; border-radius: 10px; font-size: 0.82em;
    font-family: monospace; border: 1px solid #b3e0f4;
}}
.heatmap-wrapper {{
    background: white; border-radius: 12px; padding: 18px;
    margin-top: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.06);
}}
.heatmap-grid {{ display: grid; gap: 2px; margin-top: 12px; }}
.heatmap-cell {{
    width: 28px; height: 28px; display: flex;
    align-items: center; justify-content: center;
    font-size: 0.68em; font-weight: 700; border-radius: 3px;
    overflow: hidden;
}}
.heatmap-row-header, .heatmap-col-header {{
    background: transparent; color: #64748b; font-size: 0.72em;
}}
.heatmap-row-header {{ justify-content: flex-start; padding-left: 6px; }}
.heatmap-col-header {{
    writing-mode: vertical-rl; text-orientation: mixed;
    height: 100px; width: 28px; padding: 4px 0;
}}

/* CHANGED: .dashboard-grid now has its own closing brace on the same line,
   and the media query is properly closed below */
@media (max-width: 768px) {{
    .compare-row {{ grid-template-columns: 1fr; }}
    .compare-visuals-grid {{ grid-template-columns: 1fr; }}
    .dashboard-grid {{ grid-template-columns: repeat(2, 1fr); }}
}}

.mutations-table {{
    table-layout: fixed;
}}
.mutations-table td:first-child {{
    word-wrap: break-word;
    overflow-wrap: break-word;
}}
.mutation-subtitle {{
    display: block;
    margin-top: 4px;
    font-size: 0.78em;
    font-weight: normal;
    font-style: italic;
    color: #6b7280;
    line-height: 1.35;
    word-break: break-word;
}}

@media print {{
    body * {{ visibility: hidden; }}
    .tab-content.active, .tab-content.active * {{ visibility: visible; }}
    .tab-content.active {{
        position: absolute; left: 0; top: 0;
        width: 100%; padding: 20px;
        box-shadow: none; border-radius: 0;
    }}
    .print-section-btn, .tab-navigation, .dashboard-grid,
    .search-box, .action-buttons, .grouping-controls {{ display: none !important; }}
}}
</style>"""

    # ------------------------------------------------------------------ JS
    def _get_js(self, typing_json: str) -> str:
        return f"""
<script>
var sampleTyping = {typing_json};
var originalGenomeLists = {{}};

function escapeHtml(s) {{
    return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;')
        .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}}

function switchTab(name) {{
    document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
    const content = document.getElementById(name + '-tab');
    const button = document.querySelector('.tab-button.' + name);
    if (content) content.classList.add('active');
    if (button) button.classList.add('active');
    window.location.hash = name;
    if (name === 'compare') initCompare();
}}

function searchTable(tableId, searchId) {{
    const filter = (document.getElementById(searchId)?.value || '').toUpperCase();
    const table = document.getElementById(tableId);
    if (!table) return;
    const rows = table.tBodies[0].rows;
    for (let i = 0; i < rows.length; i++) {{
        const cells = rows[i].getElementsByTagName('td');
        let found = false;
        for (let j = 0; j < cells.length; j++) {{
            if ((cells[j].textContent || '').toUpperCase().indexOf(filter) > -1) {{
                found = true; break;
            }}
        }}
        rows[i].style.display = found ? '' : 'none';
    }}
}}

function highlightGenome(tableId, searchId) {{
    const filter = (document.getElementById(searchId)?.value || '').toUpperCase().trim();
    const table = document.getElementById(tableId);
    if (!table) return;
    table.querySelectorAll('.genome-tag').forEach(tag => {{
        tag.classList.remove('highlight');
        if (filter && tag.textContent.toUpperCase().indexOf(filter) > -1) {{
            tag.classList.add('highlight');
        }}
    }});
}}

function getTypingValue(genome, groupBy) {{
    const info = sampleTyping[genome];
    if (!info) return 'Unknown';
    const st = info.ST || 'ND';
    const k = info.K || 'ND';
    const o = info.O || 'ND';
    switch (groupBy) {{
        case 'ST': return 'ST' + st;
        case 'K': return k;
        case 'O': return o;
        case 'ST-K': return 'ST' + st + '-' + k;
        case 'ST-O': return 'ST' + st + '-' + o;
        case 'ST-K:O': return 'ST' + st + '-' + k + ':' + o;
    }}
    return 'Unknown';
}}

function groupGenomesByTyping(tableId, groupBy) {{
    const table = document.getElementById(tableId);
    if (!table) return;
    const tbody = table.tBodies[0];
    if (!tbody) return;
    const rows = tbody.rows;
    for (let i = 0; i < rows.length; i++) {{
        const row = rows[i];
        const nameCell = row.cells[0];
        if (!nameCell) continue;
        const name = nameCell.textContent.trim().replace(/⚠️/g, '').trim();
        if (!originalGenomeLists[name]) {{
            let genomesCell = null;
            for (let j = 0; j < row.cells.length; j++) {{
                if (row.cells[j].querySelector('.genome-list')) {{
                    genomesCell = row.cells[j]; break;
                }}
            }}
            if (genomesCell) {{
                const tags = genomesCell.querySelectorAll('.genome-tag');
                originalGenomeLists[name] = Array.from(tags).map(t => t.textContent.trim());
            }} else {{
                originalGenomeLists[name] = [];
            }}
        }}
    }}
    for (let i = 0; i < rows.length; i++) {{
        const row = rows[i];
        const nameCell = row.cells[0];
        if (!nameCell) continue;
        const name = nameCell.textContent.trim().replace(/⚠️/g, '').trim();
        const genomes = originalGenomeLists[name] || [];
        let genomesCell = null;
        for (let j = 0; j < row.cells.length; j++) {{
            if (row.cells[j].querySelector('.genome-list')) {{
                genomesCell = row.cells[j]; break;
            }}
        }}
        if (!genomesCell) continue;
        if (genomes.length === 0) {{
            genomesCell.innerHTML = '<div class="genome-list">None</div>';
            continue;
        }}
        const groups = {{}};
        genomes.forEach(g => {{
            const key = getTypingValue(g, groupBy);
            if (!groups[key]) groups[key] = [];
            groups[key].push(g);
        }});
        let html = '<div class="genome-list">';
        for (const key in groups) {{
            const tags = groups[key].map(g => '<span class="genome-tag">' + g + '</span>').join('');
            html += '<div class="genome-group"><div class="genome-group-header">' + key +
                    '</div><div class="genome-group-tags">' + tags + '</div></div>';
        }}
        html += '</div>';
        genomesCell.innerHTML = html;
    }}
}}

function resetGenomeList(tableId) {{
    const table = document.getElementById(tableId);
    if (!table) return;
    const tbody = table.tBodies[0];
    if (!tbody) return;
    for (let i = 0; i < tbody.rows.length; i++) {{
        const row = tbody.rows[i];
        const nameCell = row.cells[0];
        if (!nameCell) continue;
        const name = nameCell.textContent.trim().replace(/⚠️/g, '').trim();
        const original = originalGenomeLists[name] || [];
        for (let j = 0; j < row.cells.length; j++) {{
            if (row.cells[j].querySelector('.genome-list')) {{
                const tags = original.map(g => '<span class="genome-tag">' + g + '</span>').join('');
                row.cells[j].innerHTML = '<div class="genome-list">' + tags + '</div>';
                break;
            }}
        }}
    }}
}}

function sortTable(tableId, colIdx, type) {{
    const table = document.getElementById(tableId);
    const tbody = table.tBodies[0];
    const rows = Array.from(tbody.rows);
    const asc = table.getAttribute('data-sort-dir') !== 'asc';
    rows.sort((a, b) => {{
        let av = a.cells[colIdx].innerText.trim();
        let bv = b.cells[colIdx].innerText.trim();
        if (type === 'number') {{
            av = parseFloat(av.replace(/,/g, '')) || 0;
            bv = parseFloat(bv.replace(/,/g, '')) || 0;
            return asc ? av - bv : bv - av;
        }}
        return asc ? av.localeCompare(bv) : bv.localeCompare(av);
    }});
    rows.forEach(r => tbody.appendChild(r));
    table.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
}}

function printSection(id) {{
    const content = document.getElementById(id);
    if (!content) return;
    const w = window.open('', '_blank');
    const style = document.querySelector('style');
    w.document.write('<html><head><title>Print</title>');
    if (style) w.document.write('<style>' + style.textContent + '</style>');
    w.document.write('</head><body>' + content.innerHTML + '</body></html>');
    w.document.close();
    w.print();
}}

function exportTableToCSV(tableId, filename) {{
    const table = document.getElementById(tableId);
    if (!table) return;
    const rows = table.querySelectorAll('tr');
    const csv = [];
    rows.forEach(row => {{
        const cols = row.querySelectorAll('td, th');
        const data = Array.from(cols).map(c => '"' + (c.innerText || '').replace(/"/g, '""') + '"');
        csv.push(data.join(','));
    }});
    const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

/* Compare */
function initCompare() {{
    const data = window._COMPARE_DATA || {{}};
    const selA = document.getElementById('compare-a');
    const selB = document.getElementById('compare-b');
    if (!selA || !selB || selA.options.length > 0) return;
    const samples = Object.keys(data).sort();
    samples.forEach((s, i) => {{
        const oA = document.createElement('option');
        oA.value = s; oA.textContent = s;
        if (i === 0) oA.selected = true;
        selA.appendChild(oA);
        const oB = document.createElement('option');
        oB.value = s; oB.textContent = s;
        if (i === 1) oB.selected = true;
        selB.appendChild(oB);
    }});
}}

function setCompareMode(mode) {{
    document.querySelectorAll('.mode-btn').forEach(b =>
        b.classList.toggle('active', b.dataset.mode === mode));
    document.querySelectorAll('.compare-mode').forEach(m => m.classList.remove('active'));
    const el = document.getElementById('mode-' + mode);
    if (el) el.classList.add('active');
    if (mode === 'cluster') runCluster();
}}

function similarityColor(pct) {{
    if (pct >= 95) return '#059669';
    if (pct >= 90) return '#10b981';
    if (pct >= 80) return '#84cc16';
    if (pct >= 70) return '#eab308';
    if (pct >= 50) return '#f59e0b';
    if (pct >= 30) return '#f97316';
    return '#dc2626';
}}

function computeSimilarity(A, B) {{
    const typingKeys = Object.keys(A.typing);
    let matches = 0;
    typingKeys.forEach(k => {{ if ((A.typing[k]||'') === (B.typing[k]||'')) matches++; }});
    const typingPct = typingKeys.length ? (matches / typingKeys.length) * 100 : 0;
    const cats = ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'];
    const catScores = [];
    cats.forEach(cat => {{
        const sa = new Set(A.genes[cat] || []);
        const sb = new Set(B.genes[cat] || []);
        const union = new Set([...sa, ...sb]);
        if (!union.size) return;
        let shared = 0;
        union.forEach(g => {{ if (sa.has(g) && sb.has(g)) shared++; }});
        catScores.push((shared / union.size) * 100);
    }});
    const genePct = catScores.length ? catScores.reduce((a, b) => a + b, 0) / catScores.length : 0;
    return {{
        typing: Math.round(typingPct),
        genes: Math.round(genePct),
        overall: Math.round(typingPct * 0.5 + genePct * 0.5),
        typingMatches: matches,
        typingTotal: typingKeys.length
    }};
}}

function drawGauge(pct) {{
    const color = similarityColor(pct);
    const radius = 70;
    const circ = 2 * Math.PI * radius;
    const offset = circ * (1 - pct / 100);
    return '<svg viewBox="0 0 180 180" class="similarity-gauge">' +
        '<circle cx="90" cy="90" r="' + radius + '" fill="none" stroke="#e2e8f0" stroke-width="14"/>' +
        '<circle cx="90" cy="90" r="' + radius + '" fill="none" stroke="' + color + '" stroke-width="14" ' +
        'stroke-dasharray="' + circ + '" stroke-dashoffset="' + offset + '" stroke-linecap="round" ' +
        'transform="rotate(-90 90 90)"/>' +
        '<text x="90" y="86" text-anchor="middle" font-size="28" font-weight="700" fill="' + color + '">' + pct + '%</text>' +
        '<text x="90" y="108" text-anchor="middle" font-size="11" fill="#64748b">overall</text>' +
        '</svg>';
}}

function typingValuePill(va, vb) {{
    const match = (va === vb);
    const color = match ? '#10b981' : '#dc2626';
    const bg = match ? '#dcfce7' : '#fee2e2';
    const icon = match ? '✓' : '✗';
    return '<span class="typing-pill" style="background:' + bg +
           ';color:' + color + ';border-color:' + color + ';">' +
           icon + ' ' + escapeHtml(va) + '</span>';
}}

function runCompare() {{
    const data = window._COMPARE_DATA || {{}};
    const a = document.getElementById('compare-a').value;
    const b = document.getElementById('compare-b').value;
    if (!a || !b || a === b) return;
    if (!data[a] || !data[b]) return;
    const A = data[a], B = data[b];
    const sim = computeSimilarity(A, B);

    document.getElementById('compare-visuals').innerHTML =
        '<div class="compare-visuals-grid">' +
            '<div class="visual-card" style="display:flex;justify-content:center;align-items:center;">' +
                drawGauge(sim.overall) + '</div>' +
            '<div class="visual-card metrics-card">' +
                '<div class="metric-block">' +
                    '<div style="font-weight:700;font-size:0.82em;color:#475569;text-transform:uppercase;">Typing match</div>' +
                    '<div class="metric-track"><div class="metric-fill" style="width:' + sim.typing +
                    '%;background:' + similarityColor(sim.typing) + ';"></div></div>' +
                    '<div style="font-size:0.85em;color:#334155;font-weight:600;">' +
                    sim.typingMatches + '/' + sim.typingTotal + ' fields (' + sim.typing + '%)</div>' +
                '</div>' +
                '<div class="metric-block">' +
                    '<div style="font-weight:700;font-size:0.82em;color:#475569;text-transform:uppercase;">Gene content match</div>' +
                    '<div class="metric-track"><div class="metric-fill" style="width:' + sim.genes +
                    '%;background:' + similarityColor(sim.genes) + ';"></div></div>' +
                    '<div style="font-size:0.85em;color:#334155;font-weight:600;">' +
                    sim.genes + '% Jaccard similarity</div>' +
                '</div>' +
            '</div>' +
        '</div>';

    const verdict = document.getElementById('compare-verdict');
    let verdictHtml = '';
    if (sim.overall >= 95) {{
        verdictHtml = '<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle fa-2x"></i>' +
            '<div><strong>Near-identical (' + sim.overall + '%)</strong> — possible transmission pair or duplicate.</div></div>';
    }} else if (sim.overall >= 85) {{
        verdictHtml = '<div class="alert-box alert-warning"><i class="fas fa-info-circle fa-2x"></i>' +
            '<div><strong>High similarity (' + sim.overall + '%)</strong> — likely same clonal lineage.</div></div>';
    }} else if (sim.overall >= 60) {{
        verdictHtml = '<div class="alert-box alert-warning"><i class="fas fa-info-circle fa-2x"></i>' +
            '<div><strong>Moderate similarity (' + sim.overall + '%)</strong> — related lineage.</div></div>';
    }} else {{
        verdictHtml = '<div class="alert-box alert-info"><i class="fas fa-info-circle fa-2x"></i>' +
            '<div><strong>Low similarity (' + sim.overall + '%)</strong> — distinct lineages.</div></div>';
    }}
    verdict.innerHTML = verdictHtml;
    verdict.style.display = 'block';

    let typingRows = '';
    Object.keys(A.typing).forEach(k => {{
        const va = A.typing[k] || 'ND';
        const vb = B.typing[k] || 'ND';
        const diff = va !== vb;
        typingRows += '<div class="compare-row ' + (diff ? 'diff' : 'match') + '">' +
            '<span class="row-label">' + k + '</span>' +
            '<span>' + typingValuePill(va, vb) + '</span>' +
            '<span>' + typingValuePill(vb, va) + '</span>' +
            '</div>';
    }});

    let geneHtml = '';
    ['AMR', 'Virulence', 'BACMET', 'Plasmids', 'Mutations'].forEach(cat => {{
        const ga = new Set(A.genes[cat] || []);
        const gb = new Set(B.genes[cat] || []);
        const union = new Set([...ga, ...gb]);
        if (!union.size) return;
        const shared = [], onlyA = [], onlyB = [];
        union.forEach(g => {{
            if (ga.has(g) && gb.has(g)) shared.push(g);
            else if (ga.has(g)) onlyA.push(g);
            else onlyB.push(g);
        }});
        const total = union.size;
        const sharedPct = shared.length / total * 100;
        const aPct = onlyA.length / total * 100;
        const bPct = onlyB.length / total * 100;
        geneHtml += '<div class="compare-section">' +
            '<h3>' + cat + ' <span style="font-size:0.85em;color:#64748b;font-weight:400;">' +
            shared.length + ' shared · ' + onlyA.length + ' only A · ' + onlyB.length + ' only B (' +
            sharedPct.toFixed(0) + '% match)</span></h3>' +
            '<div class="gene-proportion-bar">' +
                '<div class="gene-prop-segment shared" style="width:' + sharedPct + '%;"></div>' +
                '<div class="gene-prop-segment only-a" style="width:' + aPct + '%;"></div>' +
                '<div class="gene-prop-segment only-b" style="width:' + bPct + '%;"></div>' +
            '</div>' +
            '<div style="display:flex;gap:12px;flex-wrap:wrap;">' +
                '<div style="flex:1;min-width:180px;background:#f8fafc;border-radius:8px;padding:10px;border:1px solid #e2e8f0;">' +
                    '<strong style="color:#065f46;">Shared (' + shared.length + ')</strong><br>' +
                    (shared.map(g => '<span class="genome-tag" style="background:#dcfce7;color:#065f46;border-color:#86efac;">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') +
                '</div>' +
                '<div style="flex:1;min-width:180px;background:#f8fafc;border-radius:8px;padding:10px;border:1px solid #e2e8f0;">' +
                    '<strong style="color:#005a8c;">Only A (' + onlyA.length + ')</strong><br>' +
                    (onlyA.map(g => '<span class="genome-tag">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') +
                '</div>' +
                '<div style="flex:1;min-width:180px;background:#f8fafc;border-radius:8px;padding:10px;border:1px solid #e2e8f0;">' +
                    '<strong style="color:#831843;">Only B (' + onlyB.length + ')</strong><br>' +
                    (onlyB.map(g => '<span class="genome-tag" style="background:#fce7f3;color:#831843;border-color:#f9a8d4;">' + escapeHtml(g) + '</span>').join('') || '<em>None</em>') +
                '</div>' +
            '</div></div>';
    }});

    const typingTable = '<div class="compare-section">' +
        '<h3>🧬 Typing Profile</h3>' +
        '<div>' +
        '<div class="compare-row header">' +
            '<span class="row-label">Field</span>' +
            '<span><strong>' + escapeHtml(a) + '</strong></span>' +
            '<span><strong>' + escapeHtml(b) + '</strong></span>' +
        '</div>' + typingRows + '</div></div>';

    document.getElementById('compare-output').innerHTML =
        typingTable + (geneHtml || '<div class="alert-box alert-info">No gene data.</div>');
}}

function swapCompare() {{
    const a = document.getElementById('compare-a');
    const b = document.getElementById('compare-b');
    const t = a.value; a.value = b.value; b.value = t;
    runCompare();
}}

function resetCompare() {{
    const sel = document.getElementById('compare-a');
    if (sel && sel.options.length >= 2) {{
        document.getElementById('compare-a').value = sel.options[0].value;
        document.getElementById('compare-b').value = sel.options[1].value;
    }}
    document.getElementById('compare-output').innerHTML = '';
    document.getElementById('compare-visuals').innerHTML = '';
    document.getElementById('compare-verdict').style.display = 'none';
}}

function getPrecomputedSim(i, j) {{
    if (i === j) return 100;
    if (i > j) {{ const t = i; i = j; j = t; }}
    const sim = window._SIM;
    if (!sim || !sim.flat) return 0;
    const n = sim.n;
    return sim.flat[i * n - i * (i + 1) / 2 + (j - i - 1)];
}}

function runCluster() {{
    const sim = window._SIM;
    if (!sim || !sim.samples || sim.samples.length < 2) {{
        document.getElementById('cluster-clusters').innerHTML =
            '<div class="alert-box alert-warning">Need at least two samples.</div>';
        return;
    }}
    const samples = sim.samples;
    const n = sim.n;
    const threshold = parseInt(document.getElementById('cluster-threshold').value);
    const parent = {{}};
    for (let i = 0; i < n; i++) parent[i] = i;
    function find(x) {{ while (parent[x] !== x) {{ parent[x] = parent[parent[x]]; x = parent[x]; }} return x; }}
    function union(a, b) {{ const ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; }}
    for (let i = 0; i < n; i++) {{
        for (let j = i + 1; j < n; j++) {{
            if (getPrecomputedSim(i, j) >= threshold) union(i, j);
        }}
    }}
    const groups = {{}};
    for (let k = 0; k < n; k++) {{
        const r = find(k);
        if (!groups[r]) groups[r] = [];
        groups[r].push(samples[k]);
    }}
    const clusters = Object.values(groups).filter(g => g.length >= 2);
    clusters.sort((a, b) => b.length - a.length);

    let html = '';
    if (!clusters.length) {{
        html = '<div class="alert-box alert-success">No clusters at ≥' + threshold + '%.</div>';
    }} else {{
        clusters.forEach((cluster, idx) => {{
            let sum = 0, cnt = 0;
            for (let a = 0; a < cluster.length; a++) {{
                for (let b = a + 1; b < cluster.length; b++) {{
                    const ia = samples.indexOf(cluster[a]);
                    const ib = samples.indexOf(cluster[b]);
                    sum += getPrecomputedSim(ia, ib);
                    cnt++;
                }}
            }}
            const mean = cnt ? (sum / cnt).toFixed(1) : '100';
            const color = similarityColor(parseFloat(mean));
            html += '<div class="cluster-card" style="border-left-color:' + color + ';">' +
                '<div class="cluster-header">' +
                    '<span class="cluster-badge" style="background:' + color + ';">Cluster ' + (idx + 1) + '</span>' +
                    '<span class="cluster-stats">' + cluster.length + ' samples · mean ' + mean + '% similar</span>' +
                '</div>' +
                '<div class="cluster-samples">' +
                    cluster.map(s => '<span class="cluster-sample-pill">' + escapeHtml(s) + '</span>').join('') +
                '</div></div>';
        }});
    }}
    document.getElementById('cluster-clusters').innerHTML = html;
    document.getElementById('cluster-summary').textContent =
        clusters.length + ' cluster(s) at ≥' + threshold + '%';

    if (n <= 50) renderHeatmapDOM(samples, n);
    else renderHeatmapCanvas(samples, n);
}}

function renderHeatmapDOM(samples, n) {{
    const rowLabelW = 200;
    const cellSize = 26;
    const trunc = (s, m) => s.length > m ? s.substring(0, m - 1) + '…' : s;
    let html = '<div class="heatmap-wrapper"><h3>🔷 Pairwise Similarity Matrix</h3>' +
        '<div style="overflow-x:auto;padding-bottom:10px;">' +
        '<div class="heatmap-grid" style="grid-template-columns:' + rowLabelW +
        'px repeat(' + n + ',' + cellSize + 'px);">';
    html += '<div class="heatmap-cell" style="background:transparent;"></div>';
    samples.forEach(s => {{
        html += '<div class="heatmap-cell heatmap-col-header" title="' + escapeHtml(s) + '">' +
                escapeHtml(trunc(s, 15)) + '</div>';
    }});
    samples.forEach((s1, i) => {{
        html += '<div class="heatmap-cell heatmap-row-header" style="width:' + rowLabelW +
                'px;" title="' + escapeHtml(s1) + '">' + escapeHtml(trunc(s1, 22)) + '</div>';
        samples.forEach((s2, j) => {{
            const pct = (i === j) ? 100 : getPrecomputedSim(i, j);
            const color = (i === j) ? '#1f2937' : similarityColor(pct);
            html += '<div class="heatmap-cell" style="background:' + color + ';color:white;" title="' +
                    escapeHtml(s1) + ' vs ' + escapeHtml(s2) + ': ' + pct + '%">' +
                    (i === j ? '—' : Math.round(pct)) + '</div>';
        }});
    }});
    html += '</div></div></div>';
    document.getElementById('cluster-heatmap').innerHTML = html;
}}

function renderHeatmapCanvas(samples, n) {{
    const cellSize = Math.max(3, Math.min(10, Math.floor(700 / n)));
    const rowLabelW = 150;
    const headerH = 60;
    const canvas = document.createElement('canvas');
    canvas.width = rowLabelW + n * cellSize + 20;
    canvas.height = headerH + n * cellSize + 20;
    canvas.style.display = 'block';
    canvas.style.maxWidth = '100%';
    const ctx = canvas.getContext('2d');
    ctx.font = Math.max(8, Math.min(11, cellSize - 1)) + 'px sans-serif';
    ctx.textBaseline = 'middle';
    for (let i = 0; i < n; i++) {{
        for (let j = 0; j < n; j++) {{
            const pct = (i === j) ? 100 : getPrecomputedSim(i, j);
            ctx.fillStyle = (i === j) ? '#1f2937' : similarityColor(pct);
            ctx.fillRect(rowLabelW + j * cellSize, headerH + i * cellSize, cellSize, cellSize);
        }}
    }}
    ctx.fillStyle = '#333';
    ctx.textAlign = 'right';
    const step = Math.max(1, Math.floor(n / 40));
    for (let i = 0; i < n; i += step) {{
        const lb = samples[i].length > 18 ? samples[i].substring(0, 17) + '…' : samples[i];
        ctx.fillText(lb, rowLabelW - 5, headerH + i * cellSize + cellSize / 2);
    }}
    const wrapper = document.getElementById('cluster-heatmap');
    wrapper.innerHTML = '<div class="heatmap-wrapper"><h3>🔷 Pairwise Similarity Matrix</h3>' +
        '<p style="font-size:0.85em;color:#64748b;">Canvas-rendered for ' + n + ' samples.</p>' +
        '<div id="canvas-heatmap-host" style="overflow-x:auto;"></div></div>';
    wrapper.querySelector('#canvas-heatmap-host').appendChild(canvas);
}}

document.addEventListener('DOMContentLoaded', function() {{
    const hash = window.location.hash.substring(1);
    const target = hash ? document.querySelector('.tab-button.' + hash) : document.querySelector('.tab-button');
    if (target) target.click();
    document.querySelectorAll('.data-table').forEach(table => {{
        table.querySelectorAll('th').forEach((th, idx) => {{
            const type = th.getAttribute('data-sort') || 'string';
            th.style.cursor = 'pointer';
            th.addEventListener('click', () => sortTable(table.id, idx, type));
        }});
    }});
    document.querySelectorAll('.copy-btn').forEach(btn => {{
        btn.addEventListener('click', function() {{
            const cit = this.getAttribute('data-citation');
            if (cit) {{
                navigator.clipboard.writeText(cit).then(() => {{
                    const orig = this.innerHTML;
                    this.innerHTML = '✓ Copied';
                    setTimeout(() => {{ this.innerHTML = orig; }}, 2000);
                }});
            }}
        }});
    }});
    document.querySelectorAll('.accordion-header').forEach(h => {{
        h.addEventListener('click', function() {{
            const content = this.nextElementSibling;
            content.style.display = content.style.display === 'block' ? 'none' : 'block';
        }});
    }});
}});
</script>"""

    # ------------------------------------------------------------------ sections
    def _stat_card(self, value, label: str) -> str:
        return f'<div class="metric-card"><div class="metric-value">{value}</div><div class="metric-label">{label}</div></div>'

    def _generate_summary_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        patterns = data.get('patterns', {})
        gene_centric = data.get('gene_centric', {})
        total = len(samples)
        st_dist = patterns.get('st_distribution', {})
        k_dist = patterns.get('k_locus_distribution', {})
        high_risk = len(patterns.get('high_risk_combinations', []))
        most_common_st = max(st_dist.items(), key=lambda x: x[1])[0] if st_dist else 'ND'
        most_common_k = max(k_dist.items(), key=lambda x: x[1])[0] if k_dist else 'ND'
        all_genes = gene_centric.get('all_genes', [])
        carbapenemase_count = sum(
            1 for g in all_genes
            if any(k in g['gene'] for k in ('blaKPC', 'blaNDM', 'blaOXA-48', 'blaVIM', 'blaIMP'))
        )
        hv_count = len(gene_centric.get('by_category', {}).get('High-Risk Virulence', []))

        stat_cards = (
            self._stat_card(total, 'Samples') +
            self._stat_card(len(st_dist), 'Unique STs') +
            self._stat_card(len(k_dist), 'K-types') +
            self._stat_card(len(patterns.get('o_locus_distribution', {})), 'O-types') +
            self._stat_card(len(all_genes), 'Unique Genes') +
            self._stat_card(carbapenemase_count, 'Carbapenemases') +
            self._stat_card(hv_count, 'Virulence Genes') +
            self._stat_card(high_risk, 'High-Risk Combos')
        )

        alerts = ''
        if carbapenemase_count > 0:
            alerts += self._alert('danger', 'fa-biohazard',
                f'<strong>{carbapenemase_count} carbapenemase gene(s) detected.</strong> See the AMR tab.')
        if hv_count > 0:
            alerts += self._alert('warning', 'fa-virus',
                f'<strong>{hv_count} high-risk virulence gene(s) detected.</strong> See the Virulence tab.')

        return f'''
        {self._section_header('Executive Summary', 'fa-chart-pie', 'summary')}
        {self._alert('info', 'fa-info-circle',
            f'<strong>Population overview</strong> — {total} genomes analysed. '
            f'Dominant ST: <strong>ST{most_common_st}</strong>; '
            f'dominant K-type: <strong>{most_common_k}</strong>. '
            f'Use the tabs above to explore the full dataset.')}
        <div class="metrics-grid">{stat_cards}</div>
        {alerts}
        {self._alert('info', 'fa-lightbulb',
            '<strong>What you can do with this report:</strong>'
            '<ul style="margin:8px 0 0 18px;">'
            '<li><strong>Explore population structure</strong> — MLST distribution, Kaptive K/O typing, and per-sample overview.</li>'
            '<li><strong>Track resistance genes</strong> — the AMR tab shows carbapenemases (KPC, NDM, OXA-48), ESBLs (CTX-M, SHV, TEM), colistin (mcr), tigecycline (tetX), and aminoglycoside-modifying enzymes with all genomes that carry each gene.</li>'
            '<li><strong>Examine virulence markers</strong> — the Virulence tab highlights hypervirulence siderophores (ybt, iro, iuc), colibactin (clb), and hypermucoidy regulators (rmp).</li>'
            '<li><strong>Compare any two samples</strong> side-by-side or <strong>detect clusters</strong> of highly similar genomes in the Compare tab.</li>'
            '<li><strong>Inspect MGE profiles</strong> — mobileOG-db protein signatures per sample in the MGE tab.</li>'
            '<li><strong>Export any table to CSV</strong> for downstream analysis, or download the full JSON dataset.</li>'
            '</ul>')}
        '''

    def _generate_sample_overview_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        mge_per_sample = data.get('mge_data', {}).get('per_sample', {})
        rows = []
        for sample, d in sorted(samples.items()):
            st = d.get('mlst', {}).get('ST', 'ND')
            st_disp = f"ST{st}" if st != 'ND' else 'ND'
            k = d.get('kaptive', {})
            k_locus = k.get('K_Locus', 'ND')
            k_type = k.get('K_Type', 'ND')
            o_locus = k.get('O_Locus', 'ND')
            o_type = k.get('O_Type', 'ND')
            vir_genes = []
            for db in ('vfdb', 'ecoli_vf'):
                vir_genes.extend(d.get('abricate_genes', {}).get(db, []))
            vir_genes = sorted(set(vir_genes))
            vir_count = len(vir_genes)
            if vir_count > 0:
                tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in vir_genes)
                vir_cell = (
                    f'<details class="vir-details">'
                    f'<summary>{vir_count} gene(s)</summary>'
                    f'<div class="vir-gene-tags">{tags}</div>'
                    f'</details>'
                )
            else:
                vir_cell = '<span style="color:#9ca3af;">0</span>'
            mge_hits = mge_per_sample.get(sample, {}).get('mobileOG Hits', 0)
            rows.append(
                f'<tr><td><strong>{sample}</strong></td>'
                f'<td>{st_disp}</td>'
                f'<td>{k_locus}</td>'
                f'<td>{k_type}</td>'
                f'<td>{o_locus}</td>'
                f'<td>{o_type}</td>'
                f'<td>{vir_cell}</td>'
                f'<td>{mge_hits}</td></tr>'
            )
        return f'''
        {self._section_header('Sample Overview', 'fa-list', 'sample_overview')}
        {self._alert('info', 'fa-info-circle',
            '<strong>Master sample table</strong> — MLST, Kaptive K/O typing, expandable virulence gene list, and mobileOG hit counts. Click the "N gene(s)" link to reveal the full virulence list inline. Use the search boxes to filter the table or highlight specific genome tags.')}
        {self._search_pair('samples-table', 'Search by sample, ST, K, or O...')}
        {self._action_buttons('samples-table')}
        <div class="master-scrollable-container">
            <table id="samples-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Sample</th>
                    <th data-sort="string">ST</th>
                    <th data-sort="string">K Locus</th>
                    <th data-sort="string">K Type</th>
                    <th data-sort="string">O Locus</th>
                    <th data-sort="string">O Type</th>
                    <th data-sort="string">Virulence Genes</th>
                    <th data-sort="number">MobileOG Hits</th>
                </tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>'''

    def _generate_mlst_section(self, data: Dict) -> str:
        """MLST tab — distribution + combos only. No per-sample detail."""
        patterns = data.get('patterns', {})
        st_dist = patterns.get('st_distribution', Counter())
        credit = self._credit_bar(
            '🧬',
            'MLST — PubMLST typing scheme',
            'Powered by <strong>mlst</strong> (Seemann T.) and the '
            '<a href="https://pubmlst.org/" target="_blank">PubMLST database</a> '
            'for <em>K. pneumoniae</em>. '
            'Cite: Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: '
            'BIGSdb software, the PubMLST.org website and their applications. '
            '<em>Wellcome Open Res</em>. 2018;3:124. '
            '<a href="https://doi.org/10.12688/wellcomeopenres.14826.1" target="_blank">🔗 DOI</a>'
        )
        return f'''
        {self._section_header('MLST Distribution', 'fa-code-branch', 'mlst')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>MLST (Multi-Locus Sequence Typing)</strong> indexes internal fragments of seven '
            'housekeeping genes — <em>gapA</em>, <em>infB</em>, <em>mdh</em>, <em>pgi</em>, '
            '<em>phoE</em>, <em>rpoB</em>, and <em>tonB</em>. Each unique combination of alleles '
            'defines a <strong>Sequence Type (ST)</strong>, the global standard for '
            '<em>K. pneumoniae</em> population structure and epidemiological surveillance. '
            'Closely related STs are grouped into <strong>Clonal Complexes (CC)</strong> — '
            'for example, ST258 belongs to CC258, the most prevalent KPC-producing clone worldwide. '
            'Use the combination tables below to explore how STs associate with capsule types '
            'and lipopolysaccharide serotypes.')}
        <h3 class="sub-heading">📊 ST Distribution</h3>
        {self._distribution_table('mlst-dist-table', st_dist, 'ST')}
        {self._combo_table('mlst-k-table', 'ST – K Locus Combination', patterns.get('st_k_combinations', {}))}
        {self._combo_table('mlst-o-table', 'ST – O Locus Combination', patterns.get('st_o_combinations', {}))}
        {self._combo_table('mlst-ko-table', 'ST – K:O Combination', patterns.get('st_ko_combinations', {}))}
        '''

    def _generate_qc_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        all_metrics = set()
        for d in samples.values():
            all_metrics.update(d.get('qc', {}).keys())
        if not all_metrics:
            return f'{self._section_header("FASTA QC", "fa-chart-line", "qc")}{self._alert("warning", "fa-exclamation-circle", "No QC data available.")}'
        metric_list = sorted(all_metrics)
        rows = []
        for sample, d in sorted(samples.items()):
            qc = d.get('qc', {})
            cells = [f'<td><strong>{sample}</strong></td>']
            for m in metric_list:
                v = qc.get(m, 'ND')
                if isinstance(v, float):
                    v = f"{v:,.0f}" if v > 1000 else f"{v:.1f}"
                cells.append(f'<td>{v}</td>')
            rows.append('<tr>' + ''.join(cells) + '</tr>')
        headers = ''.join(f'<th data-sort="number">{m}</th>' for m in metric_list)
        credit = self._credit_bar(
            '📊',
            'Biopython + fastANI — Assembly QC and species confirmation',
            'Assembly metrics (contigs, N50, GC%, length) computed with <strong>Biopython</strong>. '
            'Species confirmed with <strong>fastANI</strong> (Jain C et al. <em>Nat Commun</em> 2018;9(1):5114). '
            '<a href="https://doi.org/10.1038/s41467-018-07641-9" target="_blank">🔗 DOI</a>'
        )
        return f'''
        {self._section_header('FASTA QC', 'fa-chart-line', 'qc')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>Assembly quality metrics</strong> — <strong>N50</strong> is the contig length at which 50% of the assembly is in larger contigs (higher is better). '
            '<strong>Number of contigs</strong> — lower is better; typical <em>K. pneumoniae</em>: &lt;200 good, &lt;100 excellent. '
            '<strong>GC%</strong> — <em>K. pneumoniae</em> is typically 55–58%. '
            '<strong>Total length</strong> — approximately 5.0–5.6 Mbp for a complete genome. '
            'Poor assemblies may miss genes or break up plasmids, affecting downstream interpretation.')}
        {self._search_pair('qc-table', 'Search sample...')}
        {self._action_buttons('qc-table')}
        <div class="master-scrollable-container">
            <table id="qc-table" class="data-table">
                <thead><tr><th data-sort="string">Sample</th>{headers}</tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>'''

    def _generate_kaptive_section(self, data: Dict) -> str:
        """Kaptive tab — distribution + combos only. No per-sample detail table."""
        patterns = data.get('patterns', {})
        credit = self._credit_bar(
            '🛡️',
            'Kaptive — K and O locus typing',
            'Powered by <strong>Kaptive</strong> (Wyres KL, Holt KE et al.) '
            'with <em>kpsc_k</em> / <em>kpsc_o</em> reference databases. '
            'Cite: Wyres KL et al. Identification of Klebsiella capsule synthesis loci from '
            'whole genome sequence data. <em>Microb Genom</em>. 2016;2(5):e000064. '
            '<a href="https://doi.org/10.1099/mgen.0.000064" target="_blank">🔗 DOI</a>'
        )
        return f'''
        {self._section_header('Kaptive K/O Typing', 'fa-shield-alt', 'kaptive')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>Kaptive</strong> identifies the <strong>K</strong> (capsule) and <strong>O</strong> '
            '(lipopolysaccharide) loci from whole-genome sequence data. '
            'The capsule protects <em>K. pneumoniae</em> from phagocytosis and complement killing, and '
            'specific K-types are strongly associated with hypervirulent clones — '
            '<strong>K1</strong> with ST23 (liver abscess), '
            '<strong>K2</strong> with ST86 and ST375, '
            '<strong>KL64</strong> with ST11 (carbapenem-resistant), and '
            '<strong>KL106</strong> with ST258 (KPC-producing). '
            'The O-locus type (O1, O2, O3, O5) influences serum resistance and immune evasion. '
            'Higher identity and coverage (≥95%) indicate reliable typing results.')}
        <h3 class="sub-heading">📊 K Locus Distribution</h3>
        {self._distribution_table('kaptive-k-dist', patterns.get('k_locus_distribution', Counter()), 'K Locus')}
        <h3 class="sub-heading">📊 O Locus Distribution</h3>
        {self._distribution_table('kaptive-o-dist', patterns.get('o_locus_distribution', Counter()), 'O Locus')}
        {self._combo_table('kaptive-st-k', 'ST – K Locus Combination', patterns.get('st_k_combinations', {}))}
        {self._combo_table('kaptive-ko', 'K : O Capsule Type', patterns.get('ko_combinations', {}))}
        '''

    def _generate_combinations_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        html_parts = [self._section_header('Combination Tables', 'fa-link', 'combinations')]
        html_parts.append(self._alert('info', 'fa-info-circle',
            '<strong>Typing combinations</strong> reveal clonal associations between Sequence Type (ST), '
            'capsule type (K locus), and lipopolysaccharide serotype (O locus). '
            'For example, <strong>ST23 is classically K1</strong> (hypervirulent liver-abscess clone), '
            '<strong>ST258 is associated with KL106 or KL107</strong> (KPC-producing epidemic), '
            '<strong>ST11 with KL64</strong> (carbapenem-resistant), and '
            '<strong>ST307 with various K-types</strong> (emerging MDR clone). '
            'These patterns help identify high-risk lineages and understand clonal expansion.'))
        html_parts.append(self._combo_table(
            'comb-st-k', 'ST – K Locus Combination',
            patterns.get('st_k_combinations', {}),
        ))
        html_parts.append(self._combo_table(
            'comb-st-o', 'ST – O Locus Combination',
            patterns.get('st_o_combinations', {}),
        ))
        html_parts.append(self._combo_table(
            'comb-ko', 'K : O Capsule Type',
            patterns.get('ko_combinations', {}),
        ))
        html_parts.append(self._combo_table(
            'comb-st-ko', 'ST – K : O Combination',
            patterns.get('st_ko_combinations', {}),
        ))
        return ''.join(html_parts)

    def _generate_amr_section(self, data: Dict) -> str:
        gene_centric = data.get('gene_centric', {})
        amr_databases = gene_centric.get('amr_databases', {})
        total_samples = len(data.get('samples', {}))
        all_genes = []
        for genes in amr_databases.values():
            all_genes.extend(genes)
        all_genes.sort(key=lambda x: x['count'], reverse=True)
        credit = self._credit_bar(
            '💊',
            'ABRicate + AMRFinderPlus — Resistance detection',
            'Powered by <strong>ABRicate</strong> (Seemann T.) with the CARD, ResFinder, ARG-ANNOT, '
            'MEGARes, and NCBI AMR databases, plus <strong>AMRFinderPlus</strong> (NCBI). '
            'Every hit preserves its database provenance — filter by database name to inspect '
            'a single source. Remember: gene presence is not the same as phenotypic resistance — '
            'AST remains the clinical gold standard.'
        )
        filter_buttons = self._action_buttons('amr-table', [
            ('KPC', 'blaKPC', 'btn-danger', 'fa-biohazard'),
            ('NDM', 'blaNDM', 'btn-danger', 'fa-biohazard'),
            ('OXA-48', 'OXA-48', 'btn-danger', 'fa-biohazard'),
            ('OXA-181', 'OXA-181', 'btn-danger', 'fa-biohazard'),
            ('OXA-232', 'OXA-232', 'btn-danger', 'fa-biohazard'),
            ('IMP', 'blaIMP', 'btn-danger', 'fa-biohazard'),
            ('VIM', 'blaVIM', 'btn-danger', 'fa-biohazard'),
            ('GES', 'blaGES', 'btn-danger', 'fa-biohazard'),
            ('CTX-M', 'CTX-M', 'btn-warning', 'fa-pills'),
            ('SHV', 'SHV', 'btn-warning', 'fa-pills'),
            ('TEM', 'TEM', 'btn-warning', 'fa-pills'),
            ('CMY', 'blaCMY', 'btn-warning', 'fa-pills'),
            ('DHA', 'blaDHA', 'btn-warning', 'fa-pills'),
            ('ACT', 'blaACT', 'btn-warning', 'fa-pills'),
            ('mcr', 'mcr', 'btn-info', 'fa-pills'),
            ('tet(X)', 'tet(X)', 'btn-info', 'fa-pills'),
            ('armA', 'armA', 'btn-info', 'fa-pills'),
            ('rmt', 'rmt', 'btn-info', 'fa-pills'),
            ('npmA', 'npmA', 'btn-info', 'fa-pills'),
            ('aac', 'aac', 'btn-secondary', 'fa-syringe'),
            ('ant', 'ant', 'btn-secondary', 'fa-syringe'),
            ('aph', 'aph', 'btn-secondary', 'fa-syringe'),
            ('aad', 'aad', 'btn-secondary', 'fa-syringe'),
            ('strA', 'strA', 'btn-secondary', 'fa-syringe'),
            ('strB', 'strB', 'btn-secondary', 'fa-syringe'),
            ('sul', 'sul', 'btn-secondary', 'fa-tablets'),
            ('dfr', 'dfr', 'btn-secondary', 'fa-tablets'),
            ('cat', 'cat', 'btn-secondary', 'fa-tablets'),
            ('floR', 'floR', 'btn-secondary', 'fa-tablets'),
            ('tetA', 'tetA', 'btn-secondary', 'fa-capsules'),
            ('tetB', 'tetB', 'btn-secondary', 'fa-capsules'),
            ('tetC', 'tetC', 'btn-secondary', 'fa-capsules'),
            ('tetD', 'tetD', 'btn-secondary', 'fa-capsules'),
            ('tetM', 'tetM', 'btn-secondary', 'fa-capsules'),
            ('fosA', 'fosA', 'btn-secondary', 'fa-tablets'),
            ('qnr', 'qnr', 'btn-secondary', 'fa-filter'),
            ('gyrA', 'gyrA_', 'btn-secondary', 'fa-filter'),
            ('parC', 'parC_', 'btn-secondary', 'fa-filter'),
            ('OXA (all)', 'blaOXA', 'btn-secondary', 'fa-biohazard'),
        ])
        return f'''
        {self._section_header('AMR Genes', 'fa-biohazard', 'amr')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>Antimicrobial resistance genes</strong> identified from all sources. '
            'Key families: <strong>carbapenemases</strong> (KPC, NDM, OXA-48, VIM, IMP, GES — these are the last-resort resistance determinants), '
            '<strong>ESBLs</strong> (CTX-M, SHV, TEM), '
            '<strong>colistin</strong> (mcr-1 through mcr-10), '
            '<strong>tigecycline</strong> (tetX variants), '
            '<strong>aminoglycoside-modifying enzymes</strong> (aac, aph, ant), '
            'and <strong>16S rRNA methyltransferases</strong> (armA, rmtA–rmtH). '
            'Use the filter buttons to quickly locate specific gene families, or use the '
            'grouping controls to see which STs carry each gene.')}
        {filter_buttons}
        {self._gene_table('amr-table', 'All AMR Genes', all_genes, total_samples)}'''

    def _generate_virulence_section(self, data: Dict) -> str:
        gene_centric = data.get('gene_centric', {})
        vir_databases = gene_centric.get('virulence_databases', {})
        total_samples = len(data.get('samples', {}))
        all_vir = []
        for genes in vir_databases.values():
            all_vir.extend(genes)
        all_vir.sort(key=lambda x: x['count'], reverse=True)
        credit = self._credit_bar(
            '🦠',
            'VFDB — Virulence factor detection',
            'Powered by <strong>ABRicate</strong> with the '
            '<a href="http://www.mgc.ac.cn/VFs/" target="_blank">VFDB database</a>. '
            'Cite: Chen L, Xiong Z, Sun L, Yang J, Jin Q. VFDB 2012 update: toward the genetic '
            'diversity and molecular evolution of bacterial virulence factors. '
            '<em>Nucleic Acids Res</em>. 2012;40(D1):D641-5. '
            '<a href="https://doi.org/10.1093/nar/gkr989" target="_blank">🔗 DOI</a>'
        )
        filter_buttons = self._action_buttons('vir-table', [
            ('ybt', 'ybt', 'btn-success', 'fa-syringe'),
            ('clb', 'clb', 'btn-success', 'fa-syringe'),
            ('iro', 'iro', 'btn-success', 'fa-syringe'),
            ('iuc', 'iuc', 'btn-success', 'fa-syringe'),
            ('rmp', 'rmp', 'btn-success', 'fa-syringe'),
            ('rmpA2', 'rmpA2', 'btn-success', 'fa-syringe'),
            ('fim', 'fim', 'btn-info', 'fa-filter'),
            ('mrk', 'mrk', 'btn-info', 'fa-filter'),
            ('ecp', 'ecp', 'btn-info', 'fa-filter'),
            ('tss', 'tss', 'btn-info', 'fa-filter'),
            ('hcp', 'hcp', 'btn-info', 'fa-filter'),
            ('vgrG', 'vgrG', 'btn-info', 'fa-filter'),
            ('ent', 'ent', 'btn-secondary', 'fa-filter'),
            ('fep', 'fep', 'btn-secondary', 'fa-filter'),
            ('fec', 'fec', 'btn-secondary', 'fa-filter'),
            ('hly', 'hly', 'btn-secondary', 'fa-filter'),
            ('cnf', 'cnf', 'btn-secondary', 'fa-filter'),
            ('cdt', 'cdt', 'btn-secondary', 'fa-filter'),
            ('sat', 'sat', 'btn-secondary', 'fa-filter'),
            ('pic', 'pic', 'btn-secondary', 'fa-filter'),
        ])
        return f'''
        {self._section_header('Virulence Genes', 'fa-virus', 'virulence')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>Virulence factors</strong> — siderophores (ybt, iro, iuc), '
            'colibactin (clb), hypermucoidy regulators (rmp), adhesins (fim, mrk, ecp), '
            'and type VI secretion systems (tss). The siderophore genes '
            '<strong>iucA–iucD</strong> (aerobactin) and <strong>iroB–iroN</strong> (salmochelin) '
            'together with <strong>rmpA/rmpA2</strong> define the classic hypervirulent hvKp '
            'signature. Use the filter buttons to quickly locate specific families, or use the '
            'grouping controls to see which STs carry each gene.')}
        {filter_buttons}
        {self._gene_table('vir-table', 'All Virulence Genes', all_vir, total_samples)}'''

    def _generate_bacmet_section(self, data: Dict) -> str:
        gene_centric = data.get('gene_centric', {})
        bac_databases = gene_centric.get('bacmet_databases', {})
        total_samples = len(data.get('samples', {}))
        all_bac = []
        for genes in bac_databases.values():
            all_bac.extend(genes)
        all_bac.sort(key=lambda x: x['count'], reverse=True)
        credit = self._credit_bar(
            '🧪',
            'BacMet — Biocide and heavy metal resistance',
            'Powered by <strong>ABRicate</strong> with the '
            '<a href="http://bacmet.biomedicine.gu.se/" target="_blank">BacMet database</a>. '
            'Cite: Pal C, Bengtsson-Palme J, Rensing C, Kristiansson E, Larsson DGJ. '
            'BacMet: antibacterial biocide and metal resistance genes database. '
            '<em>Nucleic Acids Res</em>. 2014;42(D1):D737-43. '
            '<a href="https://doi.org/10.1093/nar/gkt1252" target="_blank">🔗 DOI</a>'
        )
        filter_buttons = self._action_buttons('bac-table', [
            ('qac', 'qac', 'btn-info', 'fa-filter'),
            ('qacE', 'qacE', 'btn-info', 'fa-filter'),
            ('cep', 'cep', 'btn-info', 'fa-filter'),
            ('form', 'form', 'btn-info', 'fa-filter'),
            ('sil', 'sil', 'btn-info', 'fa-filter'),
            ('mer', 'mer', 'btn-warning', 'fa-filter'),
            ('merA', 'merA', 'btn-warning', 'fa-filter'),
            ('ars', 'ars', 'btn-warning', 'fa-filter'),
            ('arsB', 'arsB', 'btn-warning', 'fa-filter'),
            ('pco', 'pco', 'btn-warning', 'fa-filter'),
            ('cop', 'cop', 'btn-warning', 'fa-filter'),
            ('chr', 'chr', 'btn-warning', 'fa-filter'),
            ('cad', 'cad', 'btn-warning', 'fa-filter'),
            ('czc', 'czc', 'btn-warning', 'fa-filter'),
            ('znt', 'znt', 'btn-warning', 'fa-filter'),
            ('soxR', 'soxR', 'btn-secondary', 'fa-filter'),
            ('cpxR', 'cpxR', 'btn-secondary', 'fa-filter'),
            ('baeR', 'baeR', 'btn-secondary', 'fa-filter'),
            ('emr', 'emr', 'btn-secondary', 'fa-filter'),
        ])
        return f'''
        {self._section_header('Biocide & Heavy Metal Resistance', 'fa-flask', 'bacmet')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>BACMET genes</strong> confer resistance to hospital disinfectants '
            '(<strong>qac</strong> — quaternary ammonium compounds, <strong>cep</strong> — chlorhexidine) '
            'and heavy metals (<strong>mer</strong> — mercury, <strong>ars</strong> — arsenic, '
            '<strong>sil</strong> — silver, <strong>pco</strong> — copper). '
            'These are important <strong>co-selection markers</strong> — the same plasmids or transposons '
            'often carry both biocide and antibiotic resistance genes, so disinfectant use can drive '
            'antibiotic resistance in hospital settings.')}
        {filter_buttons}
        {self._gene_table('bac-table', 'All BACMET Genes', all_bac, total_samples)}'''

    def _generate_plasmids_section(self, data: Dict) -> str:
        gene_centric = data.get('gene_centric', {})
        plasmid_databases = gene_centric.get('plasmid_databases', {})
        total_samples = len(data.get('samples', {}))
        all_plas = []
        for genes in plasmid_databases.values():
            all_plas.extend(genes)
        all_plas.sort(key=lambda x: x['count'], reverse=True)
        credit = self._credit_bar(
            '🧬',
            'PlasmidFinder — Replicon typing',
            'Powered by <strong>ABRicate</strong> with the '
            '<a href="https://genepi.food.dtu.dk/PlasmidFinder/" target="_blank">PlasmidFinder database</a>. '
            'Cite: Carattoli A, Zankari E, García-Fernández A, Voldby Larsen M, Lund O, '
            'Villa L, Møller Aarestrup F, Hasman H. In silico detection and typing of plasmids '
            'using PlasmidFinder and plasmid multilocus sequence typing. '
            '<em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903. '
            '<a href="https://doi.org/10.1128/AAC.02412-14" target="_blank">🔗 DOI</a>'
        )
        return f'''
        {self._section_header('Plasmid Replicons', 'fa-dna', 'plasmids')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>Plasmid replicons</strong> reveal the horizontal gene transfer potential of your genomes. '
            'Common plasmid families in <em>K. pneumoniae</em> include <strong>IncF, IncI, IncX, IncN, '
            'IncH, IncA/C</strong>, and <strong>Col</strong>-type plasmids. '
            'The presence of specific replicons often correlates with the carriage of specific '
            'resistance genes — for example, IncX3 plasmids frequently carry <em>blaNDM</em>, '
            'IncF plasmids often carry <em>blaCTX-M</em> or <em>blaKPC</em>, and IncHI2 plasmids '
            'frequently carry <em>mcr</em> colistin resistance genes.')}
        {self._gene_table('plasmids-table', 'Plasmid Replicons', all_plas, total_samples)}'''

    def _generate_mutations_section(self, data: Dict) -> str:
        """Point mutations tab — aggregated per gene, no mutation/subclass columns.

        Mutations for each gene are shown as a small italic subtitle under the
        gene name, so the table stays at 4 columns and the Genomes column gets
        the vast majority of the row width — important when there are 1000+
        genomes to display.
        """
        mutation_data = data.get('mutation_data', {})
        mutations = mutation_data.get('mutations', [])
        total_samples = len(data.get('samples', {}))
        credit = self._credit_bar(
            '🧬',
            'AMRFinderPlus — Point mutations',
            'Powered by <strong>AMRFinderPlus</strong> (NCBI). Detects resistance-conferring '
            'mutations in <em>gyrA</em>, <em>parC</em>, <em>rpoB</em>, <em>mgrB</em>, '
            '<em>pmrAB</em>, <em>phoPQ</em>, and other targets. '
            'Cite: Feldgarden M et al. AMRFinderPlus and the Reference Gene Catalog facilitate '
            'examination of the genomic links among antimicrobial resistance, stress response, '
            'and virulence. <em>Sci Rep</em>. 2021;11(1):12728. '
            '<a href="https://doi.org/10.1038/s41598-021-91456-0" target="_blank">🔗 DOI</a>'
        )
        if not mutations:
            return f'''
            {self._section_header('Point Mutations', 'fa-dna', 'mutations')}
            {credit}
            {self._alert('warning', 'fa-exclamation-circle',
                'No mutation data found. Please ensure mutation_summary.html is present in the input directory.')}'''

        class_counts = Counter()
        for m in mutations:
            class_counts[m['class']] += 1
        summary_badges = ''.join(
            f"<span class='badge badge-info' style='margin:2px;'>{cls}: {cnt}</span>"
            for cls, cnt in class_counts.items() if cls and cls != "NA"
        )
        extra_buttons = []
        for cls in sorted(class_counts.keys()):
            if cls and cls != "NA":
                extra_buttons.append((cls, cls, 'btn-info', 'fa-tag'))
        extra_buttons.extend([
            ('gyrA', 'gyrA', 'btn-danger', 'fa-filter'),
            ('parC', 'parC', 'btn-danger', 'fa-filter'),
            ('rpoB', 'rpoB', 'btn-warning', 'fa-filter'),
            ('mgrB', 'mgrB', 'btn-info', 'fa-filter'),
            ('pmr', 'pmr', 'btn-info', 'fa-filter'),
            ('lpx', 'lpx', 'btn-info', 'fa-filter'),
            ('envZ', 'envZ', 'btn-info', 'fa-filter'),
        ])

        # Aggregate mutations by gene so each gene appears exactly once with the
        # union of its mutation names, its class labels, and all affected genomes.
        gene_groups = defaultdict(lambda: {
            'mutations': set(),
            'classes': set(),
            'genomes': set(),
        })
        for m in mutations:
            gene = m.get('gene', '')
            if not gene:
                continue
            gene_groups[gene]['mutations'].add(m.get('mutation', ''))
            if m.get('class'):
                gene_groups[gene]['classes'].add(m['class'])
            gene_groups[gene]['genomes'].update(m.get('genomes', []))

        rows = []
        for gene in sorted(gene_groups.keys()):
            grp = gene_groups[gene]
            mutation_list = ', '.join(sorted(x for x in grp['mutations'] if x))
            class_list = ', '.join(sorted(grp['classes'])) or '—'
            genome_list = sorted(grp['genomes'])
            tags = ''.join(f'<span class="genome-tag">{g}</span>' for g in genome_list)
            count = len(genome_list)
            pct = min(count / total_samples * 100 if total_samples else 0, 100)

            mutation_subtitle = (
                f'<span class="mutation-subtitle">{mutation_list}</span>'
                if mutation_list else ''
            )
            rows.append(
                f'<tr>'
                f'<td><strong>{gene}</strong>{mutation_subtitle}</td>'
                f'<td>{class_list}</td>'
                f'<td>{count} ({pct:.1f}%)</td>'
                f'<td><div class="genome-list">{tags}</div></td>'
                f'</tr>'
            )

        return f'''
        {self._section_header('Point Mutations', 'fa-dna', 'mutations')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            '<strong>Point mutations</strong> can confer resistance even without acquired genes. '
            'Common targets in <em>K. pneumoniae</em>: <strong><em>gyrA</em></strong> and '
            '<strong><em>parC</em></strong> (quinolone resistance — ciprofloxacin, levofloxacin), '
            '<strong><em>rpoB</em></strong> (rifampin resistance), '
            '<strong><em>mgrB</em></strong> (colistin resistance via PhoPQ upregulation), '
            '<strong><em>pmrA</em></strong>/<strong><em>pmrB</em></strong> (colistin resistance via lipid A modification), '
            'and <strong><em>ompK35</em></strong>/<strong><em>ompK36</em></strong> (porin loss affecting carbapenem uptake). '
            'Mutations are aggregated per gene — the specific amino-acid substitutions '
            'are shown in small italic text under each gene name.')}
        <div style="margin:12px 0;background:#f8f9fa;padding:12px;border-radius:8px;">
            <strong><i class="fas fa-chart-bar"></i> Summary by Mutation Class:</strong><br>
            {summary_badges or 'No class data.'}
        </div>
        {self._grouping_controls('mutations-table')}
        {self._search_pair('mutations-table', 'Search gene, mutation, or class...')}
        {self._action_buttons('mutations-table', extra_buttons)}
        <div class="master-scrollable-container">
            <table id="mutations-table" class="data-table mutations-table">
                <colgroup>
                    <col style="width: 240px;">
                    <col style="width: 180px;">
                    <col style="width: 120px;">
                    <col>
                </colgroup>
                <thead><tr>
                    <th data-sort="string">Gene</th>
                    <th data-sort="string">Class</th>
                    <th data-sort="string">Frequency</th>
                    <th data-sort="string">Genomes</th>
                </tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>'''

    def _generate_patterns_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        icekp = patterns.get('icekp_marker_presence', {})
        vp = patterns.get('virulence_plasmid_marker_presence', {})
        cooc = patterns.get('gene_cooccurrence', {})
        html = [self._section_header('Patterns & Associations', 'fa-project-diagram', 'patterns')]
        html.append(self._alert('info', 'fa-info-circle',
            '<strong>Cross-genome pattern discovery</strong> — Identify <strong>ICEKp-associated markers</strong> '
            '(ybt, clb, iro, rmp) that define hypervirulent lineages, '
            '<strong>virulence plasmid markers</strong> (iro, iuc, rmp, rmpA2) carried on the pK2044-like '
            'virulence plasmids, and <strong>gene co-occurrence networks</strong> '
            'that reveal physical or functional linkage between resistance and virulence determinants.'))
        if icekp:
            rows = ''.join(
                f'<tr><td><strong>{s}</strong></td><td>{", ".join(m)}</td></tr>'
                for s, m in icekp.items()
            )
            html.append(f'''
            <h3 class="sub-heading">ICEKp Associated Markers (ybt, clb, iro, rmp)</h3>
            {self._search_pair('icekp-table', 'Search sample...')}
            <div class="master-scrollable-container">
                <table id="icekp-table" class="data-table">
                    <thead><tr><th>Sample</th><th>Markers</th></tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>''')
        if vp:
            rows = ''.join(
                f'<tr><td><strong>{s}</strong></td><td>{", ".join(m)}</td></tr>'
                for s, m in vp.items()
            )
            html.append(f'''
            <h3 class="sub-heading">Virulence Plasmid Markers (iro, iuc, rmp, rmpA2)</h3>
            {self._search_pair('vp-table', 'Search sample...')}
            <div class="master-scrollable-container">
                <table id="vp-table" class="data-table">
                    <thead><tr><th>Sample</th><th>Markers</th></tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>''')
        if cooc:
            pairs = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    pairs.append((g1, g2, cnt))
            pairs.sort(key=lambda x: x[2], reverse=True)
            rows = ''.join(
                f'<tr><td>{g1}</td><td>{g2}</td><td>{c}</td></tr>'
                for g1, g2, c in pairs[:200]
            )
            html.append(f'''
            <h3 class="sub-heading">Gene Co-occurrence (Top 200)</h3>
            <div class="master-scrollable-container">
                <table class="data-table">
                    <thead><tr><th>Gene 1</th><th>Gene 2</th><th>Count</th></tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>''')
        return ''.join(html)

    def _generate_highrisk_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        high_risk = patterns.get('high_risk_combinations', [])
        header = self._section_header('High-Risk Combinations', 'fa-exclamation-triangle', 'highrisk')
        if not high_risk:
            return header + self._alert('success', 'fa-check-circle',
                'No samples carry both critical resistance and high-risk virulence genes.')
        rows = []
        for c in high_risk:
            st_disp = f"ST{c['st']}" if c['st'] != 'ND' else 'ND'
            rows.append(
                f'<tr><td><strong>{c["sample"]}</strong></td>'
                f'<td>{st_disp}</td><td>{c["k_locus"]}</td>'
                f'<td>{", ".join(c["critical_resistance"])}</td>'
                f'<td>{", ".join(c["high_risk_virulence"])}</td></tr>'
            )
        return f'''
        {header}
        {self._alert('danger', 'fa-biohazard',
            f'<strong>{len(high_risk)} high-risk combination(s)</strong> — samples that carry <em>both</em> '
            'critical resistance genes (carbapenemases, colistin resistance, tigecycline resistance) '
            '<em>and</em> high-risk virulence genes (hypervirulence siderophores, colibactin, hypermucoidy regulators). '
            'These represent the most clinically concerning isolates — the "superbug" convergence that is '
            'driving the global CR-hvKp crisis.')}
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead><tr>
                    <th>Sample</th><th>ST</th><th>K Locus</th>
                    <th>Critical Resistance</th><th>High-Risk Virulence</th>
                </tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>'''

    def _generate_mge_section(self, data: Dict) -> str:
        mge = data.get('mge_data', {})
        per_sample = mge.get('per_sample', {})
        analysis = mge.get('analysis', {})
        note = mge.get('note', '')
        credit = self._credit_bar(
            '🧬',
            'mobileOG-db — Mobile genetic elements',
            'Powered by <strong>mobileOG-db</strong> (Brown CL, Mullet J et al.). '
            'Cite: Brown CL, Mullet J, Hindi F, Stoll JE, Gupta S, Choi M, Keenum I, Vikesland P, '
            'Pruden A, Zhang L. mobileOG-db: a manually curated database of protein families '
            'mediating the life cycle of bacterial mobile genetic elements. '
            '<em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22. '
            '<a href="https://doi.org/10.1128/aem.00991-22" target="_blank">🔗 DOI</a>'
        )
        if not per_sample:
            return f'''
            {self._section_header('Mobile Genetic Elements', 'fa-mobile-alt', 'mge')}
            {credit}
            {self._alert('warning', 'fa-exclamation-circle',
                'No MGE data found. Ensure kleboscope_mge_summary.html is in the input directory.')}'''
        cat_rows = ''
        for c in analysis.get('category_totals', []):
            cat_rows += (
                f'<tr><td><strong>{c["category"]}</strong></td>'
                f'<td>{c["total"]:,}</td>'
                f'<td>{c["pct_of_all"]}%</td>'
                f'<td>{c["mean"]}</td>'
                f'<td>{c["max"]}</td>'
                f'<td>{c["max_sample"]}</td></tr>'
            )
        all_cols = ['Sample']
        seen = set()
        for s in per_sample.values():
            for k in s.keys():
                if k not in seen:
                    seen.add(k)
                    all_cols.append(k)
        header_cells = ''.join(
            f'<th data-sort="{"string" if c == "Sample" else "number"}">{c}</th>'
            for c in all_cols
        )
        sample_rows = ''
        for sample, vals in sorted(per_sample.items()):
            cells = f'<td><strong>{sample}</strong></td>'
            for c in all_cols[1:]:
                cells += f'<td>{vals.get(c, 0)}</td>'
            sample_rows += f'<tr>{cells}</tr>'
        samples_data = self._samples_data
        st_groups = defaultdict(list)
        for sample, vals in per_sample.items():
            st = samples_data.get(sample, {}).get('mlst', {}).get('ST', 'ND')
            if st and st != 'ND':
                st_groups[st].append(vals)
        agg_rows = ''
        for st, records in sorted(st_groups.items(), key=lambda x: -len(x[1])):
            means = {}
            for c in all_cols[1:]:
                nums = [r.get(c, 0) for r in records if isinstance(r.get(c, 0), (int, float))]
                means[c] = round(sum(nums) / len(nums), 1) if nums else 0
            cells = ''.join(f'<td>{means[c]}</td>' for c in all_cols[1:])
            agg_rows += f'<tr><td><strong>ST{st}</strong></td><td>{len(records)}</td>{cells}</tr>'
        agg_header = ''.join(f'<th>{c}</th>' for c in all_cols[1:])
        note_html = ''
        if note:
            note_html = self._alert('warning', 'fa-info-circle', f'<strong>Interpretation note:</strong> {note}')
        return f'''
        {self._section_header('Mobile Genetic Elements', 'fa-mobile-alt', 'mge')}
        {credit}
        {self._alert('info', 'fa-info-circle',
            'Mobile genetic elements (MGEs) are the physical vehicles of horizontal gene transfer — '
            'the mechanism by which bacteria acquire resistance and virulence genes from other bacteria. '
            'Without MGEs, most acquired resistance in <em>K. pneumoniae</em> would not exist. '
            '<strong>Plasmids</strong> carry carbapenemases (KPC, NDM, OXA-48) and colistin resistance (mcr). '
            '<strong>ICEKp elements</strong> carry siderophore systems (ybt, clb, iro) that define hypervirulent clones. '
            '<strong>Prophages</strong> carry virulence factors. '
            '<strong>IS elements and transposons</strong> move resistance genes into new genomic contexts.')}
        {note_html}
        <h3 class="sub-heading">📊 Category Totals</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead><tr>
                    <th>Category</th><th>Total hits</th><th>% of all</th>
                    <th>Mean/sample</th><th>Max/sample</th><th>Sample with max</th>
                </tr></thead>
                <tbody>{cat_rows}</tbody>
            </table>
        </div>
        <h3 class="sub-heading">🧬 Per-sample mobileOG Profile</h3>
        {self._search_pair('mge-table', 'Search sample...')}
        {self._action_buttons('mge-table')}
        <div class="master-scrollable-container">
            <table id="mge-table" class="data-table">
                <thead><tr>{header_cells}</tr></thead>
                <tbody>{sample_rows}</tbody>
            </table>
        </div>
        <h3 class="sub-heading">📈 Mean MGE Profile by ST</h3>
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead><tr><th>ST</th><th>N samples</th>{agg_header}</tr></thead>
                <tbody>{agg_rows or '<tr><td colspan="3" style="text-align:center;color:#9ca3af;">No ST-grouped data available</td></tr>'}</tbody>
            </table>
        </div>'''

    def _generate_compare_section(self, data: Dict) -> str:
        compare_data = data.get('compare_data', {})
        similarity = data.get('compare_similarity', {})
        header = self._section_header('Compare & Cluster', 'fa-balance-scale', 'compare')
        if len(compare_data) < 2:
            return header + self._alert('warning', 'fa-exclamation-circle',
                'Compare requires at least two samples.')
        compare_json = json.dumps(compare_data, default=str, ensure_ascii=False)
        sim_json = json.dumps(similarity, default=str, ensure_ascii=False)
        return f'''
        {header}
        {self._alert('info', 'fa-info-circle',
            '<strong>Pair Compare</strong> — do two isolates share MLST + K-type + O-type + resistance/virulence profile? '
            'This is the classic outbreak-confirmation question. '
            '<strong>Cluster Detection</strong> — scan the entire dataset to find groups of highly similar genomes '
            'that may represent recent transmission events or duplicate samples.')}
        <div class="mode-toggle">
            <button class="mode-btn active" data-mode="pair" onclick="setCompareMode('pair')">
                <i class="fas fa-balance-scale"></i> Pair Compare
            </button>
            <button class="mode-btn" data-mode="cluster" onclick="setCompareMode('cluster')">
                <i class="fas fa-project-diagram"></i> Cluster Detection
            </button>
        </div>
        <div id="mode-pair" class="compare-mode active">
            <div class="compare-panel">
                <div class="compare-select-row">
                    <label><span>Sample A</span><select id="compare-a"></select></label>
                    <button class="action-btn btn-primary" onclick="runCompare()">
                        <i class="fas fa-balance-scale"></i> Compare
                    </button>
                    <button class="action-btn btn-light" onclick="swapCompare()">
                        <i class="fas fa-exchange-alt"></i> Swap
                    </button>
                    <button class="action-btn btn-light" onclick="resetCompare()">
                        <i class="fas fa-sync"></i> Reset
                    </button>
                    <label><span>Sample B</span><select id="compare-b"></select></label>
                </div>
                <div id="compare-verdict" style="display:none;"></div>
                <div id="compare-visuals"></div>
                <div id="compare-output"></div>
            </div>
        </div>
        <div id="mode-cluster" class="compare-mode">
            <div class="compare-panel">
                <div class="compare-select-row">
                    <label>
                        <span>Similarity threshold</span>
                        <select id="cluster-threshold" onchange="runCluster()">
                            <option value="95">Tight (≥95%)</option>
                            <option value="90" selected>Moderate (≥90%)</option>
                            <option value="85">Loose (≥85%)</option>
                            <option value="80">Broad (≥80%)</option>
                        </select>
                    </label>
                    <button class="action-btn btn-primary" onclick="runCluster()">
                        <i class="fas fa-project-diagram"></i> Detect Clusters
                    </button>
                    <span id="cluster-summary" style="color:#64748b;font-weight:600;"></span>
                </div>
                <div id="cluster-clusters"></div>
                <div id="cluster-heatmap"></div>
            </div>
        </div>
        <script>
        window._COMPARE_DATA = {compare_json};
        window._SIM = {sim_json};
        </script>'''

    def _generate_databases_section(self, data: Dict) -> str:
        gene_centric = data.get('gene_centric', {})
        patterns = data.get('patterns', {})
        coverage = patterns.get('database_coverage', {})
        stats = gene_centric.get('database_stats', {})
        all_dbs = sorted(set(stats.keys()) | set(coverage.keys()))
        rows = []
        for db in all_dbs:
            s = stats.get(db, {})
            c = coverage.get(db, {})
            total_samples = c.get('total_samples', 0)
            samples_with_hits = c.get('samples_with_hits', s.get('total_occurrences', 0))
            cov_pct = (samples_with_hits / total_samples * 100) if total_samples else 0
            rows.append(
                f'<tr><td><strong>{db.upper()}</strong></td>'
                f'<td>{s.get("total_genes", 0)}</td>'
                f'<td>{s.get("total_occurrences", 0)}</td>'
                f'<td>{samples_with_hits} / {total_samples}</td>'
                f'<td>{cov_pct:.1f}%</td></tr>'
            )
        return f'''
        {self._section_header('Database Coverage', 'fa-database', 'databases')}
        {self._alert('info', 'fa-info-circle',
            '<strong>Database coverage</strong> — total unique genes, total occurrences, and the number of samples '
            'with at least one hit per database. This table helps evaluate the completeness of your screening. '
            'A database with zero hits may indicate: (1) genuinely absent resistance determinants in your population, '
            '(2) database not properly installed/updated, or (3) hits below the identity/coverage thresholds. '
            'Cross-check the raw TSV outputs if a database unexpectedly returns zero hits.')}
        <div class="master-scrollable-container">
            <table class="data-table">
                <thead><tr>
                    <th>Database</th><th>Unique Genes</th><th>Total Occurrences</th>
                    <th>Samples with Hits</th><th>Coverage (%)</th>
                </tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>'''

    def _generate_aiguide_section(self, data: Dict) -> str:
        return f'''
        {self._section_header('AI Assistant Guide', 'fa-robot', 'aiguide')}
        <div class="alert-box alert-info">
            <i class="fas fa-robot fa-2x"></i>
            <div>
                <h3>🤖 AI Assistant Guide — Unleash the Power of AI for Genomic Epidemiology</h3>
                <p>Use large language models (LLMs) like <strong>ChatGPT, Claude, or Gemini</strong> to interact with your <em>K. pneumoniae</em> dataset — turning static reports into dynamic conversations that uncover patterns you might otherwise miss.</p>
            </div>
        </div>
        <div style="margin:20px 0;">
            <div class="ai-guide-card">
                <h4 class="ai-guide-title"><i class="fas fa-brain"></i> Why Use AI for Genomic Data Analysis?</h4>
                <ul class="ai-guide-list">
                    <li><strong>Pattern recognition</strong> — spot epidemiological trends, clone associations, and co-occurrence networks across thousands of samples.</li>
                    <li><strong>Natural-language queries</strong> — ask in plain English, get instant answers without writing a single line of SQL or Python.</li>
                    <li><strong>Hypothesis generation</strong> — uncover unexpected correlations that merit experimental follow-up or targeted surveillance.</li>
                    <li><strong>Literature synthesis</strong> — connect your findings with published resistance mechanisms, outbreak reports, and clinical guidelines.</li>
                </ul>
                <p class="ai-guide-note-warning"><i class="fas fa-lightbulb"></i> <strong>Scientific note:</strong> AI is <em>pattern-finding</em>, not <em>causal-inferring</em>. Use it to suggest hypotheses, then verify with wet-lab or clinical correlation before drawing conclusions.</p>
            </div>

            <div class="ai-guide-card" style="margin-top:16px;">
                <h4 class="ai-guide-title"><i class="fas fa-upload"></i> How to Feed This Report to AI</h4>
                <ol class="ai-guide-list">
                    <li><strong>Upload the JSON file</strong> — <code class="ai-guide-code">kleboscope_ultimate_gene_centric_report.json</code> contains all structured data. Upload to ChatGPT (Advanced Data Analysis), Claude, or Gemini. <em>Best for quantitative queries.</em></li>
                    <li><strong>Upload the HTML report</strong> — modern AI tools parse HTML tables natively. Upload the <code class="ai-guide-code">.html</code> file directly. <em>Great for visual context.</em></li>
                    <li><strong>Copy-paste specific tables</strong> — if you only need a quick insight, paste a table directly into the chat. <em>Instant, no file upload needed.</em></li>
                </ol>
                <p class="ai-guide-note-info"><i class="fas fa-info-circle"></i> <strong>Pro tip:</strong> tell the AI: <em>"You are a bioinformatician analysing K. pneumoniae genomes. The attached data contains typing (MLST, Kaptive K/O), AMR (carbapenemases, ESBLs, colistin, tigecycline), virulence (ICEKp, hypervirulence siderophores), BACMET, mutation, and MGE information. Answer my questions with references to the data."</em></p>
            </div>

            <div class="ai-guide-card" style="margin-top:16px;">
                <h4 class="ai-guide-title"><i class="fas fa-chart-line"></i> Scientifically Relevant Questions to Ask</h4>
                <div style="display:grid;grid-template-columns:1fr 1fr;gap:14px;">
                    <div class="ai-guide-subcard">
                        <strong class="ai-guide-subtitle">🧬 Epidemiology &amp; Clonality</strong>
                        <ul class="ai-guide-sublist">
                            <li>What are the most common MLST sequence types in this dataset?</li>
                            <li>Which K-types are dominant in carbapenemase-positive vs negative isolates?</li>
                            <li>Are there any ST-K-O combinations with more than two isolates (possible clonal expansion)?</li>
                            <li>What is the K-locus distribution, and does it correlate with hypervirulence markers?</li>
                            <li>Which clonal complexes carry the highest burden of resistance genes?</li>
                        </ul>
                    </div>
                    <div class="ai-guide-subcard">
                        <strong class="ai-guide-subtitle">💊 Antimicrobial Resistance</strong>
                        <ul class="ai-guide-sublist">
                            <li>How many samples carry <em>blaKPC</em>, <em>blaNDM</em>, or <em>blaOXA-48</em>? What are their STs and K-types?</li>
                            <li>Are there any <em>mcr</em>-positive samples? What is their colistin resistance mechanism (plasmid mcr or chromosomal mgrB/pmrAB)?</li>
                            <li>Which AMR genes co-occur most frequently, and what does that suggest about plasmid linkage?</li>
                            <li>What is the distribution of tetracycline (<em>tet</em>) resistance genes across lineages?</li>
                            <li>Do any samples carry combined carbapenemase + ESBL + colistin resistance (XDR phenotype)?</li>
                        </ul>
                    </div>
                    <div class="ai-guide-subcard">
                        <strong class="ai-guide-subtitle">🦠 Virulence &amp; Hypervirulence</strong>
                        <ul class="ai-guide-sublist">
                            <li>Which samples carry <em>rmpA</em> and <em>rmpA2</em>? Are they associated with specific STs (ST23, ST86, ST375, ST380)?</li>
                            <li>List all samples with complete aerobactin (<em>iucA-iucD</em>) and salmochelin (<em>iroB-iroN</em>) systems.</li>
                            <li>Which ICEKp markers are present, and which lineages carry the ICEKp10 (ybt + clb) signature?</li>
                            <li>Is there a correlation between K-locus (K1/K2 vs others) and hypervirulence markers?</li>
                            <li>Do any isolates combine carbapenemase genes AND hypervirulence markers (CR-hvKp)?</li>
                        </ul>
                    </div>
                    <div class="ai-guide-subcard">
                        <strong class="ai-guide-subtitle">🧪 Mutations, MGE &amp; Biocides</strong>
                        <ul class="ai-guide-sublist">
                            <li>What are the most frequent point mutations in <em>gyrA</em> or <em>parC</em>, and which clones carry them?</li>
                            <li>Are there any <em>mgrB</em> or <em>pmrAB</em> mutations associated with colistin resistance?</li>
                            <li>Which samples carry <em>qac</em> genes (quaternary ammonium resistance)? Do they co-occur with antibiotic resistance?</li>
                            <li>Is <em>mer</em> (mercury) or <em>sil</em> (silver) resistance linked to specific STs or plasmid families?</li>
                            <li>Which isolates have the highest mobileOG hit counts, and what does that suggest about their mobilome?</li>
                        </ul>
                    </div>
                </div>
            </div>

            <div class="ai-guide-card" style="margin-top:16px;">
                <h4 class="ai-guide-title"><i class="fas fa-balance-scale"></i> Scientific Rigour &amp; Ethical AI Use</h4>
                <ul class="ai-guide-list">
                    <li><strong>AI is your co-pilot, not the pilot.</strong> Interpret AI-generated insights in the context of local epidemiology, clinical guidelines, and laboratory validation.</li>
                    <li><strong>Verify, verify, verify.</strong> Cross-check critical calls with primary literature, genome browsers, or secondary tools. A gene call from one database is weaker evidence than cross-database agreement.</li>
                    <li><strong>No patient-identifiable data.</strong> Only upload aggregated, de-identified genomic data. Public reference genomes are fine; clinical isolates with patient metadata are not.</li>
                    <li><strong>Transparency in publications.</strong> Mention AI-assisted pattern discovery in your methods, followed by manual curation of the key findings.</li>
                    <li><strong>AI hallucination is real.</strong> If the AI confidently tells you that <em>mecA</em> is found in <em>Klebsiella</em>, or that TSST-1 causes food poisoning in a carbapenem-resistant strain — <strong>don't believe it</strong>. Treat every AI statement as a hypothesis to be verified.</li>
                </ul>
            </div>

            <div class="ai-guide-card-warning" style="margin-top:16px;">
                <h4 class="ai-guide-title-warning"><i class="fas fa-smile-wink"></i> A (Mostly Serious) AI Survival Guide</h4>
                <ul class="ai-guide-list-warning">
                    <li><strong>If the AI says "I don't know"</strong> — trust it. It's being honest.</li>
                    <li><strong>If the AI says "It is widely known"</strong> — ask for a reference. It may have made it up.</li>
                    <li><strong>If the AI offers a completely novel evolutionary theory for Klebsiella evolution in your specific hospital</strong> — check if your coffee is spiked.</li>
                    <li><strong>Remember:</strong> AI won't take your job — but a bioinformatician who knows how to use AI might! Learn it, use it, and always keep a healthy dose of scepticism. 😉</li>
                </ul>
            </div>
        </div>'''
    def _generate_calltoaction_section(self, data: Dict) -> str:
        """Call to Action — full  narrative, Klebsiella tailored."""
        return f'''
        {self._section_header('Call to Action — AMR & the ESCAPE Initiative', 'fa-globe', 'calltoaction')}
        <div class="alert-box alert-info">
            <i class="fas fa-globe fa-2x"></i>
            <div>
                <h3>The Global Burden of AMR and Our Call to Action</h3>
                <p>Antimicrobial resistance (AMR) is one of the top global public health threats, with an estimated <strong>1.27 million direct deaths annually</strong> and nearly 5 million associated deaths. <em>Klebsiella pneumoniae</em> is a leading cause of hospital-acquired infections — bloodstream infections, ventilator-associated pneumonia, urinary tract infections — and is now the single most common carbapenem-resistant Enterobacterales (CRE) worldwide. Tracking carbapenemases, ESBLs, colistin resistance, and hypervirulence determinants is essential to inform treatment guidelines, infection control, and vaccine development.</p>
                <p>We developed <strong>Kleboscope</strong> to empower researchers and clinicians — especially in low-resource settings — to analyse their own sequencing data without requiring extensive bioinformatics expertise. Every module is freely available, openly documented, and designed to run on a laptop.</p>
            </div>
        </div>
        <div style="background:#e8f5e9;padding:20px;border-radius:12px;margin:20px 0;color:#1f2937;">
            <h3 style="color:#005a8c;"><i class="fas fa-bacterium"></i> ESCAPE AMR — Our Ongoing Project (ESKAPE Pathogens)</h3>
            <p><strong>Kleboscope</strong> is one of the first modules of a larger initiative called <strong>ESCAPE AMR</strong> (formerly ESKAPE). We target the notorious <strong>ESKAPE pathogens</strong> — the six bacterial species that "escape" the effects of most antibiotics and account for the majority of hospital-acquired infections:</p>
            <ul style="margin:10px 0 10px 20px;line-height:1.8;">
                <li><strong>E</strong>nterococcus faecium</li>
                <li><strong>S</strong>taphylococcus aureus</li>
                <li><strong>K</strong>lebsiella pneumoniae <em></em></li>
                <li><strong>A</strong>cinetobacter baumannii</li>
                <li><strong>P</strong>seudomonas aeruginosa</li>
                <li><strong>E</strong>nterobacter species</li>
            </ul>
            <p>These bacteria "escape" the effects of antibiotics — hence the name. But we also see the name as a global call to action:</p>
            <div style="background:#fff3e0;padding:15px;border-radius:8px;margin:15px 0;">
                <p><strong>🔹 E</strong>veryone must join forces — researchers, clinicians, policymakers, and citizens.<br>
                <strong>🔹 S</strong>mart surveillance is our first line of defence. No more guessing — we need genomic data.<br>
                <strong>🔹 K</strong>nowledge must be shared openly. No paywalls, no closed silos.<br>
                <strong>🔹 A</strong>frica bears a heavy AMR burden, but African solutions are already emerging.<br>
                <strong>🔹 P</strong>revention is cheaper than cure. Let's stop resistant infections before they spread.<br>
                <strong>🔹 E</strong>very day we delay, more lives are at stake. The time to act is now, not tomorrow.</p>
            </div>
            <p><i class="fas fa-laugh-squint"></i> <strong>"We didn't choose the name ESKAPE because it sounds cool (though it does). We chose it because it reminds us every single day: we must ESCAPE the AMR crisis — together, urgently, and with the best science we have."</strong><br>— Brown Beckley, lead developer (who secretly hopes this pun makes you smile, not roll your eyes 😉)</p>
        </div>
        <div style="text-align:center;margin:40px 0;color:#1f2937;">
            <i class="fas fa-star" style="font-size:3em;color:#ffc107;"></i>
            <h3 style="color:#005a8c;margin:16px 0;">🤝 We Invite You to Contribute!</h3>
            <p><strong>If you find this tool useful, please:</strong></p>
            <div class="action-buttons" style="justify-content:center;">
                <a href="https://github.com/bbeckley-hub" target="_blank" class="action-btn btn-primary" style="text-decoration:none;"><i class="fab fa-github"></i> ⭐ Star us on GitHub</a>
                <a href="mailto:brownbeckley94@gmail.com" class="action-btn btn-success" style="text-decoration:none;"><i class="fas fa-envelope"></i> Contact the team</a>
                <a href="https://github.com/bbeckley-hub" target="_blank" class="action-btn btn-warning" style="text-decoration:none;"><i class="fas fa-bug"></i> Report issues</a>
            </div>
            <p style="margin-top:20px;"><i class="fas fa-chalkboard-user"></i> <strong>We welcome collaborations</strong> to adapt this tool for other pathogens and to improve AMR surveillance in Africa and beyond.</p>
            <p><i class="fas fa-hand-holding-heart"></i> If you are a funder or organisation interested in supporting the <strong>ESCAPE AMR</strong> project, please reach out. Together we can build a free, open-source ecosystem for genomic surveillance of all ESKAPE pathogens.</p>
        </div>
        <div style="background:#f8f9fa;padding:15px;border-radius:8px;color:#1f2937;">
            <i class="fas fa-quote-left"></i> "AMR is a silent pandemic, but we have the tools to fight it — if we share them, if we teach each other, and if we act with urgency. Let's escape the era of untreatable infections."<br>
            <strong>— The ESCAPE AMR Team, University of Ghana Medical School</strong>
        </div>'''

    def _generate_citation_section(self, data: Dict) -> str:
        """Citation — full  accordion, white/black cards, suggested acknowledgement."""
        def render_citation(name, text, url, plain):
            link_btn = (
                f'<a class="citation-link" href="{url}" target="_blank">🔗 Open ↗</a>'
                if url else ''
            )
            return f'''
            <li class="citation-item">
                <div class="citation-body">
                    <div><strong class="citation-name" style="color:#005a8c;">{name}</strong> <span style="color:#333;">– {text}</span></div>
                    <div class="citation-actions">
                        {link_btn}
                        <button class="copy-btn" data-citation="{plain}">📋 Copy</button>
                    </div>
                </div>
            </li>'''

        main_html = render_citation(
            'Kleboscope',
            'Beckley B, Amarh V. Kleboscope: a species-optimized computational pipeline for rapid and accessible <em>Klebsiella pneumoniae</em> genotyping and surveillance. <em>Github</em>. 2026.',
            'https://github.com/bbeckley-hub',
            'Beckley B, Amarh V. Kleboscope: a species-optimized computational pipeline for rapid and accessible Klebsiella pneumoniae genotyping and surveillance. Bioinformatics. 2026 (in preparation).',
        )

        deps_citations = [
            ('MLST', 'Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018.',
             'https://github.com/tseemann/mlst',
             'Seemann T. MLST: Scan contig files against PubMLST typing schemes. GitHub. 2018.'),
            ('PubMLST / BIGSdb', 'Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. <em>Wellcome Open Res</em>. 2018;3:124.',
             'https://doi.org/10.12688/wellcomeopenres.14826.1',
             'Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res. 2018;3:124.'),
            ('Kaptive', 'Wyres KL, Wick RR, Gorrie C, Jenney A, Follador R, Thomson NR, Holt KE. Identification of Klebsiella capsule synthesis loci from whole genome sequence data. <em>Microb Genom</em>. 2016;2(5):e000064.',
             'https://doi.org/10.1099/mgen.0.000064',
             'Wyres KL, Wick RR, Gorrie C, Jenney A, Follador R, Thomson NR, Holt KE. Identification of Klebsiella capsule synthesis loci from whole genome sequence data. Microb Genom. 2016;2(5):e000064.'),
            ('fastANI', 'Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. <em>Nat Commun</em>. 2018;9(1):5114.',
             'https://doi.org/10.1038/s41467-018-07641-9',
             'Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun. 2018;9(1):5114.'),
            ('AMRFinderPlus', 'Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. <em>Sci Rep</em>. 2021;11(1):12728.',
             'https://doi.org/10.1038/s41598-021-91456-0',
             'Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728.'),
            ('ABRicate', 'Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024.',
             'https://github.com/tseemann/abricate',
             'Seemann T. ABRicate: mass screening of contigs for antibiotic resistance genes. GitHub. 2024.'),
            ('CARD', 'McArthur AG, Waglechner N, Nizam F, Yan A, Azad MA, Baylay AJ, Bhullar K, Canova MJ, De Pascale G, Ejim L, Kalan L, King AM, Koteva K, Morar M, Mulvey MR, O\'Brien JS, Pagliarulo NP, Piddock LJV, Spanogiannopoulos P, Sutherland AD, Tang I, Taylor PL, Thaker M, Wang W, Yan M, Yu T, Wright GD. The comprehensive antibiotic resistance database. <em>Antimicrob Agents Chemother</em>. 2013;57(7):3348-57.',
             'https://doi.org/10.1128/AAC.00419-13',
             'McArthur AG, Waglechner N, Nizam F, Yan A, Azad MA, Baylay AJ, Bhullar K, Canova MJ, De Pascale G, Ejim L, Kalan L, King AM, Koteva K, Morar M, Mulvey MR, O\'Brien JS, Pagliarulo NP, Piddock LJV, Spanogiannopoulos P, Sutherland AD, Tang I, Taylor PL, Thaker M, Wang W, Yan M, Yu T, Wright GD. The comprehensive antibiotic resistance database. Antimicrob Agents Chemother. 2013;57(7):3348-57.'),
            ('ResFinder', 'Florensa AF, Kaas RS, Clausen PTLC, Aytan-Aktug D, Aarestrup FM. ResFinder - an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes. <em>Microb Genom</em>. 2022;8(1):000748.',
             'https://doi.org/10.1099/mgen.0.000748',
             'Florensa AF, Kaas RS, Clausen PTLC, Aytan-Aktug D, Aarestrup FM. ResFinder - an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes. Microb Genom. 2022;8(1):000748.'),
            ('VFDB', 'Chen L, Xiong Z, Sun L, Yang J, Jin Q. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. <em>Nucleic Acids Res</em>. 2012;40(D1):D641-5.',
             'https://doi.org/10.1093/nar/gkr989',
             'Chen L, Xiong Z, Sun L, Yang J, Jin Q. VFDB 2012 update: toward the genetic diversity and molecular evolution of bacterial virulence factors. Nucleic Acids Res. 2012;40(D1):D641-5.'),
            ('PlasmidFinder', 'Carattoli A, Zankari E, García-Fernández A, Voldby Larsen M, Lund O, Villa L, Møller Aarestrup F, Hasman H. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. <em>Antimicrob Agents Chemother</em>. 2014;58(7):3895-903.',
             'https://doi.org/10.1128/AAC.02412-14',
             'Carattoli A, Zankari E, García-Fernández A, Voldby Larsen M, Lund O, Villa L, Møller Aarestrup F, Hasman H. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrob Agents Chemother. 2014;58(7):3895-903.'),
            ('BacMet', 'Pal C, Bengtsson-Palme J, Rensing C, Kristiansson E, Larsson DGJ. BacMet: antibacterial biocide and metal resistance genes database. <em>Nucleic Acids Res</em>. 2014;42(D1):D737-43.',
             'https://doi.org/10.1093/nar/gkt1252',
             'Pal C, Bengtsson-Palme J, Rensing C, Kristiansson E, Larsson DGJ. BacMet: antibacterial biocide and metal resistance genes database. Nucleic Acids Res. 2014;42(D1):D737-43.'),
            ('MEGARes 2.0', 'Doster E, Lakin SM, Dean CJ, Wolfe C, Young JG, Boucher C, Belk KE, Noyes NR, Morley PS. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. <em>Nucleic Acids Res</em>. 2020;48(D1):D561-D569.',
             'https://doi.org/10.1093/nar/gkz1010',
             'Doster E, Lakin SM, Dean CJ, Wolfe C, Young JG, Boucher C, Belk KE, Noyes NR, Morley PS. MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. Nucleic Acids Res. 2020;48(D1):D561-D569.'),
            ('ARG-ANNOT', 'Gupta SK, Padmanabhan BR, Diene SM, Lopez-Rojas R, Kempf M, Landraud L, Rolain JM. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. <em>Antimicrob Agents Chemother</em>. 2014;58(1):212-20.',
             'https://doi.org/10.1128/AAC.01310-13',
             'Gupta SK, Padmanabhan BR, Diene SM, Lopez-Rojas R, Kempf M, Landraud L, Rolain JM. ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrob Agents Chemother. 2014;58(1):212-20.'),
            ('mobileOG-db', 'Brown CL, Mullet J, Hindi F, Stoll JE, Gupta S, Choi M, Keenum I, Vikesland P, Pruden A, Zhang L. mobileOG-db: a manually curated database of protein families mediating the life cycle of bacterial mobile genetic elements. <em>Appl Environ Microbiol</em>. 2022;88(18):e00991-22.',
             'https://doi.org/10.1128/aem.00991-22',
             'Brown CL, Mullet J, Hindi F, Stoll JE, Gupta S, Choi M, Keenum I, Vikesland P, Pruden A, Zhang L. mobileOG-db: a manually curated database of protein families mediating the life cycle of bacterial mobile genetic elements. Appl Environ Microbiol. 2022;88(18):e00991-22.'),
            ('Biopython', 'Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL. Biopython: freely available Python tools for computational molecular biology and bioinformatics. <em>Bioinformatics</em>. 2009;25(11):1422-3.',
             'https://doi.org/10.1093/bioinformatics/btp163',
             'Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3.'),
        ]
        deps_html = ''.join(render_citation(name, text, url, plain) for name, text, url, plain in deps_citations)

        return f'''
        {self._section_header('Citations', 'fa-book', 'citation')}
        <div class="alert-box alert-info">
            <i class="fas fa-quote-right fa-2x"></i>
            <div>
                <h3>📚 How to Cite Kleboscope and Its Dependencies</h3>
                <p>If you use Kleboscope in your research, please cite the main tool and the relevant third-party tools and databases. Each entry below links directly to its DOI or source repository, and the copy button gives you the exact citation text ready to paste into your manuscript.</p>
            </div>
        </div>
        <div class="accordion">
            <div class="accordion-item">
                <div class="accordion-header">
                    <span>📄 From the ESCAPE AMR Platform</span>
                    <i class="fas fa-chevron-down"></i>
                </div>
                <div class="accordion-content" style="display:block;">
                    <ul style="padding-left:0;margin:0;">{main_html}</ul>
                </div>
            </div>
            <div class="accordion-item">
                <div class="accordion-header">
                    <span>🔧 Key Databases &amp; Methods</span>
                    <i class="fas fa-chevron-down"></i>
                </div>
                <div class="accordion-content" style="display:block;">
                    <ul style="padding-left:0;margin:0;">{deps_html}</ul>
                </div>
            </div>
        </div>
        <div class="alert-box alert-success">
            <i class="fas fa-hand-peace fa-2x"></i>
            <div>
                <strong>Suggested acknowledgement:</strong><br>
                "Genomic analysis was performed using <strong>Kleboscope</strong> [Beckley &amp; Amarh, 2026], which integrates MLST [Seemann, 2018] using the PubMLST database [Jolley et al., 2018], Kaptive [Wyres et al., 2016] for K and O locus typing, ABRicate [Seemann, 2024], and AMRFinderPlus [Feldgarden et al., 2021] for comprehensive <em>Klebsiella pneumoniae</em> characterization. Species confirmation was performed with fastANI [Jain et al., 2018]. Antimicrobial resistance genes were identified using the CARD [McArthur et al., 2013], ResFinder [Florensa et al., 2022], MEGARes [Doster et al., 2020], and ARG-ANNOT [Gupta et al., 2014] databases. For biocide and heavy metal resistance genes, BacMet [Pal et al., 2014] was used. Virulence and plasmid screening were performed with ABRicate using the VFDB [Chen et al., 2012] and PlasmidFinder [Carattoli et al., 2014] databases. Mutation detection was performed using AMRFinderPlus. Mobile genetic element profiling used mobileOG-db [Brown et al., 2022]. FASTA QC and assembly statistics were computed with Biopython [Cock et al., 2009]."
            </div>
        </div>'''

    def _generate_funding_section(self, data: Dict) -> str:
        """Funding statement — full  text, Klebsiella tailored."""
        return f'''
        {self._section_header('Funding & Support', 'fa-coffee', 'funding')}
        <div class="alert-box alert-info">
            <i class="fas fa-coffee fa-2x"></i>
            <div>
                <h3>☕ Funding &amp; Support — Keeping the Lights On (with code and caffeine)</h3>
                <p>Kleboscope is an <strong>independent, unfunded project</strong> born out of passion for genomic surveillance and AMR research at the University of Ghana Medical School.</p>
                <p>No grants, no sponsors, no institutional backing — just a laptop, a lot of coffee, and a burning desire to help researchers fight antimicrobial resistance.</p>
            </div>
        </div>
        <div class="alert-box alert-warning">
            <i class="fas fa-heart fa-2x"></i>
            <div>
                <h3>💡 How You Can Help (Without Opening Your Wallet)</h3>
                <ul>
                    <li><strong>⭐ Star us on GitHub</strong> – It takes two seconds and makes us feel like rockstars.</li>
                    <li><strong>🐛 Report bugs</strong> – If something breaks, let us know. We'll fix it with joy.</li>
                    <li><strong>💡 Suggest features</strong> – Have an idea? We're all ears (and we actually implement them).</li>
                    <li><strong>🧬 Share your data</strong> – If you've used Kleboscope and want to collaborate, we'd love to hear your story.</li>
                    <li><strong>📢 Spread the word</strong> – Tell your colleagues, tweet about it, or mention it in your next Zoom call.</li>
                    <li><strong>👋 Say hello!</strong> – Just drop an email to <strong>brownbeckley94@gmail.com</strong>. It makes our day.</li>
                </ul>
                <p><i class="fas fa-microbe"></i> <strong>Fun fact:</strong> This project runs on 100% volunteer tears, 0% grant money. But we're not bitter — we're just caffeinated.</p>
            </div>
        </div>
        <div class="alert-box alert-success">
            <i class="fas fa-hand-holding-heart fa-2x"></i>
            <div>
                <h3>🤝 Contribute to the ESCAPE AMR Platform</h3>
                <p>We also maintain pipelines for other ESKAPE pathogens (AcinetoScope, , Pseudoscope, etc.). If you're a developer, bioinformatician, or just someone who loves clean code and bacteria, we welcome pull requests, issues, documentation improvements, and ideas for new databases.</p>
                <p>Visit our GitHub: <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a> — star, fork, and let's fight AMR together!</p>
                <p><strong>Brown Beckley</strong> — <i class="fas fa-envelope"></i> brownbeckley94@gmail.com</p>
                <p><i class="fas fa-laugh-beam"></i> <strong>P.S.</strong> If you ever meet Brown in person, buy him a coffee. He'll probably talk your ear off about carbapenemase variants, but it's worth it.</p>
            </div>
        </div>'''

    def _generate_export_section(self, data: Dict) -> str:
        exports = [
            ('samples-table', 'Sample Overview', 'sample_overview.csv', 'fa-table'),
            ('amr-table', 'AMR Genes', 'amr_genes.csv', 'fa-biohazard'),
            ('vir-table', 'Virulence Genes', 'virulence_genes.csv', 'fa-virus'),
            ('bac-table', 'BACMET Genes', 'bacmet_genes.csv', 'fa-flask'),
            ('plasmids-table', 'Plasmids', 'plasmid_replicons.csv', 'fa-dna'),
            ('mutations-table', 'Mutations', 'mutations.csv', 'fa-dna'),
            ('qc-table', 'FASTA QC', 'fasta_qc.csv', 'fa-chart-line'),
            ('mge-table', 'MGE Profile', 'mge_profile.csv', 'fa-mobile-alt'),
        ]
        cards = ''
        for tid, label, fname, icon in exports:
            cards += f'''
            <div class="dashboard-card" onclick="exportTableToCSV('{tid}', '{fname}')">
                <div class="card-icon"><i class="fas {icon}"></i></div>
                <div class="card-label">{label}</div>
            </div>'''
        cards += '''
            <div class="dashboard-card" onclick="window.location.href='kleboscope_ultimate_gene_centric_report.json'">
                <div class="card-icon"><i class="fas fa-file-code"></i></div>
                <div class="card-label">Complete JSON</div>
            </div>'''
        return f'''
        {self._section_header('Export Data', 'fa-download', 'export')}
        {self._alert('info', 'fa-info-circle',
            'Click any card to download that table as a CSV file. The complete JSON dataset contains '
            'every detail from the report and is the best option for feeding data into downstream '
            'analysis tools, R/Python scripts, or AI assistants.')}
        <div class="dashboard-grid" style="margin-top:16px;">{cards}</div>'''

    # ------------------------------------------------------------------ assembly
    def generate_main_report(self, integrated_data: Dict, output_dir: Path) -> str:
        print("\nGenerating Kleboscope Ultimate HTML report...")
        self._samples_data = integrated_data.get('samples', {})
        html = self._build_html(integrated_data)
        output_file = output_dir / "kleboscope_ultimate_gene_centric_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    HTML report saved: {output_file}")
        return str(output_file)

    def _build_html(self, data: Dict) -> str:
        metadata = data.get('metadata', {})
        samples = data.get('samples', {})
        patterns = data.get('patterns', {})
        gene_centric = data.get('gene_centric', {})
        mge_data = data.get('mge_data', {})
        total_samples = len(samples)

        sample_typing = {}
        for sample, d in samples.items():
            sample_typing[sample] = {
                'ST': d.get('mlst', {}).get('ST', 'ND'),
                'K': d.get('kaptive', {}).get('K_Locus', 'ND'),
                'O': d.get('kaptive', {}).get('O_Locus', 'ND'),
            }

        css = self._get_css()
        js = self._get_js(json.dumps(sample_typing))

        tabs = [
            ('summary', 'Summary', 'fa-chart-pie'),
            ('sample_overview', 'Samples', 'fa-list'),
            ('mlst', 'MLST', 'fa-code-branch'),
            ('qc', 'QC', 'fa-chart-line'),
            ('kaptive', 'Kaptive', 'fa-shield-alt'),
            ('combinations', 'Combinations', 'fa-link'),
            ('amr', 'AMR', 'fa-biohazard'),
            ('virulence', 'Virulence', 'fa-virus'),
            ('bacmet', 'Bacmet', 'fa-flask'),
            ('plasmids', 'Plasmids', 'fa-dna'),
            ('mutations', 'Mutations', 'fa-dna'),
            ('patterns', 'Patterns', 'fa-project-diagram'),
            ('highrisk', 'High Risk', 'fa-exclamation-triangle'),
            ('mge', 'MGE', 'fa-mobile-alt'),
            ('compare', 'Compare', 'fa-balance-scale'),
            ('databases', 'Databases', 'fa-database'),
            ('aiguide', 'AI Guide', 'fa-robot'),
            ('calltoaction', 'Call to Action', 'fa-globe'),
            ('citation', 'Citation', 'fa-book'),
            ('funding', 'Funding', 'fa-coffee'),
            ('export', 'Export', 'fa-download'),
        ]

        nav_html = ''
        for i, (name, label, icon) in enumerate(tabs):
            active = ' active' if i == 0 else ''
            nav_html += (
                f'<button class="tab-button {name}{active}" '
                f'onclick="switchTab(\'{name}\')">'
                f'<i class="fas {icon}"></i> {label}</button>'
            )

        sections = {
            'summary': self._generate_summary_section,
            'sample_overview': self._generate_sample_overview_section,
            'mlst': self._generate_mlst_section,
            'qc': self._generate_qc_section,
            'kaptive': self._generate_kaptive_section,
            'combinations': self._generate_combinations_section,
            'amr': self._generate_amr_section,
            'virulence': self._generate_virulence_section,
            'bacmet': self._generate_bacmet_section,
            'plasmids': self._generate_plasmids_section,
            'mutations': self._generate_mutations_section,
            'patterns': self._generate_patterns_section,
            'highrisk': self._generate_highrisk_section,
            'mge': self._generate_mge_section,
            'compare': self._generate_compare_section,
            'databases': self._generate_databases_section,
            'aiguide': self._generate_aiguide_section,
            'calltoaction': self._generate_calltoaction_section,
            'citation': self._generate_citation_section,
            'funding': self._generate_funding_section,
            'export': self._generate_export_section,
        }

        tabs_html = ''
        for i, (name, label, icon) in enumerate(tabs):
            active = ' active' if i == 0 else ''
            content = sections[name](data)
            tabs_html += f'<div id="{name}-tab" class="tab-content{active}">{content}</div>'

        all_genes = gene_centric.get('all_genes', [])
        mge_n = len(mge_data.get('per_sample', {}))
        carbapenemase_count = sum(
            1 for g in all_genes
            if any(k in g['gene'] for k in ('blaKPC', 'blaNDM', 'blaOXA-48', 'blaVIM', 'blaIMP'))
        )
        dashboard_cards = [
            ('summary', 'fa-vial', total_samples, 'Samples'),
            ('mlst', 'fa-code-branch', len(patterns.get('st_distribution', {})), 'Unique STs'),
            ('kaptive', 'fa-shield-alt', len(patterns.get('k_locus_distribution', {})), 'K-types'),
            ('kaptive', 'fa-dna', len(patterns.get('o_locus_distribution', {})), 'O-types'),
            ('amr', 'fa-biohazard', len(all_genes), 'Genes'),
            ('amr', 'fa-biohazard', carbapenemase_count, 'Carbapenemases'),
            ('mge', 'fa-mobile-alt', mge_n, 'MGE Samples'),
            ('highrisk', 'fa-exclamation-triangle', len(patterns.get('high_risk_combinations', [])), 'High-Risk Combos'),
        ]
        dash_html = '<div class="dashboard-grid">'
        for tab, icon, value, label in dashboard_cards:
            dash_html += (
                f'<div class="dashboard-card card-{tab}" onclick="switchTab(\'{tab}\')">'
                f'<div class="card-icon"><i class="fas {icon}"></i></div>'
                f'<div class="card-number">{value}</div>'
                f'<div class="card-label">{label}</div></div>'
            )
        dash_html += '</div>'

        html = f'''<!DOCTYPE html>
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
        <h1><i class="fas fa-bacterium"></i> Kleboscope Ultimate Report</h1>
        <p>Gene-centric cross-genome analysis for <em>Klebsiella pneumoniae</em></p>s
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar-alt"></i> {metadata.get('analysis_date', 'Unknown')}</div>
            <div class="metadata-item"><i class="fas fa-database"></i> {total_samples} Samples</div>
            <div class="metadata-item"><i class="fas fa-university"></i> University of Ghana Medical School</div>
        </div>
    </div>
    {dash_html}
    <div class="tab-navigation">{nav_html}</div>
    {tabs_html}
    <div class="footer">
        <h3>Kleboscope Ultimate Reporter v3.0.0</h3>
        <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
        <p>Generated on {metadata.get('analysis_date', 'Unknown')}</p>
        <p><strong>Critical Genes Tracked:</strong> Carbapenemases (KPC, NDM, OXA-48, VIM, IMP, GES) • Colistin (mcr-1 to mcr-10, mgrB, pmrAB) • Tigecycline (tetX) • ICEKp Markers (ybt, clb, iro, rmp) • Virulence Plasmid Markers (iro, iuc, rmp, rmpA2) • Biocides & Heavy Metals (qac, sil, mer, ars, pco) • Adhesins (fim, mrk, ecp) • Secretion Systems (tss) • Siderophores (aerobactin, salmochelin, yersiniabactin) • Toxins (hly, cnf, cdt, sat, pic)</p>
        <p>If you find this useful, please <a href="https://github.com/bbeckley-hub" target="_blank">⭐ star us on GitHub</a>.</p>
    </div>
</div>
</body>
</html>'''
        return html


# =============================================================================
# REPORTER
# =============================================================================
class KleboscopeUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "KLEBOSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = KleboHTMLParser()
        self.analyzer = KleboDataAnalyzer()
        self.generator = KleboHTMLGenerator(self.analyzer)
        self.metadata = {
            'tool_name': 'Kleboscope Ultimate Reporter',
            'version': '2.3.0',
            'author': 'Brown Beckley',
            'affiliation': 'University of Ghana Medical School',
            'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'input_directory': str(self.input_dir),
        }

    def find_files(self) -> Dict[str, Any]:
        found = {
            'mlst_tsv': None, 'kaptive_tsv': None, 'qc_tsv': None,
            'mge_html': None, 'amrfinder_html': None,
            'abricate_html': defaultdict(list), 'mutation_html': None,
        }
        for f in self.input_dir.glob("**/*.tsv"):
            name = f.name.lower()
            if 'mlst' in name and not found['mlst_tsv']:
                found['mlst_tsv'] = f
            elif 'kaptive' in name and not found['kaptive_tsv']:
                found['kaptive_tsv'] = f
            elif ('fasta_qc' in name or 'qc_summary' in name) and not found['qc_tsv']:
                found['qc_tsv'] = f
        for f in self.input_dir.glob("**/*.html"):
            name = f.name.lower()
            if 'mge' in name and 'summary' in name and not found['mge_html']:
                found['mge_html'] = f
            elif 'amrfinder_summary_report' in name and not found['amrfinder_html']:
                found['amrfinder_html'] = f
            elif 'mutation_summary' in name and not found['mutation_html']:
                found['mutation_html'] = f
            else:
                for db in self.parser.abricate_databases:
                    if f'klebo_{db}_summary_report' in name:
                        found['abricate_html'][db].append(f)
                        break
        print(f"  MLST TSV:     {found['mlst_tsv'].name if found['mlst_tsv'] else 'NOT FOUND'}")
        print(f"  Kaptive TSV:  {found['kaptive_tsv'].name if found['kaptive_tsv'] else 'NOT FOUND'}")
        print(f"  QC TSV:       {found['qc_tsv'].name if found['qc_tsv'] else 'NOT FOUND'}")
        print(f"  AMRfinder:    {found['amrfinder_html'].name if found['amrfinder_html'] else 'NOT FOUND'}")
        print(f"  MGE HTML:     {found['mge_html'].name if found['mge_html'] else 'NOT FOUND'}")
        print(f"  Mutation:     {found['mutation_html'].name if found['mutation_html'] else 'NOT FOUND'}")
        print(f"  Abricate DBs: {len(found['abricate_html'])}")
        return found

    def integrate(self, files: Dict) -> Dict:
        integrated = {'metadata': self.metadata, 'samples': {}}
        all_samples = set()

        if files['mlst_tsv']:
            for s, d in self.parser.parse_mlst_tsv(files['mlst_tsv']).items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['mlst'] = d
        if files['qc_tsv']:
            for s, d in self.parser.parse_qc_tsv(files['qc_tsv']).items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['qc'] = d
        if files['kaptive_tsv']:
            for s, d in self.parser.parse_kaptive_tsv(files['kaptive_tsv']).items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['kaptive'] = d

        total = len(all_samples)
        amr_gene_freq = {}
        if files['amrfinder_html']:
            genes_by_sample, amr_gene_freq = self.parser.parse_amrfinder_report(files['amrfinder_html'], total)
            for s, genes in genes_by_sample.items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['amr_genes'] = genes

        abricate_gene_freq = {}
        for db, fs in files['abricate_html'].items():
            if fs:
                by_sample, freq = self.parser.parse_abricate_database_report(fs[0], total)
                for s, genes in by_sample.items():
                    all_samples.add(s)
                    integrated['samples'].setdefault(s, {}).setdefault('abricate_genes', {})[db] = genes
                for gene, data in freq.items():
                    abricate_gene_freq.setdefault(db, {})[gene] = data

        for s in all_samples:
            integrated['samples'].setdefault(s, {})

        if files['mutation_html']:
            integrated['mutation_data'] = self.parser.parse_mutation_summary_html(files['mutation_html'])
        else:
            integrated['mutation_data'] = {}

        if files['mge_html']:
            mge_data = self.parser.parse_mge_summary_html(files['mge_html'])
            mge_data['analysis'] = self.analyzer.compute_mge_stats(mge_data)
            integrated['mge_data'] = mge_data
        else:
            integrated['mge_data'] = {'per_sample': {}, 'note': '', 'analysis': {'category_totals': [], 'grand_total': 0}}

        gene_freqs = {}
        if amr_gene_freq:
            gene_freqs['amrfinder'] = amr_gene_freq
        for db, freq in abricate_gene_freq.items():
            gene_freqs[db] = freq
        integrated['gene_frequencies'] = gene_freqs
        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(gene_freqs, len(all_samples))
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(
            integrated['samples'], integrated['gene_centric']
        )

        coverage = {}
        for db, genes in gene_freqs.items():
            with_hits = set()
            for gd in genes.values():
                with_hits.update(gd['genomes'])
            coverage[db] = {
                'samples_with_hits': len(with_hits),
                'total_samples': len(all_samples),
                'coverage_percentage': (len(with_hits) / len(all_samples) * 100) if all_samples else 0,
            }
        integrated['patterns']['database_coverage'] = coverage

        compare_data, compare_sim = self.analyzer.build_compare_payload(
            integrated['samples'], integrated.get('mutation_data', {})
        )
        integrated['compare_data'] = compare_data
        integrated['compare_similarity'] = compare_sim

        return integrated

    def write_json(self, integrated: Dict) -> Path:
        out = self.output_dir / "kleboscope_ultimate_gene_centric_report.json"
        with open(out, 'w', encoding='utf-8') as f:
            json.dump(integrated, f, indent=2, default=str)
        return out

    def write_csvs(self, integrated: Dict):
        samples = integrated['samples']
        rows = []
        for sample, d in samples.items():
            vir_count = sum(
                len(d.get('abricate_genes', {}).get(db, []))
                for db in ('vfdb', 'ecoli_vf')
            )
            mge_hits = integrated.get('mge_data', {}).get('per_sample', {}).get(sample, {}).get('mobileOG Hits', 0)
            rows.append({
                'Sample': sample,
                'ST': d.get('mlst', {}).get('ST', 'ND'),
                'K_Locus': d.get('kaptive', {}).get('K_Locus', 'ND'),
                'K_Type': d.get('kaptive', {}).get('K_Type', 'ND'),
                'O_Locus': d.get('kaptive', {}).get('O_Locus', 'ND'),
                'O_Type': d.get('kaptive', {}).get('O_Type', 'ND'),
                'Virulence_Genes': vir_count,
                'MobileOG_Hits': mge_hits,
            })
        pd.DataFrame(rows).to_csv(self.output_dir / "sample_overview.csv", index=False)

        gc = integrated['gene_centric']
        for cat_key, fname in [
            ('amr_databases', 'amr_genes.csv'),
            ('virulence_databases', 'virulence_genes.csv'),
            ('bacmet_databases', 'bacmet_genes.csv'),
            ('plasmid_databases', 'plasmid_replicons.csv'),
        ]:
            out_rows = []
            for db, genes in gc.get(cat_key, {}).items():
                for g in genes:
                    out_rows.append({
                        'Gene': g['gene'],
                        'Database': g['database'],
                        'Count': g['count'],
                        'Frequency': g['frequency_display'],
                        'Genomes': ';'.join(g['genomes']),
                    })
            if out_rows:
                pd.DataFrame(out_rows).to_csv(self.output_dir / fname, index=False)

        muts = integrated.get('mutation_data', {}).get('mutations', [])
        if muts:
            pd.DataFrame([{
                'Gene': m['gene'],
                'Mutation': m['mutation'],
                'Class': m['class'],
                'Subclass': m['subclass'],
                'Count': m['count'],
                'Genomes': ';'.join(m['genomes']),
            } for m in muts]).to_csv(self.output_dir / "mutations.csv", index=False)

        qc_rows = []
        for sample, d in samples.items():
            if d.get('qc'):
                row = {'Sample': sample}
                row.update(d['qc'])
                qc_rows.append(row)
        if qc_rows:
            pd.DataFrame(qc_rows).to_csv(self.output_dir / "fasta_qc.csv", index=False)

        mge_ps = integrated.get('mge_data', {}).get('per_sample', {})
        if mge_ps:
            pd.DataFrame([{'Sample': s, **v} for s, v in mge_ps.items()]).to_csv(
                self.output_dir / "mge_profile.csv", index=False
            )

    def run(self) -> bool:
        print("=" * 80)
        print("Kleboscope Ultimate Reporter v3.0.0")
        print("=" * 80)
        print(f"Input:  {self.input_dir}")
        print(f"Output: {self.output_dir}")
        files = self.find_files()
        if not any([files['mlst_tsv'], files['kaptive_tsv'], files['qc_tsv'],
                    files['amrfinder_html'], files['abricate_html'], files['mge_html']]):
            print("No recognized input files found.")
            return False
        integrated = self.integrate(files)
        if not integrated['samples']:
            print("No samples integrated.")
            return False
        self.write_json(integrated)
        self.write_csvs(integrated)
        self.generator.generate_main_report(integrated, self.output_dir)
        print()
        print("=" * 80)
        print("REPORT COMPLETE")
        print("=" * 80)
        print(f"Samples:  {len(integrated['samples'])}")
        print(f"Output:   {self.output_dir}")
        print("=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(
        description='Kleboscope Ultimate Reporter - Gene-centric K. pneumoniae analysis'
    )
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory containing Kleboscope reports')
    args = parser.parse_args()
    reporter = KleboscopeUltimateReporter(Path(args.input_dir))
    if not reporter.run():
        sys.exit(1)


if __name__ == "__main__":
    main()