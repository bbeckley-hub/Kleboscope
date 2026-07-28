#!/usr/bin/env python3
"""
Kleboscope Ultimate Reporter – Hybrid Gene-Centric & Sample-Centric
====================================================================
Version 1.0.0

Gene-centric for MLST/Kaptive/QC/Combinations/Patterns.
Sample-centric interactive boxes for AMR, Virulence, BACMET, Plasmids, and Mutations.

Each sample-centric tab displays per-isolate boxes with:
- Sample name, total count, and typing badges (ST, K-locus, O-locus, hypervirulence).
- Horizontally scrollable tables with full details from TSV files.
- Filtering by sample name and by database (where applicable).

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Date: 2026-07-28
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


# -----------------------------------------------------------------------------
# PARSER CLASS – extends with TSV loaders for sample-centric data
# -----------------------------------------------------------------------------
class KleboHTMLParser:
    """Parse all Kleboscope reports, including TSV details for sample-centric boxes."""

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
        self.abricate_tsv_files = {
            'card': 'klebo_card_abricate_summary.tsv',
            'resfinder': 'klebo_resfinder_abricate_summary.tsv',
            'argannot': 'klebo_argannot_abricate_summary.tsv',
            'vfdb': 'klebo_vfdb_abricate_summary.tsv',
            'plasmidfinder': 'klebo_plasmidfinder_abricate_summary.tsv',
            'megares': 'klebo_megares_abricate_summary.tsv',
            'ncbi': 'klebo_ncbi_abricate_summary.tsv',
            'ecoh': 'klebo_ecoh_abricate_summary.tsv',
            'ecoli_vf': 'klebo_ecoli_vf_abricate_summary.tsv',
            'bacmet2': 'klebo_bacmet2_abricate_summary.tsv',
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id).strip()
        for ext in ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample

    # --------------------------------------------------------------------------
    # TSV LOADERS (for sample-centric details)
    # --------------------------------------------------------------------------
    def load_amrfinder_details_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, List[Dict]], Dict[str, int]]:
        tsv_path = input_dir / 'klebo_amrfinder_summary.tsv'
        if not tsv_path.exists():
            return {}, {}
        df = pd.read_csv(tsv_path, sep='\t')
        details = defaultdict(list)
        gene_counts = Counter()
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(str(row['Genome']))
            gene_dict = row.to_dict()
            gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
            # Rename columns to standard keys (optional)
            if 'Gene_Symbol' in gene_dict:
                gene_dict['gene'] = gene_dict.pop('Gene_Symbol')
            details[sample].append(gene_dict)
            gene_counts[gene_dict.get('gene', 'unknown')] += 1
        return dict(details), dict(gene_counts)

    def load_abricate_details_from_tsv(self, input_dir: Path) -> Tuple[Dict[str, Dict[str, List[Dict]]], Dict[str, Dict[str, int]]]:
        details = defaultdict(lambda: defaultdict(list))
        gene_counts = defaultdict(lambda: defaultdict(int))
        for db, fname in self.abricate_tsv_files.items():
            path = input_dir / fname
            if not path.exists():
                continue
            df = pd.read_csv(path, sep='\t')
            sample_col = 'genome' if 'genome' in df.columns else 'file'
            for _, row in df.iterrows():
                sample = self.normalize_sample_id(str(row[sample_col]))
                gene_dict = row.to_dict()
                gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
                details[sample][db].append(gene_dict)
                gene = gene_dict.get('gene', 'unknown')
                gene_counts[db][gene] += 1
        return dict(details), dict(gene_counts)

    def load_mutation_details_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        tsv_path = input_dir / 'mutation_summary.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        details = defaultdict(list)
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(str(row['genome']))
            gene_dict = row.to_dict()
            gene_dict = {k: (None if pd.isna(v) else v) for k, v in gene_dict.items()}
            # Rename for consistency
            if 'gene_symbol' in gene_dict:
                gene_dict['gene'] = gene_dict.pop('gene_symbol')
            if 'element_name' in gene_dict:
                gene_dict['mutation'] = gene_dict.pop('element_name')
            details[sample].append(gene_dict)
        return dict(details)

    def load_mlst_from_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        tsv_path = input_dir / 'mlst_summary.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        results = {}
        sample_col = None
        for col in df.columns:
            if 'sample' in col.lower():
                sample_col = col
                break
        if not sample_col:
            sample_col = df.columns[0]
        st_col = None
        for col in df.columns:
            if col.lower() == 'st' or ('st' in col.lower() and 'sample' not in col.lower()):
                st_col = col
                break
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row[sample_col])
            st = row[st_col] if st_col else 'ND'
            if pd.isna(st) or str(st).lower() in ['', 'nan', 'none', 'nd', 'unknown']:
                st = 'ND'
            else:
                st_val = str(st).strip()
                if st_val.startswith('ST'):
                    st = st_val[2:]
                else:
                    st = st_val
            results[sample] = {'ST': st}
        return results

    def load_kaptive_from_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        tsv_path = input_dir / 'klebo_kaptive_summary.tsv'
        if not tsv_path.exists():
            return {}
        df = pd.read_csv(tsv_path, sep='\t')
        results = {}
        genome_col = None
        for col in df.columns:
            if 'genome' in col.lower():
                genome_col = col
                break
        if not genome_col and len(df.columns) > 0:
            genome_col = df.columns[0]
        for _, row in df.iterrows():
            sample = self.normalize_sample_id(row[genome_col])
            k_locus = row.get('K Locus', 'ND')
            o_locus = row.get('O Locus', 'ND')
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
                'K_Identity': row.get('K Identity', 'ND'),
                'K_Coverage': row.get('K Coverage', 'ND'),
            }
        return results

    # --------------------------------------------------------------------------
    # HTML FALLBACK PARSERS 
    # --------------------------------------------------------------------------
    def parse_html_table(self, html_content: str, table_index: int = 0) -> pd.DataFrame:
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
            print(f"  ⚠️ Table parsing error: {e}")
            return pd.DataFrame()

    def parse_mlst_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing MLST HTML: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}
            sample_col = None
            for col in df.columns:
                if 'sample' in col.lower():
                    sample_col = col
                    break
            if not sample_col and len(df.columns) > 0:
                sample_col = df.columns[0]
            st_col = None
            for col in df.columns:
                if col.lower() == 'st' or ('st' in col.lower() and 'sample' not in col.lower()):
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

    def parse_qc_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            df = self.parse_html_table(html, 0)
            if df.empty:
                return {}
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
                qc_data = {}
                for col in df.columns:
                    if col == sample_col:
                        continue
                    val = row[col]
                    if pd.isna(val) or val == '' or val == 'ND':
                        qc_data[col] = 'ND'
                    else:
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

    def parse_kaptive_report(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing Kaptive HTML: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')

            # Find the table containing K and O locus columns
            tables = soup.find_all('table')
            target_table = None
            for table in tables:
                headers = [th.get_text().strip() for th in table.find_all('th')]
                header_lower = [h.lower() for h in headers]
                if any('k locus' in h for h in header_lower) and any('o locus' in h for h in header_lower):
                    target_table = table
                    break
            if not target_table:
                print("    ⚠️ Could not find Kaptive table with K/O columns")
                return {}

            rows = target_table.find_all('tr')
            if len(rows) < 2:
                return {}

            # Parse headers
            header_cells = rows[0].find_all('th')
            if not header_cells:
                header_cells = rows[0].find_all('td')
            headers = [cell.get_text().strip() for cell in header_cells]

            # Map column indices
            col_map = {}
            for idx, h in enumerate(headers):
                h_lower = h.lower()
                if 'genome' in h_lower or 'sample' in h_lower:
                    col_map['sample'] = idx
                elif 'k locus' in h_lower:
                    col_map['k_locus'] = idx
                elif 'o locus' in h_lower:
                    col_map['o_locus'] = idx
                elif 'k identity' in h_lower:
                    col_map['k_identity'] = idx
                elif 'k coverage' in h_lower:
                    col_map['k_coverage'] = idx

            if 'sample' not in col_map or 'k_locus' not in col_map or 'o_locus' not in col_map:
                print(f"    ⚠️ Missing required columns. Found: {headers}")
                return {}

            results = {}
            for row in rows[1:]:
                cells = row.find_all('td')
                if len(cells) <= max(col_map.values()):
                    continue
                sample_raw = cells[col_map['sample']].get_text().strip()
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                k_locus = cells[col_map['k_locus']].get_text().strip()
                o_locus = cells[col_map['o_locus']].get_text().strip()

                # Clean "unknown" KL/OC patterns
                if 'unknown' in k_locus.lower():
                    k_match = re.search(r'KL(\d+)', k_locus, re.I)
                    if k_match:
                        k_locus = f"KL{k_match.group(1)}"
                if 'unknown' in o_locus.lower():
                    o_match = re.search(r'OCL?(\d+)', o_locus, re.I)
                    if o_match:
                        o_locus = f"OC{o_match.group(1)}"

                results[sample] = {
                    'K_Locus': k_locus if k_locus else 'ND',
                    'O_Locus': o_locus if o_locus else 'ND',
                    'K_Identity': cells[col_map.get('k_identity', -1)].get_text().strip() if 'k_identity' in col_map else 'ND',
                    'K_Coverage': cells[col_map.get('k_coverage', -1)].get_text().strip() if 'k_coverage' in col_map else 'ND',
                }
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error parsing Kaptive: {e}")
            return {}

    def parse_amrfinder_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        # Keep existing HTML fallback; but we now prefer TSV
        print(f"  🧬 Parsing AMRfinder HTML (fallback): {file_path.name}")
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
            print(f"    ❌ Error parsing AMRfinder HTML: {e}")
            return {}, {}

    def parse_abricate_database_report(self, file_path: Path, total_samples: int = 0) -> Tuple[Dict[str, List], Dict[str, Dict]]:
        # HTML fallback
        print(f"  🧬 Parsing ABRicate HTML (fallback): {file_path.name}")
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
            print(f"    ❌ Error parsing ABRicate HTML: {e}")
            return {}, {}

    def parse_mutation_summary_html(self, file_path: Path) -> Dict[str, Any]:
        # Keep existing HTML mutation parser (used if TSV missing)
        print(f"  🧬 Parsing mutation summary HTML (fallback): {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            soup = BeautifulSoup(html_content, 'html.parser')
            mutation_table = None
            for table in soup.find_all('table'):
                if table.find(string=re.compile(r'Gene', re.I)) and table.find(string=re.compile(r'Mutation', re.I)):
                    mutation_table = table
                    break
            if not mutation_table:
                print("    ⚠️ Could not find mutation table")
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
                h_lower = h.lower()
                if 'gene' in h_lower:
                    col_idx['gene'] = idx
                elif 'mutation' in h_lower:
                    col_idx['mutation'] = idx
                elif 'count' in h_lower:
                    col_idx['count'] = idx
                elif 'genome' in h_lower:
                    col_idx['genomes'] = idx
                elif 'class' in h_lower:
                    col_idx['class'] = idx
                elif 'subclass' in h_lower:
                    col_idx['subclass'] = idx
            required = ['gene', 'mutation', 'count', 'genomes']
            for req in required:
                if req not in col_idx:
                    print(f"    ⚠️ Missing required column: {req}. Found: {headers}")
                    return {}
            tbody = mutation_table.find('tbody')
            if tbody:
                rows = tbody.find_all('tr')
            else:
                rows = mutation_table.find_all('tr')[1:]
            mutations_list = []
            genome_counts = defaultdict(int)
            for row in rows:
                cells = row.find_all('td')
                if len(cells) <= max(col_idx.values()):
                    continue
                gene = cells[col_idx['gene']].get_text().strip()
                mutation = cells[col_idx['mutation']].get_text().strip()
                count_str = cells[col_idx['count']].get_text().strip()
                count_match = re.search(r'(\d+)', count_str)
                count = int(count_match.group(1)) if count_match else 0
                genomes_str = cells[col_idx['genomes']].get_text().strip()
                genomes = [g.strip() for g in genomes_str.split(',') if g.strip()]
                if not genomes:
                    continue
                for g in genomes:
                    genome_counts[g] += 1
                class_name = cells[col_idx['class']].get_text().strip() if 'class' in col_idx else ''
                subclass = cells[col_idx['subclass']].get_text().strip() if 'subclass' in col_idx else ''
                mutations_list.append({
                    'gene': gene,
                    'mutation': mutation,
                    'class': class_name,
                    'subclass': subclass,
                    'count': count,
                    'genomes': genomes
                })
            mutations_list.sort(key=lambda x: x['count'], reverse=True)
            print(f"    ✓ Parsed {len(mutations_list)} unique mutations across {len(genome_counts)} genomes")
            return {
                'mutations': mutations_list,
                'genome_mutation_counts': dict(genome_counts)
            }
        except Exception as e:
            print(f"    ❌ Error parsing mutation summary HTML: {e}")
            return {}


# -----------------------------------------------------------------------------
# DATA ANALYZER 
# -----------------------------------------------------------------------------
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
        if any(crit.lower() in g for crit in self.critical_resistance_genes):
            return 'Critical Resistance'
        elif any(vir.lower() in g for vir in self.high_risk_virulence_genes):
            return 'High-Risk Virulence'
        elif any(bla.lower() in g for bla in self.beta_lactamase_genes):
            return 'Beta-Lactamase'
        else:
            return 'Other'

    def create_gene_centric_tables(self, gene_freqs: Dict[str, Dict], total_samples: int) -> Dict:
        gene_centric = {
            'all_genes': [],
            'by_category': defaultdict(list),
            'database_stats': {},
            'amr_databases': {},
            'virulence_databases': {},
            'plasmid_databases': {},
            'bacmet_databases': {}
        }
        all_genes = {}
        for db, genes in gene_freqs.items():
            if db == 'bacmet2':
                target_dict = gene_centric['bacmet_databases']
            elif db in ['vfdb', 'ecoli_vf']:
                target_dict = gene_centric['virulence_databases']
            elif db == 'plasmidfinder':
                target_dict = gene_centric['plasmid_databases']
            else:
                target_dict = gene_centric['amr_databases']

            gene_list = []
            for gene, data in genes.items():
                if gene not in all_genes:
                    all_genes[gene] = {
                        'count': 0, 'percentage': 0, 'frequency_display': '',
                        'genomes': set(), 'databases': set(), 'risk_level': 'Standard'
                    }
                all_genes[gene]['count'] += data['count']
                all_genes[gene]['percentage'] += data.get('percentage', 0)
                all_genes[gene]['frequency_display'] = data.get('frequency_display', f"{data['count']} ({data['count']/total_samples*100:.1f}%)")
                all_genes[gene]['genomes'].update(data['genomes'])
                all_genes[gene]['databases'].add(data.get('database', db))
                if data.get('risk_level') and data['risk_level'] != 'Standard':
                    all_genes[gene]['risk_level'] = data['risk_level']

                gene_list.append({
                    'gene': gene,
                    'database': db.upper(),
                    'frequency_display': data.get('frequency_display', f"{data['count']} ({data.get('percentage',0):.1f}%)"),
                    'count': data['count'],
                    'percentage': data.get('percentage', 0),
                    'genomes': list(data['genomes']),
                    'risk_level': data.get('risk_level', 'Standard')
                })

            if gene_list:
                gene_list.sort(key=lambda x: x['count'], reverse=True)
                target_dict[db] = gene_list

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
            'gene_cooccurrence': defaultdict(Counter)
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
            crit_res = [g for g in genes if self.categorize_gene(g) == 'Critical Resistance']
            hv = [g for g in genes if self.categorize_gene(g) == 'High-Risk Virulence']
            if crit_res and hv:
                patterns['high_risk_combinations'].append({
                    'sample': sample, 'st': st, 'k_locus': k,
                    'critical_resistance': crit_res,
                    'high_risk_virulence': hv
                })
            icekp_found = [g for g in genes if any(marker in g.lower() for marker in self.icekp_markers)]
            if icekp_found:
                patterns['icekp_marker_presence'][sample] = icekp_found
            vp_found = [g for g in genes if any(marker in g.lower() for marker in self.virulence_plasmid_markers)]
            if vp_found:
                patterns['virulence_plasmid_marker_presence'][sample] = vp_found
            genes_list = list(genes)
            for i, g1 in enumerate(genes_list):
                for g2 in genes_list[i+1:]:
                    patterns['gene_cooccurrence'][g1][g2] += 1
        return patterns


# -----------------------------------------------------------------------------
# HTML GENERATOR – now with sample‑centric boxes
# -----------------------------------------------------------------------------
class KleboHTMLGenerator:
    def __init__(self, analyzer: KleboDataAnalyzer):
        self.analyzer = analyzer
        self.tab_colors = {
            'summary': '#4CAF50',
            'sample_overview': '#2196F3',
            'mlst': '#FF9800',
            'qc': '#607D8B',
            'kaptive': '#9C27B0',
            'combinations': '#009688',
            'amr': '#F44336',
            'virulence': '#E91E63',
            'bacmet': '#FF5722',
            'plasmids': '#673AB7',
            'mutations': '#00BCD4',
            'patterns': '#3F51B5',
            'highrisk': '#DC143C',
            'databases': '#795548',
            'credit': '#8BC34A',
            'aiguide': '#00BCD4',
            'citation': '#FFC107',
            'funding': '#FFB300',
            'export': '#9E9E9E'
        }

    def generate_main_report(self, integrated_data: Dict, output_dir: Path) -> str:
        print("\n🎨 Generating Kleboscope Ultimate HTML report v1.0.0 (Hybrid) ...")
        html = self._create_ultimate_html(integrated_data)
        output_file = output_dir / "kleboscope_ultimate_sample_centric_report.html"
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

        # Build sample typing JSON for badges (used in sample‑centric boxes)
        sample_typing = {}
        hv_markers_set = set(self.analyzer.high_risk_virulence_genes)
        for sample, d in samples.items():
            st = d.get('mlst', {}).get('ST', 'ND')
            k = d.get('kaptive', {}).get('K_Locus', 'ND')
            o = d.get('kaptive', {}).get('O_Locus', 'ND')
            # Hypervirulence if any high‑risk virulence gene is present (from abricate)
            vir_genes = []
            for db, genes in d.get('abricate_genes', {}).items():
                if db in ['vfdb', 'ecoli_vf']:
                    vir_genes.extend(genes)
            hv_status = 'Hypervirulent' if any(g in hv_markers_set or any(marker in g for marker in ['ybt','clb','iro','iuc','rmp']) for g in vir_genes) else 'No'
            sample_typing[sample] = {
                "ST": st,
                "K": k,
                "O": o,
                "HV": hv_status
            }
        sample_typing_json = json.dumps(sample_typing)

        # CSS and JS
        css = self._get_css()
        js = self._get_js(sample_typing_json)

        # Generate all sections
        summary_html = self._generate_summary_section(data)
        samples_html = self._generate_sample_overview_section(data)
        mlst_html = self._generate_mlst_section(data)
        qc_html = self._generate_qc_section(data)
        kaptive_html = self._generate_kaptive_section(data)
        combinations_html = self._generate_combinations_section(data)

        # SAMPLE-CENTRIC TABS
        amr_html = self._generate_sample_centric_boxes(data, 'amr', 'AMR',
                    ['amrfinder', 'card', 'resfinder', 'argannot', 'megares', 'ncbi'])
        virulence_html = self._generate_sample_centric_boxes(data, 'vir', 'Virulence',
                    ['vfdb', 'ecoli_vf'])
        bacmet_html = self._generate_sample_centric_boxes(data, 'bac', 'BACMET',
                    ['bacmet2'])
        plasmids_html = self._generate_sample_centric_boxes(data, 'plasmids', 'Plasmids',
                    ['plasmidfinder'])
        mutations_html = self._generate_mutation_boxes(data)

        # Gene‑centric pattern and high‑risk tabs (unchanged)
        patterns_html = self._generate_patterns_section(data)
        highrisk_html = self._generate_highrisk_section(data)
        databases_html = self._generate_databases_section(data)
        credit_html = self._generate_credit_section(data)
        aiguide_html = self._generate_aiguide_section(data)
        citation_html = self._generate_citation_section(data)
        funding_html = self._generate_funding_section(data)
        export_html = self._generate_export_section(data)

        # Build tab buttons
        tab_buttons = []
        tab_order = [
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
            ('databases', 'Databases', 'fa-database'),
            ('credit', 'Credit', 'fa-thumbs-up'),
            ('aiguide', 'AI Guide', 'fa-robot'),
            ('citation', 'Citation', 'fa-book'),
            ('funding', 'Funding', 'fa-coffee'),
            ('export', 'Export', 'fa-download')
        ]
        for name, label, icon in tab_order:
            color = self.tab_colors.get(name, '#6c757d')
            btn = f'<button class="tab-button {name}" onclick="switchTab(\'{name}\')" style="background-color:{color};color:white;"><i class="fas {icon}"></i> {label}</button>'
            tab_buttons.append(btn)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Kleboscope Ultimate Report – Hybrid</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    {css}
    {js}
</head>
<body>
<div class="container">
    <div class="main-header">
        <h1><i class="fas fa-bacterium"></i> Kleboscope Ultimate Report</h1>
        <p>Hybrid: Gene‑Centric for Typing / Sample‑Centric for AMR, Virulence, Bacmet, Plasmids, Mutations</p>
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar-alt"></i> {metadata.get('analysis_date', 'Unknown')}</div>
            <div class="metadata-item"><i class="fas fa-database"></i> {total_samples} Samples</div>
            <div class="metadata-item"><i class="fas fa-university"></i> University of Ghana Medical School</div>
        </div>
    </div>

    <div class="dashboard-grid">
        <div class="dashboard-card card-summary" onclick="switchTab('summary')"><div class="card-number">{total_samples}</div><div class="card-label">Samples</div></div>
        <div class="dashboard-card card-mlst" onclick="switchTab('mlst')"><div class="card-number">{len(patterns.get('st_distribution', {}))}</div><div class="card-label">STs</div></div>
        <div class="dashboard-card card-kaptive" onclick="switchTab('kaptive')"><div class="card-number">{len(patterns.get('k_locus_distribution', {}))}</div><div class="card-label">Capsule Types</div></div>
        <div class="dashboard-card card-amr" onclick="switchTab('amr')"><div class="card-number">{len(gene_centric.get('all_genes', []))}</div><div class="card-label">Total Genes</div></div>
        <div class="dashboard-card card-highrisk" onclick="switchTab('highrisk')"><div class="card-number">{len(patterns.get('high_risk_combinations', []))}</div><div class="card-label">High‑Risk Combos</div></div>
    </div>

    <div class="tab-navigation">
        {''.join(tab_buttons)}
    </div>

    <div id="summary-tab" class="tab-content active">{summary_html}</div>
    <div id="sample_overview-tab" class="tab-content">{samples_html}</div>
    <div id="mlst-tab" class="tab-content">{mlst_html}</div>
    <div id="qc-tab" class="tab-content">{qc_html}</div>
    <div id="kaptive-tab" class="tab-content">{kaptive_html}</div>
    <div id="combinations-tab" class="tab-content">{combinations_html}</div>
    <div id="amr-tab" class="tab-content">{amr_html}</div>
    <div id="virulence-tab" class="tab-content">{virulence_html}</div>
    <div id="bacmet-tab" class="tab-content">{bacmet_html}</div>
    <div id="plasmids-tab" class="tab-content">{plasmids_html}</div>
    <div id="mutations-tab" class="tab-content">{mutations_html}</div>
    <div id="patterns-tab" class="tab-content">{patterns_html}</div>
    <div id="highrisk-tab" class="tab-content">{highrisk_html}</div>
    <div id="databases-tab" class="tab-content">{databases_html}</div>
    <div id="credit-tab" class="tab-content">{credit_html}</div>
    <div id="aiguide-tab" class="tab-content">{aiguide_html}</div>
    <div id="citation-tab" class="tab-content">{citation_html}</div>
    <div id="funding-tab" class="tab-content">{funding_html}</div>
    <div id="export-tab" class="tab-content">{export_html}</div>

    <div class="footer">
        <h3>Kleboscope Ultimate Reporter v1.0.0</h3>
        <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
        <p>Generated on {metadata.get('analysis_date', 'Unknown')}</p>
        <p><strong>Critical Genes Tracked:</strong> Carbapenemases (KPC, NDM, OXA-48) • Colistin (mcr) • Tigecycline (tetX) • ICEKp Markers (ybt, clb, iro, rmp) • Virulence Plasmid Markers (iro, iuc, rmp, rmpA2) • Biocides & Heavy Metals (qac, sil, mer, ars, pco) • Adhesins (fim, mrk, ecp) • Secretion Systems (tss) • Siderophores • Toxins</p>
        <p>If you find this useful, please <a href="https://github.com/bbeckley-hub" target="_blank">⭐ star us on GitHub</a> and share with your network.</p>
    </div>
</div>
</body>
</html>
"""
        return html

    # --------------------------------------------------------------------------
    # CSS and JavaScript
    # --------------------------------------------------------------------------
    def _get_css(self) -> str:
        return """
        <style>
        :root {
            --summary-color: #4CAF50;
            --sample_overview-color: #2196F3;
            --qc-color: #607D8B;
            --mlst-color: #FF9800;
            --kaptive-color: #9C27B0;
            --combinations-color: #009688;
            --amr-color: #F44336;
            --virulence-color: #E91E63;
            --bacmet-color: #FF5722;
            --plasmids-color: #673AB7;
            --mutations-color: #00BCD4;
            --patterns-color: #3F51B5;
            --highrisk-color: #DC143C;
            --databases-color: #795548;
            --credit-color: #8BC34A;
            --aiguide-color: #00BCD4;
            --citation-color: #FFC107;
            --funding-color: #FFB300;
            --export-color: #9E9E9E;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; background: #f5f5f5; min-width: 1200px; }
        .container { max-width: none; margin: 0 auto; padding: 20px; width: 100%; overflow-x: auto; }
        .main-header { background: linear-gradient(135deg, #1e5a3a 0%, #2c7a4d 100%); color: white; padding: 30px; border-radius: 15px; box-shadow: 0 10px 30px rgba(0,0,0,0.2); margin-bottom: 30px; text-align: center; }
        .main-header h1 { font-size: 2.8em; margin-bottom: 10px; color: white; }
        .metadata-bar { background: rgba(255,255,255,0.1); padding: 15px; border-radius: 10px; margin: 20px 0; display: flex; justify-content: space-around; flex-wrap: wrap; gap: 15px; backdrop-filter: blur(10px); }
        .metadata-item { display: flex; align-items: center; gap: 8px; font-size: 0.95em; }
        .dashboard-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin-bottom: 30px; }
        .dashboard-card { background: white; padding: 25px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); text-align: center; transition: all 0.3s ease; cursor: pointer; border-left: 5px solid; position: relative; overflow: hidden; }
        .dashboard-card::before { content: ''; position: absolute; top: 0; left: 0; right: 0; height: 3px; background: linear-gradient(90deg, transparent, rgba(255,255,255,0.8), transparent); }
        .dashboard-card:hover { transform: translateY(-10px); box-shadow: 0 15px 30px rgba(0,0,0,0.2); }
        .card-summary { border-left-color: var(--summary-color); }
        .card-sample_overview { border-left-color: var(--sample_overview-color); }
        .card-mlst { border-left-color: var(--mlst-color); }
        .card-qc { border-left-color: var(--qc-color); }
        .card-kaptive { border-left-color: var(--kaptive-color); }
        .card-combinations { border-left-color: var(--combinations-color); }
        .card-amr { border-left-color: var(--amr-color); }
        .card-virulence { border-left-color: var(--virulence-color); }
        .card-bacmet { border-left-color: var(--bacmet-color); }
        .card-plasmids { border-left-color: var(--plasmids-color); }
        .card-mutations { border-left-color: var(--mutations-color); }
        .card-patterns { border-left-color: var(--patterns-color); }
        .card-highrisk { border-left-color: var(--highrisk-color); }
        .card-databases { border-left-color: var(--databases-color); }
        .card-credit { border-left-color: var(--credit-color); }
        .card-aiguide { border-left-color: var(--aiguide-color); }
        .card-citation { border-left-color: var(--citation-color); }
        .card-funding { border-left-color: var(--funding-color); }
        .card-export { border-left-color: var(--export-color); }
        .card-number { font-size: 3em; font-weight: bold; margin: 15px 0; background: linear-gradient(90deg, #1e5a3a, #2c7a4d); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
        .tab-navigation { display: flex; gap: 5px; margin-bottom: 20px; flex-wrap: wrap; background: white; padding: 15px; border-radius: 12px; box-shadow: 0 5px 20px rgba(0,0,0,0.1); position: sticky; top: 10px; z-index: 100; }
        .tab-button { padding: 12px 20px; background: #f5f5f5; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; color: #666; transition: all 0.3s ease; display: flex; align-items: center; gap: 8px; position: relative; overflow: hidden; font-size: 0.9em; }
        .tab-button::after { content: ''; position: absolute; bottom: 0; left: 50%; right: 50%; height: 3px; background: currentColor; transition: all 0.3s ease; }
        .tab-button:hover::after { left: 10%; right: 10%; }
        .tab-button.active { color: white; }
        .tab-button.active::after { left: 10%; right: 10%; }
        .tab-button.summary.active { background: var(--summary-color); }
        .tab-button.sample_overview.active { background: var(--sample_overview-color); }
        .tab-button.mlst.active { background: var(--mlst-color); }
        .tab-button.qc.active { background: var(--qc-color); }
        .tab-button.kaptive.active { background: var(--kaptive-color); }
        .tab-button.combinations.active { background: var(--combinations-color); }
        .tab-button.amr.active { background: var(--amr-color); }
        .tab-button.virulence.active { background: var(--virulence-color); }
        .tab-button.bacmet.active { background: var(--bacmet-color); }
        .tab-button.plasmids.active { background: var(--plasmids-color); }
        .tab-button.mutations.active { background: var(--mutations-color); }
        .tab-button.patterns.active { background: var(--patterns-color); }
        .tab-button.highrisk.active { background: var(--highrisk-color); }
        .tab-button.databases.active { background: var(--databases-color); }
        .tab-button.credit.active { background: var(--credit-color); }
        .tab-button.aiguide.active { background: var(--aiguide-color); }
        .tab-button.citation.active { background: var(--citation-color); }
        .tab-button.funding.active { background: var(--funding-color); }
        .tab-button.export.active { background: var(--export-color); }
        .tab-content { display: none; background: white; padding: 30px; border-radius: 15px; box-shadow: 0 10px 30px rgba(0,0,0,0.1); margin-bottom: 30px; animation: fadeIn 0.5s ease; width: 100%; overflow-x: auto; }
        .tab-content.active { display: block; }
        @keyframes fadeIn { from { opacity: 0; transform: translateY(20px); } to { opacity: 1; transform: translateY(0); } }
        .section-header { color: #2c3e50; margin-bottom: 25px; padding-bottom: 15px; border-bottom: 3px solid; font-size: 1.8em; display: flex; align-items: center; justify-content: space-between; }
        .summary-header { border-color: var(--summary-color); }
        .sample_overview-header { border-color: var(--sample_overview-color); }
        .mlst-header { border-color: var(--mlst-color); }
        .qc-header { border-color: var(--qc-color); }
        .kaptive-header { border-color: var(--kaptive-color); }
        .combinations-header { border-color: var(--combinations-color); }
        .amr-header { border-color: var(--amr-color); }
        .virulence-header { border-color: var(--virulence-color); }
        .bacmet-header { border-color: var(--bacmet-color); }
        .plasmids-header { border-color: var(--plasmids-color); }
        .mutations-header { border-color: var(--mutations-color); }
        .patterns-header { border-color: var(--patterns-color); }
        .highrisk-header { border-color: var(--highrisk-color); }
        .databases-header { border-color: var(--databases-color); }
        .credit-header { border-color: var(--credit-color); }
        .aiguide-header { border-color: var(--aiguide-color); }
        .citation-header { border-color: var(--citation-color); }
        .funding-header { border-color: var(--funding-color); }
        .export-header { border-color: var(--export-color); }
        .data-table { width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.95em; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; table-layout: auto; }
        .data-table th { background: #2c3e50; color: white; padding: 15px; text-align: left; font-weight: 600; position: sticky; top: 0; white-space: nowrap; cursor: pointer; }
        .data-table th:hover { background: #1a252f; }
        .data-table td { padding: 12px; border-bottom: 1px solid #e0e0e0; vertical-align: top; word-wrap: break-word; white-space: nowrap; }
        .data-table tr:hover { background: #f8f9fa; }
        .scrollable-table { max-height: none; overflow-y: auto; border: 1px solid #e0e0e0; border-radius: 8px; margin: 20px 0; width: 100%; }
        .master-scrollable-container { width: 100%; overflow-x: auto; border: 1px solid #e0e0e0; border-radius: 8px; margin: 20px 0; }
        .genome-list { display: flex; flex-wrap: wrap; gap: 5px; max-height: 200px; overflow-y: auto; padding: 5px; background: #f8f9fa; border-radius: 5px; }
        .genome-group { margin-bottom: 10px; width: 100%; }
        .genome-group-header { font-weight: bold; background: #e0e0e0; padding: 4px 8px; border-radius: 4px; margin: 5px 0; font-size: 0.85em; display: inline-block; }
        .genome-group-tags { display: flex; flex-wrap: wrap; gap: 5px; margin-left: 10px; }
        .genome-tag { display: inline-block; background: #e6ffe6; color: #006400; padding: 3px 10px; border-radius: 12px; font-size: 0.85em; border: 1px solid #b3ffb3; white-space: nowrap; margin: 2px; }
        .genome-tag.highlight { background-color: #ffff99 !important; color: #000 !important; border: 1px solid #ffc107; }
        .search-box { width: 100%; padding: 12px; margin-bottom: 20px; border: 2px solid #e0e0e0; border-radius: 8px; font-size: 1em; transition: all 0.3s ease; }
        .search-box:focus { outline: none; border-color: #006400; box-shadow: 0 0 0 3px rgba(0, 100, 0, 0.1); }
        .badge { display: inline-block; padding: 5px 15px; border-radius: 20px; font-size: 0.85em; font-weight: 600; margin: 2px; }
        .badge-critical { background: #DC143C; color: white; }
        .badge-high { background: #FF4500; color: white; }
        .badge-medium { background: #FF8C00; color: black; }
        .badge-low { background: #32CD32; color: white; }
        .alert-box { padding: 20px; border-radius: 10px; margin: 20px 0; display: flex; align-items: center; gap: 20px; border-left: 5px solid; }
        .alert-success { background: #d4edda; color: #155724; border-left-color: #28a745; }
        .alert-warning { background: #fff3cd; color: #856404; border-left-color: #ffc107; }
        .alert-danger { background: #f8d7da; color: #721c24; border-left-color: #dc3545; }
        .alert-info { background: #d1ecf1; color: #0c5460; border-left-color: #17a2b8; }
        .action-buttons { display: flex; gap: 10px; margin: 20px 0; flex-wrap: wrap; }
        .action-btn { padding: 10px 20px; border: none; border-radius: 8px; cursor: pointer; font-weight: 600; display: flex; align-items: center; gap: 8px; transition: all 0.3s ease; }
        .action-btn:hover { transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0,0,0,0.2); }
        .btn-primary { background: #006400; color: white; }
        .btn-success { background: #28a745; color: white; }
        .btn-danger { background: #dc3545; color: white; }
        .btn-warning { background: #ffc107; color: black; }
        .btn-info { background: #17a2b8; color: white; }
        .btn-secondary { background: #6c757d; color: white; }
        .btn-light { background: #f8f9fa; color: #212529; border: 1px solid #dee2e6; }
        .database-section { margin: 30px 0; padding: 25px; border-radius: 12px; background: #f8f9fa; box-shadow: 0 3px 15px rgba(0,0,0,0.08); }
        .database-header { font-size: 1.4em; color: #2c3e50; margin-bottom: 20px; padding-bottom: 10px; border-bottom: 2px solid #006400; display: flex; align-items: center; justify-content: space-between; }
        .print-section-btn { background: #006400; color: white; border: none; border-radius: 5px; padding: 8px 15px; cursor: pointer; display: flex; align-items: center; gap: 5px; font-size: 0.9em; }
        .print-section-btn:hover { background: #228B22; }
        .footer { text-align: center; padding: 30px; color: white; margin-top: 40px; border-radius: 15px; background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%); }
        .footer a { color: #ffc107; text-decoration: none; }
        .footer a:hover { text-decoration: underline; }
        .sort-icon { margin-left: 5px; font-size: 0.8em; opacity: 0.6; }
        .grouping-controls { background: #f0f7f0; padding: 12px; border-radius: 8px; margin: 15px 0; display: flex; flex-wrap: wrap; gap: 10px; align-items: center; border-left: 4px solid #006400; }
        .grouping-controls label { font-weight: bold; margin-right: 5px; }
        .group-btn { background: white; border: 1px solid #006400; color: #006400; padding: 6px 12px; border-radius: 20px; cursor: pointer; font-size: 0.85em; transition: all 0.2s; }
        .group-btn:hover { background: #006400; color: white; }
        .group-btn.active { background: #006400; color: white; }
        .citation-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 15px; margin: 15px 0;
        }
        .citation-card {
            padding: 16px; border-radius: 12px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.06);
            border-left: 5px solid var(--primary-light);
            background: #fafafa;
        }
        .citation-card .title { font-weight: 700; color: var(--primary-dark); }
        .citation-card .authors { font-style: italic; color: #2c3e50; }
        .citation-card .journal { color: #2c3e50; }
        .citation-card .doi { font-size: 0.85em; color: #007bff; }
        .citation-card .copy-btn {
            background: var(--primary-light); color: white; border: none;
            padding: 3px 12px; border-radius: 30px; cursor: pointer;
            font-size: 0.75em; margin-top: 6px;
        }
        .citation-card .copy-btn:hover { background: var(--primary-dark); }
        .credit-grid {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: 15px; margin: 15px 0;
        }
        .credit-card {
            padding: 16px; border-radius: 12px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.06);
            text-align: center;
            border-top: 4px solid var(--primary-light);
        }
        .credit-card .tool-name { font-weight: 700; font-size: 1.2em; }
        .credit-card .tool-desc { font-size: 0.9em; color: #555; }
        .credit-card .tool-link a { color: var(--primary-light); text-decoration: none; font-weight: 600; }
        .feature-cards {
            display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 20px; margin: 30px 0;
        }
        .feature-card {
            background: white; padding: 20px; border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.06);
            border-top: 4px solid var(--primary-light);
            text-align: center;
        }
        .feature-card i { font-size: 2.5em; color: var(--primary-light); margin-bottom: 10px; }
        .feature-card h4 { color: var(--primary-dark); margin: 10px 0; }

        /* SAMPLE-CENTRIC BOX STYLES (from StaphScope) */
        .isolate-box {
            border: 1px solid #ddd;
            border-radius: 12px;
            margin-bottom: 30px;
            padding: 20px;
            background: #fafafa;
            box-shadow: 0 2px 8px rgba(0,0,0,0.06);
            border-bottom: 2px dashed #ccc;
        }
        .isolate-box .sample-header { display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 15px; border-bottom: 2px solid #e0e0e0; padding-bottom: 10px; }
        .isolate-box .sample-header h3 { font-size: 1.4em; margin: 0; }
        .isolate-box .sample-header .total-badge { background: #006400; color: white; padding: 4px 16px; border-radius: 20px; font-weight: bold; font-size: 0.9em; }
        .isolate-box .sample-header .typing-info { display: flex; flex-wrap: wrap; gap: 8px; margin-left: auto; }
        .isolate-box .sample-header .typing-badge { display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 0.8em; font-weight: 600; background: #e0e0e0; color: #333; border: 1px solid #ccc; }
        .isolate-box .sample-header .typing-badge.badge-ST { background: #FF9800; color: white; border-color: #FF9800; }
        .isolate-box .sample-header .typing-badge.badge-K { background: #9C27B0; color: white; border-color: #9C27B0; }
        .isolate-box .sample-header .typing-badge.badge-O { background: #009688; color: white; border-color: #009688; }
        .isolate-box .sample-header .typing-badge.badge-HV { background: #E91E63; color: white; border-color: #E91E63; }
        .isolate-box .database-table-wrapper { margin: 15px 0; overflow-x: auto; }
        .isolate-box .database-table-wrapper table { width: 100%; border-collapse: collapse; font-size: 0.85em; min-width: 800px; }
        .isolate-box .database-table-wrapper table th { background: #2c3e50; color: white; padding: 8px 12px; text-align: left; white-space: nowrap; }
        .isolate-box .database-table-wrapper table td { padding: 8px 12px; border-bottom: 1px solid #e0e0e0; white-space: nowrap; }
        .isolate-box .database-table-wrapper table tr:hover { background: #f1f1f1; }
        .isolate-box .db-title { font-weight: bold; color: #006400; margin: 10px 0 5px 0; font-size: 1.1em; border-left: 4px solid #006400; padding-left: 10px; }
        .filter-controls { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; }
        .filter-controls select { padding: 10px; border-radius: 8px; border: 2px solid #ddd; background: white; min-width: 150px; }

        @media print { body * { visibility: hidden; } .tab-content.active, .tab-content.active * { visibility: visible; } .tab-content.active { position: absolute; left: 0; top: 0; width: 100%; padding: 20px; box-shadow: none; border-radius: 0; } .print-section-btn, .tab-navigation, .dashboard-grid, .search-box, .action-buttons, .grouping-controls, .filter-controls { display: none !important; } .data-table { page-break-inside: auto; } .data-table tr { page-break-inside: avoid; page-break-after: auto; } }
        @media (max-width: 768px) { .container { padding: 10px; } .main-header { padding: 20px; } .main-header h1 { font-size: 2em; } .tab-button { padding: 8px 12px; font-size: 0.8em; } .dashboard-grid { grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); } .data-table { font-size: 0.8em; } body { min-width: auto; overflow-x: auto; } }
        </style>
        """

    def _get_js(self, typing_json: str) -> str:
        return f"""
        <script>
        var sampleTyping = {typing_json};
        var originalGenomeLists = {{}};

        function switchTab(tabName) {{
            document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(button => button.classList.remove('active'));
            document.getElementById(tabName + '-tab').classList.add('active');
            event.currentTarget.classList.add('active');
            window.location.hash = tabName;
        }}

        function searchTable(tableId, searchId) {{
            const input = document.getElementById(searchId);
            const filter = input.value.toUpperCase();
            const table = document.getElementById(tableId);
            const rows = table.getElementsByTagName('tr');
            for (let i = 1; i < rows.length; i++) {{
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {{
                    const cell = cells[j];
                    if (cell) {{
                        const txtValue = cell.textContent || cell.innerText;
                        if (txtValue.toUpperCase().indexOf(filter) > -1) {{
                            found = true;
                            break;
                        }}
                    }}
                }}
                rows[i].style.display = found ? '' : 'none';
            }}
        }}

        function highlightGenome(tableId, searchId) {{
            const filter = document.getElementById(searchId).value.toUpperCase().trim();
            const table = document.getElementById(tableId);
            const allTags = table.querySelectorAll('.genome-tag');
            allTags.forEach(tag => tag.classList.remove('highlight'));
            if (filter === '') return;
            allTags.forEach(tag => {{
                if (tag.textContent.toUpperCase().indexOf(filter) > -1) {{
                    tag.classList.add('highlight');
                }}
            }});
        }}

        function getTypingValue(genome, groupBy) {{
            var info = sampleTyping[genome];
            if (!info) return "Unknown";
            if (groupBy === "ST") return "ST" + info.ST;
            if (groupBy === "K") return info.K;
            if (groupBy === "O") return info.O;
            if (groupBy === "ST-K") return "ST" + info.ST + "-" + info.K;
            if (groupBy === "ST-O") return "ST" + info.ST + "-" + info.O;
            if (groupBy === "ST-K:O") return "ST" + info.ST + "-" + info.K + ":" + info.O;
            return "Unknown";
        }}

        function groupRowGenomes(row, groupBy, originalList) {{
            let genomesCell = null;
            for (let i = 0; i < row.cells.length; i++) {{
                if (row.cells[i].querySelector('.genome-list')) {{
                    genomesCell = row.cells[i];
                    break;
                }}
            }}
            if (!genomesCell) {{
                console.warn("Could not find genomes cell in row");
                return;
            }}
            var genomes = originalList.slice();
            if (genomes.length === 0) {{
                genomesCell.innerHTML = '<div class="genome-list">None</div>';
                return;
            }}
            var groups = {{}};
            genomes.forEach(function(genome) {{
                var key = getTypingValue(genome, groupBy);
                if (!groups[key]) groups[key] = [];
                groups[key].push(genome);
            }});
            var html = '<div class="genome-list">';
            for (var key in groups) {{
                var tags = groups[key].map(g => `<span class="genome-tag">${{g}}</span>`).join('');
                html += `<div class="genome-group"><div class="genome-group-header">${{key}}</div><div class="genome-group-tags">${{tags}}</div></div>`;
            }}
            html += '</div>';
            genomesCell.innerHTML = html;
        }}

        function groupGenomesByTyping(tableId, groupBy) {{
            var table = document.getElementById(tableId);
            if (!table) {{
                console.error("Table not found:", tableId);
                return;
            }}
            var tbody = table.tBodies[0];
            if (!tbody) {{
                console.error("No tbody found in table", tableId);
                return;
            }}
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                if (!originalGenomeLists[geneName]) {{
                    var genomesCell = null;
                    for (var j = 0; j < row.cells.length; j++) {{
                        if (row.cells[j].querySelector('.genome-list')) {{
                            genomesCell = row.cells[j];
                            break;
                        }}
                    }}
                    if (genomesCell) {{
                        var tags = genomesCell.querySelectorAll('.genome-tag');
                        var genomes = Array.from(tags).map(tag => tag.textContent.trim());
                        originalGenomeLists[geneName] = genomes;
                    }} else {{
                        originalGenomeLists[geneName] = [];
                    }}
                }}
            }}
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                groupRowGenomes(row, groupBy, original);
            }}
            var container = table.closest('.tab-content');
            if (container) {{
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
                var activeBtn = container.querySelector(`.group-btn[data-group="${{groupBy}}"]`);
                if (activeBtn) activeBtn.classList.add('active');
            }}
        }}

        function resetGenomeList(tableId) {{
            var table = document.getElementById(tableId);
            if (!table) return;
            var tbody = table.tBodies[0];
            if (!tbody) return;
            var rows = tbody.rows;
            for (var i = 0; i < rows.length; i++) {{
                var row = rows[i];
                var geneNameCell = row.cells[0];
                if (!geneNameCell) continue;
                var geneName = geneNameCell.textContent.trim().replace(/⚠️/g, '').trim();
                var original = originalGenomeLists[geneName] || [];
                var genomesCell = null;
                for (var j = 0; j < row.cells.length; j++) {{
                    if (row.cells[j].querySelector('.genome-list')) {{
                        genomesCell = row.cells[j];
                        break;
                    }}
                }}
                if (genomesCell) {{
                    var tags = original.map(g => `<span class="genome-tag">${{g}}</span>`).join('');
                    genomesCell.innerHTML = `<div class="genome-list">${{tags}}</div>`;
                }}
            }}
            var container = table.closest('.tab-content');
            if (container) {{
                var btns = container.querySelectorAll('.group-btn');
                btns.forEach(btn => btn.classList.remove('active'));
            }}
        }}

        function sortTable(tableId, colIndex, type) {{
            const table = document.getElementById(tableId);
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.rows);
            const isAscending = table.getAttribute('data-sort-dir') !== 'asc';
            rows.sort((a, b) => {{
                let aVal = a.cells[colIndex].innerText.trim();
                let bVal = b.cells[colIndex].innerText.trim();
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
            const headers = table.querySelectorAll('th');
            headers.forEach((th, idx) => {{
                const icon = th.querySelector('.sort-icon');
                if (icon) icon.innerHTML = '⇅';
            }});
            const currentHeader = headers[colIndex];
            const icon = currentHeader.querySelector('.sort-icon');
            if (icon) icon.innerHTML = isAscending ? '↑' : '↓';
        }}

        function printSection(sectionId) {{
            const content = document.getElementById(sectionId);
            const printWindow = window.open('', '_blank');
            printWindow.document.write('<html><head><title>Print Section</title>');
            printWindow.document.write('<style>' + document.querySelector('style').textContent + '</style>');
            printWindow.document.write('</head><body>');
            printWindow.document.write(content.innerHTML);
            printWindow.document.write('</body></html>');
            printWindow.document.close();
            printWindow.print();
        }}

        function exportTableToCSV(tableId, filename) {{
            const table = document.getElementById(tableId);
            const rows = table.querySelectorAll('tr');
            const csv = [];
            for (let i = 0; i < rows.length; i++) {{
                const row = [], cols = rows[i].querySelectorAll('td, th');
                for (let j = 0; j < cols.length; j++) {{
                    row.push('"' + (cols[j].innerText || '').replace(/"/g, '""') + '"');
                }}
                csv.push(row.join(','));
            }}
            const csvFile = new Blob([csv.join('\\n')], {{type: 'text/csv'}});
            const downloadLink = document.createElement('a');
            downloadLink.download = filename;
            downloadLink.href = window.URL.createObjectURL(csvFile);
            downloadLink.style.display = 'none';
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }}

        // Sample-centric box filtering
        function filterBoxes(tabId) {{
            var search = document.getElementById('search-' + tabId).value.toUpperCase();
            var dbFilter = document.getElementById('dbFilter-' + tabId).value;
            var boxes = document.querySelectorAll('#' + tabId + '-tab .isolate-box');
            boxes.forEach(function(box) {{
                var sample = box.getAttribute('data-sample') || '';
                var show = true;
                if (search && sample.toUpperCase().indexOf(search) === -1) show = false;
                var dbWrappers = box.querySelectorAll('.database-table-wrapper');
                if (dbFilter !== 'all') {{
                    dbWrappers.forEach(function(wrapper) {{
                        var dbName = wrapper.getAttribute('data-db') || '';
                        if (dbName === dbFilter) {{
                            wrapper.style.display = '';
                        }} else {{
                            wrapper.style.display = 'none';
                        }}
                    }});
                    var dbTitles = box.querySelectorAll('.db-title');
                    dbTitles.forEach(function(title) {{
                        var wrapper = title.nextElementSibling;
                        if (wrapper && wrapper.style.display === 'none') {{
                            title.style.display = 'none';
                        }} else {{
                            title.style.display = '';
                        }}
                    }});
                }} else {{
                    dbWrappers.forEach(function(wrapper) {{ wrapper.style.display = ''; }});
                    box.querySelectorAll('.db-title').forEach(function(title) {{ title.style.display = ''; }});
                }}
                box.style.display = show ? '' : 'none';
            }});
        }}

        function resetBoxFilters(tabId) {{
            document.getElementById('search-' + tabId).value = '';
            document.getElementById('dbFilter-' + tabId).value = 'all';
            filterBoxes(tabId);
        }}

        document.addEventListener('DOMContentLoaded', function() {{
            const hash = window.location.hash.substring(1);
            if (hash) {{
                const tabButton = document.querySelector(`.tab-button.${{hash}}`);
                if (tabButton) tabButton.click();
            }} else {{
                document.querySelector('.tab-button').click();
            }}
            document.querySelectorAll('.data-table').forEach(table => {{
                const headers = table.querySelectorAll('th');
                headers.forEach((header, idx) => {{
                    const type = header.getAttribute('data-sort') || 'string';
                    header.style.cursor = 'pointer';
                    header.addEventListener('click', () => sortTable(table.id, idx, type));
                    const icon = document.createElement('span');
                    icon.className = 'sort-icon';
                    icon.innerHTML = '⇅';
                    header.appendChild(icon);
                }});
            }});
            document.querySelectorAll('.copy-btn').forEach(btn => {{
                btn.addEventListener('click', function() {{
                    const citation = this.getAttribute('data-citation');
                    if (citation) {{
                        navigator.clipboard.writeText(citation).then(() => {{
                            const originalText = this.innerHTML;
                            this.innerHTML = '✓ Copied!';
                            setTimeout(() => {{ this.innerHTML = originalText; }}, 2000);
                        }});
                    }}
                }});
            }});
        }});
        </script>
        """

    # --------------------------------------------------------------------------
    # SAMPLE-CENTRIC BOX GENERATORS (core new methods)
    # --------------------------------------------------------------------------
    def _generate_sample_centric_boxes(self, data: Dict, tab_id: str, title: str, db_list: List[str]) -> str:
        samples_data = data.get('samples', {})
        amr_details = data.get('amrfinder_details', {})
        abricate_details = data.get('abricate_details', {})

        relevant_samples = []
        for sample in samples_data:
            has_data = False
            if 'amrfinder' in db_list and sample in amr_details and amr_details[sample]:
                has_data = True
            for db in db_list:
                if db != 'amrfinder' and sample in abricate_details and db in abricate_details[sample] and abricate_details[sample][db]:
                    has_data = True
            if has_data:
                relevant_samples.append(sample)
        relevant_samples.sort()

        if not relevant_samples:
            return f"""
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No {title} Data Available</h3><p>No samples with {title} genes were found.</p></div>
            </div>
            """

        html = f"""
        <div class="section-header {tab_id}-header"><h2><i class="fas fa-{self._icon_for_tab(tab_id)}"></i> {title} – Interactive Isolate Boxes</h2><button class="print-section-btn" onclick="printSection('{tab_id}-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 {title} – Sample‑Centric View</h3>
                <p>Each box represents one isolate. Inside, you will find separate tables for each database with full gene details – horizontally scrollable.</p>
                <p>Use the filters below to search by sample name or to show only a specific database.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}" onkeyup="filterBoxes('{tab_id}')" placeholder="🔍 Search sample...">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">
                <option value="all">All Databases</option>
        """
        for db in db_list:
            display = db.upper()
            if db == 'amrfinder':
                display = 'AMRfinder'
            html += f'<option value="{db}">{display}</option>'
        html += f"""
            </select>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-{tab_id}">
        """

        for sample in relevant_samples:
            # Gather typing info from samples_data
            st = samples_data.get(sample, {}).get('mlst', {}).get('ST', 'ND')
            k = samples_data.get(sample, {}).get('kaptive', {}).get('K_Locus', 'ND')
            o = samples_data.get(sample, {}).get('kaptive', {}).get('O_Locus', 'ND')
            # Hypervirulence status from abricate genes
            vir_genes = []
            for db, genes in samples_data.get(sample, {}).get('abricate_genes', {}).items():
                if db in ['vfdb', 'ecoli_vf']:
                    vir_genes.extend(genes)
            hv_status = 'Hypervirulent' if any(g in self.analyzer.high_risk_virulence_genes or any(marker in g for marker in ['ybt','clb','iro','iuc','rmp']) for g in vir_genes) else ''

            html += f'<div class="isolate-box" data-sample="{sample}">'
            total_genes = 0
            if 'amrfinder' in db_list and sample in amr_details:
                total_genes += len(amr_details[sample])
            for db in db_list:
                if db != 'amrfinder' and sample in abricate_details and db in abricate_details[sample]:
                    total_genes += len(abricate_details[sample][db])

            html += f"""
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Genes: {total_genes}</span>
                    <div class="typing-info">
                        <span class="typing-badge badge-ST">ST: {st}</span>
                        <span class="typing-badge badge-K">K: {k}</span>
                        <span class="typing-badge badge-O">O: {o}</span>
                        {f'<span class="typing-badge badge-HV">Hypervirulent</span>' if hv_status else ''}
                    </div>
                </div>
            """

            if 'amrfinder' in db_list and sample in amr_details and amr_details[sample]:
                html += self._make_database_table('AMRfinder', amr_details[sample], 'amrfinder')

            for db in db_list:
                if db == 'amrfinder':
                    continue
                if sample in abricate_details and db in abricate_details[sample] and abricate_details[sample][db]:
                    html += self._make_database_table(db.upper(), abricate_details[sample][db], db)

            html += '</div>'

        html += """
        </div>
        """
        return html

    def _make_database_table(self, db_name: str, genes: List[Dict], db_key: str) -> str:
        if not genes:
            return ''
        # Determine columns from first gene dict
        keys = list(genes[0].keys()) if genes else []
        # Prefer certain order
        priority = ['gene', 'product', 'coverage_percent', 'identity_percent', 'accession', 'contig', 'start', 'end', 'strand', 'class', 'subclass', 'scope', 'resistance']
        ordered = []
        for p in priority:
            if p in keys:
                ordered.append(p)
        for k in keys:
            if k not in ordered:
                ordered.append(k)

        html = f"""
                <div class="database-table-wrapper" data-db="{db_key}">
                    <div class="db-title">{db_name}</div>
                    <table>
                        <thead>
                            <tr>
        """
        for col in ordered:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for gene_dict in genes:
            html += '<tr>'
            for col in ordered:
                val = gene_dict.get(col, '')
                if val is None:
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                        </tbody>
                    </table>
                </div>
        """
        return html

    def _generate_mutation_boxes(self, data: Dict) -> str:
        samples_data = data.get('samples', {})
        mutation_details = data.get('mutation_details', {})
        relevant_samples = [s for s in samples_data if s in mutation_details and mutation_details[s]]
        relevant_samples.sort()

        if not relevant_samples:
            return """
            <div class="section-header mutations-header"><h2><i class="fas fa-dna"></i> Mutations – Interactive Isolate Boxes</h2></div>
            <div class="alert-box alert-warning">
                <i class="fas fa-exclamation-circle fa-2x"></i>
                <div><h3>No Mutation Data Available</h3><p>No mutations found for any sample.</p></div>
            </div>
            """

        html = f"""
        <div class="section-header mutations-header"><h2><i class="fas fa-dna"></i> Mutations – Interactive Isolate Boxes</h2><button class="print-section-btn" onclick="printSection('mutations-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info">
            <i class="fas fa-info-circle fa-2x"></i>
            <div>
                <h3>🧬 Point Mutations – Sample‑Centric View</h3>
                <p>Each box represents one isolate. Inside, you will find a table of all detected point mutations with full details – horizontally scrollable.</p>
                <p>Use the search bar below to filter boxes by sample name.</p>
            </div>
        </div>
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-mutation" onkeyup="filterBoxes('mutation')" placeholder="🔍 Search sample...">
            <button class="action-btn btn-success" onclick="resetBoxFilters('mutation')"><i class="fas fa-sync"></i> Clear Filters</button>
        </div>
        <div id="box-container-mutation">
        """

        for sample in relevant_samples:
            st = samples_data.get(sample, {}).get('mlst', {}).get('ST', 'ND')
            k = samples_data.get(sample, {}).get('kaptive', {}).get('K_Locus', 'ND')
            o = samples_data.get(sample, {}).get('kaptive', {}).get('O_Locus', 'ND')
            mutations = mutation_details.get(sample, [])
            total_muts = len(mutations)

            html += f'<div class="isolate-box" data-sample="{sample}">'
            html += f"""
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">Total Mutations: {total_muts}</span>
                    <div class="typing-info">
                        <span class="typing-badge badge-ST">ST: {st}</span>
                        <span class="typing-badge badge-K">K: {k}</span>
                        <span class="typing-badge badge-O">O: {o}</span>
                    </div>
                </div>
            """
            if mutations:
                html += self._make_mutation_table(mutations)
            html += '</div>'

        html += """
        </div>
        """
        return html

    def _make_mutation_table(self, mutations: List[Dict]) -> str:
        if not mutations:
            return ''
        cols = ['gene', 'mutation', 'class', 'subclass', 'contig', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
        html = """
        <div class="database-table-wrapper" data-db="mutations">
            <div class="db-title">Mutations</div>
            <table>
                <thead>
                    <tr>
        """
        for col in cols:
            display = col.replace('_', ' ').title()
            html += f'<th>{display}</th>'
        html += '</tr></thead><tbody>'
        for mut in mutations:
            html += '<tr>'
            for col in cols:
                val = mut.get(col, '')
                if val is None or (isinstance(val, float) and val != val):
                    val = ''
                html += f'<td>{val}</td>'
            html += '</tr>'
        html += """
                </tbody>
            </table>
        </div>
        """
        return html

    def _icon_for_tab(self, tab_id: str) -> str:
        mapping = {
            'amr': 'biohazard',
            'vir': 'virus',
            'bac': 'flask',
            'plasmids': 'dna',
            'mutations': 'dna'
        }
        return mapping.get(tab_id, 'box')

    # --------------------------------------------------------------------------
    # ALL OTHER SECTION GENERATORS 
    # --------------------------------------------------------------------------
    def _generate_summary_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        patterns = data.get('patterns', {})
        gene_centric = data.get('gene_centric', {})
        total = len(samples)
        st_dist = patterns.get('st_distribution', {})
        k_dist = patterns.get('k_locus_distribution', {})
        high_risk = len(patterns.get('high_risk_combinations', []))
        most_common_st = max(st_dist.items(), key=lambda x: x[1])[0] if st_dist else 'None'
        most_common_k = max(k_dist.items(), key=lambda x: x[1])[0] if k_dist else 'None'
        carbapenemase_count = sum(1 for g in gene_centric.get('all_genes', []) if 'blaKPC' in g['gene'] or 'blaNDM' in g['gene'] or 'OXA-48' in g['gene'])
        hv_count = len(gene_centric.get('by_category', {}).get('High-Risk Virulence', []))

        about = f"""
        <div class="section-header" style="border-bottom-color: #2c7a4d; margin-top: 40px;">
            <h2><i class="fas fa-info-circle"></i> About This Report</h2>
        </div>
        <div class="alert-box alert-info" style="border-left-color: #2c7a4d;">
            <i class="fas fa-lightbulb" style="font-size:2em;"></i>
            <div>
                <h3>What you can do with this report</h3>
                <ul style="margin:10px 0 0 20px;">
                    <li><strong>Explore the population structure</strong> – MLST, capsule types (K/O), and sample overview.</li>
                    <li><strong>Track resistance and virulence per isolate</strong> – AMR, Virulence, Bacmet, Plasmids, and Mutations tabs now show each sample in its own interactive box with detailed tables.</li>
                    <li><strong>Discover patterns</strong> – Combinations tab reveals ST–K:O associations; Patterns tab shows high‑risk combos and co‑occurrences.</li>
                    <li><strong>Export and share</strong> – Each table can be exported as CSV, and the full data as JSON.</li>
                </ul>
            </div>
        </div>
        <div class="feature-cards">
            <div class="feature-card"><i class="fas fa-box"></i><h4>Sample‑Centric Boxes</h4><p>Each isolate gets a detailed box with all its genes/mutations and typing badges – perfect for clinical reports.</p></div>
            <div class="feature-card"><i class="fas fa-layer-group"></i><h4>Gene‑Centric Views</h4><p>MLST, Kaptive, and Combination tabs remain gene‑centric for cross‑genome comparisons.</p></div>
            <div class="feature-card"><i class="fas fa-flask"></i><h4>Bacmet &amp; Plasmids</h4><p>Track biocide/heavy metal resistance and plasmid replicons – key for hospital hygiene and horizontal gene transfer.</p></div>
        </div>
        """

        return f"""
        <div class="section-header summary-header"><h2><i class="fas fa-chart-pie"></i> Executive Summary</h2><button class="print-section-btn" onclick="printSection('summary-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><strong>Population overview</strong> – {total} genomes analysed. The most common ST is <strong>ST{most_common_st}</strong> and the most common K‑locus is <strong>{most_common_k}</strong>. Use the tabs to explore in depth.</div></div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin:15px 0;">
            <div style="background:var(--summary-color);color:white;padding:12px;border-radius:10px;"><div style="font-size:1.8em;font-weight:bold;">{total}</div><div>Total Samples</div></div>
            <div style="background:#FF9800;color:white;padding:12px;border-radius:10px;"><div style="font-size:1.8em;font-weight:bold;">ST{most_common_st}</div><div>Most Common ST</div></div>
            <div style="background:#E91E63;color:white;padding:12px;border-radius:10px;"><div style="font-size:1.8em;font-weight:bold;">{most_common_k}</div><div>Most Common K Locus</div></div>
            <div style="background:#F44336;color:white;padding:12px;border-radius:10px;"><div style="font-size:1.8em;font-weight:bold;">{carbapenemase_count}</div><div>Carbapenemase Genes</div></div>
            <div style="background:#E91E63;color:white;padding:12px;border-radius:10px;"><div style="font-size:1.8em;font-weight:bold;">{hv_count}</div><div>High‑Risk Virulence Genes</div></div>
            <div style="background:#FF5722;color:white;padding:12px;border-radius:10px;"><div style="font-size:1.8em;font-weight:bold;">{high_risk}</div><div>High‑Risk Combos</div></div>
        </div>
        {f'<div class="alert-box alert-danger"><i class="fas fa-exclamation-triangle"></i><div><strong>⚠️ Carbapenemase genes detected!</strong> Check the AMR tab for details.</div></div>' if carbapenemase_count > 0 else ''}
        {f'<div class="alert-box alert-warning"><i class="fas fa-virus"></i><div><strong>High‑risk virulence genes present.</strong> See Virulence tab for details.</div></div>' if hv_count > 0 else ''}
        {about}
        """

    def _generate_sample_overview_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        rows = []
        for sample, d in sorted(samples.items()):
            st = d.get('mlst', {}).get('ST', 'ND')
            st_display = f"ST{st}" if st != 'ND' else 'ND'
            k = d.get('kaptive', {}).get('K_Locus', 'ND')
            o = d.get('kaptive', {}).get('O_Locus', 'ND')
            vir_count = sum(len(genes) for db, genes in d.get('abricate_genes', {}).items() if db in ['vfdb', 'ecoli_vf'])
            rows.append(f"<tr><td><span class='genome-tag'>{sample}</span></td><td>{st_display}</td><td>{k}</td><td>{o}</td><td>{vir_count}</td></tr>")
        return f"""
        <div class="section-header sample_overview-header"><h2><i class="fas fa-list"></i> Sample Overview</h2><button class="print-section-btn" onclick="printSection('sample_overview-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><strong>What this tab shows:</strong> Each sample with its MLST, K‑locus, O‑locus, and total number of virulence genes detected. MLST uses seven housekeeping genes (<em>gapA, infB, mdh, pgi, phoE, rpoB, tonB</em>) and is the gold standard for global epidemiology. K and O loci define the capsule and lipopolysaccharide antigens, which are key for immune evasion and virulence. Click column headers to sort.</div></div>
        <input type="text" class="search-box" id="search-samples" onkeyup="searchTable('samples-table','search-samples')" placeholder="🔍 Search by sample, ST, K, or O...">
        <input type="text" class="search-box" id="highlight-samples" onkeyup="highlightGenome('samples-table','highlight-samples')" placeholder="🔍 Highlight genome tags (e.g., sample name)">
        <div class="action-buttons"><button class="action-btn btn-primary" onclick="exportTableToCSV('samples-table','sample_overview.csv')"><i class="fas fa-download"></i> Export CSV</button></div>
        <div class="master-scrollable-container"><table id="samples-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">ST</th><th data-sort="string">K Locus</th><th data-sort="string">O Locus</th><th data-sort="number">Virulence Count</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_mlst_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        st_dist = patterns.get('st_distribution', {})
        if not st_dist:
            return "<div class='alert-box alert-warning'>No MLST data available.</div>"
        rows = []
        total = sum(st_dist.values())
        for st, cnt in sorted(st_dist.items(), key=lambda x: x[1], reverse=True):
            pct = cnt / total * 100
            st_display = f"ST{st}" if st != 'ND' else 'ND'
            rows.append(f"<tr><td><strong>{st_display}</strong></td><td>{cnt}</td><td>{pct:.1f}%</td></tr>")
        return f"""
        <div class="section-header mlst-header"><h2><i class="fas fa-code-branch"></i> MLST Distribution</h2><button class="print-section-btn" onclick="printSection('mlst-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><strong>MLST (Multi‑Locus Sequence Typing)</strong> is based on internal fragments of seven housekeeping genes: <em>gapA, infB, mdh, pgi, phoE, rpoB, tonB</em>. Each unique combination of alleles defines a Sequence Type (ST). Closely related STs belong to the same clonal complex (CC). This is the standard for global surveillance of <em>K. pneumoniae</em>.</div></div>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th data-sort="string">ST</th><th data-sort="number">Count</th><th data-sort="number">Percentage</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_qc_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        all_metrics = set()
        for d in samples.values():
            qc = d.get('qc', {})
            all_metrics.update(qc.keys())
        if not all_metrics:
            return "<div class='alert-box alert-warning'>No QC data available.</div>"
        metric_list = sorted(all_metrics)
        rows = []
        for sample, d in sorted(samples.items()):
            qc = d.get('qc', {})
            row = [f"<td><span class='genome-tag'>{sample}</span>"]
            for m in metric_list:
                val = qc.get(m, 'ND')
                if isinstance(val, float):
                    val = f"{val:,.0f}" if val > 1000 else f"{val:.1f}"
                row.append(f"<td>{val}")
            rows.append("<tr>" + "".join(row) + "</tr>")
        header = "<thead><tr><th data-sort='string'>Sample</th>" + "".join([f"<th data-sort='number'>{m}</th>" for m in metric_list]) + "</tr></thead>"
        return f"""
        <div class="section-header qc-header"><h2><i class="fas fa-chart-line"></i> FASTA QC</h2><button class="print-section-btn" onclick="printSection('qc-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><strong>Assembly quality metrics</strong> – N50 (contig length at 50% of assembly), total length, GC%, and contig count. Higher N50 and fewer contigs indicate better assembly quality. Poor assemblies can miss genes or break up plasmids.</div></div>
        <input type="text" class="search-box" id="search-qc" onkeyup="searchTable('qc-table','search-qc')" placeholder="🔍 Search sample...">
        <input type="text" class="search-box" id="highlight-qc" onkeyup="highlightGenome('qc-table','highlight-qc')" placeholder="🔍 Highlight samples...">
        <div class="master-scrollable-container"><table id="qc-table" class="data-table">{header}<tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_kaptive_section(self, data: Dict) -> str:
        samples = data.get('samples', {})
        rows = []
        for sample, d in sorted(samples.items()):
            k = d.get('kaptive', {})
            if k:
                rows.append(f"<tr><td><span class='genome-tag'>{sample}</span></td><td>{k.get('K_Locus', 'ND')}</td><td>{k.get('O_Locus', 'ND')}</td><td>{k.get('K_Identity', 'ND')}</td><td>{k.get('K_Coverage', 'ND')}</td></tr>")
        if not rows:
            return "<div class='alert-box alert-warning'>No Kaptive data available.</div>"
        return f"""
        <div class="section-header kaptive-header"><h2><i class="fas fa-shield-alt"></i> Kaptive Capsule Typing</h2><button class="print-section-btn" onclick="printSection('kaptive-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><strong>Kaptive</strong> identifies the K (capsule) and O (lipopolysaccharide) loci. The capsule protects <em>K. pneumoniae</em> from phagocytosis and complement killing. Specific K types (e.g., K1, K2, K5) are associated with hypervirulence. O‐locus types influence serum resistance. High identity/coverage (≥95%) indicate reliable typing.</div></div>
        <input type="text" class="search-box" id="search-kaptive" onkeyup="searchTable('kaptive-table','search-kaptive')" placeholder="🔍 Search...">
        <input type="text" class="search-box" id="highlight-kaptive" onkeyup="highlightGenome('kaptive-table','highlight-kaptive')" placeholder="🔍 Highlight samples...">
        <div class="master-scrollable-container"><table id="kaptive-table" class="data-table"><thead><tr><th data-sort="string">Sample</th><th data-sort="string">K Locus</th><th data-sort="string">O Locus</th><th data-sort="string">K Identity</th><th data-sort="number">K Coverage</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_combinations_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        combos = [
            ('st_k_combinations', 'ST – K Locus', 'ST-K'),
            ('st_o_combinations', 'ST – O Locus', 'ST-O'),
            ('ko_combinations', 'K : O Capsule Type', 'K:O'),
            ('st_ko_combinations', 'ST – K : O', 'ST-K:O')
        ]
        html_parts = []
        for key, title, label in combos:
            combo_dict = patterns.get(key, {})
            if not combo_dict:
                continue
            rows = []
            for combo, samples in sorted(combo_dict.items(), key=lambda x: len(x[1]), reverse=True):
                sample_tags = ''.join([f'<span class="genome-tag">{s}</span>' for s in samples])
                rows.append(f"<tr><td><strong>{combo}</strong></td><td>{len(samples)}</td><td><div class='genome-list'>{sample_tags}</div></td></tr>")
            html_parts.append(f"""
            <h3>{title}</h3>
            <input type="text" class="search-box" id="search-{key}" onkeyup="searchTable('{key}-table','search-{key}')" placeholder="🔍 Search {label}...">
            <input type="text" class="search-box" id="highlight-{key}" onkeyup="highlightGenome('{key}-table','highlight-{key}')" placeholder="🔍 Highlight genomes...">
            <div class="master-scrollable-container"><table id="{key}-table" class="data-table"><thead><tr><th data-sort="string">{label}</th><th data-sort="number">Count</th><th data-sort="string">Samples</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
            """)
        if not html_parts:
            return "<div class='alert-box alert-warning'>No combination data available.</div>"
        return f"""
        <div class="section-header combinations-header"><h2><i class="fas fa-link"></i> Combination Tables</h2><button class="print-section-btn" onclick="printSection('combinations-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><strong>Combinations</strong> reveal associations between ST, K‑locus, and O‑locus. For example, ST23 often carries K1, and ST258 is associated with K64. These associations can indicate hypervirulent clones or successful epidemic lineages.</div></div>
        {"".join(html_parts)}
        """

    def _generate_patterns_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        icekp = patterns.get('icekp_marker_presence', {})
        vp = patterns.get('virulence_plasmid_marker_presence', {})
        cooc = patterns.get('gene_cooccurrence', {})
        html = "<div class='section-header patterns-header'><h2><i class='fas fa-project-diagram'></i> Patterns & Associations</h2><button class='print-section-btn' onclick=\"printSection('patterns-tab')\"><i class='fas fa-print'></i> Print</button></div>"
        html += "<div class='alert-box alert-info'><i class='fas fa-info-circle'></i><div><strong>This tab reveals associations</strong> between typing and gene content, including ICEKp markers, virulence plasmid markers, and co‑occurrence patterns. For high‑risk combinations (critical resistance + virulence), see the dedicated <strong>High Risk</strong> tab.</div></div>"
        if icekp:
            rows = []
            for sample, markers in icekp.items():
                rows.append(f"<tr><td><strong>{sample}</strong></td><td>{', '.join(markers)}</td></tr>")
            html += f"""
            <h3>ICEKp Associated Markers (ybt, clb, iro, rmp)</h3>
            <input type="text" class="search-box" id="search-icekp" onkeyup="searchTable('icekp-table','search-icekp')" placeholder="🔍 Search...">
            <input type="text" class="search-box" id="highlight-icekp" onkeyup="highlightGenome('icekp-table','highlight-icekp')" placeholder="🔍 Highlight samples...">
            <div class='master-scrollable-container'><table id='icekp-table' class='data-table'><thead><tr><th>Sample</th><th>Markers</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
            """
        if vp:
            rows = []
            for sample, markers in vp.items():
                rows.append(f"<tr><td><strong>{sample}</strong></td><td>{', '.join(markers)}</td></tr>")
            html += f"""
            <h3>Virulence Plasmid Markers (iro, iuc, rmp, rmpA2)</h3>
            <input type="text" class="search-box" id="search-vp" onkeyup="searchTable('vp-table','search-vp')" placeholder="🔍 Search...">
            <input type="text" class="search-box" id="highlight-vp" onkeyup="highlightGenome('vp-table','highlight-vp')" placeholder="🔍 Highlight samples...">
            <div class='master-scrollable-container'><table id='vp-table' class='data-table'><thead><tr><th>Sample</th><th>Markers</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
            """
        if cooc:
            cooc_list = []
            for g1, partners in cooc.items():
                for g2, cnt in partners.items():
                    cooc_list.append((g1, g2, cnt))
            cooc_list.sort(key=lambda x: x[2], reverse=True)
            if cooc_list:
                rows = []
                for g1, g2, cnt in cooc_list[:200]:
                    rows.append(f"<tr><td>{g1}</td><td>{g2}</td><td>{cnt}</td></tr>")
                html += f"""
                <h3>Gene Co‑occurrence (Top 200)</h3>
                <div class='master-scrollable-container'><table class='data-table'><thead><tr><th data-sort='string'>Gene 1</th><th data-sort='string'>Gene 2</th><th data-sort='number'>Count</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
                """
        return html

    def _generate_highrisk_section(self, data: Dict) -> str:
        patterns = data.get('patterns', {})
        high_risk = patterns.get('high_risk_combinations', [])
        if not high_risk:
            return """
            <div class="section-header highrisk-header"><h2><i class="fas fa-exclamation-triangle"></i> High‑Risk Combinations</h2></div>
            <div class="alert-box alert-success"><i class="fas fa-check-circle"></i><div>No high‑risk combinations (critical resistance + high‑risk virulence) detected.</div></div>
            """
        rows = []
        for c in high_risk:
            st_display = f"ST{c['st']}" if c['st'] != 'ND' else 'ND'
            rows.append(f"<tr><td><strong>{c['sample']}</strong></td><td>{st_display}</td><td>{c['k_locus']}</td><td>{', '.join(c['critical_resistance'])}</td><td>{', '.join(c['high_risk_virulence'])}</td></tr>")
        return f"""
        <div class="section-header highrisk-header"><h2><i class="fas fa-exclamation-triangle"></i> High‑Risk Combinations</h2><button class="print-section-btn" onclick="printSection('highrisk-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-danger"><i class='fas fa-radiation'></i><div><strong>{len(high_risk)} high‑risk combinations</strong> – samples carrying both critical resistance and high‑risk virulence genes.</div></div>
        <input type="text" class="search-box" id="search-highrisk" onkeyup="searchTable('highrisk-table','search-highrisk')" placeholder="🔍 Search...">
        <input type="text" class="search-box" id="highlight-highrisk" onkeyup="highlightGenome('highrisk-table','highlight-highrisk')" placeholder="🔍 Highlight samples...">
        <div class='master-scrollable-container'><table id='highrisk-table' class='data-table'><thead><tr><th>Sample</th><th>ST</th><th>K Locus</th><th>Critical Resistance</th><th>High‑Risk Virulence</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_databases_section(self, data: Dict) -> str:
        gene_centric = data.get('gene_centric', {})
        patterns = data.get('patterns', {})
        coverage = patterns.get('database_coverage', {})
        stats = gene_centric.get('database_stats', {})
        if not stats and not coverage:
            return "<div class='alert-box alert-warning'>No database statistics available.</div>"
        all_dbs = sorted(set(stats.keys()) | set(coverage.keys()))
        rows = []
        for db in all_dbs:
            s = stats.get(db, {})
            c = coverage.get(db, {})
            total_samples = c.get('total_samples', 0)
            samples_with_hits = c.get('samples_with_hits', s.get('total_occurrences', 0))
            cov_pct = (samples_with_hits / total_samples * 100) if total_samples else 0
            rows.append(f"""
            <tr>
                <td><strong>{db.upper()}</strong></td>
                <td>{s.get('total_genes', 0)}</td>
                <td>{s.get('total_occurrences', 0)}</td>
                <td>{samples_with_hits} / {total_samples}</td>
                <td>{cov_pct:.1f}%</td>
            </tr>
            """)
        return f"""
        <div class="section-header databases-header"><h2><i class="fas fa-database"></i> Database Coverage & Statistics</h2><button class="print-section-btn" onclick="printSection('databases-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-info-circle"></i><div><strong>Database coverage</strong> – shows how many samples have hits per database, total unique genes, and total occurrences. This helps evaluate the completeness of your screening.</div></div>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Database</th><th>Unique Genes</th><th>Total Occurrences</th><th>Samples with Hits</th><th>Coverage (%)</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
        """

    def _generate_credit_section(self, data: Dict) -> str:
        return f"""
        <div class="section-header credit-header"><h2><i class="fas fa-thumbs-up"></i> Credits & Acknowledgments</h2><button class="print-section-btn" onclick="printSection('credit-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-info"><i class="fas fa-heart"></i><div><strong>Kleboscope</strong> is built on the shoulders of many excellent open‑source tools and databases. We gratefully acknowledge the following:</div></div>
        <div class="credit-grid">
            <div class="credit-card" style="border-top-color:#1e5a3a;"><div class="tool-name">🧬 Kleboscope</div><div class="tool-desc">This pipeline – integrating everything into a user‑friendly report.</div><div class="tool-link"><a href="https://github.com/bbeckley-hub/kleboscope" target="_blank">GitHub</a></div></div>
            <div class="credit-card" style="border-top-color:#006400;"><div class="tool-name">📦 ABRicate</div><div class="tool-desc">Mass screening of contigs for resistance and virulence genes.</div><div class="tool-link"><a href="https://github.com/tseemann/abricate" target="_blank">GitHub</a></div></div>
            <div class="credit-card" style="border-top-color:#FF9800;"><div class="tool-name">🔬 MLST (PubMLST)</div><div class="tool-desc">The PubMLST database for <em>K. pneumoniae</em> typing schemes.</div><div class="tool-link"><a href="https://pubmlst.org/" target="_blank">PubMLST</a></div></div>
            <div class="credit-card" style="border-top-color:#9C27B0;"><div class="tool-name">🛡️ Kaptive</div><div class="tool-desc">Capsule and O‑antigen typing for <em>K. pneumoniae</em>.</div><div class="tool-link"><a href="https://github.com/katholt/Kaptive" target="_blank">GitHub</a></div></div>
            <div class="credit-card" style="border-top-color:#F44336;"><div class="tool-name">🧪 AMRFinderPlus</div><div class="tool-desc">Antimicrobial resistance gene detection.</div><div class="tool-link"><a href="https://www.ncbi.nlm.nih.gov/pathogens/amr/" target="_blank">NCBI</a></div></div>
            <div class="credit-card" style="border-top-color:#17a2b8;"><div class="tool-name">📚 CARD</div><div class="tool-desc">Comprehensive Antibiotic Resistance Database.</div><div class="tool-link"><a href="https://card.mcmaster.ca/" target="_blank">CARD</a></div></div>
            <div class="credit-card" style="border-top-color:#28a745;"><div class="tool-name">🔍 ResFinder</div><div class="tool-desc">Resistance gene database.</div><div class="tool-link"><a href="https://genepi.food.dtu.dk/resfinder" target="_blank">ResFinder</a></div></div>
            <div class="credit-card" style="border-top-color:#E91E63;"><div class="tool-name">🦠 VFDB</div><div class="tool-desc">Virulence Factors Database.</div><div class="tool-link"><a href="http://www.mgc.ac.cn/VFs/" target="_blank">VFDB</a></div></div>
            <div class="credit-card" style="border-top-color:#673AB7;"><div class="tool-name">🧬 PlasmidFinder</div><div class="tool-desc">Plasmid replicon detection.</div><div class="tool-link"><a href="https://genepi.food.dtu.dk/plasmidfinder" target="_blank">PlasmidFinder</a></div></div>
            <div class="credit-card" style="border-top-color:#FF5722;"><div class="tool-name">🧪 BacMet</div><div class="tool-desc">Biocide and metal resistance genes.</div><div class="tool-link"><a href="http://bacmet.biomedicine.gu.se/" target="_blank">BacMet</a></div></div>
            <div class="credit-card" style="border-top-color:#3F51B5;"><div class="tool-name">🌊 MEGARes</div><div class="tool-desc">Antimicrobial, biocide, and metal resistance.</div><div class="tool-link"><a href="https://megares.meglab.org/" target="_blank">MEGARes</a></div></div>
            <div class="credit-card" style="border-top-color:#795548;"><div class="tool-name">🔬 ARG-ANNOT</div><div class="tool-desc">Antibiotic resistance gene database.</div><div class="tool-link"><a href="https://www.mediterranee-infection.com/arg-annot/" target="_blank">ARG-ANNOT</a></div></div>
            <div class="credit-card" style="border-top-color:#00BCD4;"><div class="tool-name">🐍 Biopython</div><div class="tool-desc">Python tools for computational biology.</div><div class="tool-link"><a href="https://biopython.org/" target="_blank">Biopython</a></div></div>
        </div>
        <div class="alert-box alert-success"><i class="fas fa-lightbulb"></i><div><strong>All tools are part of the ESCAPE AMR initiative</strong> – open‑source genomic surveillance for ESKAPE pathogens. Ideas and code are shared between tools to improve each other. We thank the open‑source community and all contributors.</div></div>
        """

    def _generate_aiguide_section(self, data: Dict) -> str:
        return """
        <div class="section-header aiguide-header"><h2><i class="fas fa-robot"></i> AI Assistant Guide</h2></div>
        <div class="alert-box alert-info" style="background:#e3f2fd;border-left-color:#0d47a1;"><i class="fas fa-lightbulb" style="color:#0d47a1;"></i><div><strong>🤖 Upload the JSON file</strong> to ChatGPT, Claude, or Gemini and ask questions about your K. pneumoniae dataset. <span style="font-size:1.2em;">😉</span></div></div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:15px;margin:15px 0;">
            <div style="background:#fce4ec;padding:14px;border-radius:10px;"><i class="fas fa-chart-bar" style="color:#c62828;"></i> <strong>Population structure</strong><br>"What are the most common STs? Which K types dominate?"</div>
            <div style="background:#e8f5e9;padding:14px;border-radius:10px;"><i class="fas fa-biohazard" style="color:#2e7d32;"></i> <strong>Resistance patterns</strong><br>"Which samples carry blaKPC? Are there any with mcr?"</div>
            <div style="background:#fff3e0;padding:14px;border-radius:10px;"><i class="fas fa-virus" style="color:#e65100;"></i> <strong>Virulence</strong><br>"Which STs have ICEKp markers? List samples with rmpA and iuc."</div>
            <div style="background:#e0f7fa;padding:14px;border-radius:10px;"><i class="fas fa-project-diagram" style="color:#00695c;"></i> <strong>Combinations</strong><br>"Show me ST‑K:O combos with carbapenemases."</div>
            <div style="background:#f3e5f5;padding:14px;border-radius:10px;"><i class="fas fa-dna" style="color:#4a148c;"></i> <strong>Mutations</strong><br>"List all gyrA mutations and their STs."</div>
            <div style="background:#ffebee;padding:14px;border-radius:10px;"><i class="fas fa-exclamation-triangle" style="color:#b71c1c;"></i> <strong>High‑risk</strong><br>"Identify samples with both carbapenemase and colistin resistance."</div>
        </div>
        <div class="alert-box alert-warning" style="background:#fff8e1;border-left-color:#ffa000;"><i class="fas fa-laugh-beam" style="color:#ff6f00;"></i><div><strong>⚠️ Warning:</strong> AI is like a hyper‑caffeinated grad student – great at spotting patterns but may hallucinate. Always verify critical findings!</div></div>
        """

    def _generate_citation_section(self, data: Dict) -> str:
        citations = [
            {"title": "Kleboscope", "authors": "Beckley B et al.", "journal": "GitHub", "url": "https://github.com/bbeckley-hub/kleboscope", "year": "2026"},
            {"title": "MLST (PubMLST / BIGSdb)", "authors": "Jolley KA, Bray JE, Maiden MCJ", "journal": "Wellcome Open Res, 2018;3:124", "doi": "10.12688/wellcomeopenres.14826.1"},
            {"title": "Kaptive (capsule typing)", "authors": "Wyres KL, et al.", "journal": "Microb Genom, 2016;2(5):e000064", "doi": "10.1099/mgen.0.000064"},
            {"title": "AMRFinderPlus", "authors": "Feldgarden M, et al.", "journal": "Sci Rep, 2021;11(1):12728", "doi": "10.1038/s41598-021-91456-0"},
            {"title": "ABRicate (mass screening)", "authors": "Seemann T", "journal": "GitHub, 2024", "url": "https://github.com/tseemann/abricate"},
            {"title": "CARD", "authors": "McArthur AG, et al.", "journal": "Antimicrob Agents Chemother, 2013;57(7):3348-57", "doi": "10.1128/AAC.00419-13"},
            {"title": "ResFinder", "authors": "Florensa AF, et al.", "journal": "Microb Genom, 2022;8(1):000748", "doi": "10.1099/mgen.0.000748"},
            {"title": "VFDB", "authors": "Chen L, et al.", "journal": "Nucleic Acids Res, 2012;40(D1):D641-5", "doi": "10.1093/nar/gkr989"},
            {"title": "PlasmidFinder", "authors": "Carattoli A, et al.", "journal": "Antimicrob Agents Chemother, 2014;58(7):3895-903", "doi": "10.1128/AAC.02412-14"},
            {"title": "BacMet", "authors": "Pal C, et al.", "journal": "Nucleic Acids Res, 2014;42(D1):D737-43", "doi": "10.1093/nar/gkt1252"},
            {"title": "MEGARes 2.0", "authors": "Doster E, et al.", "journal": "Nucleic Acids Res, 2020;48(D1):D561-D569", "doi": "10.1093/nar/gkz1010"},
            {"title": "ARG-ANNOT", "authors": "Gupta SK, et al.", "journal": "Antimicrob Agents Chemother, 2014;58(1):212-20", "doi": "10.1128/AAC.01310-13"},
            {"title": "Biopython", "authors": "Cock PJ, et al.", "journal": "Bioinformatics, 2009;25(11):1422-3", "doi": "10.1093/bioinformatics/btp163"},
        ]
        colors = ['#e3f2fd', '#e8f5e9', '#fff3e0', '#fce4ec', '#f3e5f5', '#e0f7fa', '#fff8e1', '#fbe9e7', '#e0f2f1', '#f1f8e9', '#ede7f6', '#f9fbe7', '#efebe9']
        cards = []
        for i, c in enumerate(citations):
            col = colors[i % len(colors)]
            if 'doi' in c and c['doi']:
                link = f"https://doi.org/{c['doi']}"
                label = c['doi']
            elif 'url' in c and c['url']:
                link = c['url']
                label = c['url']
            else:
                link = '#'
                label = 'No link available'
            citation_text = f"{c['title']} – {c['authors']}, {c['journal']}"
            if 'doi' in c and c['doi']:
                citation_text += f" (doi:{c['doi']})"
            elif 'url' in c and c['url']:
                citation_text += f" ({c['url']})"
            cards.append(f"""
            <div class="citation-card" style="background:{col};border-left-color:{col};">
                <div class="title">{c['title']}</div>
                <div class="authors">{c['authors']}</div>
                <div class="journal">{c['journal']}</div>
                <div class="doi"><a href="{link}" target="_blank">{label}</a></div>
                <button class="copy-btn" data-citation="{citation_text}">📋 Copy citation</button>
            </div>
            """)
        return f"""
        <div class="section-header citation-header"><h2><i class="fas fa-book"></i> Citations & References</h2><button class="print-section-btn" onclick="printSection('citation-tab')"><i class="fas fa-print"></i> Print</button></div>
        <div class="alert-box alert-success"><i class="fas fa-quote-right"></i><div>If you use Kleboscope, please cite the following tools and databases.</div></div>
        <div class="citation-grid">{"".join(cards)}</div>
        """

    def _generate_funding_section(self, data: Dict) -> str:
        return """
        <div class="section-header funding-header"><h2><i class="fas fa-coffee"></i> Funding & Support</h2></div>
        <div class="alert-box alert-info" style="background:#e8f5e9;border-left-color:#2e7d32;"><i class="fas fa-heart" style="color:#2e7d32;"></i><div><strong>No grants, just passion.</strong> Kleboscope is built with caffeine and curiosity at the University of Ghana Medical School.</div></div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:15px;margin:15px 0;">
            <div style="background:#fff3e0;padding:14px;border-radius:10px;text-align:center;"><i class="fas fa-star" style="color:#ff6f00;font-size:2em;"></i><div><strong>⭐ Star us on GitHub</strong><br>It takes 2 seconds and makes our day!</div></div>
            <div style="background:#fce4ec;padding:14px;border-radius:10px;text-align:center;"><i class="fas fa-bug" style="color:#c62828;font-size:2em;"></i><div><strong>🐛 Report bugs</strong><br>Help us improve.</div></div>
            <div style="background:#e0f7fa;padding:14px;border-radius:10px;text-align:center;"><i class="fas fa-lightbulb" style="color:#00695c;font-size:2em;"></i><div><strong>💡 Suggest features</strong><br>We listen and implement.</div></div>
            <div style="background:#f3e5f5;padding:14px;border-radius:10px;text-align:center;"><i class="fas fa-share-alt" style="color:#4a148c;font-size:2em;"></i><div><strong>📢 Spread the word</strong><br>Share with your network.</div></div>
        </div>
        <div class="alert-box alert-warning" style="background:#fff8e1;border-left-color:#ffa000;"><i class="fas fa-hand-holding-heart"></i><div><strong>🤝 Contribute to the ESCAPE AMR Platform</strong><br>We also have pipelines for Acinetobacter, E. coli, Staphylococcus, and more. Join us at <a href="https://github.com/bbeckley-hub" target="_blank">github.com/bbeckley-hub</a>.</div></div>
        <div style="background:#f5f5f5;padding:12px;border-radius:8px;margin-top:10px;font-style:italic;text-align:center;">
            <i class="fas fa-laugh-beam"></i> “This project runs on 0% grant money, 100% volunteer tears – but we’re not bitter, we’re caffeinated!” ☕😉
        </div>
        """

    def _generate_export_section(self, data: Dict) -> str:
        return f"""
        <div class="section-header export-header"><h2><i class="fas fa-download"></i> Export Data</h2></div>
        <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:15px 0;">
            <div class="dashboard-card" onclick="exportTableToCSV('samples-table','sample_overview.csv')"><i class="fas fa-file-csv fa-2x"></i><div>Sample Overview</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('amr-table','amr_genes.csv')"><i class="fas fa-biohazard fa-2x"></i><div>AMR Genes</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('vir-table','virulence_genes.csv')"><i class="fas fa-virus fa-2x"></i><div>Virulence Genes</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('bac-table','bacmet_genes.csv')"><i class="fas fa-flask fa-2x"></i><div>Bacmet Genes</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('plasmids-table','plasmid_replicons.csv')"><i class="fas fa-dna fa-2x"></i><div>Plasmids</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('mutations-table','mutations.csv')"><i class="fas fa-dna fa-2x"></i><div>Mutations</div></div>
            <div class="dashboard-card" onclick="exportTableToCSV('qc-table','fasta_qc.csv')"><i class="fas fa-chart-line fa-2x"></i><div>FASTA QC</div></div>
            <div class="dashboard-card" onclick="window.location.href='kleboscope_ultimate_sample_centric_report.json'"><i class="fas fa-file-code fa-2x"></i><div>Complete JSON</div></div>
        </div>
        """


# -----------------------------------------------------------------------------
# MAIN REPORTER
# -----------------------------------------------------------------------------
class KleboscopeUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "KLEBOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = KleboHTMLParser()
        self.analyzer = KleboDataAnalyzer()
        self.generator = KleboHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "Kleboscope Ultimate Reporter (Hybrid)",
            "version": "1.0.0",
            "author": "Brown Beckley",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def find_html_files(self) -> Dict[str, List[Path]]:
        html_files = {
            'mlst': [], 'qc': [], 'kaptive': [], 'amrfinder': [],
            'abricate': defaultdict(list)
        }
        all_html = list(self.input_dir.glob("**/*.html"))
        if not all_html:
            print("  ⚠️ No HTML files found.")
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
                for db in self.parser.abricate_databases:
                    if db in name:
                        html_files['abricate'][db].append(f)
                        break
        return html_files

    def integrate_all_data(self, html_files: Dict) -> Dict:
        integrated = {'metadata': self.metadata, 'samples': {}}
        all_samples = set()

        # 1. MLST: TSV first, then HTML
        mlst_data = self.parser.load_mlst_from_tsv(self.input_dir)
        if not mlst_data and html_files.get('mlst'):
            mlst_data = self.parser.parse_mlst_report(html_files['mlst'][0])
        if mlst_data:
            for s, d in mlst_data.items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['mlst'] = d

        # 2. QC from HTML (TSV fallback not needed, but keep)
        qc_data = {}
        if html_files.get('qc'):
            qc_data = self.parser.parse_qc_report(html_files['qc'][0])
        # Optionally try TSV if HTML missing
        if not qc_data:
            qc_tsv = self.input_dir / 'klebo_fasta_qc_summary.tsv'
            if qc_tsv.exists():
                # parse TSV if needed
                pass
        for s, d in qc_data.items():
            all_samples.add(s)
            integrated['samples'].setdefault(s, {})['qc'] = d

        # 3. Kaptive: HTML first, then TSV
        kaptive_data = {}
        if html_files.get('kaptive'):
            kaptive_data = self.parser.parse_kaptive_report(html_files['kaptive'][0])
        if not kaptive_data:
            kaptive_data = self.parser.load_kaptive_from_tsv(self.input_dir)
        if kaptive_data:
            for s, d in kaptive_data.items():
                all_samples.add(s)
                integrated['samples'].setdefault(s, {})['kaptive'] = d

        # 4. AMRfinder details from TSV
        amr_details, _ = self.parser.load_amrfinder_details_from_tsv(self.input_dir)
        if not amr_details and html_files.get('amrfinder'):
            amr_by_sample, _ = self.parser.parse_amrfinder_report(html_files['amrfinder'][0], len(all_samples))
            amr_details = {}
            for sample, genes in amr_by_sample.items():
                amr_details[sample] = [{'gene': g} for g in genes]
        integrated['amrfinder_details'] = amr_details
        for s, genes in amr_details.items():
            all_samples.add(s)
            integrated['samples'].setdefault(s, {})['amr_genes'] = [g.get('gene', '') for g in genes]

        # 5. ABRicate details from TSV
        abricate_details, _ = self.parser.load_abricate_details_from_tsv(self.input_dir)
        if not abricate_details and html_files.get('abricate'):
            for db, files in html_files['abricate'].items():
                if files:
                    genes_by_sample, _ = self.parser.parse_abricate_database_report(files[0], len(all_samples))
                    for sample, genes in genes_by_sample.items():
                        abricate_details[sample][db] = [{'gene': g} for g in genes]
        integrated['abricate_details'] = abricate_details
        for s, dbs in abricate_details.items():
            all_samples.add(s)
            integrated['samples'].setdefault(s, {})['abricate_genes'] = {
                db: [g.get('gene', '') for g in genes] for db, genes in dbs.items()
            }

        # 6. Mutation details from TSV
        mutation_details = self.parser.load_mutation_details_from_tsv(self.input_dir)
        if not mutation_details:
            mutation_html = self.input_dir / "mutation_summary.html"
            if mutation_html.exists():
                mutation_data = self.parser.parse_mutation_summary_html(mutation_html)
                if mutation_data and 'mutations' in mutation_data:
                    for mut in mutation_data['mutations']:
                        for genome in mut.get('genomes', []):
                            mutation_details.setdefault(genome, []).append({
                                'gene': mut['gene'],
                                'mutation': mut['mutation'],
                                'class': mut['class'],
                                'subclass': mut['subclass'],
                                'contig': '',
                                'start': '',
                                'stop': '',
                                'strand': '',
                                'coverage': '',
                                'identity': '',
                                'accession': ''
                            })
        integrated['mutation_details'] = mutation_details
        for s in mutation_details.keys():
            all_samples.add(s)

        # Ensure all samples have an entry
        for s in all_samples:
            integrated['samples'].setdefault(s, {})

      
        total_samples = len(all_samples)

        # Build gene_freqs from per‑sample details (for gene‑centric tabs)
        gene_freqs = {}
        if amr_details:
            gene_freqs['amrfinder'] = {}
            for sample, genes in amr_details.items():
                for g in genes:
                    gene = g.get('gene')
                    if gene:
                        if gene not in gene_freqs['amrfinder']:
                            gene_freqs['amrfinder'][gene] = {'count': 0, 'genomes': set(), 'database': 'amrfinder', 'risk_level': 'Standard'}
                        gene_freqs['amrfinder'][gene]['count'] += 1
                        gene_freqs['amrfinder'][gene]['genomes'].add(sample)
            for gene, data in gene_freqs['amrfinder'].items():
                data['genomes'] = list(data['genomes'])

        if abricate_details:
            db_gene_samples = defaultdict(lambda: defaultdict(set))
            for sample, dbs in abricate_details.items():
                for db, genes in dbs.items():
                    for g in genes:
                        gene = g.get('gene')
                        if gene:
                            db_gene_samples[db][gene].add(sample)
            for db, gene_samples in db_gene_samples.items():
                gene_freqs[db] = {}
                for gene, samples in gene_samples.items():
                    gene_freqs[db][gene] = {
                        'count': len(samples),
                        'genomes': list(samples),
                        'database': db,
                        'risk_level': 'Standard'
                    }

        integrated['gene_frequencies'] = gene_freqs
        integrated['gene_centric'] = self.analyzer.create_gene_centric_tables(gene_freqs, total_samples)
        integrated['patterns'] = self.analyzer.create_cross_genome_patterns(
            integrated['samples'], integrated['gene_centric']
        )

        if integrated['samples']:
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
        out = self.output_dir / "kleboscope_ultimate_sample_centric_report.json"
        with open(out, 'w', encoding='utf-8') as f:
            json.dump(integrated, f, indent=2, default=str)
        return out

    def generate_csv_reports(self, integrated: Dict):
        samples = integrated['samples']
        rows = []
        for sample, d in samples.items():
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

        # AMR, Virulence, Bacmet, Plasmids, Mutations CSV exports (from gene‑centric aggregated data)
        gene_centric = integrated.get('gene_centric', {})
        all_genes = gene_centric.get('all_genes', [])
        if all_genes:
            amr_rows = [g for g in all_genes if 'plasmidfinder' not in g['databases'] and 'bacmet2' not in g['databases'] and not any(db in ['vfdb','ecoli_vf'] for db in g['databases'])]
            if amr_rows:
                pd.DataFrame(amr_rows).to_csv(self.output_dir / "amr_genes.csv", index=False)
            vir_rows = [g for g in all_genes if g['category'] == 'High-Risk Virulence' or any(db in ['vfdb','ecoli_vf'] for db in g['databases'])]
            if vir_rows:
                pd.DataFrame(vir_rows).to_csv(self.output_dir / "virulence_genes.csv", index=False)
            bac_rows = [g for g in all_genes if 'bacmet2' in g['databases']]
            if bac_rows:
                pd.DataFrame(bac_rows).to_csv(self.output_dir / "bacmet_genes.csv", index=False)
            plasmid_rows = [g for g in all_genes if 'plasmidfinder' in g['databases']]
            if plasmid_rows:
                pd.DataFrame(plasmid_rows).to_csv(self.output_dir / "plasmid_replicons.csv", index=False)

        mutation_details = integrated.get('mutation_details', {})
        if mutation_details:
            mut_rows = []
            for sample, muts in mutation_details.items():
                for m in muts:
                    mut_rows.append({
                        'Sample': sample,
                        'Gene': m.get('gene', ''),
                        'Mutation': m.get('mutation', ''),
                        'Class': m.get('class', ''),
                        'Subclass': m.get('subclass', ''),
                        'Contig': m.get('contig', ''),
                        'Start': m.get('start', ''),
                        'Stop': m.get('stop', ''),
                        'Strand': m.get('strand', ''),
                        'Coverage': m.get('coverage', ''),
                        'Identity': m.get('identity', ''),
                        'Accession': m.get('accession', '')
                    })
            if mut_rows:
                pd.DataFrame(mut_rows).to_csv(self.output_dir / "mutations.csv", index=False)

    def run(self):
        print("=" * 80)
        print("🧬 Kleboscope Ultimate Reporter v1.0.0 (Hybrid)")
        print("   Gene‑centric for typing / Sample‑centric for AMR, Virulence, Bacmet, Plasmids, Mutations")
        print("=" * 80)
        html_files = self.find_html_files()
        if not any(html_files.values()):
            print("❌ No HTML files found.")
            return False
        integrated = self.integrate_all_data(html_files)
        if not integrated['samples']:
            print("❌ No data integrated.")
            return False
        self.generate_json_report(integrated)
        self.generate_csv_reports(integrated)
        self.generator.generate_main_report(integrated, self.output_dir)
        print("✅ All reports generated successfully.")
        print(f"📂 Output directory: {self.output_dir}")
        print("📄 Open kleboscope_ultimate_sample_centric_report.html in your browser.")
        return True


def main():
    parser = argparse.ArgumentParser(description='Kleboscope Ultimate Reporter – Hybrid Gene‑Centric & Sample‑Centric v1.0.0')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing Kleboscope TSV and HTML reports')
    args = parser.parse_args()
    reporter = KleboscopeUltimateReporter(Path(args.input_dir))
    reporter.run()

if __name__ == "__main__":
    main()