#!/usr/bin/env python3
"""
Kleboscope Ultimate Reporter – Sample-Centric
==============================================
Lazy-loaded isolate boxes for AMR, Virulence, Bacmet, Plasmids, Mutations, MGE.
Single typing tab (Sample Overview). No QC, no aggregation, no compare, no AI/CTA/funding/patterns tabs.

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.2.0
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
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# PARSER
# =============================================================================
class KleboHTMLParser:
    """Load typing, per-sample gene, mutation, and MGE data from TSVs."""

    def __init__(self):
        self.abricate_tsv_files = {
            'card':         'klebo_card_abricate_summary.tsv',
            'resfinder':    'klebo_resfinder_abricate_summary.tsv',
            'argannot':     'klebo_argannot_abricate_summary.tsv',
            'vfdb':         'klebo_vfdb_abricate_summary.tsv',
            'plasmidfinder':'klebo_plasmidfinder_abricate_summary.tsv',
            'megares':      'klebo_megares_abricate_summary.tsv',
            'ncbi':         'klebo_ncbi_abricate_summary.tsv',
            'ecoh':         'klebo_ecoh_abricate_summary.tsv',
            'ecoli_vf':     'klebo_ecoli_vf_abricate_summary.tsv',
            'bacmet2':      'klebo_bacmet2_abricate_summary.tsv',
        }

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id).strip()
        for ext in ['.fna', '.fasta', '.fa', '.gb', '.gbk', '.gbff', '.txt', '.tsv', '.csv']:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        if '/' in sample or '\\' in sample:
            sample = Path(sample).name
        return sample

    @staticmethod
    def _clean(val):
        return '' if pd.isna(val) else val

    def load_mlst_from_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        p = input_dir / 'kleboscope_mlst_summary.tsv'
        if not p.exists():
            return {}
        df = pd.read_csv(p, sep='\t', dtype=str).fillna('ND')
        out = {}
        for _, row in df.iterrows():
            s = self.normalize_sample_id(row.get('sample', ''))
            if not s:
                continue
            st = str(row.get('st', 'ND')).strip()
            if st.startswith('ST'):
                st = st[2:]
            out[s] = {'ST': st}
        return out

    def load_kaptive_from_tsv(self, input_dir: Path) -> Dict[str, Dict]:
        p = input_dir / 'kleboscope_kaptive_summary.tsv'
        if not p.exists():
            return {}
        df = pd.read_csv(p, sep='\t', dtype=str).fillna('ND')
        out = {}
        for _, row in df.iterrows():
            s = self.normalize_sample_id(row.get('sample', ''))
            if not s:
                continue
            out[s] = {
                'K_Locus': str(row.get('k_locus', 'ND')).strip(),
                'O_Locus': str(row.get('o_locus', 'ND')).strip(),
            }
        return out

    def load_amrfinder_details_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        p = input_dir / 'klebo_amrfinder_summary.tsv'
        if not p.exists():
            return {}
        df = pd.read_csv(p, sep='\t')
        details = defaultdict(list)
        for _, row in df.iterrows():
            s = self.normalize_sample_id(str(row.get('Genome', '')))
            if not s:
                continue
            d = {k: self._clean(v) for k, v in row.to_dict().items()}
            if 'Gene_Symbol' in d:
                d['gene'] = d.pop('Gene_Symbol')
            details[s].append(d)
        return dict(details)

    def load_abricate_details_from_tsv(self, input_dir: Path) -> Dict[str, Dict[str, List[Dict]]]:
        details = defaultdict(lambda: defaultdict(list))
        for db, fname in self.abricate_tsv_files.items():
            p = input_dir / fname
            if not p.exists():
                continue
            df = pd.read_csv(p, sep='\t')
            sample_col = 'genome' if 'genome' in df.columns else 'file'
            for _, row in df.iterrows():
                s = self.normalize_sample_id(str(row.get(sample_col, '')))
                if not s:
                    continue
                d = {k: self._clean(v) for k, v in row.to_dict().items()}
                details[s][db].append(d)
        return {s: dict(v) for s, v in details.items()}

    def load_mutations_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        p = input_dir / 'mutation_summary.tsv'
        if not p.exists():
            return {}
        df = pd.read_csv(p, sep='\t')
        details = defaultdict(list)
        for _, row in df.iterrows():
            s = self.normalize_sample_id(str(row.get('genome', '')))
            if not s:
                continue
            d = {k: self._clean(v) for k, v in row.to_dict().items()}
            if 'gene_symbol' in d:
                d['gene'] = d.pop('gene_symbol')
            if 'sequence_name' in d:
                d['mutation'] = d.pop('sequence_name')
            if 'contig_id' in d:
                d['contig'] = d.pop('contig_id')
            details[s].append(d)
        return dict(details)

    def load_mge_from_tsv(self, input_dir: Path) -> Dict[str, List[Dict]]:
        p = input_dir / 'kleboscope_mge_all_hits.tsv'
        if not p.exists():
            return {}
        df = pd.read_csv(p, sep='\t')
        details = defaultdict(list)
        for _, row in df.iterrows():
            s = self.normalize_sample_id(str(row.get('Sample', '')))
            if not s:
                continue
            d = {k: self._clean(v) for k, v in row.to_dict().items()}
            d.pop('Sample', None)
            details[s].append(d)
        return dict(details)


# =============================================================================
# HTML GENERATOR
# =============================================================================
class KleboHTMLGenerator:
    """Render the interactive report with lazy-loaded isolate boxes."""

    TAB_COLORS = {
        'sample_overview': '#2196F3',
        'amr':             '#F44336',
        'virulence':       '#E91E63',
        'bacmet':          '#FF5722',
        'plasmids':        '#673AB7',
        'mutations':       '#00BCD4',
        'mge':             '#16A085',
        'citation':        '#FFC107',
        'export':          '#9E9E9E',
    }

    def __init__(self):
        pass

    # ---------- helpers ----------
    def _section_header(self, title: str, icon: str, tab_id: str) -> str:
        return (
            f'<div class="section-header {tab_id}-header">'
            f'<h2><i class="fas fa-{icon}"></i> {title}</h2>'
            f'<button class="print-section-btn" onclick="printSection(\'{tab_id}-tab\')">'
            f'<i class="fas fa-print"></i> Print</button></div>'
        )

    def _credit_bar(self, color: str, icon: str, title: str, body: str) -> str:
        return f'''
        <div class="scientific-note" style="border-left-color:{color};">
            <div style="font-size:1.4em;">{icon}</div>
            <div>
                <strong style="color:{color};font-size:1.05em;">{title}</strong><br>
                <span style="font-size:0.92em;color:#4b5563;">{body}</span>
            </div>
        </div>'''

    def _alert(self, kind: str, icon: str, body: str) -> str:
        return f'<div class="alert-box alert-{kind}"><i class="fas fa-{icon} fa-2x"></i><div>{body}</div></div>'

    # ---------- sample-centric lazy boxes ----------
    def _generate_sample_centric_tab(
        self, tab_id: str, title: str, icon: str,
        samples_data: Dict, box_data: Dict, db_options: Dict[str, str],
        credit_html: str, count_label: str
    ) -> str:
        """Generate a tab with one lazy-loaded box per sample."""
        relevant = sorted(box_data.keys())

        if not relevant:
            return (
                self._section_header(title, icon, tab_id) + credit_html +
                self._alert('warning', 'exclamation-circle',
                            f'<h3>No {title} Data</h3><p>No samples with {title} hits.</p>')
            )

        box_json = json.dumps(box_data, default=str, ensure_ascii=False)

        options_html = '<option value="all">All Sources</option>'
        for key, label in db_options.items():
            options_html += f'<option value="{key}">{label}</option>'

        filter_html = f'''
        <div class="filter-controls">
            <input type="text" class="search-box" id="search-{tab_id}"
                   onkeyup="filterBoxes('{tab_id}')"
                   placeholder="🔍 Search sample..." style="max-width:320px;">
            <select id="dbFilter-{tab_id}" onchange="filterBoxes('{tab_id}')">{options_html}</select>
            <button class="action-btn btn-info" onclick="expandAllBoxes('{tab_id}')">
                <i class="fas fa-expand-alt"></i> Expand All (visible)</button>
            <button class="action-btn btn-light" onclick="collapseAllBoxes('{tab_id}')">
                <i class="fas fa-compress-alt"></i> Collapse All</button>
            <button class="action-btn btn-success" onclick="resetBoxFilters('{tab_id}')">
                <i class="fas fa-sync"></i> Clear</button>
            <span class="results-counter" id="counter-{tab_id}">{len(relevant)} shown</span>
        </div>
        <div id="box-container-{tab_id}">'''

        boxes_html = ''
        for sample in relevant:
            sd = samples_data.get(sample, {})
            st = sd.get('mlst', {}).get('ST', 'ND')
            k = sd.get('kaptive', {}).get('K_Locus', 'ND')
            o = sd.get('kaptive', {}).get('O_Locus', 'ND')
            st_display = f"ST{st}" if st != 'ND' else '—'

            entry = box_data.get(sample, {})
            n_items = self._count_items(tab_id, entry)

            boxes_html += f'''
            <div class="isolate-box" data-sample="{sample}" data-tab="{tab_id}">
                <div class="sample-header">
                    <h3><i class="fas fa-microbe"></i> {sample}</h3>
                    <span class="total-badge">{count_label}: {n_items}</span>
                    <div class="typing-info">
                        <span class="typing-badge badge-ST">ST: {st_display}</span>
                        <span class="typing-badge badge-K">K: {k}</span>
                        <span class="typing-badge badge-O">O: {o}</span>
                    </div>
                    <button class="toggle-btn" onclick="toggleBoxDetails('{tab_id}', this.closest('.isolate-box').getAttribute('data-sample'))">
                        <i class="fas fa-chevron-down"></i> Show Details
                    </button>
                </div>
                <div class="box-details" style="display:none;"></div>
            </div>'''

        script_block = f'''
        <script>
        window.KLEBOSCOPE_BOX_DATA = window.KLEBOSCOPE_BOX_DATA || {{}};
        window.KLEBOSCOPE_BOX_DATA["{tab_id}"] = {box_json};
        </script>'''

        return (
            self._section_header(title, icon, tab_id) + credit_html +
            filter_html + boxes_html + '</div>' + script_block
        )

    @staticmethod
    def _count_items(tab_id: str, entry: Dict) -> int:
        if tab_id in ('mutation', 'mutations'):
            return len(entry.get('mutations', []))
        if tab_id == 'mge':
            return len(entry.get('mge', []))
        total = 0
        if entry.get('amrfinder'):
            total += len(entry['amrfinder'])
        for db, hits in entry.get('abricate', {}).items():
            total += len(hits)
        return total

    # ---------- sample overview ----------
    def _generate_sample_overview_section(self, data: Dict) -> str:
        samples = data['samples']
        mge_details = data.get('mge_details', {})

        rows = ''
        for sample, sd in sorted(samples.items()):
            st = sd.get('mlst', {}).get('ST', 'ND')
            st_disp = f"ST{st}" if st != 'ND' else 'ND'
            k = sd.get('kaptive', {}).get('K_Locus', 'ND')
            o = sd.get('kaptive', {}).get('O_Locus', 'ND')
            mge_count = len(mge_details.get(sample, []))

            rows += (
                f'<tr>'
                f'<td><strong>{sample}</strong></td>'
                f'<td><span class="typing-badge badge-ST">{st_disp}</span></td>'
                f'<td><span class="typing-badge badge-K">{k}</span></td>'
                f'<td><span class="typing-badge badge-O">{o}</span></td>'
                f'<td>{mge_count}</td>'
                f'</tr>'
            )

        return f'''
        {self._section_header('Sample Overview', 'list', 'sample_overview')}
        {self._alert('info', 'info-circle',
            '<strong>One row per isolate</strong> with ST, K-locus, O-locus, and total MGE hits. '
            'Gene and mutation details live in the AMR / Virulence / Bacmet / Plasmids / '
            'Mutations / MGE tabs as lazy-loaded isolate boxes.')}
        <input type="text" class="search-box" id="search-samples"
               onkeyup="searchTable('samples-table','search-samples')"
               placeholder="🔍 Search sample, ST, K, or O...">
        <div class="action-buttons">
            <button class="action-btn btn-primary"
                    onclick="exportTableToCSV('samples-table','sample_overview.csv')">
                <i class="fas fa-download"></i> Export CSV</button>
            <button class="action-btn btn-light"
                    onclick="document.getElementById('search-samples').value='';searchTable('samples-table','search-samples')">
                <i class="fas fa-sync"></i> Clear</button>
        </div>
        <div class="master-scrollable-container">
            <table id="samples-table" class="data-table">
                <thead><tr>
                    <th data-sort="string">Sample</th>
                    <th data-sort="string">ST</th>
                    <th data-sort="string">K Locus</th>
                    <th data-sort="string">O Locus</th>
                    <th data-sort="number">MGE Hits</th>
                </tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </div>'''

    # ---------- AMR ----------
    def _generate_amr_section(self, data: Dict) -> str:
        amr_details = data.get('amrfinder_details', {})
        abricate = data.get('abricate_details', {})
        amr_dbs = ('amrfinder', 'card', 'resfinder', 'argannot', 'megares', 'ncbi')

        box_data = {}
        for sample in data['samples']:
            entry = {'amrfinder': [], 'abricate': {}}
            if sample in amr_details and amr_details[sample]:
                entry['amrfinder'] = amr_details[sample]
            if sample in abricate:
                for db in amr_dbs:
                    if db == 'amrfinder':
                        continue
                    hits = abricate[sample].get(db, [])
                    if hits:
                        entry['abricate'][db] = hits
            if entry['amrfinder'] or entry['abricate']:
                box_data[sample] = entry

        credit = self._credit_bar('#F44336', '💊',
            'ABRicate + AMRFinderPlus — Resistance Detection',
            'Powered by <strong>ABRicate</strong> (Seemann T.) with CARD, ResFinder, ARG-ANNOT, '
            'MEGARes and NCBI AMR databases, plus <strong>AMRFinderPlus</strong> (NCBI). '
            'Every hit preserves database provenance. '
            '<em>Gene presence ≠ phenotypic resistance — AST remains the gold standard.</em>')

        db_options = {
            'amrfinder': 'AMRfinderPlus',
            'card': 'CARD', 'resfinder': 'ResFinder', 'argannot': 'ARG-ANNOT',
            'megares': 'MEGARes', 'ncbi': 'NCBI AMR',
        }
        return self._generate_sample_centric_tab(
            'amr', 'AMR Genes', 'biohazard',
            data['samples'], box_data, db_options, credit, 'Total Gene Hits'
        )

    # ---------- Virulence ----------
    def _generate_virulence_section(self, data: Dict) -> str:
        abricate = data.get('abricate_details', {})
        box_data = {}
        for sample in data['samples']:
            entry = {'abricate': {}}
            if sample in abricate:
                for db in ('vfdb', 'ecoli_vf'):
                    hits = abricate[sample].get(db, [])
                    if hits:
                        entry['abricate'][db] = hits
            if entry['abricate']:
                box_data[sample] = entry

        credit = self._credit_bar('#E91E63', '🦠',
            'VFDB — Virulence Factor Detection',
            'Powered by <strong>ABRicate</strong> with the '
            '<a href="http://www.mgc.ac.cn/VFs/" target="_blank">VFDB database</a>. '
            'Cite: Chen L et al. <em>Nucleic Acids Res</em> 2012;40(D1):D641-5.')

        db_options = {'vfdb': 'VFDB', 'ecoli_vf': 'Ecoli VF'}
        return self._generate_sample_centric_tab(
            'virulence', 'Virulence Genes', 'virus',
            data['samples'], box_data, db_options, credit, 'Total Virulence Hits'
        )

    # ---------- Bacmet ----------
    def _generate_bacmet_section(self, data: Dict) -> str:
        abricate = data.get('abricate_details', {})
        box_data = {}
        for sample in data['samples']:
            entry = {'abricate': {}}
            if sample in abricate:
                hits = abricate[sample].get('bacmet2', [])
                if hits:
                    entry['abricate']['bacmet2'] = hits
            if entry['abricate']:
                box_data[sample] = entry

        credit = self._credit_bar('#FF5722', '🧪',
            'BacMet — Biocide & Heavy Metal Resistance',
            'Powered by <strong>ABRicate</strong> with the '
            '<a href="http://bacmet.biomedicine.gu.se/" target="_blank">BacMet database</a>. '
            'Cite: Pal C et al. <em>Nucleic Acids Res</em> 2014;42(D1):D737-43.')

        db_options = {'bacmet2': 'BacMet'}
        return self._generate_sample_centric_tab(
            'bacmet', 'Biocide & Metal Resistance', 'flask',
            data['samples'], box_data, db_options, credit, 'Total BACMET Hits'
        )

    # ---------- Plasmids ----------
    def _generate_plasmids_section(self, data: Dict) -> str:
        abricate = data.get('abricate_details', {})
        box_data = {}
        for sample in data['samples']:
            entry = {'abricate': {}}
            if sample in abricate:
                hits = abricate[sample].get('plasmidfinder', [])
                if hits:
                    entry['abricate']['plasmidfinder'] = hits
            if entry['abricate']:
                box_data[sample] = entry

        credit = self._credit_bar('#673AB7', '🔌',
            'PlasmidFinder — Replicon Typing',
            'Powered by <strong>ABRicate</strong> with the '
            '<a href="https://genepi.food.dtu.dk/PlasmidFinder/" target="_blank">PlasmidFinder database</a>. '
            'Cite: Carattoli A et al. <em>Antimicrob Agents Chemother</em> 2014;58(7):3895-903.')

        db_options = {'plasmidfinder': 'PlasmidFinder'}
        return self._generate_sample_centric_tab(
            'plasmids', 'Plasmid Replicons', 'dna',
            data['samples'], box_data, db_options, credit, 'Total Replicon Hits'
        )

    # ---------- Mutations ----------
    def _generate_mutations_section(self, data: Dict) -> str:
        mutation_details = data.get('mutation_details', {})
        box_data = {}
        for sample in data['samples']:
            hits = mutation_details.get(sample, [])
            if hits:
                box_data[sample] = {'mutations': hits}

        credit = self._credit_bar('#00BCD4', '🧬',
            'AMRFinderPlus — Point Mutations',
            'Powered by <strong>AMRFinderPlus</strong> (NCBI). '
            'Cite: Feldgarden M et al. <em>Sci Rep</em> 2021;11(1):12728.')

        db_options = {'mutations': 'All Mutations'}
        return self._generate_sample_centric_tab(
            'mutations', 'Point Mutations', 'dna',
            data['samples'], box_data, db_options, credit, 'Total Mutations'
        )

    # ---------- MGE ----------
    def _generate_mge_section(self, data: Dict) -> str:
        mge_details = data.get('mge_details', {})
        box_data = {}
        for sample in data['samples']:
            hits = mge_details.get(sample, [])
            if hits:
                box_data[sample] = {'mge': hits}

        credit = self._credit_bar('#16A085', '🧬',
            'mobileOG-db — Mobile Genetic Elements',
            'Powered by <strong>mobileOG-db</strong> (Brown CL, Mullet J et al.). '
            'Cite: Brown CL et al. <em>Appl Environ Microbiol</em> 2022;88(18):e00991-22.')

        db_options = {'mge': 'mobileOG Hits'}
        return self._generate_sample_centric_tab(
            'mge', 'Mobile Genetic Elements', 'mobile-alt',
            data['samples'], box_data, db_options, credit, 'Total MGE Hits'
        )

    # ---------- Citation ----------
    def _generate_citation_section(self, data: Dict) -> str:
        citations = [
            ('Kleboscope', 'Beckley B et al.', 'GitHub', 'https://github.com/bbeckley-hub/kleboscope'),
            ('MLST (PubMLST)', 'Jolley KA, Bray JE, Maiden MCJ', 'Wellcome Open Res 2018;3:124', 'https://doi.org/10.12688/wellcomeopenres.14826.1'),
            ('Kaptive', 'Wyres KL et al.', 'Microb Genom 2016;2(5):e000064', 'https://doi.org/10.1099/mgen.0.000064'),
            ('AMRFinderPlus', 'Feldgarden M et al.', 'Sci Rep 2021;11(1):12728', 'https://doi.org/10.1038/s41598-021-91456-0'),
            ('ABRicate', 'Seemann T', 'GitHub 2024', 'https://github.com/tseemann/abricate'),
            ('CARD', 'McArthur AG et al.', 'Antimicrob Agents Chemother 2013;57(7):3348-57', 'https://doi.org/10.1128/AAC.00419-13'),
            ('ResFinder', 'Florensa AF et al.', 'Microb Genom 2022;8(1):000748', 'https://doi.org/10.1099/mgen.0.000748'),
            ('VFDB', 'Chen L et al.', 'Nucleic Acids Res 2012;40(D1):D641-5', 'https://doi.org/10.1093/nar/gkr989'),
            ('PlasmidFinder', 'Carattoli A et al.', 'Antimicrob Agents Chemother 2014;58(7):3895-903', 'https://doi.org/10.1128/AAC.02412-14'),
            ('BacMet', 'Pal C et al.', 'Nucleic Acids Res 2014;42(D1):D737-43', 'https://doi.org/10.1093/nar/gkt1252'),
            ('MEGARes 2.0', 'Doster E et al.', 'Nucleic Acids Res 2020;48(D1):D561-D569', 'https://doi.org/10.1093/nar/gkz1010'),
            ('ARG-ANNOT', 'Gupta SK et al.', 'Antimicrob Agents Chemother 2014;58(1):212-20', 'https://doi.org/10.1128/AAC.01310-13'),
            ('mobileOG-db', 'Brown CL et al.', 'Appl Environ Microbiol 2022;88(18):e00991-22', 'https://doi.org/10.1128/aem.00991-22'),
            ('Biopython', 'Cock PJ et al.', 'Bioinformatics 2009;25(11):1422-3', 'https://doi.org/10.1093/bioinformatics/btp163'),
        ]
        items = ''
        for name, authors, journal, url in citations:
            plain = f"{name} – {authors}. {journal}."
            items += f'''
            <li class="citation-item">
                <div class="citation-body">
                    <div><strong style="color:#005a8c;">{name}</strong> – <span style="color:#333;">{authors}. <em>{journal}</em></span></div>
                    <div class="citation-actions">
                        <a class="citation-link" href="{url}" target="_blank">🔗 Open</a>
                        <button class="copy-btn" data-citation="{plain}">📋 Copy</button>
                    </div>
                </div>
            </li>'''

        return f'''
        {self._section_header('Citations', 'book', 'citation')}
        {self._alert('success', 'quote-right',
            'If you use Kleboscope, please cite the underlying tools and databases.')}
        <ul style="padding:0;margin:0;">{items}</ul>'''

    # ---------- Export ----------
    def _generate_export_section(self, data: Dict) -> str:
        exports = [
            ('samples-table', 'Sample Overview', 'sample_overview.csv', 'table'),
        ]
        cards = ''
        for tid, label, fname, icon in exports:
            cards += f'''
            <div class="dashboard-card" onclick="exportTableToCSV('{tid}', '{fname}')">
                <div class="card-icon"><i class="fas fa-{icon}"></i></div>
                <div class="card-label">{label}</div>
            </div>'''
        cards += '''
            <div class="dashboard-card" onclick="window.location.href='kleboscope_ultimate_sample_centric_report.json'">
                <div class="card-icon"><i class="fas fa-file-code"></i></div>
                <div class="card-label">Complete JSON</div>
            </div>'''
        return f'''
        {self._section_header('Export Data', 'download', 'export')}
        {self._alert('info', 'info-circle',
            'Download tables as CSV, or grab the full JSON for AI-assisted analysis. '
            'Gene and MGE tables are available in the JSON (lazy-loaded on screen).')}
        <div class="dashboard-grid" style="margin-top:16px;">{cards}</div>'''

    # ---------- Main assembly ----------
    def generate_main_report(self, data: Dict, output_dir: Path) -> str:
        print("\n🎨 Generating Kleboscope Ultimate HTML report...")
        html = self._build_html(data)
        out = output_dir / "kleboscope_ultimate_sample_centric_report.html"
        with open(out, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {out}")
        return str(out)

    def _build_html(self, data: Dict) -> str:
        metadata = data.get('metadata', {})
        samples = data.get('samples', {})
        total = len(samples)

        # ---- Precompute dashboard numbers (avoids f-string comprehension bugs) ----
        total_sts = len({
            s.get('mlst', {}).get('ST')
            for s in samples.values()
            if s.get('mlst', {}).get('ST') not in (None, '', 'ND')
        })
        total_amr = sum(len(v) for v in data.get('amrfinder_details', {}).values())
        total_vir = 0
        for sample_dbs in data.get('abricate_details', {}).values():
            if isinstance(sample_dbs, dict):
                for db in ('vfdb', 'ecoli_vf'):
                    total_vir += len(sample_dbs.get(db, []) or [])
        total_mge = sum(len(v) for v in data.get('mge_details', {}).values())

        tabs = [
            ('sample_overview', 'Sample Overview', 'list'),
            ('amr',             'AMR',             'biohazard'),
            ('virulence',       'Virulence',       'virus'),
            ('bacmet',          'Bacmet',          'flask'),
            ('plasmids',        'Plasmids',        'dna'),
            ('mutations',       'Mutations',       'dna'),
            ('mge',             'MGE',             'mobile-alt'),
            ('citation',        'Citation',        'book'),
            ('export',          'Export',          'download'),
        ]

        section_builders = {
            'sample_overview': self._generate_sample_overview_section,
            'amr':             self._generate_amr_section,
            'virulence':       self._generate_virulence_section,
            'bacmet':          self._generate_bacmet_section,
            'plasmids':        self._generate_plasmids_section,
            'mutations':       self._generate_mutations_section,
            'mge':             self._generate_mge_section,
            'citation':        self._generate_citation_section,
            'export':          self._generate_export_section,
        }

        nav_html = ''
        for i, (tid, label, icon) in enumerate(tabs):
            active = ' active' if i == 0 else ''
            nav_html += (
                f'<button class="tab-button {tid}{active}" '
                f'onclick="switchTab(\'{tid}\')">'
                f'<i class="fas fa-{icon}"></i> {label}</button>'
            )

        tabs_html = ''
        for i, (tid, label, icon) in enumerate(tabs):
            active = ' active' if i == 0 else ''
            content = section_builders[tid](data)
            tabs_html += f'<div id="{tid}-tab" class="tab-content{active}">{content}</div>'

        dash = f'''
        <div class="dashboard-grid">
            <div class="dashboard-card card-sample_overview" onclick="switchTab('sample_overview')">
                <div class="card-icon"><i class="fas fa-vial"></i></div>
                <div class="card-number">{total}</div>
                <div class="card-label">Samples</div></div>
            <div class="dashboard-card card-sample_overview" onclick="switchTab('sample_overview')">
                <div class="card-icon"><i class="fas fa-list"></i></div>
                <div class="card-number">{total_sts}</div>
                <div class="card-label">Unique STs</div></div>
            <div class="dashboard-card card-amr" onclick="switchTab('amr')">
                <div class="card-icon"><i class="fas fa-biohazard"></i></div>
                <div class="card-number">{total_amr}</div>
                <div class="card-label">AMR Hits</div></div>
            <div class="dashboard-card card-virulence" onclick="switchTab('virulence')">
                <div class="card-icon"><i class="fas fa-virus"></i></div>
                <div class="card-number">{total_vir}</div>
                <div class="card-label">Virulence Hits</div></div>
            <div class="dashboard-card card-mge" onclick="switchTab('mge')">
                <div class="card-icon"><i class="fas fa-mobile-alt"></i></div>
                <div class="card-number">{total_mge}</div>
                <div class="card-label">MGE Hits</div></div>
        </div>'''

        return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Kleboscope Ultimate Report</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
{self._get_css()}
{self._get_js()}
</head>
<body>
<div class="container">
    <div class="main-header">
        <h1><i class="fas fa-bacterium"></i> Kleboscope Ultimate Report</h1>
        <p>Lazy-Loaded Isolate Boxes with Typing Overview</p>
        <div class="metadata-bar">
            <div class="metadata-item"><i class="fas fa-calendar-alt"></i> {metadata.get('analysis_date', 'Unknown')}</div>
            <div class="metadata-item"><i class="fas fa-database"></i> {total} Samples</div>
            <div class="metadata-item"><i class="fas fa-university"></i> University of Ghana Medical School</div>
        </div>
    </div>
    {dash}
    <div class="tab-navigation">{nav_html}</div>
    {tabs_html}
    <div class="footer">
        <h3>Kleboscope Ultimate Reporterv2.0.0</h3>
        <p>University of Ghana Medical School | Brown Beckley &lt;brownbeckley94@gmail.com&gt;</p>
        <p>Generated on {metadata.get('analysis_date', 'Unknown')}</p>
        <p>If you find this useful, please <a href="https://github.com/bbeckley-hub" target="_blank">⭐ star us on GitHub</a>.</p>
    </div>
</div>
</body>
</html>'''

    # ---------- CSS ----------
    def _get_css(self) -> str:
        colors_css = '\n'.join(f'    --{t}-color: {c};' for t, c in self.TAB_COLORS.items())
        active_css = '\n'.join(
            f'.tab-button.active.{t} {{ background: var(--{t}-color); color: white; border-color: var(--{t}-color); }}'
            for t in self.TAB_COLORS
        )
        header_css = '\n'.join(
            f'.{t}-header {{ border-color: var(--{t}-color); }}'
            for t in self.TAB_COLORS
        )
        hover_css = '\n'.join(
            f'.tab-button.{t}:hover:not(.active) {{ color: var(--{t}-color); }}'
            for t in self.TAB_COLORS
        )
        dash_css = '\n'.join(
            f'.dashboard-card.card-{t} {{ border-left-color: var(--{t}-color); }}'
            for t in self.TAB_COLORS
        )
        return f"""
<style>
:root {{
    --primary: #0077be;
    --primary-dark: #005a8c;
    --success: #10b981;
    --warning: #f59e0b;
    --danger: #dc2626;
    --bg-light: #f8f9fa;
{colors_css}
}}
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    color: white; padding:20px; min-height:100vh;
}}
.container {{ max-width:1800px; margin:0 auto; }}

.main-header {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 100%);
    color: white; padding:30px; border-radius:15px;
    box-shadow: 0 10px 30px rgba(0,0,0,0.3); margin-bottom:24px;
    text-align:center; border:2px solid rgba(0,119,190,0.3);
}}
.main-header h1 {{ font-size:2.4em; margin-bottom:8px; }}
.main-header p {{ opacity:0.9; margin-bottom:12px; }}
.metadata-bar {{
    background: rgba(255,255,255,0.1); padding:12px; border-radius:10px;
    margin-top:16px; display:flex; justify-content:center;
    flex-wrap:wrap; gap:22px; backdrop-filter:blur(10px);
}}
.metadata-item {{ display:flex; align-items:center; gap:8px; font-size:0.92em; }}

.dashboard-grid {{
    display: grid; grid-template-columns:repeat(auto-fit, minmax(180px, 1fr));
    gap:15px; margin-bottom:24px;
}}
.dashboard-card {{
    background: white; color:#1f2937; padding:20px;
    border-radius:12px; text-align:center; cursor:pointer;
    box-shadow: 0 4px 15px rgba(0,0,0,0.15);
    transition: transform 0.2s, box-shadow 0.2s;
    border-left: 5px solid #0077be;
}}
.dashboard-card:hover {{ transform:translateY(-4px); box-shadow:0 8px 25px rgba(0,0,0,0.25); }}
.dashboard-card .card-icon {{ font-size:1.6em; opacity:0.8; margin-bottom:6px; }}
.dashboard-card .card-number {{ font-size:2.2em; font-weight:bold; margin:4px 0; color:#1f2937; }}
.dashboard-card .card-label {{ font-size:0.82em; text-transform:uppercase; letter-spacing:0.5px; color:#6b7280; }}
{dash_css}

.tab-navigation {{
    display:flex; gap:5px; flex-wrap:wrap;
    background: rgba(255,255,255,0.95); padding:12px; border-radius:12px;
    box-shadow: 0 5px 20px rgba(0,0,0,0.2); margin-bottom:20px;
    position:sticky; top:10px; z-index:100;
}}
.tab-button {{
    padding:10px 16px; background:#e5e7eb; border:2px solid #e5e7eb;
    border-radius:8px; cursor:pointer; font-weight:600;
    color:#4b5563; font-size:0.85em;
    display:flex; align-items:center; gap:6px;
    transition: background 0.2s, color 0.2s, border-color 0.2s;
}}
{hover_css}
{active_css}

.tab-content {{
    display:none; background:rgba(255,255,255,0.97); color:#1f2937;
    padding:28px; border-radius:15px;
    box-shadow: 0 10px 30px rgba(0,0,0,0.2); margin-bottom:24px;
}}
.tab-content.active {{ display:block; animation:fadeIn 0.3s ease; }}
@keyframes fadeIn {{ from {{ opacity:0; transform:translateY(10px); }} to {{ opacity:1; transform:translateY(0); }} }}

.section-header {{
    display:flex; justify-content:space-between; align-items:center;
    border-bottom:3px solid #0077be; padding-bottom:12px; margin-bottom:22px;
}}
.section-header h2 {{ color:#1f2937; font-size:1.6em; display:flex; align-items:center; gap:10px; }}
{header_css}

.data-table {{
    width:100%; border-collapse:collapse; font-size:0.9em;
    background:white; border-radius:8px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.05); margin:12px 0;
}}
.data-table th {{
    background:#005a8c; color:white; padding:12px 10px;
    text-align:left; font-weight:600; position:sticky; top:0;
    cursor:pointer; white-space:nowrap;
}}
.data-table th:hover {{ background:#003f63; }}
.data-table td {{ padding:10px; border-bottom:1px solid #e5e7eb; vertical-align:top; }}
.data-table tr:hover {{ background:#f3f4f6; }}
.master-scrollable-container {{
    width:100%; overflow-x:auto; border:1px solid #e5e7eb;
    border-radius:8px; margin:14px 0;
}}

.alert-box {{
    padding:18px 22px; border-radius:10px; margin:16px 0;
    display:flex; align-items:flex-start; gap:16px;
    border-left:5px solid; font-size:0.94em;
}}
.alert-info    {{ background:#e0f2fe; color:#0c4a6e; border-color:#0077be; }}
.alert-success {{ background:#d1fae5; color:#065f46; border-color:#10b981; }}
.alert-warning {{ background:#fef3c7; color:#92400e; border-color:#f59e0b; }}
.alert-danger  {{ background:#fee2e2; color:#7f1d1d; border-color:#dc2626; }}

.scientific-note {{
    background:linear-gradient(135deg,#f0f9ff 0%,#e0f2fe 100%);
    border-left:5px solid #0077be; padding:15px 20px;
    border-radius:8px; margin-bottom:18px;
    display:flex; align-items:flex-start; gap:14px;
    color:#1f2937; font-size:0.92em;
}}

.search-box {{
    width:100%; padding:10px 14px; margin-bottom:10px;
    border:2px solid #e5e7eb; border-radius:8px; font-size:0.95em;
}}
.search-box:focus {{
    outline:none; border-color:#0077be;
    box-shadow: 0 0 0 3px rgba(0,119,190,0.15);
}}
.action-buttons {{ display:flex; gap:8px; margin:12px 0; flex-wrap:wrap; }}
.action-btn {{
    padding:8px 16px; border:none; border-radius:8px;
    cursor:pointer; font-weight:600; font-size:0.85em;
    display:flex; align-items:center; gap:6px;
}}
.action-btn:hover {{ transform:translateY(-1px); }}
.btn-primary   {{ background:#0077be; color:white; }}
.btn-success   {{ background:#10b981; color:white; }}
.btn-info      {{ background:#3b82f6; color:white; }}
.btn-light     {{ background:#f8f9fa; color:#1f2937; border:1px solid #d1d5db; }}

.print-section-btn {{
    background:#0077be; color:white; border:none;
    border-radius:6px; padding:7px 14px; cursor:pointer;
    font-size:0.85em; display:flex; align-items:center; gap:5px;
}}
.print-section-btn:hover {{ background:#005a8c; }}

.typing-badge {{
    display:inline-block; padding:3px 10px; border-radius:12px;
    font-size:0.8em; font-weight:600; background:#e0e0e0;
    color:#333; border:1px solid #ccc; white-space:nowrap;
}}
.typing-badge.badge-ST  {{ background:#FF9800; color:white; border-color:#FF9800; }}
.typing-badge.badge-K   {{ background:#9C27B0; color:white; border-color:#9C27B0; }}
.typing-badge.badge-O   {{ background:#009688; color:white; border-color:#009688; }}

.isolate-box {{
    border:1px solid #ddd; border-radius:12px;
    margin-bottom:20px; padding:20px;
    background:#fafafa; box-shadow:0 2px 8px rgba(0,0,0,0.06);
    transition:box-shadow 0.2s;
}}
.isolate-box:hover {{ box-shadow:0 4px 14px rgba(0,0,0,0.10); }}
.isolate-box .sample-header {{
    display:flex; align-items:center; gap:15px;
    flex-wrap:wrap;
}}
.isolate-box .sample-header h3 {{ font-size:1.4em; margin:0; }}
.isolate-box .sample-header .total-badge {{
    background:#0077be; color:white; padding:4px 16px;
    border-radius:20px; font-weight:bold; font-size:0.9em;
}}
.isolate-box .sample-header .typing-info {{
    display:flex; flex-wrap:wrap; gap:8px;
}}
.isolate-box .sample-header .toggle-btn {{
    background:#0077be; color:white; border:none;
    padding:8px 16px; border-radius:8px; cursor:pointer;
    font-weight:600; font-size:0.85em;
    display:inline-flex; align-items:center; gap:6px;
    margin-left:auto; transition:background 0.2s;
}}
.isolate-box .sample-header .toggle-btn:hover {{ background:#005a8c; }}
.isolate-box .box-details {{ margin-top:15px; }}

.database-table-wrapper {{
    margin:15px 0; overflow-x:auto;
    border:1px solid #e0e0e0; border-radius:8px;
}}
.database-table-wrapper table {{
    width:100%; border-collapse:collapse;
    font-size:0.85em; min-width:800px;
}}
.database-table-wrapper table th {{
    background:#2c3e50; color:white;
    padding:8px 12px; text-align:left; white-space:nowrap;
}}
.database-table-wrapper table td {{
    padding:8px 12px; border-bottom:1px solid #e0e0e0;
    white-space:nowrap;
}}
.database-table-wrapper table tr:hover {{ background:#f1f1f1; }}
.db-title {{
    font-weight:bold; color:#0077be;
    margin:10px 0 5px 0; font-size:1.05em;
    border-left:4px solid #0077be; padding-left:10px;
}}

.filter-controls {{
    display:flex; flex-wrap:wrap; gap:10px;
    align-items:center; background:#f8f9fa;
    padding:15px; border-radius:8px; margin-bottom:20px;
}}
.filter-controls select {{
    padding:10px; border-radius:8px;
    border:2px solid #ddd; background:white; min-width:150px;
}}
.results-counter {{
    font-size:0.9em; color:#555; font-weight:600;
    padding:8px 12px; background:#f0f0f0;
    border-radius:6px; white-space:nowrap;
}}

.citation-item {{
    background:white; color:#1f2937;
    border:1px solid #e5e7eb; border-left:4px solid #0077be;
    border-radius:8px; padding:14px 18px; margin-bottom:12px;
    list-style:none; box-shadow:0 1px 3px rgba(0,0,0,0.05);
}}
.citation-item:hover {{ box-shadow:0 4px 12px rgba(0,0,0,0.08); }}
.citation-body {{ display:flex; flex-direction:column; gap:8px; font-size:0.92em; }}
.citation-actions {{ display:flex; gap:8px; flex-wrap:wrap; }}
.citation-link {{
    display:inline-flex; align-items:center; gap:4px;
    padding:4px 14px; border-radius:14px; font-size:0.8em;
    font-weight:600; color:white; text-decoration:none;
    background:#0077be;
}}
.citation-link:hover {{ background:#005a8c; }}
.copy-btn {{
    background:#6b7280; color:white; border:none;
    padding:4px 14px; border-radius:14px; font-size:0.8em;
    font-weight:600; cursor:pointer;
}}
.copy-btn:hover {{ background:#4b5563; }}

.footer {{
    text-align:center; padding:30px; color:white;
    margin-top:40px; border-radius:15px;
    background:linear-gradient(135deg,#2c3e50 0%,#34495e 100%);
}}
.footer a {{ color:#fbbf24; text-decoration:none; }}
.footer a:hover {{ text-decoration:underline; }}

@media print {{
    body * {{ visibility:hidden; }}
    .tab-content.active, .tab-content.active * {{ visibility:visible; }}
    .tab-content.active {{
        position:absolute; left:0; top:0; width:100%;
        padding:20px; box-shadow:none; border-radius:0;
    }}
    .print-section-btn, .tab-navigation, .dashboard-grid,
    .search-box, .action-buttons, .filter-controls, .toggle-btn {{ display:none !important; }}
}}
</style>"""

    # ---------- JS ----------
    def _get_js(self) -> str:
        return """
<script>
window.KLEBOSCOPE_BOX_DATA = window.KLEBOSCOPE_BOX_DATA || {};

function escapeHtml(s) {
    return String(s)
        .replace(/&/g,'&amp;').replace(/</g,'&lt;')
        .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

function switchTab(name) {
    document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
    const content = document.getElementById(name + '-tab');
    const button = document.querySelector('.tab-button.' + name);
    if (content) content.classList.add('active');
    if (button) button.classList.add('active');
    window.location.hash = name;
}

function searchTable(tableId, searchId) {
    const filter = (document.getElementById(searchId)?.value || '').toUpperCase();
    const table = document.getElementById(tableId);
    if (!table || !table.tBodies[0]) return;
    const rows = table.tBodies[0].rows;
    for (let i = 0; i < rows.length; i++) {
        let found = false;
        for (const cell of rows[i].getElementsByTagName('td')) {
            if ((cell.textContent || '').toUpperCase().indexOf(filter) > -1) { found = true; break; }
        }
        rows[i].style.display = found ? '' : 'none';
    }
}

function sortTable(tableId, colIdx, type) {
    const table = document.getElementById(tableId);
    if (!table || !table.tBodies[0]) return;
    const tbody = table.tBodies[0];
    const rows = Array.from(tbody.rows);
    const asc = table.getAttribute('data-sort-dir') !== 'asc';
    rows.sort((a, b) => {
        let av = a.cells[colIdx].innerText.trim();
        let bv = b.cells[colIdx].innerText.trim();
        if (type === 'number') {
            av = parseFloat(av.replace(/,/g, '')) || 0;
            bv = parseFloat(bv.replace(/,/g, '')) || 0;
            return asc ? av - bv : bv - av;
        }
        return asc ? av.localeCompare(bv) : bv.localeCompare(av);
    });
    rows.forEach(r => tbody.appendChild(r));
    table.setAttribute('data-sort-dir', asc ? 'asc' : 'desc');
}

function exportTableToCSV(tableId, filename) {
    const table = document.getElementById(tableId);
    if (!table) return;
    const rows = table.querySelectorAll('tr');
    const csv = [];
    rows.forEach(row => {
        const cols = row.querySelectorAll('td, th');
        csv.push(Array.from(cols).map(c => '"' + (c.innerText || '').replace(/"/g, '""') + '"').join(','));
    });
    const blob = new Blob([csv.join('\\n')], { type: 'text/csv' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
}

function printSection(sectionId) {
    const content = document.getElementById(sectionId);
    if (!content) return;
    content.querySelectorAll('.isolate-box').forEach(box => {
        const tabId = box.getAttribute('data-tab');
        const sample = box.getAttribute('data-sample');
        const details = box.querySelector('.box-details');
        if (details && !details.innerHTML.trim()) {
            details.innerHTML = buildBoxContent(tabId, sample);
        }
        if (details) details.style.display = 'block';
    });
    const w = window.open('', '_blank');
    const style = document.querySelector('style');
    w.document.write('<html><head><title>Print</title>');
    if (style) w.document.write('<style>' + style.textContent + '</style>');
    w.document.write('</head><body>' + content.innerHTML + '</body></html>');
    w.document.close();
    w.print();
}

/* ---------------- LAZY BOX RENDERING ---------------- */

function toggleBoxDetails(tabId, sample) {
    const tab = document.getElementById(tabId + '-tab');
    if (!tab) return;
    let box = null;
    for (const candidate of tab.querySelectorAll('.isolate-box')) {
        if (candidate.getAttribute('data-sample') === sample) { box = candidate; break; }
    }
    if (!box) return;
    const details = box.querySelector('.box-details');
    const btn = box.querySelector('.toggle-btn');
    const isOpen = details.style.display === 'block';
    if (isOpen) {
        details.style.display = 'none';
        if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
        return;
    }
    if (!details.innerHTML.trim()) {
        details.innerHTML = buildBoxContent(tabId, sample);
    }
    details.style.display = 'block';
    if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
}

function buildBoxContent(tabId, sample) {
    const tabData = window.KLEBOSCOPE_BOX_DATA[tabId];
    if (!tabData) return '<p>No data</p>';
    const entry = tabData[sample];
    if (!entry) return '<p>No data</p>';

    if (tabId === 'mutations') {
        return buildMutationTable(entry.mutations || []);
    }
    if (tabId === 'mge') {
        return buildMgeTable(entry.mge || []);
    }

    let html = '';
    if (entry.amrfinder && entry.amrfinder.length) {
        html += buildGenericTable('AMRfinderPlus', entry.amrfinder, 'amrfinder');
    }
    if (entry.abricate) {
        for (const db in entry.abricate) {
            if (entry.abricate[db] && entry.abricate[db].length) {
                html += buildGenericTable(db.toUpperCase(), entry.abricate[db], db);
            }
        }
    }
    return html || '<p>No gene details</p>';
}

function buildGenericTable(dbName, rows, dbKey) {
    if (!rows || !rows.length) return '';
    const keys = Object.keys(rows[0]);
    const priority = ['gene', 'product', 'coverage_percent', 'identity_percent',
                      'coverage', 'identity', 'class', 'subclass', 'resistance',
                      'accession', 'contig', 'sequence', 'start', 'end', 'stop',
                      'strand', 'scope', 'element_type', 'database'];
    const ordered = [];
    priority.forEach(p => { if (keys.indexOf(p) !== -1) ordered.push(p); });
    keys.forEach(k => { if (ordered.indexOf(k) === -1) ordered.push(k); });

    let html = '<div class="database-table-wrapper" data-db="' + escapeHtml(dbKey) + '">';
    html += '<div class="db-title">' + escapeHtml(dbName) + '</div>';
    html += '<table><thead><tr>';
    ordered.forEach(c => {
        const display = c.replace(/_/g, ' ').replace(/\\b\\w/g, x => x.toUpperCase());
        html += '<th>' + escapeHtml(display) + '</th>';
    });
    html += '</tr></thead><tbody>';
    rows.forEach(r => {
        html += '<tr>';
        ordered.forEach(c => {
            let v = r[c];
            if (v == null) v = '';
            html += '<td>' + escapeHtml(v) + '</td>';
        });
        html += '</tr>';
    });
    html += '</tbody></table></div>';
    return html;
}

function buildMutationTable(mutations) {
    if (!mutations || !mutations.length) return '';
    const cols = ['gene', 'mutation', 'class', 'subclass', 'contig', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession'];
    let html = '<div class="database-table-wrapper" data-db="mutations">';
    html += '<div class="db-title">Mutations</div>';
    html += '<table><thead><tr>';
    cols.forEach(c => {
        html += '<th>' + c.replace(/_/g, ' ').replace(/\\b\\w/g, x => x.toUpperCase()) + '</th>';
    });
    html += '</tr></thead><tbody>';
    mutations.forEach(m => {
        html += '<tr>';
        cols.forEach(c => {
            let v = m[c];
            if (v == null) v = '';
            html += '<td>' + escapeHtml(v) + '</td>';
        });
        html += '</tr>';
    });
    html += '</tbody></table></div>';
    return html;
}

function buildMgeTable(mgeHits) {
    if (!mgeHits || !mgeHits.length) return '';
    const cols = ['Rank', 'mobileOG_ID', 'Gene', 'Contig', 'Start', 'Stop', 'Strand',
                  'Major_Category', 'Minor_Category', 'Percent_Identity',
                  'Align_Over_Qlen', 'Coverage_Percent'];
    const headers = ['Rank', 'mobileOG ID', 'Gene', 'Contig', 'Start', 'Stop', 'Strand',
                     'Major Category', 'Minor Category', '% Identity', 'Align / Qlen', 'Coverage %'];
    let html = '<div class="database-table-wrapper" data-db="mge">';
    html += '<div class="db-title">MGE (mobileOG-db)</div>';
    html += '<table><thead><tr>';
    headers.forEach(h => { html += '<th>' + h + '</th>'; });
    html += '</tr></thead><tbody>';
    mgeHits.forEach(m => {
        html += '<tr>';
        cols.forEach(c => {
            let v = m[c];
            if (v == null) v = '';
            html += '<td>' + escapeHtml(v) + '</td>';
        });
        html += '</tr>';
    });
    html += '</tbody></table></div>';
    return html;
}

/* ---------------- BULK EXPAND / COLLAPSE ---------------- */

function expandAllBoxes(tabId) {
    const tab = document.getElementById(tabId + '-tab');
    if (!tab) return;
    const boxes = tab.querySelectorAll('.isolate-box');
    const visible = [];
    boxes.forEach(b => { if (b.style.display !== 'none') visible.push(b); });
    if (visible.length > 100 && !confirm('You have ' + visible.length + ' visible samples. Expanding all may take a moment. Continue?')) return;
    let i = 0;
    const batch = 20;
    function step() {
        const end = Math.min(i + batch, visible.length);
        for (; i < end; i++) {
            const box = visible[i];
            const sample = box.getAttribute('data-sample');
            const details = box.querySelector('.box-details');
            const btn = box.querySelector('.toggle-btn');
            if (details && details.style.display !== 'block') {
                if (!details.innerHTML.trim()) details.innerHTML = buildBoxContent(tabId, sample);
                details.style.display = 'block';
                if (btn) btn.innerHTML = '<i class="fas fa-chevron-up"></i> Hide Details';
            }
        }
        if (i < visible.length) setTimeout(step, 30);
    }
    step();
}

function collapseAllBoxes(tabId) {
    const tab = document.getElementById(tabId + '-tab');
    if (!tab) return;
    tab.querySelectorAll('.isolate-box').forEach(box => {
        const details = box.querySelector('.box-details');
        if (details) details.style.display = 'none';
        const btn = box.querySelector('.toggle-btn');
        if (btn) btn.innerHTML = '<i class="fas fa-chevron-down"></i> Show Details';
    });
}

function filterBoxes(tabId) {
    const search = (document.getElementById('search-' + tabId)?.value || '').toUpperCase();
    const dbFilter = (document.getElementById('dbFilter-' + tabId)?.value || 'all');
    const tab = document.getElementById(tabId + '-tab');
    if (!tab) return;
    const boxes = tab.querySelectorAll('.isolate-box');
    let shown = 0;
    boxes.forEach(box => {
        const sample = box.getAttribute('data-sample') || '';
        const show = !search || sample.toUpperCase().indexOf(search) !== -1;
        box.style.display = show ? '' : 'none';
        if (show) shown++;
        if (show && dbFilter !== 'all') {
            box.querySelectorAll('.database-table-wrapper').forEach(w => {
                w.style.display = (w.getAttribute('data-db') === dbFilter) ? '' : 'none';
            });
        } else if (show) {
            box.querySelectorAll('.database-table-wrapper').forEach(w => { w.style.display = ''; });
        }
    });
    const counter = document.getElementById('counter-' + tabId);
    if (counter) counter.textContent = shown + ' shown';
}

function resetBoxFilters(tabId) {
    const s = document.getElementById('search-' + tabId);
    const d = document.getElementById('dbFilter-' + tabId);
    if (s) s.value = '';
    if (d) d.value = 'all';
    filterBoxes(tabId);
}

document.addEventListener('DOMContentLoaded', function() {
    const hash = window.location.hash.substring(1);
    const target = hash ? document.querySelector('.tab-button.' + hash) : document.querySelector('.tab-button');
    if (target) target.click();
    document.querySelectorAll('.data-table').forEach(table => {
        table.querySelectorAll('th').forEach((th, idx) => {
            const type = th.getAttribute('data-sort') || 'string';
            th.style.cursor = 'pointer';
            th.addEventListener('click', () => sortTable(table.id, idx, type));
        });
    });
    document.querySelectorAll('.copy-btn').forEach(btn => {
        btn.addEventListener('click', function() {
            const c = this.getAttribute('data-citation');
            if (c) navigator.clipboard.writeText(c).then(() => {
                const orig = this.innerHTML;
                this.innerHTML = '✓ Copied';
                setTimeout(() => { this.innerHTML = orig; }, 2000);
            });
        });
    });
});
</script>"""


# =============================================================================
# REPORTER
# =============================================================================
class KleboscopeUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "KLEBOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = KleboHTMLParser()
        self.generator = KleboHTMLGenerator()
        self.metadata = {
            'tool_name': 'Kleboscope Ultimate Reporter',
            'version': '1.2.0',
            'author': 'Brown Beckley',
            'affiliation': 'University of Ghana Medical School',
            'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'input_directory': str(self.input_dir),
        }

    def integrate(self) -> Dict:
        integrated = {'metadata': self.metadata, 'samples': {}}
        print("\nLoading data from TSV sources...")

        mlst = self.parser.load_mlst_from_tsv(self.input_dir)
        print(f"  MLST: {len(mlst)} samples")
        kaptive = self.parser.load_kaptive_from_tsv(self.input_dir)
        print(f"  Kaptive: {len(kaptive)} samples")
        amr_details = self.parser.load_amrfinder_details_from_tsv(self.input_dir)
        print(f"  AMRfinder: {len(amr_details)} samples")
        abricate_details = self.parser.load_abricate_details_from_tsv(self.input_dir)
        print(f"  ABRicate: {len(abricate_details)} samples")
        mutation_details = self.parser.load_mutations_from_tsv(self.input_dir)
        print(f"  Mutations: {len(mutation_details)} samples")
        mge_details = self.parser.load_mge_from_tsv(self.input_dir)
        print(f"  MGE: {len(mge_details)} samples")

        all_samples = (
            set(mlst) | set(kaptive) | set(amr_details) |
            set(abricate_details) | set(mutation_details) | set(mge_details)
        )
        for s in sorted(all_samples):
            entry = {
                'mlst': mlst.get(s, {}),
                'kaptive': kaptive.get(s, {}),
                'abricate_genes': {},
            }
            if s in abricate_details:
                for db, hits in abricate_details[s].items():
                    entry['abricate_genes'][db] = [
                        h.get('gene', '') for h in hits if h.get('gene')
                    ]
            integrated['samples'][s] = entry

        integrated['amrfinder_details'] = amr_details
        integrated['abricate_details'] = abricate_details
        integrated['mutation_details'] = mutation_details
        integrated['mge_details'] = mge_details

        print(f"\nTotal unique samples: {len(integrated['samples'])}")
        return integrated

    def write_json(self, data: Dict) -> Path:
        out = self.output_dir / "kleboscope_ultimate_sample_centric_report.json"
        with open(out, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, default=str, ensure_ascii=False)
        print(f"  ✅ JSON saved: {out}")
        return out

    def write_csvs(self, data: Dict):
        samples = data['samples']
        rows = []
        for s, sd in samples.items():
            rows.append({
                'Sample': s,
                'ST': sd.get('mlst', {}).get('ST', 'ND'),
                'K_Locus': sd.get('kaptive', {}).get('K_Locus', 'ND'),
                'O_Locus': sd.get('kaptive', {}).get('O_Locus', 'ND'),
            })
        pd.DataFrame(rows).to_csv(self.output_dir / 'sample_overview.csv', index=False)
        print(f"  ✅ CSVs saved to {self.output_dir}")

    def run(self) -> bool:
        print("=" * 80)
        print("Kleboscope Ultimate Reporterv2.0.0")
        print("=" * 80)
        print(f"Input:  {self.input_dir}")
        print(f"Output: {self.output_dir}")

        integrated = self.integrate()
        if not integrated['samples']:
            print("❌ No samples integrated.")
            return False

        self.write_json(integrated)
        self.write_csvs(integrated)
        self.generator.generate_main_report(integrated, self.output_dir)

        print()
        print("=" * 80)
        print("✅ REPORT COMPLETE")
        print("=" * 80)
        print(f"  Samples: {len(integrated['samples'])}")
        print(f"  Output:  {self.output_dir}")
        print("=" * 80)
        return True


def main():
    parser = argparse.ArgumentParser(
        description='Kleboscope Ultimate Reporter – Sample-Centric'
    )
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Directory with Kleboscope TSV outputs')
    args = parser.parse_args()

    inp = Path(args.input_dir)
    if not inp.exists():
        print(f"❌ Input directory not found: {inp}")
        sys.exit(1)

    reporter = KleboscopeUltimateReporter(inp)
    if not reporter.run():
        sys.exit(1)


if __name__ == '__main__':
    main()