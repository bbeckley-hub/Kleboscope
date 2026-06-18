#!/usr/bin/env python3
"""
Kleboscope AMRfinderPlus - K. pneumoniae Antimicrobial Resistance
Comprehensive AMR analysis with interactive per‑genome reports and clean batch summary
Author: Brown Beckley
Affiliation: University of Ghana Medical School
Version: 1.1.0
Uses bundled AMRfinderPlus with dynamic latest database version.
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import json
import random
import shutil
from collections import defaultdict

class KleboAMRfinderPlus:
    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_update = os.path.join(self.module_dir, "bin", "amrfinder_update")
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)

        self.critical_carbapenemases = {
            'blaKPC', 'blaNDM', 'blaIMP', 'blaVIM', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaOXA-244', 'blaGES-2', 'blaGES-4', 'blaGES-5', 'blaGES-6', 'blaGES-7', 'blaGES-8',
            'blaGES-11', 'blaGES-14', 'blaGES-15', 'blaGES-16', 'blaGES-17', 'blaGES-18', 'blaGES-19',
            'blaGES-20', 'blaGES-21', 'blaGES-22', 'blaGES-23', 'blaGES-24', 'blaGES-25', 'blaGES-26',
            'blaGES-27', 'blaGES-28', 'blaIMI', 'blaSME', 'blaNMC', 'blaCcrA', 'blaBIC', 'blaDIM',
            'blaTMB', 'blaSPM', 'blaGIM', 'blaSIM', 'blaAIM', 'blaFRI', 'blaFRI-1', 'blaFRI-2'
        }
        self.critical_esbls = {
            'blaCTX-M-15', 'blaCTX-M-14', 'blaCTX-M-1', 'blaCTX-M-2', 'blaCTX-M-3', 'blaCTX-M-4',
            'blaCTX-M-5', 'blaCTX-M-6', 'blaCTX-M-7', 'blaCTX-M-8', 'blaCTX-M-9', 'blaCTX-M-10',
            'blaSHV', 'blaTEM', 'blaPER', 'blaVEB', 'blaGES-1', 'blaGES-3', 'blaBEL', 'blaBES'
        }
        self.critical_colistin = {
            'mcr-1', 'mcr-1.1', 'mcr-2', 'mcr-3', 'mcr-3.1', 'mcr-3.2', 'mcr-3.3', 'mcr-3.4', 'mcr-3.5',
            'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 'mcr-8', 'mcr-9', 'mcr-10',
            'pmrA', 'pmrB', 'pmrC', 'lpxA', 'lpxC', 'lpxD', 'arnA', 'arnB', 'arnC', 'arnD', 'eptA'
        }
        self.critical_aminoglycoside = {
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH', 'npmA'
        }
        self.high_risk_resistance = {
            'blaCTX-M-31', 'blaCTX-M-32', 'blaCTX-M-33', 'blaCTX-M-34', 'blaCTX-M-35', 'blaCTX-M-36',
            'blaCTX-M-37', 'blaCTX-M-38', 'blaCTX-M-39', 'blaCTX-M-40', 'blaCTX-M-41', 'blaCTX-M-42',
            'blaCTX-M-43', 'blaCTX-M-44', 'blaCTX-M-45', 'blaCTX-M-46', 'blaCTX-M-47', 'blaCTX-M-48',
            'blaCTX-M-49', 'blaCTX-M-50', 'blaCTX-M-51', 'blaCTX-M-52', 'blaCTX-M-53', 'blaCTX-M-54',
            'blaCTX-M-55', 'blaCTX-M-56', 'blaCTX-M-57', 'blaCTX-M-58', 'blaCTX-M-59', 'blaCTX-M-60',
            'blaCMY', 'blaDHA', 'blaACC', 'blaMIR', 'blaACT', 'blaFOX', 'blaMOX', 'blaCIT', 'blaEBC',
            'aac(3)-I', 'aac(3)-II', 'aac(3)-III', 'aac(3)-IV', 'aac(6\')-I', 'aac(6\')-II', 'aadA', 'aadB',
            'aph(3\')-I', 'aph(3\')-II', 'aph(6)-I', 'aph(6)-II',
            'qnrA', 'qnrB', 'qnrC', 'qnrD', 'qnrS', 'qnrVC', 'aac(6\')-Ib-cr', 'qepA',
            'tetA', 'tetB', 'tetC', 'tetD', 'tetE', 'tetG', 'tetH', 'tetK', 'tetL', 'tetM', 'tetO', 'tetQ', 'tetS', 'tetX',
            'sul1', 'sul2', 'sul3', 'sul4',
            'dfrA', 'dfrA1', 'dfrA5', 'dfrA7', 'dfrA8', 'dfrA12', 'dfrA14', 'dfrA17', 'dfrA19', 'dfrA20', 'dfrA21',
            'dfrB1', 'dfrB2', 'dfrB3', 'dfrB4', 'dfrB5', 'dfrB6', 'dfrB7',
            'cat', 'catA1', 'catA2', 'catB2', 'catB3', 'cmlA', 'cmlA1', 'cmlA5', 'floR',
            'ermA', 'ermB', 'ermC', 'ermF', 'ermX', 'mphA', 'mphB', 'msrA',
            'fosA', 'fosB', 'fosC', 'fosX',
            'acrA', 'acrB', 'acrD', 'acrF', 'acrS', 'acrZ',
            'adeA', 'adeB', 'adeC', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeM', 'adeN',
            'mexA', 'mexB', 'mexC', 'mexD', 'mexE', 'mexF', 'mexG', 'mexH', 'mexI',
            'mexJ', 'mexK', 'mexL', 'mexM', 'mexN', 'mexO', 'mexP', 'mexQ',
            'mdtA', 'mdtB', 'mdtC', 'mdtD', 'mdtE', 'mdtF', 'mdtG', 'mdtH', 'mdtI', 'mdtJ',
            'mdtK', 'mdtL', 'mdtM', 'mdtN', 'mdtO', 'mdtP',
            'emrA', 'emrB', 'emrD', 'emrE', 'emrK', 'emrY'
        }
        self.high_risk_genes = set.union(
            self.critical_carbapenemases,
            self.critical_esbls,
            self.critical_colistin,
            self.critical_aminoglycoside,
            self.high_risk_resistance
        )
        self.critical_risk_genes = set.union(
            self.critical_carbapenemases,
            self.critical_colistin
        )

        self.ascii_art = r"""
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
"""
        self.metadata = {
            "tool_name": "Kleboscope AMRfinderPlus",
            "version": "1.1.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7",
            "database_version": "auto-detected"
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
            "Kleboscope turns genomic complexity into actionable insights for AMR surveillance. - Brown Beckley"
        ]

        self.bundled_database = self._get_latest_database_path()
        self.metadata["database_version"] = os.path.basename(self.bundled_database) if self.bundled_database else "None"
        self.logger.info(f"Using database version: {self.metadata['database_version']}")

    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

    def _get_available_ram(self) -> int:
        try:
            return psutil.virtual_memory().available / (1024 ** 3)
        except Exception:
            return 8

    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total <= 4:
                optimal = total
            elif total <= 8:
                optimal = total - 1
            elif total <= 16:
                optimal = max(8, total - 2)
            elif total <= 32:
                optimal = max(16, total - 3)
            else:
                optimal = min(32, int(total * 0.95))
            optimal = max(1, min(optimal, total))
            self._log_resource_info(optimal, total)
            return optimal
        except Exception:
            return os.cpu_count() or 4

    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            self.logger.info(f"Using CPU cores: {cpus} ({cpus/total_cores*100:.1f}% of available)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus <= 4:
            self.logger.info("💡 Performance: Multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("💡 Performance: High-speed mode")
        else:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")

    def _run_amrfinder_update(self, force: bool = False):
        if not os.path.exists(self.bundled_update):
            raise FileNotFoundError(f"amrfinder_update not found at {self.bundled_update}. Run database.sh first.")
        db_parent = os.path.join(self.module_dir, "data", "amrfinder_db")
        if force and os.path.exists(db_parent):
            self.logger.info("Forcing full database update – removing old folders...")
            for item in os.listdir(db_parent):
                full = os.path.join(db_parent, item)
                if os.path.isdir(full) and re.match(r'^\d{4}-\d{2}-\d{2}\.\d+$', item):
                    self.logger.info(f"  Removing {item}")
                    shutil.rmtree(full)
        os.makedirs(db_parent, exist_ok=True)
        self.logger.info("Downloading latest AMR database via amrfinder_update (this may take several minutes)...")
        try:
            subprocess.run([self.bundled_update, "--database", db_parent], check=True, capture_output=True, text=True)
            self.logger.info("AMR database download completed successfully.")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"amrfinder_update failed: {e.stderr}")
            raise RuntimeError("Failed to download AMR database. Check internet connection and permissions.")

    def _get_latest_database_path(self) -> str:
        db_parent = os.path.join(self.module_dir, "data", "amrfinder_db")
        if not os.path.isdir(db_parent):
            self.logger.warning(f"Database parent directory not found: {db_parent}")
            self._run_amrfinder_update(force=False)
            db_parent = os.path.join(self.module_dir, "data", "amrfinder_db")
        version_dirs = []
        for entry in os.listdir(db_parent):
            full = os.path.join(db_parent, entry)
            if os.path.isdir(full) and re.match(r'^\d{4}-\d{2}-\d{2}\.\d+$', entry):
                version_dirs.append((entry, full))
        if not version_dirs:
            self.logger.warning(f"No version directory found in {db_parent}. Running update...")
            self._run_amrfinder_update(force=False)
            version_dirs = []
            for entry in os.listdir(db_parent):
                full = os.path.join(db_parent, entry)
                if os.path.isdir(full) and re.match(r'^\d{4}-\d{2}-\d{2}\.\d+$', entry):
                    version_dirs.append((entry, full))
            if not version_dirs:
                raise RuntimeError(f"Still no version directory found after update in {db_parent}")
        version_dirs.sort(key=lambda x: x[0], reverse=True)
        latest = version_dirs[0][1]
        self.logger.info(f"Latest AMR database version: {os.path.basename(latest)}")
        return latest

    def update_database(self, force: bool = False) -> bool:
        try:
            self._run_amrfinder_update(force=force)
            self.bundled_database = self._get_latest_database_path()
            self.metadata["database_version"] = os.path.basename(self.bundled_database)
            self.logger.info(f"New database version: {self.metadata['database_version']}")
            return True
        except Exception as e:
            self.logger.error(f"Database update failed: {e}")
            return False

    def check_amrfinder_installed(self) -> bool:
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
                return False
            if not os.access(self.bundled_amrfinder, os.X_OK):
                self.logger.warning("Fixing permissions on bundled AMRfinderPlus")
                os.chmod(self.bundled_amrfinder, 0o755)
            result = subprocess.run([self.bundled_amrfinder, '--version'],
                                    capture_output=True, text=True, check=True)
            self.logger.info(f"Bundled AMRfinderPlus version: {result.stdout.strip()}")
            if self.bundled_database and os.path.exists(self.bundled_database):
                self.logger.info(f"✅ Bundled database found: {self.bundled_database}")
                self.metadata["database_version"] = os.path.basename(self.bundled_database)
                self.logger.info(f"✅ Database version: {self.metadata['database_version']}")
                return True
            else:
                self.logger.warning(f"⚠️ Bundled database not found at: {self.bundled_database}")
                return False
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Bundled AMRfinderPlus check failed: {e}")
            return False

    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str,
                                     min_identity: Optional[float] = None,
                                     min_coverage: Optional[float] = None,
                                     report_mutations: bool = True) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        cmd = [
            self.bundled_amrfinder,
            '-n', genome_file,
            '--output', output_file,
            '--threads', str(self.cpus),
            '--plus',
            '--organism', 'Klebsiella_pneumoniae'
        ]
        if self.bundled_database and os.path.exists(self.bundled_database):
            cmd.extend(['--database', self.bundled_database])
            self.logger.info(f"Using bundled database: {self.bundled_database}")
        else:
            self.logger.warning("No database found; using AMRfinderPlus default (may fail).")
        if min_identity is not None:
            cmd.extend(['--ident_min', str(min_identity)])
            self.logger.info(f"Minimum identity: {min_identity}")
        if min_coverage is not None:
            cmd.extend(['--coverage_min', str(min_coverage)])
            self.logger.info(f"Minimum coverage: {min_coverage}")

        mut_file = None
        if report_mutations:
            mut_file = os.path.join(output_dir, f"{genome_name}_mutations.tsv")
            cmd.extend(['--mutation_all', mut_file])
            self.logger.info(f"Will report mutations to {mut_file}")

        self.logger.info(f"Running BUNDLED AMRfinderPlus: {genome_name} (using {self.cpus} CPU cores)")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = self._parse_amrfinder_output(output_file)
            self._create_amrfinder_html_report(genome_name, hits, output_dir)
            if report_mutations and mut_file and os.path.exists(mut_file):
                self._create_mutation_html_report(genome_name, mut_file, output_dir)
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'mutations_file': mut_file if (report_mutations and mut_file and os.path.exists(mut_file)) else None,
                'status': 'success'
            }
        except subprocess.CalledProcessError as e:
            self.logger.error(f"AMRfinderPlus failed for {genome_name}: {e.stderr}")
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'mutations_file': None,
                'status': 'failed'
            }

    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
            if len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            for line in lines[1:]:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = dict(zip(headers, parts))
                    processed = {
                        'protein_id': hit.get('Protein id', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': str(hit.get('Element symbol', '')),
                        'sequence_name': hit.get('Element name', ''),
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Type', ''),
                        'element_subtype': hit.get('Subtype', ''),
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference', '').replace('%', ''),
                        'identity': hit.get('% Identity to reference', '').replace('%', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Closest reference accession', ''),
                        'closest_name': hit.get('Closest reference name', ''),
                        'hmm_id': hit.get('HMM accession', ''),
                        'hmm_description': hit.get('HMM description', '')
                    }
                    hits.append(processed)
        except Exception as e:
            self.logger.error(f"Error parsing {amrfinder_file}: {e}")
        self.logger.info(f"Parsed {len(hits)} AMR hits from {amrfinder_file}")
        return hits

    def _parse_mutations_file(self, mut_file: str) -> List[Dict]:
        return self._parse_amrfinder_output(mut_file)

    def _create_mutation_html_report(self, genome_name: str, mutations_file: str, output_dir: str):
        mutations = self._parse_mutations_file(mutations_file)
        if not mutations:
            self.logger.info(f"No mutations found for {genome_name}, skipping mutation HTML.")
            return
        random_quote = random.choice(self.science_quotes)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Kleboscope - Mutation Report: {genome_name}</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', sans-serif;
            color: #fff;
            padding:20px;
        }}
        .container {{ max-width:1400px; margin:auto; }}
        .header {{ text-align:center; margin-bottom:30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding:20px;
            border-radius:15px;
            margin-bottom:20px;
            border:2px solid rgba(0,119,190,0.3);
        }}
        .ascii-art {{
            font-family: monospace;
            font-size:12px;
            white-space:pre;
            color:#0077be;
            text-shadow:0 0 10px rgba(0,119,190,0.5);
            overflow-x:auto;
        }}
        .card {{
            background: rgba(255,255,255,0.95);
            color:#1f2937;
            padding:25px;
            margin:20px 0;
            border-radius:12px;
        }}
        .card h2 {{
            color:#1e3a8a;
            border-bottom:3px solid #3b82f6;
            padding-bottom:10px;
            margin-bottom:20px;
        }}
        .gene-table {{
            width:100%;
            border-collapse:collapse;
            margin:20px 0;
            background:white;
        }}
        .gene-table th, .gene-table td {{
            padding:12px;
            text-align:left;
            border-bottom:1px solid #e5e7eb;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color:white;
        }}
        .gene-table tr:nth-child(even) {{ background-color:#f8fafc; }}
        .gene-table tr:hover {{ background-color:#e0f2fe; }}
        .footer {{
            background: rgba(0,0,0,0.8);
            color:white;
            padding:20px;
            border-radius:12px;
            margin-top:40px;
            text-align:center;
        }}
        .timestamp {{ color:#fbbf24; }}
        .authorship {{ margin-top:15px; font-size:0.9em; }}
        .quote-container {{
            background: rgba(255,255,255,0.1);
            color:white;
            padding:20px;
            border-radius:12px;
            margin:20px 0;
            text-align:center;
            font-style:italic;
            border-left:4px solid #fff;
        }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="card">
            <h1>🧬 Kleboscope – Point Mutation Report</h1>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {current_time}</p>
        </div>
    </div>
    <div class="quote-container">
        <div style="font-size:1.1em;">"{random_quote}"</div>
    </div>
    <div class="card">
        <h2>All Point Mutations Detected</h2>
        <table class="gene-table">
            <thead><tr><th>Gene Symbol</th><th>Mutation</th><th>Class</th><th>Subclass</th>
                <th>Contig</th><th>Start</th><th>Stop</th><th>Strand</th>
                <th>Coverage (%)</th><th>Identity (%)</th><th>Accession</th></tr></thead>
            <tbody>
"""
        for m in mutations:
            html += f"""
                <tr>
                    <td>{m.get('gene_symbol', '')}</td>
                    <td>{m.get('sequence_name', '')}</td>
                    <td>{m.get('class', '')}</td>
                    <td>{m.get('subclass', '')}</td>
                    <td>{m.get('contig_id', '')}</td>
                    <td>{m.get('start', '')}</td>
                    <td>{m.get('stop', '')}</td>
                    <td>{m.get('strand', '')}</td>
                    <td>{m.get('coverage', '')}</td>
                    <td>{m.get('identity', '')}</td>
                    <td>{m.get('accession', '')}</td>
                </tr>
"""
        html += f"""
            </tbody>
        </table>
    </div>
    <div class="footer">
        <p><strong>KLEBOSCOPE</strong> – Mutation Analysis Module</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>University of Ghana Medical School – Department of Medical Biochemistry</p>
        </div>
    </div>
</div>
</body>
</html>"""
        out_file = os.path.join(output_dir, f"{genome_name}_mutations.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Mutation HTML report: {out_file}")

    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        analysis = self._analyze_kpneumoniae_amr_results(hits)
        interactive_js = f"""
        <script>
            function searchTable(tableId, searchTerm) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let visibleCount = 0;
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const text = row.textContent.toLowerCase();
                    if (text.includes(searchTerm.toLowerCase())) {{
                        row.style.display = '';
                        visibleCount++;
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
                const counter = document.getElementById('result-counter-' + tableId);
                if (counter) counter.textContent = visibleCount + ' results found';
            }}
            function exportToCSV(tableId, filename) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let csv = [];
                const headerCells = rows[0].getElementsByTagName('th');
                const headerRow = [];
                for (let cell of headerCells) headerRow.push(cell.textContent);
                csv.push(headerRow.join(','));
                for (let i = 1; i < rows.length; i++) {{
                    if (rows[i].style.display !== 'none') {{
                        const cells = rows[i].getElementsByTagName('td');
                        const row = [];
                        for (let cell of cells) {{
                            let text = cell.textContent.trim();
                            text = text.replace(/\\n/g, ' ').replace(/,/g, ';');
                            row.push(text);
                        }}
                        csv.push(row.join(','));
                    }}
                }}
                const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function exportToJSON(dataVar, filename) {{
                const data = window[dataVar];
                const jsonStr = JSON.stringify(data, null, 2);
                const blob = new Blob([jsonStr], {{ type: 'application/json' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function printReport() {{
                const printWindow = window.open('', '_blank');
                printWindow.document.write('<html><head><title>{genome_name} - AMR Report</title><style>body {{ font-family: Arial; margin: 20px; }} .no-print {{ display: none; }} table {{ border-collapse: collapse; width: 100%; }} th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }} .critical-row {{ background-color: #f8d7da; }} .high-risk-row {{ background-color: #fff3cd; }}</style></head><body>');
                const content = document.querySelector('.container').cloneNode(true);
                const noPrint = content.querySelectorAll('.no-print');
                noPrint.forEach(el => el.remove());
                printWindow.document.write(content.innerHTML);
                printWindow.document.write('</body></html>');
                printWindow.document.close();
                printWindow.print();
            }}
            function quickSearch() {{
                const term = document.getElementById('quick-search').value.toLowerCase();
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    const text = section.textContent.toLowerCase();
                    if (text.includes(term)) {{
                        section.style.border = '2px solid #3b82f6';
                        section.style.backgroundColor = 'rgba(59,130,246,0.1)';
                        section.scrollIntoView({{ behavior: 'smooth', block: 'nearest' }});
                    }} else {{
                        section.style.border = '';
                        section.style.backgroundColor = '';
                    }}
                }});
            }}
            function clearSearch() {{
                document.getElementById('quick-search').value = '';
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    section.style.border = '';
                    section.style.backgroundColor = '';
                }});
            }}
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                const quoteDiv = document.getElementById('science-quote');
                if (quoteDiv) quoteDiv.innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});
        </script>
        """
        export_data_js = f"""
        <script>
            window.reportData = {{
                metadata: {{
                    genome: '{genome_name}',
                    date: '{self.metadata['analysis_date']}',
                    tool: '{self.metadata['tool_name']}',
                    version: '{self.metadata['version']}'
                }},
                summary: {json.dumps(analysis)},
                hits: {json.dumps(hits)}
            }};
        </script>
        """
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Kleboscope AMRfinderPlus Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding:20px;
            min-height:100vh;
        }}
        .container {{ max-width:1400px; margin:0 auto; }}
        .header {{ text-align:center; margin-bottom:30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding:20px;
            border-radius:15px;
            margin-bottom:20px;
            box-shadow:0 8px 32px rgba(0,0,0,0.4);
            border:2px solid rgba(0,119,190,0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size:12px;
            line-height:1.2;
            white-space:pre;
            color:#0077be;
            text-shadow:0 0 10px rgba(0,119,190,0.5);
            overflow-x:auto;
        }}
        .card {{
            background: rgba(255,255,255,0.95);
            color:#1f2937;
            padding:25px;
            margin:20px 0;
            border-radius:12px;
            box-shadow:0 4px 20px rgba(0,0,0,0.2);
        }}
        .card h2 {{
            color:#1e3a8a;
            border-bottom:3px solid #3b82f6;
            padding-bottom:10px;
            margin-bottom:20px;
        }}
        .gene-table, .class-table {{
            width:100%;
            border-collapse:collapse;
            margin:20px 0;
            background:white;
            border-radius:8px;
            overflow:hidden;
        }}
        .gene-table th, .gene-table td, .class-table th, .class-table td {{
            padding:15px;
            text-align:left;
            border-bottom:1px solid #e0e0e0;
        }}
        .gene-table th, .class-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color:white;
            font-weight:600;
        }}
        tr:hover {{ background-color:#f8f9fa; }}
        .summary-stats {{
            display:flex;
            justify-content:space-around;
            margin:20px 0;
            flex-wrap:wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color:white;
            padding:20px;
            border-radius:12px;
            text-align:center;
            box-shadow:0 4px 15px rgba(0,0,0,0.2);
            margin:10px;
            flex:1;
            min-width:200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc3545 0%, #c82333 100%);
            color:white;
            padding:20px;
            border-radius:12px;
            text-align:center;
            box-shadow:0 4px 15px rgba(0,0,0,0.2);
            margin:10px;
            flex:1;
            min-width:200px;
        }}
        .quote-container {{
            background: rgba(255,255,255,0.1);
            color:white;
            padding:20px;
            border-radius:12px;
            margin:20px 0;
            text-align:center;
            font-style:italic;
            border-left:4px solid #fff;
        }}
        .footer {{
            background: rgba(0,0,0,0.8);
            color:white;
            padding:30px;
            border-radius:12px;
            margin-top:40px;
        }}
        .resistance-badge {{
            display:inline-block;
            background:#dc3545;
            color:white;
            padding:5px 10px;
            border-radius:15px;
            margin:2px;
            font-size:0.9em;
        }}
        .critical-risk-badge {{
            display:inline-block;
            background:#8b0000;
            color:white;
            padding:5px 10px;
            border-radius:15px;
            margin:2px;
            font-size:0.9em;
            font-weight:bold;
        }}
        .warning-badge {{
            display:inline-block;
            background:#ffc107;
            color:black;
            padding:5px 10px;
            border-radius:15px;
            margin:2px;
            font-size:0.9em;
        }}
        .success-badge {{
            display:inline-block;
            background:#28a745;
            color:white;
            padding:5px 10px;
            border-radius:15px;
            margin:2px;
            font-size:0.9em;
        }}
        .present {{ background-color:#d4edda; }}
        .critical-row {{ background-color:#f8d7da; font-weight:bold; border-left:4px solid #dc3545; }}
        .high-risk-row {{ background-color:#fff3cd; border-left:4px solid #ffc107; }}
        .interactive-controls {{
            background:#f8f9fa;
            padding:15px;
            border-radius:8px;
            margin:15px 0;
            display:flex;
            flex-wrap:wrap;
            gap:10px;
            align-items:center;
        }}
        .search-box input {{
            width:100%;
            padding:8px 12px;
            border:1px solid #ddd;
            border-radius:4px;
            font-size:14px;
        }}
        .export-buttons button {{
            padding:8px 16px;
            background:#3b82f6;
            color:white;
            border:none;
            border-radius:4px;
            cursor:pointer;
        }}
        .export-buttons button.print {{ background:#10b981; }}
        .result-counter {{ font-size:0.9em; color:#666; font-style:italic; }}
        .sequence-name {{ max-width:300px; overflow-wrap:break-word; }}
        .quick-search-bar {{
            background: rgba(255,255,255,0.95);
            padding:15px;
            border-radius:8px;
            margin:10px 0;
            display:flex;
            gap:10px;
        }}
        .quick-search-bar input {{ flex:1; padding:8px 12px; border:2px solid #3b82f6; border-radius:4px; }}
        .quick-search-bar button {{ padding:8px 20px; background:#3b82f6; color:white; border:none; border-radius:4px; cursor:pointer; }}
        @media print {{
            .no-print, .interactive-controls, .quick-search-bar {{ display:none; }}
            body {{ background:white; color:black; }}
            .card {{ background:white; box-shadow:none; }}
        }}
    </style>
    {interactive_js}
    {export_data_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="card">
                <h1 style="color:#333; font-size:2.5em;">🧬 Kleboscope AMRfinderPlus Analysis Report</h1>
                <p style="color:#666; font-size:1.2em;">Comprehensive K. pneumoniae Antimicrobial Resistance Analysis</p>
                <div class="quick-search-bar no-print">
                    <input type="text" id="quick-search" placeholder="🔍 Quick search across entire report...">
                    <button onclick="quickSearch()">Search</button>
                    <button onclick="clearSearch()" style="background:#6c757d;">Clear</button>
                </div>
                <div class="export-buttons no-print" style="margin-top:15px;">
                    <button onclick="exportToJSON('reportData', '{genome_name}_amr_report.json')">📥 Export JSON</button>
                    <button onclick="printReport()" class="print">🖨️ Print Report</button>
                </div>
            </div>
        </div>
        <div class="quote-container">
            <div id="science-quote" style="font-size:1.1em;"></div>
        </div>
"""
        if analysis['critical_risk_genes'] > 0:
            html_content += f"""
        <div class="card" style="border-left:4px solid #dc3545; background:#f8d7da;">
            <h2 style="color:#dc3545;">🚨 CRITICAL RISK AMR GENES DETECTED</h2>
            <p><strong>{analysis['critical_risk_genes']} CRITICAL RISK antimicrobial resistance genes found:</strong></p>
            <div style="margin:10px 0;">
                <p style="color:#721c24; font-weight:bold;">
                    ⚠️ These genes confer resistance to last-resort antibiotics and represent a serious public health concern.
                </p>
"""
            for gene in analysis['critical_risk_list']:
                html_content += f'<span class="critical-risk-badge">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        html_content += f"""
        <div class="card">
            <h2>📊 K. pneumoniae AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total AMR Genes</h3><p style="font-size:2em;">{analysis['total_genes']}</p></div>
                <div class="stat-card"><h3>High Risk Genes</h3><p style="font-size:2em;">{analysis['high_risk_genes']}</p></div>
                <div class="critical-stat-card"><h3>Critical Risk</h3><p style="font-size:2em;">{analysis['critical_risk_genes']}</p></div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus Version:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database Version:</strong> {self.metadata['database_version']}</p>
        </div>
"""
        if analysis['high_risk_genes'] > 0 and analysis['critical_risk_genes'] == 0:
            html_content += f"""
        <div class="card" style="border-left:4px solid #ffc107;">
            <h2 style="color:#856404;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{analysis['high_risk_genes']} high-risk antimicrobial resistance genes found:</strong></p>
            <div style="margin:10px 0;">
"""
            for gene in analysis['high_risk_list']:
                html_content += f'<span class="resistance-badge">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        if any(analysis['resistance_mechanisms'].values()):
            html_content += """
        <div class="card">
            <h2>🔬 Resistance Mechanism Breakdown</h2>
"""
            mech = analysis['resistance_mechanisms']
            if mech['carbapenemase']:
                html_content += f"""
            <div style="margin:10px 0; padding:10px; background:#f8d7da; border-radius:5px;">
                <strong>Carbapenemase Genes (CRITICAL):</strong> {', '.join(mech['carbapenemase'])}
            </div>
"""
            if mech['colistin_resistance']:
                html_content += f"""
            <div style="margin:10px 0; padding:10px; background:#f8d7da; border-radius:5px;">
                <strong>Colistin Resistance (CRITICAL):</strong> {', '.join(mech['colistin_resistance'])}
            </div>
"""
            if mech['esbl']:
                html_content += f"""
            <div style="margin:10px 0; padding:10px; background:#fff3cd; border-radius:5px;">
                <strong>ESBL Genes:</strong> {', '.join(mech['esbl'])}
            </div>
"""
            if mech['aminoglycoside_resistance']:
                html_content += f"""
            <div style="margin:10px 0; padding:10px; background:#d1ecf1; border-radius:5px;">
                <strong>Aminoglycoside Resistance:</strong> {', '.join(mech['aminoglycoside_resistance'])}
            </div>
"""
            if mech['fluoroquinolone_resistance']:
                html_content += f"""
            <div style="margin:10px 0; padding:10px; background:#d1ecf1; border-radius:5px;">
                <strong>Fluoroquinolone Resistance:</strong> {', '.join(mech['fluoroquinolone_resistance'])}
            </div>
"""
            if mech['efflux_pumps']:
                html_content += f"""
            <div style="margin:10px 0; padding:10px; background:#e2e3e5; border-radius:5px;">
                <strong>Efflux Pumps:</strong> {', '.join(mech['efflux_pumps'])}
            </div>
"""
            if mech['other_amr']:
                html_content += f"""
            <div style="margin:10px 0; padding:10px; background:#f8f9fa; border-radius:5px;">
                <strong>Other AMR Genes:</strong> {', '.join(mech['other_amr'])}
            </div>
"""
            html_content += "</div>"
        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2>🧪 Resistance Classes Detected</h2>
            <table class="class-table">
                <thead><tr><th>Resistance Class</th><th>Gene Count</th><th>Genes</th></tr></thead>
                <tbody>
"""
            for class_name, genes in analysis['resistance_classes'].items():
                html_content += f"""
                    <tr><td><strong>{class_name}</strong></td><td>{len(genes)}</td><td>{', '.join(genes)}</td></tr>
"""
            html_content += "</tbody></table></div>"
        if hits:
            html_content += f"""
        <div class="card">
            <h2>🔬 Detailed AMR Genes Detected</h2>
            <div class="interactive-controls no-print">
                <div class="search-box"><input type="text" id="search-detailed-amr" placeholder="Search genes, classes, or sequence names..." onkeyup="searchTable('detailed-amr-table', this.value)"></div>
                <div class="export-buttons"><button onclick="exportToCSV('detailed-amr-table', '{genome_name}_amr_genes.csv')">📥 Export CSV</button></div>
                <div class="result-counter" id="result-counter-detailed-amr-table">{len(hits)} results found</div>
            </div>
            <table class="gene-table" id="detailed-amr-table">
                <thead><tr><th>Gene Symbol</th><th>Sequence Name</th><th>Class</th><th>Subclass</th><th>Coverage</th><th>Identity</th><th>Scope</th></tr></thead>
                <tbody>
"""
            for hit in hits:
                gene = hit.get('gene_symbol', '')
                row_class = "present"
                if gene in analysis['critical_risk_list']:
                    row_class = "critical-row"
                elif gene in analysis['high_risk_list']:
                    row_class = "high-risk-row"
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{gene}</strong></td>
                        <td class="sequence-name">{hit.get('sequence_name', '')}</td>
                        <td>{hit.get('class', '')}</td>
                        <td>{hit.get('subclass', '')}</td>
                        <td>{hit.get('coverage', '')}%</td>
                        <td>{hit.get('identity', '')}%</td>
                        <td>{hit.get('scope', '')}</td>
                    </tr>
"""
            html_content += "</tbody></table></div>"
        else:
            html_content += """
        <div class="card"><h2>✅ No AMR Genes Detected</h2><p>No antimicrobial resistance genes found in this K. pneumoniae genome.</p></div>
"""
        html_content += f"""
        <div class="footer">
            <h3 style="color:#fff;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank" style="color:#3b82f6;">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
        </div>
    </div>
</body>
</html>"""
        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info(f"K. pneumoniae AMRfinderPlus HTML report generated: {html_file}")

    def _analyze_kpneumoniae_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'total_classes': 0,
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'carbapenemase': [],
                'colistin_resistance': [],
                'esbl': [],
                'aminoglycoside_resistance': [],
                'fluoroquinolone_resistance': [],
                'efflux_pumps': [],
                'other_amr': []
            }
        }
        for hit in hits:
            gene = str(hit.get('gene_symbol', ''))
            if not gene:
                continue
            self._categorize_resistance_mechanism(gene, analysis)
            if gene in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene)
            if gene in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene)
            cls = hit.get('class', '')
            if cls:
                if cls not in analysis['resistance_classes']:
                    analysis['resistance_classes'][cls] = []
                if gene not in analysis['resistance_classes'][cls]:
                    analysis['resistance_classes'][cls].append(gene)
        analysis['total_classes'] = len(analysis['resistance_classes'])
        return analysis

    def _categorize_resistance_mechanism(self, gene: str, analysis: Dict[str, Any]):
        gene = str(gene)
        if any(kw in gene for kw in ('KPC', 'NDM', 'IMP', 'VIM', 'OXA-48', 'GES-2', 'GES-4', 'GES-5',
                                     'GES-6', 'GES-7', 'GES-8', 'GES-11', 'IMI', 'SME', 'NMC', 'FRI')):
            analysis['resistance_mechanisms']['carbapenemase'].append(gene)
        elif any(kw in gene for kw in ('mcr-', 'MCR-')):
            analysis['resistance_mechanisms']['colistin_resistance'].append(gene)
        elif any(kw in gene for kw in ('CTX-M', 'SHV', 'TEM', 'PER', 'VEB', 'GES-1', 'GES-3', 'BEL', 'BES')):
            analysis['resistance_mechanisms']['esbl'].append(gene)
        elif any(kw in gene for kw in ('aac', 'aad', 'aph', 'armA', 'rmt', 'npmA')):
            analysis['resistance_mechanisms']['aminoglycoside_resistance'].append(gene)
        elif any(kw in gene for kw in ('qnr', 'aac(6\')-Ib-cr', 'qepA')):
            analysis['resistance_mechanisms']['fluoroquinolone_resistance'].append(gene)
        elif any(kw in gene for kw in ('acr', 'ade', 'mex', 'mdt', 'emr')):
            analysis['resistance_mechanisms']['efflux_pumps'].append(gene)
        else:
            analysis['resistance_mechanisms']['other_amr'].append(gene)

    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating K. pneumoniae AMR summary files...")
        summary_file = os.path.join(output_base, "klebo_amrfinder_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("Genome\tGene_Symbol\tSequence_Name\tClass\tSubclass\tCoverage\tIdentity\tScope\tElement_Type\tAccession\tContig\tStart\tStop\n")
            for genome_name, result in all_results.items():
                for hit in result['hits']:
                    row = [
                        genome_name,
                        hit.get('gene_symbol', ''),
                        hit.get('sequence_name', ''),
                        hit.get('class', ''),
                        hit.get('subclass', ''),
                        hit.get('coverage', ''),
                        hit.get('identity', ''),
                        hit.get('scope', ''),
                        hit.get('element_type', ''),
                        hit.get('accession', ''),
                        hit.get('contig_id', ''),
                        hit.get('start', ''),
                        hit.get('stop', '')
                    ]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"✓ Created summary: {summary_file}")

        stats_file = os.path.join(output_base, "klebo_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\tResistance_Classes\tGene_List\n")
            for genome_name, result in all_results.items():
                genes = list(set(hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol')))
                gene_list = ",".join(genes)
                high_risk_count = sum(1 for g in genes if g in self.high_risk_genes)
                critical_risk_count = sum(1 for g in genes if g in self.critical_risk_genes)
                classes = list(set(hit.get('class', '') for hit in result['hits'] if hit.get('class')))
                class_list = ",".join(classes)
                f.write(f"{genome_name}\t{result['hit_count']}\t{high_risk_count}\t{critical_risk_count}\t{class_list}\t{gene_list}\n")
        self.logger.info(f"✓ Created statistics: {stats_file}")

        master_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results)
            },
            'genome_summaries': {},
            'cross_genome_patterns': {}
        }
        all_genes_by_gene = defaultdict(lambda: {'count': 0, 'genomes': set()})
        for genome_name, result in all_results.items():
            genes = [hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol')]
            unique_genes = set(genes)
            critical_genes = [g for g in unique_genes if g in self.critical_risk_genes]
            high_risk_genes = [g for g in unique_genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            master_summary['genome_summaries'][genome_name] = {
                'total_hits': result['hit_count'],
                'unique_genes': len(unique_genes),
                'critical_genes': critical_genes,
                'high_risk_genes': high_risk_genes,
                'genes': list(unique_genes),
                'status': result['status']
            }
            for g in unique_genes:
                all_genes_by_gene[g]['count'] += 1
                all_genes_by_gene[g]['genomes'].add(genome_name)
        cross = {}
        for g, data in all_genes_by_gene.items():
            cross[g] = {
                'frequency': data['count'],
                'genomes': list(data['genomes']),
                'risk_level': 'CRITICAL' if g in self.critical_risk_genes else 'HIGH' if g in self.high_risk_genes else 'STANDARD'
            }
        master_summary['cross_genome_patterns'] = {
            'total_unique_genes': len(all_genes_by_gene),
            'genomes_with_critical': sum(1 for gs in master_summary['genome_summaries'].values() if gs['critical_genes']),
            'gene_frequency': cross
        }
        master_json = os.path.join(output_base, "klebo_amrfinder_master_summary.json")
        with open(master_json, 'w') as f:
            json.dump(master_summary, f, indent=2)
        self.logger.info(f"✓ Created master JSON: {master_json}")

        self._create_summary_html_report(all_results, output_base)

        # Mutation summary if any mutations reported
        self._create_mutation_summary(all_results, output_base)

    def _create_mutation_summary(self, all_results: Dict[str, Any], output_base: str):
        all_mutations = []
        genome_mutation_counts = {}
        for genome_name, result in all_results.items():
            mut_file = result.get('mutations_file')
            if mut_file and os.path.exists(mut_file):
                muts = self._parse_mutations_file(mut_file)
                if muts:
                    genome_mutation_counts[genome_name] = len(muts)
                    for m in muts:
                        m_copy = m.copy()
                        m_copy['genome'] = genome_name
                        all_mutations.append(m_copy)
                else:
                    genome_mutation_counts[genome_name] = 0
            else:
                genome_mutation_counts[genome_name] = 0

        if not all_mutations:
            self.logger.info("No mutations found in any genome; skipping mutation summaries.")
            return

        tsv_file = os.path.join(output_base, "mutation_summary.tsv")
        with open(tsv_file, 'w') as f:
            fieldnames = ['genome', 'gene_symbol', 'element_name', 'class', 'subclass',
                          'contig_id', 'start', 'stop', 'strand', 'coverage', 'identity', 'accession']
            f.write('\t'.join(fieldnames) + '\n')
            for m in all_mutations:
                row = [m.get('genome', ''),
                       m.get('gene_symbol', ''),
                       m.get('sequence_name', ''),
                       m.get('class', ''),
                       m.get('subclass', ''),
                       m.get('contig_id', ''),
                       m.get('start', ''),
                       m.get('stop', ''),
                       m.get('strand', ''),
                       m.get('coverage', ''),
                       m.get('identity', ''),
                       m.get('accession', '')]
                f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"✓ Mutation TSV: {tsv_file}")

        # Mutation HTML summary
        self._create_mutation_summary_html(all_mutations, genome_mutation_counts, output_base)

        # Mutation JSON
        self._create_mutation_json_summaries(all_mutations, genome_mutation_counts, output_base)

    def _create_mutation_summary_html(self, all_mutations: List[Dict], genome_counts: Dict[str, int], output_base: str):
        gene_freq = {}
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mutation = m.get('sequence_name', '')
            key = f"{gene}_{mutation}" if mutation else gene
            if key not in gene_freq:
                gene_freq[key] = {'count': 0, 'genomes': set(), 'gene': gene, 'mutation': mutation,
                                  'class': m.get('class',''), 'subclass': m.get('subclass','')}
            gene_freq[key]['count'] += 1
            gene_freq[key]['genomes'].add(m.get('genome',''))
        for k in gene_freq:
            gene_freq[k]['genomes'] = ', '.join(sorted(gene_freq[k]['genomes']))
        sorted_freq = sorted(gene_freq.values(), key=lambda x: x['count'], reverse=True)

        random_quote = random.choice(self.science_quotes)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Kleboscope - Mutation Summary</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', sans-serif;
            color: #fff;
            padding:20px;
        }}
        .container {{ max-width:1400px; margin:auto; }}
        .header {{ text-align:center; margin-bottom:30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding:20px;
            border-radius:15px;
            margin-bottom:20px;
            border:2px solid rgba(0,119,190,0.3);
        }}
        .ascii-art {{
            font-family: monospace;
            font-size:12px;
            white-space:pre;
            color:#0077be;
            text-shadow:0 0 10px rgba(0,119,190,0.5);
            overflow-x:auto;
        }}
        .card {{
            background: rgba(255,255,255,0.95);
            color:#1f2937;
            padding:25px;
            margin:20px 0;
            border-radius:12px;
        }}
        .card h2 {{
            color:#1e3a8a;
            border-bottom:3px solid #3b82f6;
            padding-bottom:10px;
            margin-bottom:20px;
        }}
        .gene-table {{
            width:100%;
            border-collapse:collapse;
            margin:20px 0;
            background:white;
        }}
        .gene-table th, .gene-table td {{
            padding:12px;
            text-align:left;
            border-bottom:1px solid #e5e7eb;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color:white;
        }}
        .gene-table tr:nth-child(even) {{ background-color:#f8fafc; }}
        .gene-table tr:hover {{ background-color:#e0f2fe; }}
        .footer {{
            background: rgba(0,0,0,0.8);
            color:white;
            padding:20px;
            border-radius:12px;
            margin-top:40px;
            text-align:center;
        }}
        .timestamp {{ color:#fbbf24; }}
        .authorship {{ margin-top:15px; font-size:0.9em; }}
        .quote-container {{
            background: rgba(255,255,255,0.1);
            color:white;
            padding:20px;
            border-radius:12px;
            margin:20px 0;
            text-align:center;
            font-style:italic;
            border-left:4px solid #fff;
        }}
        .sequence-cell {{ word-wrap:break-word; max-width:400px; min-width:200px; }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="card">
            <h1>🧬 Kleboscope – Mutation Summary Across All Genomes</h1>
            <p>Total genomes with mutations: {len([c for c in genome_counts.values() if c > 0])} / {len(genome_counts)}<br>
            Total mutation events: {len(all_mutations)}</p>
            <p><strong>Date:</strong> {current_time}</p>
        </div>
    </div>
    <div class="quote-container">
        <div style="font-size:1.1em;">"{random_quote}"</div>
    </div>
    <div class="card">
        <h2>📊 Mutation Frequency by Gene/Mutation</h2>
        <table class="gene-table">
            <thead><tr><th>Gene</th><th>Mutation</th><th>Count</th><th>Genomes</th><th>Class</th><th>Subclass</th></tr></thead>
            <tbody>
"""
        for item in sorted_freq:
            html += f"""
                <tr>
                    <td><strong>{item['gene']}</strong></td>
                    <td>{item['mutation']}</td>
                    <td>{item['count']}</td>
                    <td class="sequence-cell">{item['genomes']}</td>
                    <td>{item['class']}</td>
                    <td>{item['subclass']}</td>
                </tr>
"""
        html += f"""
            </tbody>
        </table>
    </div>
    <div class="footer">
        <p><strong>KLEBOSCOPE</strong> – Mutation Batch Summary</p>
        <p class="timestamp">Generated: {current_time}</p>
        <div class="authorship">
            <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>University of Ghana Medical School – Department of Medical Biochemistry</p>
        </div>
    </div>
</div>
</body>
</html>"""
        out_file = os.path.join(output_base, "mutation_summary.html")
        with open(out_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Mutation HTML summary: {out_file}")

    def _create_mutation_json_summaries(self, all_mutations: List[Dict], genome_counts: Dict[str, int], output_base: str):
        genome_summary = {genome: {'total_mutations': count} for genome, count in genome_counts.items()}
        gene_mutation_map = defaultdict(lambda: {'count': 0, 'genomes': set(), 'details': []})
        for m in all_mutations:
            gene = m.get('gene_symbol', 'unknown')
            mut_name = m.get('sequence_name', '')
            key = f"{gene}_{mut_name}"
            gene_mutation_map[key]['count'] += 1
            gene_mutation_map[key]['genomes'].add(m.get('genome',''))
            gene_mutation_map[key]['details'].append({
                'genome': m.get('genome'),
                'gene': gene,
                'mutation': mut_name,
                'class': m.get('class'),
                'subclass': m.get('subclass'),
                'contig': m.get('contig_id'),
                'start': m.get('start'),
                'stop': m.get('stop')
            })
        for v in gene_mutation_map.values():
            v['genomes'] = list(v['genomes'])
        master_json = {
            'metadata': {
                'tool': self.metadata['tool_name'] + ' Mutation Module',
                'version': self.metadata['version'],
                'analysis_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'total_genomes_analyzed': len(genome_counts),
                'total_mutations_detected': len(all_mutations)
            },
            'genome_summary': genome_summary,
            'mutation_frequency': {k: {'count': v['count'], 'genomes': v['genomes']} for k, v in gene_mutation_map.items()},
            'all_mutations': all_mutations
        }
        json_file = os.path.join(output_base, "mutation_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master_json, f, indent=2, default=str)
        self.logger.info(f"✓ Mutation master JSON: {json_file}")

    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        total_genomes = len(all_results)
        total_hits = sum(r['hit_count'] for r in all_results.values())
        genes_per_genome = {}
        gene_frequency = defaultdict(set)
        for gn, res in all_results.items():
            genes = [h.get('gene_symbol', '') for h in res['hits'] if h.get('gene_symbol')]
            genes_per_genome[gn] = genes
            for g in genes:
                gene_frequency[g].add(gn)

        critical_genes_found = set()
        high_risk_genes_found = set()
        genomes_with_critical = 0
        for gn, res in all_results.items():
            genes = set(genes_per_genome[gn])
            if any(g in self.critical_risk_genes for g in genes):
                genomes_with_critical += 1
                critical_genes_found.update(genes.intersection(self.critical_risk_genes))
            high_risk_genes_found.update(genes.intersection(self.high_risk_genes) - self.critical_risk_genes)

        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Kleboscope AMRfinderPlus - Summary Report</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', sans-serif;
            color: #fff;
            padding:20px;
        }}
        .container {{ max-width:1400px; margin:auto; }}
        .header {{ text-align:center; margin-bottom:30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding:20px;
            border-radius:15px;
            margin-bottom:20px;
            border:2px solid rgba(0,119,190,0.3);
        }}
        .ascii-art {{
            font-family: monospace;
            font-size:12px;
            white-space:pre;
            color:#0077be;
            text-shadow:0 0 10px rgba(0,119,190,0.5);
            overflow-x:auto;
        }}
        .card {{
            background: rgba(255,255,255,0.95);
            color:#1f2937;
            padding:25px;
            margin:20px 0;
            border-radius:12px;
        }}
        .card h2 {{
            color:#1e3a8a;
            border-bottom:3px solid #3b82f6;
            padding-bottom:10px;
            margin-bottom:20px;
        }}
        .summary-stats {{
            display:flex;
            gap:20px;
            justify-content:center;
            flex-wrap:wrap;
            margin:20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color:white;
            padding:20px;
            border-radius:8px;
            text-align:center;
            flex:1;
            min-width:200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%);
        }}
        .gene-table {{
            width:100%;
            border-collapse:collapse;
            margin-top:20px;
            background:white;
        }}
        .gene-table th, .gene-table td {{
            padding:12px;
            text-align:left;
            border-bottom:1px solid #e5e7eb;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color:white;
        }}
        .gene-table tr:nth-child(even) {{ background-color:#f8fafc; }}
        .gene-table tr:hover {{ background-color:#e0f2fe; }}
        .critical {{ background-color:#fee2e2; font-weight:bold; }}
        .high-risk {{ background-color:#fef3c7; }}
        .risk-badge {{
            display:inline-block;
            background:#dc2626;
            color:white;
            padding:4px 8px;
            border-radius:12px;
            margin:2px;
            font-size:0.8em;
        }}
        .warning-badge {{
            background:#f59e0b;
            color:black;
        }}
        .footer {{
            background: rgba(0,0,0,0.8);
            color:white;
            padding:20px;
            border-radius:12px;
            margin-top:40px;
            text-align:center;
        }}
        .timestamp {{ color:#fbbf24; }}
        .authorship {{ margin-top:15px; font-size:0.9em; }}
        .table-responsive {{ overflow-x:auto; }}
        .sequence-cell {{ word-wrap:break-word; max-width:400px; min-width:200px; }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="card">
            <h1>🧬 Kleboscope AMRfinderPlus – Batch Summary</h1>
            <p>K. pneumoniae Antimicrobial Resistance Analysis</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']} | <strong>Tool:</strong> {self.metadata['tool_name']} v{self.metadata['version']}</p>
        </div>
    </div>

    <div class="card">
        <h2>📊 Overall Summary</h2>
        <div class="summary-stats">
            <div class="stat-card"><div style="font-size:1.2em;">Total Genomes</div><div style="font-size:2em;">{total_genomes}</div></div>
            <div class="stat-card"><div style="font-size:1.2em;">Total AMR Genes</div><div style="font-size:2em;">{total_hits}</div></div>
            <div class="stat-card critical-stat-card"><div style="font-size:1.2em;">Genomes with Critical Risk</div><div style="font-size:2em;">{genomes_with_critical}</div></div>
        </div>
    </div>
"""
        if critical_genes_found:
            html += f"""
    <div class="card" style="border-left:4px solid #dc2626;">
        <h2 style="color:#dc2626;">🚨 CRITICAL RISK GENES ACROSS ALL GENOMES</h2>
        <div>
            <p><strong>{len(critical_genes_found)} unique critical risk genes found in {genomes_with_critical} genomes:</strong></p>
            <div>"""
            for g in sorted(critical_genes_found):
                html += f'<span class="risk-badge">🚨 {g}</span>'
            html += "</div></div></div>"
        if high_risk_genes_found and not critical_genes_found:
            html += f"""
    <div class="card" style="border-left:4px solid #f59e0b;">
        <h2 style="color:#f59e0b;">⚠️ High‑Risk Genes Detected</h2>
        <div>
            <p><strong>{len(high_risk_genes_found)} unique high‑risk genes found across {genomes_with_critical} genomes:</strong></p>
            <div>"""
            for g in sorted(high_risk_genes_found):
                html += f'<span class="risk-badge warning-badge">{g}</span>'
            html += "</div></div></div>"

        html += """
    <div class="card">
        <h2>🔍 Genes by Genome</h2>
        <div class="table-responsive">
            <table class="gene-table">
                <thead><tr><th>Genome</th><th>Gene Count</th><th>Genes Detected</th></tr></thead>
                <tbody>
"""
        for gn in sorted(genes_per_genome.keys()):
            genes = genes_per_genome[gn]
            row_class = "critical" if any(g in self.critical_risk_genes for g in genes) else "high-risk" if any(g in self.high_risk_genes for g in genes) else ""
            html += f"<tr class='{row_class}'> <td><strong>{gn}</strong></td> <td>{len(genes)}</td> <td>{', '.join(sorted(genes))}</td> </tr>"
        html += """
                </tbody>
            </table>
        </div>
    </div>

    <div class="card">
        <h2>📈 Gene Frequency</h2>
        <div class="table-responsive">
            <table class="gene-table">
                <thead> <tr><th>Gene</th><th>Frequency</th><th>Prevalence</th><th>Risk Level</th><th>Genomes</th> </tr> </thead>
                <tbody>
"""
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            freq = len(genomes)
            pct = (freq / total_genomes) * 100 if total_genomes else 0
            if gene in self.critical_risk_genes:
                risk = '<span class="risk-badge">CRITICAL</span>'
                row_class = "critical"
            elif gene in self.high_risk_genes:
                risk = '<span class="risk-badge warning-badge">HIGH</span>'
                row_class = "high-risk"
            else:
                risk = '<span class="risk-badge">Standard</span>'
                row_class = ""
            if pct >= 75:
                prev = '<span class="risk-badge">Very High</span>'
            elif pct >= 50:
                prev = '<span class="risk-badge warning-badge">High</span>'
            elif pct >= 25:
                prev = '<span class="risk-badge warning-badge">Medium</span>'
            elif pct >= 10:
                prev = '<span class="risk-badge">Low</span>'
            else:
                prev = '<span class="risk-badge">Rare</span>'
            html += f"<tr class='{row_class}'> <td><strong>{gene}</strong></td> <td>{freq} ({pct:.1f}%)</td> <td>{prev}</td> <td>{risk}</td> <td class='sequence-cell'>{', '.join(sorted(genomes))}</td> </tr>"
        html += f"""
                </tbody>
            </table>
        </div>
    </div>

    <div class="footer">
        <p><strong>KLEBOSCOPE</strong> – AMRfinderPlus Batch Analysis</p>
        <p class="timestamp">Generated: {self.metadata['analysis_date']}</p>
        <div class="authorship">
            <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>University of Ghana Medical School – Department of Medical Biochemistry</p>
        </div>
    </div>
</div>
</body>
</html>"""
        html_file = os.path.join(output_base, "klebo_amrfinder_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Created summary HTML report: {html_file}")

    def process_single_genome(self, genome_file: str, output_base: str = "klebo_amrfinder_results",
                              min_identity: Optional[float] = None,
                              min_coverage: Optional[float] = None,
                              report_mutations: bool = True) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        os.makedirs(results_dir, exist_ok=True)
        self.logger.info(f"=== PROCESSING GENOME: {genome_name} ===")
        result = self.run_amrfinder_single_genome(genome_file, results_dir,
                                                  min_identity=min_identity,
                                                  min_coverage=min_coverage,
                                                  report_mutations=report_mutations)
        status = "✓" if result['status'] == 'success' else "✗"
        self.logger.info(f"{status} {genome_name}: {result['hit_count']} AMR hits")
        return result

    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "klebo_amrfinder_results",
                                 min_identity: Optional[float] = None,
                                 min_coverage: Optional[float] = None,
                                 report_mutations: bool = True) -> Dict[str, Any]:
        if not self.check_amrfinder_installed():
            raise RuntimeError("Bundled AMRfinderPlus not properly installed")
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa",
                          f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        genome_files = []
        for pat in fasta_patterns:
            genome_files.extend(glob.glob(pat))
        genome_files = list(set(genome_files))
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching: {genome_pattern}")
        self.logger.info(f"Found {len(genome_files)} genomes")
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        max_concurrent = max(1, min(self.cpus, len(genome_files), int(self.available_ram / 2.5)))
        self.logger.info(f"🚀 MAXIMUM SPEED: Using {max_concurrent} concurrent jobs")
        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            future_to_genome = {
                executor.submit(self.process_single_genome, g, output_base,
                                min_identity, min_coverage, report_mutations): g
                for g in genome_files
            }
            for future in as_completed(future_to_genome):
                genome = future_to_genome[future]
                try:
                    result = future.result()
                    all_results[result['genome']] = result
                    self.logger.info(f"✓ COMPLETED: {result['genome']} ({result['hit_count']} hits)")
                except Exception as e:
                    self.logger.error(f"✗ FAILED: {genome} - {e}")
                    all_results[Path(genome).stem] = {'genome': Path(genome).stem, 'hits': [], 'hit_count': 0, 'status': 'failed'}
        self.create_amr_summary(all_results, output_base)
        self.logger.info("=== K. PNEUMONIAE AMR ANALYSIS COMPLETE ===")
        return all_results

def main():
    parser = argparse.ArgumentParser(
        description='Kleboscope AMRfinderPlus - K. pneumoniae Antimicrobial Resistance - MAXIMUM SPEED',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python klebo_amrfinder.py "*.fasta"
  python klebo_amrfinder.py "*.fna" --output results --cpus 8
  python klebo_amrfinder.py --update-db          # Incremental database update
  python klebo_amrfinder.py --force-update       # Full database overwrite
  python klebo_amrfinder.py --db-version         # Show current database version
  python klebo_amrfinder.py "*.fasta" --min-identity 0.95 --min-coverage 0.9
        """
    )
    parser.add_argument('pattern', nargs='?', help='File pattern for genomes (e.g., "*.fasta")')
    parser.add_argument('--cpus', '-c', type=int, default=None, help='CPU cores (auto-detect if not set)')
    parser.add_argument('--output', '-o', default='klebo_amrfinder_results', help='Output directory')
    parser.add_argument('--min-identity', type=float, default=None, help='Minimum identity (0..1)')
    parser.add_argument('--min-coverage', type=float, default=None, help='Minimum coverage (0..1)')
    parser.add_argument('--skip-mutations', action='store_true', help='Skip mutation reporting (mutations reported by default)')
    parser.add_argument('--update-db', action='store_true', help='Update AMRfinderPlus database incrementally and exit')
    parser.add_argument('--force-update', action='store_true', help='Force full database update (overwrites old) and exit')
    parser.add_argument('--db-version', action='store_true', help='Show current database version and exit')
    args = parser.parse_args()

    executor = KleboAMRfinderPlus(cpus=args.cpus)

    if args.db_version:
        print(f"Database version: {executor.metadata['database_version']}")
        print(f"Database path: {executor.bundled_database}")
        sys.exit(0)

    if args.update_db or args.force_update:
        success = executor.update_database(force=args.force_update)
        sys.exit(0 if success else 1)

    if not args.pattern:
        parser.error("Please provide a file pattern for genomes (or use --update-db / --force-update / --db-version)")

    try:
        results = executor.process_multiple_genomes(
            args.pattern,
            args.output,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            report_mutations=not args.skip_mutations
        )
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧫 Kleboscope AMRfinderPlus FINAL SUMMARY")
        executor.logger.info("="*50)
        total_hits = sum(r['hit_count'] for r in results.values())
        total_crit = 0
        total_high = 0
        for r in results.values():
            genes = [h.get('gene_symbol') for h in r['hits'] if h.get('gene_symbol')]
            total_crit += sum(1 for g in genes if g in executor.critical_risk_genes)
            total_high += sum(1 for g in genes if g in executor.high_risk_genes)
        executor.logger.info(f"Genomes processed: {len(results)}")
        executor.logger.info(f"Total AMR hits: {total_hits}")
        executor.logger.info(f"High‑risk genes detected: {total_high}")
        executor.logger.info(f"CRITICAL RISK genes detected: {total_crit}")
        executor.logger.info(f"Results saved to: {args.output}")
        executor.logger.info(f"CPU cores used: {executor.cpus}")
        executor.logger.info(f"Bundled AMRfinderPlus: {executor.metadata['amrfinder_version']}")
        executor.logger.info(f"Bundled database: {executor.metadata['database_version']}")
        if not args.skip_mutations:
            executor.logger.info("Mutation reporting: ENABLED")
        else:
            executor.logger.info("Mutation reporting: SKIPPED")
        if args.min_identity:
            executor.logger.info(f"Minimum identity: {args.min_identity}")
        if args.min_coverage:
            executor.logger.info(f"Minimum coverage: {args.min_coverage}")
    except Exception as e:
        executor.logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()