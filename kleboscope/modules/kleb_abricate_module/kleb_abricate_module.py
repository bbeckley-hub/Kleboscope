#!/usr/bin/env python3
"""
Kleboscope ABRicate Standalone Module
Comprehensive ABRicate analysis for Klebsiella pneumoniae with HTML, TSV, and JSON reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-05-24
Version: 1.0.0 (Klebsiella-optimized, complete)
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Set, Tuple
import argparse
import re
from datetime import datetime
import psutil
import json
import random
from collections import defaultdict, Counter

class AbricateExecutor:
    """ABRicate executor for Klebsiella with comprehensive reporting - MAXIMUM SPEED"""
    
    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # Databases to be used
        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'ecoh', 'ecoli_vf', 'bacmet2'
        ]
        
        # ==================== COMPREHENSIVE KLEBSIELLA GENE DICTIONARIES ====================
        # Critical resistance genes (carbapenemases, mcr, tet(X), 16S rRNA methyltransferases, etc.)
        self.critical_resistance_genes = {
            # Carbapenemases (class A, B, D)
            'blaKPC', 'blaNDM', 'blaIMP', 'blaVIM', 'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
            'blaOXA-244', 'blaOXA-162', 'blaOXA-163', 'blaOXA-247', 'blaGES', 'blaIMI', 'blaSME', 
            'blaNMC', 'blaCcrA', 'blaBIC', 'blaDIM', 'blaTMB', 'blaSPM', 'blaGIM', 'blaSIM', 
            'blaAIM', 'blaFRI', 'blaAFM', 'blaFLC', 'blaFMC', 'blaFIM', 'blaFPH', 'blaFTU', 
            # Colistin resistance (plasmid-mediated)
            'mcr-1', 'mcr-2', 'mcr-3', 'mcr-3-like', 'mcr-4', 'mcr-5', 'mcr-6', 'mcr-7', 
            'mcr-8', 'mcr-9', 'mcr-10', 'mcr-3.1', 'mcr-3.2', 'mcr-3.3', 'mcr-3.4', 'mcr-3.5',
            # Tigecycline resistance
            'tet(X)', 'tet(X1)', 'tet(X2)', 'tet(X3)', 'tet(X4)', 'tet(X5)', 'tet(X6)',
            'tet(X7)', 'tet(X8)', 'tet(X9)', 'tet(X10)', 'tet(X11)', 'tet(X12)', 'tet(X13)', 'tet(X14)',
            # 16S rRNA methyltransferases (high-level aminoglycoside resistance)
            'armA', 'rmtA', 'rmtB', 'rmtC', 'rmtD', 'rmtE', 'rmtF', 'rmtG', 'rmtH', 'npmA',
            # Other critical (plasmid-mediated fosfomycin resistance)
            'fosA3', 'fosA4', 'fosA5', 'fosA6', 'fosA7', 'fosC2', 'fosK', 'fosX', 'fosL', 'fosB',
            'fosB1', 'fosB2', 'fosB3', 'fosC', 'fosE', 'fosF', 'fosG', 'fosH', 'fosI', 'fosJ',
            # OXA-type carbapenemases (others)
            'blaOXA-23', 'blaOXA-24', 'blaOXA-40', 'blaOXA-51', 'blaOXA-58', 'blaOXA-143', 'blaOXA-235',
            # Rare carbapenemases
            'blaFRI-1', 'blaFRI-2', 'blaFRI-3', 'blaFRI-4', 'blaFRI-5', 'blaFRI-6', 'blaFRI-7',
            'blaFRI-8', 'blaFRI-9', 'blaFRI-10', 'blaFRI-11', 'blaFRI-12', 'blaFRI-13', 'blaFRI-14',
            'blaIMI-2', 'blaIMI-3', 'blaIMI-4', 'blaIMI-5', 'blaIMI-6', 'blaIMI-7', 'blaIMI-8',
            'blaSME-2', 'blaSME-3', 'blaSME-4', 'blaSME-5',
            # Beta-lactamases with carbapenemase activity (GES variants)
            'blaGES-2', 'blaGES-4', 'blaGES-5', 'blaGES-6', 'blaGES-7', 'blaGES-8', 'blaGES-9',
            'blaGES-11', 'blaGES-14', 'blaGES-15', 'blaGES-16', 'blaGES-17', 'blaGES-18', 'blaGES-19',
            'blaGES-20', 'blaGES-21', 'blaGES-22', 'blaGES-23', 'blaGES-24', 'blaGES-25', 'blaGES-26',
            'blaGES-27', 'blaGES-28'
        }
        
        # High-risk virulence genes (hypervirulence-associated and other key virulence factors)
        self.high_risk_virulence_genes = {
            # Regulators of hypermucoviscosity
            'rmpA', 'rmpA2', 'rmpC', 'rmpD', 'rmp', 'rmpA1', 'rmpA_2',
            # Aerobactin siderophore system
            'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
            # Salmochelin siderophore system
            'iroB', 'iroC', 'iroD', 'iroE', 'iroN',
            # Colibactin genotoxin cluster
            'clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH',
            'clbI', 'clbJ', 'clbK', 'clbL', 'clbM', 'clbN', 'clbO', 'clbP', 'clbQ', 'clbR', 'clbS',
            'clb1', 'clb2', 'clb3', 'clb4', 'clb5', 'clb6', 'clb7', 'clb8', 'clb9', 'clb10',
            # Additional markers of hypervirulent clones
            'peg-344', 'allS', 'kfuA', 'kfuB', 'kfuC', 'mrkD', 'mrkD', 'mrkABC', 'mrkA', 'mrkB', 'mrkC', 
            'mrkF', 'mrkI', 'wcaG', 'uge', 'wabG', 'ureA', 'ureB', 'ureC', 'ureD', 'ureE', 'ureF', 'ureG',
            'fimA', 'fimH', 'fimH', 'fimC', 'fimD', 'fimF', 'fimG', 'fimI',
            # Iron uptake systems (yersiniabactin, enterobactin, etc.)
            'ybtA', 'ybtE', 'ybtP', 'ybtQ', 'ybtS', 'ybtT', 'ybtU', 'ybtX', 'irp1', 'irp2',
            'entA', 'entB', 'entC', 'entD', 'entE', 'entF', 'entS', 'fepA', 'fepB', 'fepC', 'fepD', 'fepG',
            'fes', 'fhuA', 'fhuB', 'fhuC', 'fhuD', 'fhuE', 'fhuF',
            # LPS and capsule synthesis
            'wzi', 'wzy', 'wzc', 'wzx', 'wbaP', 'gnd', 'galF', 'manB', 'ugd',
            'wbbM', 'wbbN', 'wbbO', 'rfb', 'rfc', 'rmlA', 'rmlB', 'rmlC', 'rmlD',
            # Other virulence factors
            'hly', 'hlyA', 'hlyB', 'hlyC', 'hlyD', 'cnf1', 'cnf2', 'cdtB', 'sat', 'pic', 'vat',
            'tia', 'toxA', 'traT', 'ompA', 'ompC', 'ompF', 'ompX', 'ompT', 'ompP',
            'cvaA', 'cvaB', 'cvaC', 'cvi', 'cma', 'cmaA', 'cmaB', 'cmaC',
            'papA', 'papB', 'papC', 'papD', 'papE', 'papF', 'papG', 'papH', 'papI', 'papJ', 'papK',
            'sfaA', 'sfaB', 'sfaC', 'sfaD', 'sfaE', 'sfaF', 'sfaG', 'sfaH',
            'afaA', 'afaB', 'afaC', 'afaD', 'afaE', 'afaF', 'afaG',
            'focA', 'focB', 'focC', 'focD', 'focE', 'focF', 'focG', 'focH', 'focI', 'focJ',
            'hra', 'hra1', 'hra2', 'ibcA', 'ibcB', 'ibcC',
            'kpsE', 'kpsM', 'kpsT', 'kpsD', 'kpsF', 'kpsS', 'neuA', 'neuB', 'neuC', 'neuD', 'neuE'
        }
        
        # Beta-lactam resistance genes (ESBLs, AmpC, narrow-spectrum) - not critical but important
        self.beta_lactamase_genes = {
            'blaTEM', 'blaSHV', 'blaCTX-M', 'blaOXA-1', 'blaOXA-2', 'blaOXA-9', 'blaOXA-10',
            'blaOXA-30', 'blaOXA-32', 'blaCMY', 'blaDHA', 'blaACC', 'blaMIR', 'blaACT', 'blaFOX',
            'blaMOX', 'blaCIT', 'blaEBC', 'blaCMH', 'blaBIL', 'blaCepA', 'blaCblA', 'blaCKO',
            'blaBEL', 'blaPER', 'blaVEB', 'blaGES-1', 'blaGES-3', 'blaIBC', 'blaSFO', 'blaTLA',
            'blaZ', 'blaC', 'blaB', 'blaA', 'blaP', 'blaO', 'blaU', 'blaQ', 'blaR', 'blaS',
            'blaCTX-M-1', 'blaCTX-M-2', 'blaCTX-M-3', 'blaCTX-M-4', 'blaCTX-M-5', 'blaCTX-M-6',
            'blaCTX-M-7', 'blaCTX-M-8', 'blaCTX-M-9', 'blaCTX-M-10', 'blaCTX-M-11', 'blaCTX-M-12',
            'blaCTX-M-13', 'blaCTX-M-14', 'blaCTX-M-15', 'blaCTX-M-16', 'blaCTX-M-17', 'blaCTX-M-18',
            'blaCTX-M-19', 'blaCTX-M-20', 'blaCTX-M-21', 'blaCTX-M-22', 'blaCTX-M-23', 'blaCTX-M-24',
            'blaCTX-M-25', 'blaCTX-M-26', 'blaCTX-M-27', 'blaCTX-M-28', 'blaCTX-M-29', 'blaCTX-M-30',
            'blaCTX-M-31', 'blaCTX-M-32', 'blaCTX-M-33', 'blaCTX-M-34', 'blaCTX-M-35', 'blaCTX-M-36',
            'blaCTX-M-37', 'blaCTX-M-38', 'blaCTX-M-39', 'blaCTX-M-40', 'blaCTX-M-41', 'blaCTX-M-42',
            'blaCTX-M-43', 'blaCTX-M-44', 'blaCTX-M-45', 'blaCTX-M-46', 'blaCTX-M-47', 'blaCTX-M-48',
            'blaCTX-M-49', 'blaCTX-M-50', 'blaCTX-M-51', 'blaCTX-M-52', 'blaCTX-M-53', 'blaCTX-M-54',
            'blaCTX-M-55', 'blaCTX-M-56', 'blaCTX-M-57', 'blaCTX-M-58', 'blaCTX-M-59', 'blaCTX-M-60',
            'blaCTX-M-61', 'blaCTX-M-62', 'blaCTX-M-63', 'blaCTX-M-64', 'blaCTX-M-65', 'blaCTX-M-66',
            'blaCTX-M-67', 'blaCTX-M-68', 'blaCTX-M-69', 'blaCTX-M-70', 'blaCTX-M-71', 'blaCTX-M-72',
            'blaCTX-M-73', 'blaCTX-M-74', 'blaCTX-M-75', 'blaCTX-M-76', 'blaCTX-M-77', 'blaCTX-M-78',
            'blaCTX-M-79', 'blaCTX-M-80', 'blaCTX-M-81', 'blaCTX-M-82', 'blaCTX-M-83', 'blaCTX-M-84',
            'blaCTX-M-85', 'blaCTX-M-86', 'blaCTX-M-87', 'blaCTX-M-88', 'blaCTX-M-89', 'blaCTX-M-90',
            'blaCTX-M-91', 'blaCTX-M-92', 'blaCTX-M-93', 'blaCTX-M-94', 'blaCTX-M-95', 'blaCTX-M-96',
            'blaCTX-M-97', 'blaCTX-M-98', 'blaCTX-M-99', 'blaCTX-M-100', 'blaCTX-M-101', 'blaCTX-M-102',
            'blaCTX-M-103', 'blaCTX-M-104', 'blaCTX-M-105', 'blaCTX-M-106', 'blaCTX-M-107', 'blaCTX-M-108',
            'blaCTX-M-109', 'blaCTX-M-110', 'blaCTX-M-111', 'blaCTX-M-112', 'blaCTX-M-113', 'blaCTX-M-114',
            'blaCTX-M-115', 'blaCTX-M-116', 'blaCTX-M-117', 'blaCTX-M-118', 'blaCTX-M-119', 'blaCTX-M-120',
            'blaCTX-M-121', 'blaCTX-M-122', 'blaCTX-M-123', 'blaCTX-M-124', 'blaCTX-M-125', 'blaCTX-M-126',
            'blaCTX-M-127', 'blaCTX-M-128', 'blaCTX-M-129', 'blaCTX-M-130', 'blaCTX-M-131', 'blaCTX-M-132',
            'blaCTX-M-133', 'blaCTX-M-134', 'blaCTX-M-135', 'blaCTX-M-136', 'blaCTX-M-137', 'blaCTX-M-138',
            'blaCTX-M-139', 'blaCTX-M-140', 'blaCTX-M-141', 'blaCTX-M-142', 'blaCTX-M-143', 'blaCTX-M-144',
            'blaCTX-M-145', 'blaCTX-M-146', 'blaCTX-M-147', 'blaCTX-M-148', 'blaCTX-M-149', 'blaCTX-M-150',
            'blaCTX-M-151', 'blaCTX-M-152', 'blaCTX-M-153', 'blaCTX-M-154', 'blaCTX-M-155', 'blaCTX-M-156',
            'blaCTX-M-157', 'blaCTX-M-158', 'blaCTX-M-159', 'blaCTX-M-160', 'blaCTX-M-161', 'blaCTX-M-162',
            'blaCTX-M-163', 'blaCTX-M-164', 'blaCTX-M-165', 'blaCTX-M-166', 'blaCTX-M-167', 'blaCTX-M-168',
            'blaCTX-M-169', 'blaCTX-M-170', 'blaCTX-M-171', 'blaCTX-M-172', 'blaCTX-M-173', 'blaCTX-M-174',
            'blaCTX-M-175', 'blaCTX-M-176', 'blaCTX-M-177', 'blaCTX-M-178', 'blaCTX-M-179', 'blaCTX-M-180',
            'blaCTX-M-181', 'blaCTX-M-182', 'blaCTX-M-183', 'blaCTX-M-184', 'blaCTX-M-185', 'blaCTX-M-186',
            'blaCTX-M-187', 'blaCTX-M-188', 'blaCTX-M-189', 'blaCTX-M-190', 'blaCTX-M-191', 'blaCTX-M-192',
            'blaCTX-M-193', 'blaCTX-M-194', 'blaCTX-M-195', 'blaCTX-M-196', 'blaCTX-M-197', 'blaCTX-M-198',
            'blaCTX-M-199', 'blaCTX-M-200', 'blaCTX-M-201', 'blaCTX-M-202', 'blaCTX-M-203', 'blaCTX-M-204',
            'blaCTX-M-205', 'blaCTX-M-206', 'blaCTX-M-207', 'blaCTX-M-208', 'blaCTX-M-209', 'blaCTX-M-210',
            'blaCTX-M-211', 'blaCTX-M-212', 'blaCTX-M-213', 'blaCTX-M-214', 'blaCTX-M-215', 'blaCTX-M-216',
            'blaCTX-M-217', 'blaCTX-M-218', 'blaCTX-M-219', 'blaCTX-M-220', 'blaCTX-M-221', 'blaCTX-M-222',
            'blaCTX-M-223', 'blaCTX-M-224', 'blaCTX-M-225', 'blaCTX-M-226', 'blaCTX-M-227', 'blaCTX-M-228',
            'blaCTX-M-229', 'blaCTX-M-230', 'blaCTX-M-231', 'blaCTX-M-232', 'blaCTX-M-233', 'blaCTX-M-234',
            'blaCTX-M-235', 'blaCTX-M-236', 'blaCTX-M-237', 'blaCTX-M-238', 'blaCTX-M-239', 'blaCTX-M-240',
            'blaCTX-M-241', 'blaCTX-M-242', 'blaCTX-M-243', 'blaCTX-M-244', 'blaCTX-M-245', 'blaCTX-M-246',
            'blaCTX-M-247', 'blaCTX-M-248', 'blaCTX-M-249', 'blaCTX-M-250', 'blaCTX-M-251', 'blaCTX-M-252',
            'blaCTX-M-253', 'blaCTX-M-254', 'blaCTX-M-255', 'blaCTX-M-256', 'blaCTX-M-257', 'blaCTX-M-258',
            'blaCTX-M-259', 'blaCTX-M-260', 'blaCTX-M-261', 'blaCTX-M-262', 'blaCTX-M-263', 'blaCTX-M-264',
            'blaCTX-M-265', 'blaCTX-M-266', 'blaCTX-M-267', 'blaCTX-M-268', 'blaCTX-M-269', 'blaCTX-M-270',
            'blaCTX-M-271', 'blaCTX-M-272', 'blaCTX-M-273', 'blaCTX-M-274', 'blaCTX-M-275', 'blaCTX-M-276',
            'blaCTX-M-277', 'blaCTX-M-278', 'blaCTX-M-279', 'blaCTX-M-280', 'blaCTX-M-281', 'blaCTX-M-282',
            'blaCTX-M-283', 'blaCTX-M-284', 'blaCTX-M-285', 'blaCTX-M-286', 'blaCTX-M-287', 'blaCTX-M-288',
            'blaCTX-M-289', 'blaCTX-M-290', 'blaCTX-M-291', 'blaCTX-M-292', 'blaCTX-M-293', 'blaCTX-M-294',
            'blaCTX-M-295', 'blaCTX-M-296', 'blaCTX-M-297', 'blaCTX-M-298', 'blaCTX-M-299', 'blaCTX-M-300'
        }
        
        self.metadata = {
            "tool_name": "Kleboscope ABRicate",
            "version": "1.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        # Science quotes
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "Kleboscope turns genomic complexity into actionable insights for AMR surveillance.", "author": "Brown Beckley"}
        ]
    
    def get_random_quote(self):
        return random.choice(self.science_quotes)
    
    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> int:
        try:
            ram_gb = psutil.virtual_memory().available / (1024 ** 3)
            return ram_gb
        except Exception:
            return 8
    
    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 2)
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 3)
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
        except Exception:
            return os.cpu_count() or 4
    
    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus <= 4:
            self.logger.info("💡 Performance: Multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("💡 Performance: High-speed mode")
        else:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")
    
    def check_abricate_installed(self) -> bool:
        try:
            result = subprocess.run(['abricate', '--version'],
                                    capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match and version_match.group(1) >= "1.2.0":
                self.logger.info("✓ ABRicate version meets requirement (>=1.2.0)")
                return True
            self.logger.info("✓ ABRicate installed")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found. Please install with: conda install -c bioconda abricate")
            return False
    
    def setup_abricate_databases(self):
        self.logger.info("Setting up ABRicate databases...")
        available_dbs = []
        missing_dbs = []
        try:
            check_result = subprocess.run(['abricate', '--list'],
                                          capture_output=True, text=True, check=True)
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("✓ Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    subprocess.run(['abricate', '--setupdb', '--db', db],
                                   capture_output=True, text=True, check=True)
                    self.logger.info("✓ Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
        except Exception as e:
            self.logger.error("Error setting up databases: %s", e)
    
    def run_abricate_single_db(self, genome_file: str, database: str, output_dir: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"abricate_{database}.txt")
        cmd = [
            'abricate',
            genome_file,
            '--db', database,
            '--minid', '80',
            '--mincov', '80'
        ]
        self.logger.info("Running ABRicate: %s --db %s", genome_name, database)
        try:
            with open(output_file, 'w') as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            hits = self._parse_abricate_output(output_file)
            self._create_database_html_report(genome_name, database, hits, output_dir)
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s", database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
            if not lines:
                return hits
            headers = []
            data_lines = []
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            if not headers:
                return hits
            expected_columns = len(headers)
            for line_num, line in enumerate(data_lines, 1):
                parts = line.split('\t')
                if len(parts) > expected_columns:
                    combined_parts = parts[:expected_columns-1]
                    combined_parts.append('\t'.join(parts[expected_columns-1:]))
                    parts = combined_parts
                elif len(parts) < expected_columns:
                    parts.extend([''] * (expected_columns - len(parts)))
                if len(parts) == expected_columns:
                    hit = {}
                    for i, header in enumerate(headers):
                        hit[header] = parts[i] if i < len(parts) else ''
                    processed_hit = {
                        'file': hit.get('FILE', ''),
                        'sequence': hit.get('SEQUENCE', ''),
                        'start': hit.get('START', ''),
                        'end': hit.get('END', ''),
                        'strand': hit.get('STRAND', ''),
                        'gene': hit.get('GENE', ''),
                        'coverage': hit.get('COVERAGE', ''),
                        'coverage_map': hit.get('COVERAGE_MAP', ''),
                        'gaps': hit.get('GAPS', ''),
                        'coverage_percent': hit.get('%COVERAGE', ''),
                        'identity_percent': hit.get('%IDENTITY', ''),
                        'database': hit.get('DATABASE', ''),
                        'accession': hit.get('ACCESSION', ''),
                        'product': hit.get('PRODUCT', ''),
                        'resistance': hit.get('RESISTANCE', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d", line_num, len(parts), expected_columns)
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits
    
    def _create_database_html_report(self, genome_name: str, database: str, hits: List[Dict], output_dir: str):
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KLEBOSCOPE - ABRicate {database.upper()} Database Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 119, 190, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #0077be;
            text-shadow: 0 0 10px rgba(0, 119, 190, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
        }}
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }}
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        .summary-table td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        .summary-table tr:hover {{
            background-color: #e0f2fe;
        }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .high-risk {{ background-color: #fef3c7; }}
        .present {{ background-color: #d1fae5; }}
        .product-cell {{
            white-space: normal !important;
            word-wrap: break-word;
            max-width: 400px;
            min-width: 200px;
        }}
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 20px 0;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        .confidence-badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: bold;
            text-transform: uppercase;
        }}
        .confidence-high {{
            background-color: #16a34a;
            color: white;
        }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .metrics-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
                </div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Database Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Database Name</div>
                    <div class="metric-value">{database.upper()}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Genome</div>
                    <div class="metric-value">{genome_name}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Total Hits</div>
                    <div class="metric-value">{len(hits)}</div>
                </div>
            </div>
        </div>
"""
        if hits:
            html_content += """
        <div class="report-section">
            <h2>🔍 Genes Detected</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Product</th>
                            <th>Coverage</th>
                            <th>Identity</th>
                            <th>Accession</th>
                        </tr>
                    </thead>
                    <tbody>
"""
            for hit in hits:
                gene = hit['gene']
                row_class = "present"
                # Determine risk class based on dictionaries
                if any(crit in gene for crit in self.critical_resistance_genes):
                    row_class = "critical"
                elif any(vir in gene for vir in self.high_risk_virulence_genes):
                    row_class = "high-risk"
                html_content += f"""
                    <tr class="{row_class}">
                        <td><strong>{hit['gene']}</strong></td>
                        <td class="product-cell">{hit['product']}</td>
                        <td>{hit['coverage_percent']}%</td>
                        <td>{hit['identity_percent']}%</td>
                        <td>{hit['accession']}</td>
                    </tr>
"""
            html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""
        else:
            html_content += f"""
        <div class="report-section">
            <h2>✅ No Genes Detected</h2>
            <p>No significant hits found in the {database.upper()} database.</p>
        </div>
"""
        html_content += f"""
        <div class="footer">
            <p><strong>KLEBOSCOPE</strong> - ABRicate Analysis Module</p>
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
        const quotes = {json.dumps(self.science_quotes)};
        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');
        function getRandomQuote() {{ return quotes[Math.floor(Math.random() * quotes.length)]; }}
        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>"""
        html_file = os.path.join(output_dir, f"abricate_{database}_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Individual database report: %s", html_file)
    
    def analyze_klebsiella_genes(self, all_hits: List[Dict]) -> Dict[str, Any]:
        """Classify hits into critical resistance, high-risk virulence, and other."""
        analysis = {
            'critical_resistance_genes': [],
            'high_risk_virulence_genes': [],
            'beta_lactamase_genes': [],  # non-critical beta-lactamases
            'other_genes': [],
            'resistance_classes': {},
            'total_critical_resistance': 0,
            'total_high_risk_virulence': 0,
            'total_beta_lactamase': 0,
            'total_other': 0,
            'total_hits': len(all_hits)
        }
        for hit in all_hits:
            gene = hit['gene']
            # Check critical resistance
            if any(crit in gene for crit in self.critical_resistance_genes):
                analysis['critical_resistance_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'CRITICAL-RES'
                })
            # Check high-risk virulence
            elif any(vir in gene for vir in self.high_risk_virulence_genes):
                analysis['high_risk_virulence_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'HIGH-VIRULENCE'
                })
            # Check beta-lactamases (non-critical)
            elif any(bla in gene for bla in self.beta_lactamase_genes):
                analysis['beta_lactamase_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'BETA-LACTAMASE'
                })
            else:
                analysis['other_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent']
                })
            # Classify resistance class based on product
            res_class = self._classify_resistance(hit['product'])
            if res_class:
                if res_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][res_class] = []
                if gene not in [g['gene'] for g in analysis['resistance_classes'][res_class]]:
                    analysis['resistance_classes'][res_class].append({'gene': gene, 'product': hit['product']})
        analysis['total_critical_resistance'] = len(analysis['critical_resistance_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        analysis['total_beta_lactamase'] = len(analysis['beta_lactamase_genes'])
        analysis['total_other'] = len(analysis['other_genes'])
        return analysis
    
    def _classify_resistance(self, product: str) -> str:
        product_lower = product.lower()
        if any(term in product_lower for term in ['carbapenem', 'kpc', 'ndm', 'imp', 'vim', 'oxa']):
            return 'Carbapenem resistance'
        elif any(term in product_lower for term in ['esbl', 'ctx-m', 'shv', 'tem', 'cephem']):
            return 'ESBL/AmpC resistance'
        elif any(term in product_lower for term in ['colistin', 'mcr']):
            return 'Colistin resistance'
        elif any(term in product_lower for term in ['tigecycline', 'tet(x)']):
            return 'Tigecycline resistance'
        elif any(term in product_lower for term in ['aminoglycoside', 'aac', 'aad', 'aph', 'rmt', 'arm']):
            return 'Aminoglycoside resistance'
        elif any(term in product_lower for term in ['fluoroquinolone', 'qnr', 'gyr']):
            return 'Fluoroquinolone resistance'
        elif any(term in product_lower for term in ['macrolide', 'erm', 'mph']):
            return 'Macrolide resistance'
        elif any(term in product_lower for term in ['tetracycline', 'tet']):
            return 'Tetracycline resistance'
        elif any(term in product_lower for term in ['sulfonamide', 'sul']):
            return 'Sulfonamide resistance'
        elif any(term in product_lower for term in ['trimethoprim', 'dfr']):
            return 'Trimethoprim resistance'
        elif any(term in product_lower for term in ['phenicol', 'cat', 'cml', 'flo']):
            return 'Phenicol resistance'
        elif any(term in product_lower for term in ['fosfomycin', 'fos']):
            return 'Fosfomycin resistance'
        elif any(term in product_lower for term in ['efflux', 'multidrug']):
            return 'Efflux pumps'
        else:
            return 'Other resistance'
    
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])
        analysis = self.analyze_klebsiella_genes(all_hits)
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KLEBOSCOPE - ABRicate Analysis Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 119, 190, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #0077be;
            text-shadow: 0 0 10px rgba(0, 119, 190, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .metric-label {{ font-size: 14px; opacity: 0.9; margin-bottom: 5px; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .risk-badge {{
            display: inline-block;
            background: #dc2626;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .warning-badge {{
            display: inline-block;
            background: #f59e0b;
            color: black;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
        }}
        .summary-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #e0f2fe; }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .high-risk {{ background-color: #fef3c7; }}
        .present {{ background-color: #d1fae5; }}
        .product-cell {{ white-space: normal; word-wrap: break-word; max-width: 500px; }}
        .table-responsive {{ width: 100%; overflow-x: auto; margin: 20px 0; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .metrics-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
                </div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Klebsiella AMR/Virulence Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Total Genes</div>
                    <div class="metric-value">{analysis['total_hits']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Critical Resistance</div>
                    <div class="metric-value">{analysis['total_critical_resistance']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">High-Risk Virulence</div>
                    <div class="metric-value">{analysis['total_high_risk_virulence']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Beta-Lactamases</div>
                    <div class="metric-value">{analysis['total_beta_lactamase']}</div>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        </div>
"""
        # Critical resistance alerts
        if analysis['critical_resistance_genes']:
            html_content += """
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">⚠️ CRITICAL RESISTANCE GENES DETECTED</h2>
            <div style="margin: 10px 0;">
"""
            for gene_info in analysis['critical_resistance_genes']:
                html_content += f'<span class="risk-badge">{gene_info["gene"]}</span> '
            html_content += """
            </div>
        </div>
"""
        # High-risk virulence alerts
        if analysis['high_risk_virulence_genes']:
            html_content += """
        <div class="report-section" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #f59e0b;">🟡 HIGH-RISK VIRULENCE GENES DETECTED</h2>
            <div style="margin: 10px 0;">
"""
            for gene_info in analysis['high_risk_virulence_genes']:
                html_content += f'<span class="warning-badge">{gene_info["gene"]}</span> '
            html_content += """
            </div>
        </div>
"""
        # Resistance classes
        if analysis['resistance_classes']:
            html_content += """
        <div class="report-section">
            <h2>🧪 Resistance Classes Detected</h2>
"""
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join([g['gene'] for g in genes])
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 8px;">
                <strong style="color: #667eea;">{class_name}</strong> ({len(genes)} genes)
                <br><span style="color: #666; font-size: 0.9em;">{gene_list}</span>
            </div>
"""
            html_content += "</div>"
        
        # Critical resistance table
        if analysis['critical_resistance_genes']:
            html_content += """
        <div class="report-section">
            <h2>🔴 Critical Resistance Genes</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['critical_resistance_genes']:
                html_content += f"""
                    <tr class="critical">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            html_content += "</tbody></table></div></div>"
        
        # High-risk virulence table
        if analysis['high_risk_virulence_genes']:
            html_content += """
        <div class="report-section">
            <h2>🟡 High-Risk Virulence Genes</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['high_risk_virulence_genes']:
                html_content += f"""
                    <tr class="high-risk">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            html_content += "</tbody></table></div></div>"
        
        # Beta-lactamase table
        if analysis['beta_lactamase_genes']:
            html_content += """
        <div class="report-section">
            <h2>🔵 Beta-Lactamase Genes (non-critical)</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['beta_lactamase_genes'][:1000]:  # limit to 1000
                html_content += f"""
                    <tr class="present">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            if len(analysis['beta_lactamase_genes']) > 1000:
                html_content += f"<tr><td colspan='5'>... and {len(analysis['beta_lactamase_genes'])-1000} more</td></tr>"
            html_content += "</tbody></table></div></div>"
        
        # Other genes table (if any)
        if analysis['other_genes']:
            html_content += """
        <div class="report-section">
            <h2>🔵 Other Detected Genes</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['other_genes'][:5000]:  # limit to 5000
                html_content += f"""
                    <tr class="present">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            if len(analysis['other_genes']) > 5000:
                html_content += f"<tr><td colspan='5'>... and {len(analysis['other_genes'])-5000} more</td></tr>"
            html_content += "</tbody></table></div></div>"
        
        # Database summary
        html_content += """
        <div class="report-section">
            <h2>🗃️ Database Results Summary</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead><tr><th>Database</th><th>Hits</th><th>Status</th></tr></thead>
                    <tbody>
"""
        for db, result in results.items():
            status_icon = "✅" if result['status'] == 'success' else "❌"
            html_content += f"""
                    <tr><td>{db}</td><td>{result['hit_count']}</td><td>{status_icon} {result['status']}</td></tr>
"""
        html_content += f"""
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>KLEBOSCOPE</strong> - ABRicate Analysis Module</p>
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
        const quotes = {json.dumps(self.science_quotes)};
        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');
        function getRandomQuote() {{ return quotes[Math.floor(Math.random() * quotes.length)]; }}
        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>"""
        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Comprehensive Klebsiella HTML report generated: %s", html_file)
    
    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating database summary files and HTML reports...")
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = []
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db].append(hit_with_genome)
        for db, hits in db_results.items():
            if hits:
                summary_file = os.path.join(output_base, f"klebo_{db}_abricate_summary.tsv")
                headers = list(hits[0].keys())
                with open(summary_file, 'w') as f:
                    f.write('\t'.join(headers) + '\n')
                    for hit in hits:
                        row = [str(hit.get(header, '')) for header in headers]
                        f.write('\t'.join(row) + '\n')
                self.logger.info("✓ Created %s summary: %s (%d hits)", db, summary_file, len(hits))
                self._create_database_summary_html(db, hits, output_base)
            else:
                self.logger.info("No hits for database %s, skipping summary", db)
    
    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        """Create HTML summary for a specific database across all genomes."""
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        unique_genomes = list(set(hit['genome'] for hit in hits))
        genes_per_genome = {}
        for hit in hits:
            genome = hit['genome']
            if genome not in genes_per_genome:
                genes_per_genome[genome] = set()
            genes_per_genome[genome].add(hit['gene'])
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KLEBOSCOPE - Batch ABRicate Summary ({database.upper()})</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(0, 119, 190, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #0077be;
            text-shadow: 0 0 10px rgba(0, 119, 190, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
        }}
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }}
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        .summary-table td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{
            background-color: #f8fafc;
        }}
        .summary-table tr:hover {{
            background-color: #e0f2fe;
        }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .high-risk {{ background-color: #fef3c7; }}
        .present {{ background-color: #d1fae5; }}
        .product-cell {{
            white-space: normal !important;
            word-wrap: break-word;
            max-width: 400px;
            min-width: 200px;
        }}
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 20px 0;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .metrics-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
                </div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Batch ABRicate Summary - {database.upper()} Database</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{len(hits)}</div>
                    <div class="metric-label">Total Hits</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{len(unique_genomes)}</div>
                    <div class="metric-label">Genomes</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{len(set(h['gene'] for h in hits))}</div>
                    <div class="metric-label">Unique Genes</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{current_time.split()[0]}</div>
                    <div class="metric-label">Date</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🔍 Genes by Genome</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Genome</th>
                            <th>Count</th>
                            <th>Genes</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        for genome in sorted(unique_genomes):
            genes = genes_per_genome.get(genome, set())
            gene_list = ", ".join(sorted(genes))
            html_content += f"""
                    <tr>
                        <td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{gene_list}</td>
                    </tr>
"""
        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📈 Gene Frequency</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Frequency</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        gene_frequency = {}
        for hit in hits:
            gene = hit['gene']
            if gene not in gene_frequency:
                gene_frequency[gene] = set()
            gene_frequency[gene].add(hit['genome'])
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            genome_list = ", ".join(sorted(genomes))
            html_content += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{len(genomes)}</td>
                        <td>{genome_list}</td>
                    </tr>
"""
        html_content += f"""
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>KLEBOSCOPE</strong> - Batch ABRicate Summary</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p><strong>Technical Support:</strong> Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>University of Ghana Medical School</p>
            </div>
        </div>
    </div>
    <script>
        const quotes = {json.dumps(self.science_quotes)};
        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');
        function getRandomQuote() {{ return quotes[Math.floor(Math.random() * quotes.length)]; }}
        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>"""
        html_file = os.path.join(output_base, f"klebo_{database}_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Database summary HTML report: %s", html_file)
    
    def create_database_json_summaries(self, all_results: Dict[str, Any], output_base: str):
        """Create JSON summary for each database across all genomes."""
        self.logger.info("Creating JSON database summaries...")
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = {
                        'hits': [],
                        'genomes': [],
                        'gene_frequency': {}
                    }
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db]['hits'].append(hit_with_genome)
                if genome_name not in db_results[db]['genomes']:
                    db_results[db]['genomes'].append(genome_name)
        for db, data in db_results.items():
            if not data['hits']:
                continue
            gene_frequency = {}
            for hit in data['hits']:
                gene = hit['gene']
                if gene not in gene_frequency:
                    gene_frequency[gene] = {
                        'count': 0,
                        'genomes': set(),
                        'details': []
                    }
                gene_frequency[gene]['count'] += 1
                gene_frequency[gene]['genomes'].add(hit['genome'])
                gene_frequency[gene]['details'].append({
                    'genome': hit['genome'],
                    'product': hit['product'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'accession': hit['accession']
                })
            for gene in gene_frequency:
                gene_frequency[gene]['genomes'] = list(gene_frequency[gene]['genomes'])
            json_summary = {
                'metadata': {
                    'database': db,
                    'analysis_date': self.metadata['analysis_date'],
                    'tool': self.metadata['tool_name'],
                    'version': self.metadata['version'],
                    'total_hits': len(data['hits']),
                    'total_genomes': len(data['genomes']),
                    'unique_genes': len(gene_frequency)
                },
                'gene_frequency': gene_frequency,
                'hits': data['hits'][:100000]  # limit to first 100000 hits
            }
            json_file = os.path.join(output_base, f"klebo_{db}_summary.json")
            with open(json_file, 'w') as f:
                json.dump(json_summary, f, indent=2, default=str)
            self.logger.info("✓ Created JSON summary: %s", json_file)
    
    def create_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        """Create master JSON summary combining all databases."""
        self.logger.info("Creating master JSON summary...")
        master_summary = {
            'metadata': {
                'tool': 'Kleboscope ABRicate Module',
                'version': self.metadata['version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results),
                'databases_used': self.required_databases
            },
            'genome_summaries': {},
            'critical_findings': {},
            'cross_database_patterns': {}
        }
        all_hits_by_gene = {}
        for genome_name, genome_result in all_results.items():
            all_genome_hits = []
            for db_result in genome_result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            analysis = self.analyze_klebsiella_genes(all_genome_hits)
            master_summary['genome_summaries'][genome_name] = {
                'total_hits': genome_result['total_hits'],
                'critical_resistance': len(analysis['critical_resistance_genes']),
                'high_risk_virulence': len(analysis['high_risk_virulence_genes']),
                'beta_lactamases': len(analysis['beta_lactamase_genes']),
                'other_genes': len(analysis['other_genes'])
            }
            if analysis['critical_resistance_genes']:
                if 'critical_genomes' not in master_summary['critical_findings']:
                    master_summary['critical_findings']['critical_genomes'] = []
                master_summary['critical_findings']['critical_genomes'].append(genome_name)
            if analysis['high_risk_virulence_genes']:
                if 'hv_genomes' not in master_summary['critical_findings']:
                    master_summary['critical_findings']['hv_genomes'] = []
                master_summary['critical_findings']['hv_genomes'].append(genome_name)
            for hit in all_genome_hits:
                gene = hit['gene']
                if gene not in all_hits_by_gene:
                    all_hits_by_gene[gene] = {
                        'count': 0,
                        'genomes': set(),
                        'products': set(),
                        'databases': set()
                    }
                all_hits_by_gene[gene]['count'] += 1
                all_hits_by_gene[gene]['genomes'].add(genome_name)
                all_hits_by_gene[gene]['products'].add(hit['product'])
                all_hits_by_gene[gene]['databases'].add(hit['database'])
        for gene in all_hits_by_gene:
            all_hits_by_gene[gene]['genomes'] = list(all_hits_by_gene[gene]['genomes'])
            all_hits_by_gene[gene]['products'] = list(all_hits_by_gene[gene]['products'])
            all_hits_by_gene[gene]['databases'] = list(all_hits_by_gene[gene]['databases'])
        master_summary['cross_database_patterns'] = {
            'total_genes_found': len(all_hits_by_gene),
            'common_genes': {g: d for g, d in all_hits_by_gene.items() if d['count'] > 1},
            'top_genes': sorted(
                [(g, d) for g, d in all_hits_by_gene.items()],
                key=lambda x: x[1]['count'],
                reverse=True
            )[:50]
        }
        json_file = os.path.join(output_base, "klebo_abricate_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master_summary, f, indent=2, default=str)
        self.logger.info("✓ Created master JSON summary: %s", json_file)
    
    def process_single_genome(self, genome_file: str, output_base: str = "abricate_results") -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        os.makedirs(results_dir, exist_ok=True)
        databases = self.required_databases
        results = {}
        for db in databases:
            result = self.run_abricate_single_db(genome_file, db, results_dir)
            results[db] = result
        self.create_comprehensive_html_report(genome_name, results, results_dir)
        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values())
        }
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "abricate_results") -> Dict[str, Any]:
        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")
        self.setup_abricate_databases()
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa",
                          f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        self.logger.info("Found %d genomes", len(genome_files))
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        if len(genome_files) > 1 and self.cpus > 1:
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                future_to_genome = {executor.submit(self.process_single_genome, g, output_base): g for g in genome_files}
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("✓ Completed: %s (%d hits)", result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("✗ Failed: %s - %s", genome, e)
        else:
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("✓ Completed: %s (%d hits)", result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("✗ Failed: %s - %s", genome, e)
        self.create_database_summaries(all_results, output_base)
        self.create_database_json_summaries(all_results, output_base)
        self.create_master_json_summary(all_results, output_base)
        self.logger.info("=== ANALYSIS COMPLETE ===")
        return all_results

def main():
    parser = argparse.ArgumentParser(
        description='Kleboscope ABRicate Analysis - MAXIMUM SPEED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python kleb_abricate_module.py "*.fasta"
  python kleb_abricate_module.py "*.fna" --output my_results --cpus 8"""
    )
    parser.add_argument('pattern', help='File pattern for genomes (e.g., "*.fasta")')
    parser.add_argument('--cpus', '-c', type=int, default=None, help='Number of CPU cores (default: auto-detect)')
    parser.add_argument('--output', '-o', default='abricate_results', help='Output directory')
    args = parser.parse_args()
    executor = AbricateExecutor(cpus=args.cpus)
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        executor.logger.info("\n" + "="*50)
        executor.logger.info("📊 FINAL SUMMARY")
        executor.logger.info("="*50)
        total_crit = 0
        total_high = 0
        for genome_name, result in results.items():
            all_hits = []
            for db_result in result['results'].values():
                all_hits.extend(db_result['hits'])
            analysis = executor.analyze_klebsiella_genes(all_hits)
            executor.logger.info("✓ %s: %d hits (critical: %d, high-risk virulence: %d, beta-lactamases: %d)",
                               genome_name, result['total_hits'],
                               analysis['total_critical_resistance'],
                               analysis['total_high_risk_virulence'],
                               analysis['total_beta_lactamase'])
            total_crit += analysis['total_critical_resistance']
            total_high += analysis['total_high_risk_virulence']
        executor.logger.info("\n📁 Results saved to: %s", args.output)
        executor.logger.info("⚡ MAXIMUM SPEED: %d cores utilized", executor.cpus)
    except Exception as e:
        executor.logger.error("Analysis failed: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()