#!/usr/bin/env python3
"""
Kleboscope Main Orchestrator – Parallel Execution with Scientific Quotes
Complete K. pneumoniae typing & resistance pipeline
Author: Brown Beckley <brownbeckley94@gmail.com>
Version: 1.2.0
Date: 2026-07-28
Affiliation: University of Ghana Medical School – Department of Medical Biochemistry
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import tempfile
import random
import logging
import traceback
import signal
import atexit
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Set
from concurrent.futures import ThreadPoolExecutor, as_completed

__version__ = "1.2.0"

# Global flag for graceful shutdown
_should_exit = False

def signal_handler(sig, frame):
    global _should_exit
    print("\n" + "=" * 60)
    print("⚠️  Interrupt signal received. Cleaning up and exiting gracefully...")
    print("=" * 60)
    _should_exit = True
    # Allow some time for cleanup
    sys.exit(1)

# Register signal handlers
signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)


class Color:
    """ANSI color codes for coloured output."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'


class KleboscopeOrchestrator:
    """Kleboscope orchestrator with temporary directories and resource control."""

    def __init__(self):
        self.base_dir = Path(__file__).parent
        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]
        self.keep_temp = False
        self.logger = None
        self.user_output_dir = None
        self.log_file = None
        self.temp_dirs: Set[str] = set()  # Track temp dirs for cleanup
        atexit.register(self._cleanup_temp_dirs)

        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'mlst': 'mlst_results',
            'kaptive': 'kaptive_results',
            'abricate': 'abricate_results',
            'amr': 'klebo_amrfinder_results',
            'summary': 'KLEBOSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS',
            'sample_centric': 'KLEBOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS'
        }

        self.summary_files = {
            'mlst_summary.html': ('mlst_results', 'mlst_summary.html'),
            'mlst_summary.tsv': ('mlst_results', 'mlst_summary.tsv'),
            'mlst_summary.json': ('mlst_results', 'mlst_summary.json'),
            'klebo_fasta_qc_summary.html': ('fasta_qc_results', 'klebo_fasta_qc_summary.html'),
            'klebo_kaptive_summary.html': ('kaptive_results', 'klebo_kaptive_summary.html'),
            'klebo_amrfinder_summary_report.html': ('klebo_amrfinder_results', 'klebo_amrfinder_summary_report.html'),
            'mutation_summary.html': ('klebo_amrfinder_results', 'mutation_summary.html'),
            'mutation_summary.tsv': ('klebo_amrfinder_results', 'mutation_summary.tsv'),
            'mutation_master_summary.json': ('klebo_amrfinder_results', 'mutation_master_summary.json'),
            'klebo_card_summary_report.html': ('abricate_results', 'klebo_card_summary_report.html'),
            'klebo_ncbi_summary_report.html': ('abricate_results', 'klebo_ncbi_summary_report.html'),
            'klebo_resfinder_summary_report.html': ('abricate_results', 'klebo_resfinder_summary_report.html'),
            'klebo_vfdb_summary_report.html': ('abricate_results', 'klebo_vfdb_summary_report.html'),
            'klebo_argannot_summary_report.html': ('abricate_results', 'klebo_argannot_summary_report.html'),
            'klebo_megares_summary_report.html': ('abricate_results', 'klebo_megares_summary_report.html'),
            'klebo_ecoli_vf_summary_report.html': ('abricate_results', 'klebo_ecoli_vf_summary_report.html'),
            'klebo_bacmet2_summary_report.html': ('abricate_results', 'klebo_bacmet2_summary_report.html'),
            'klebo_plasmidfinder_summary_report.html': ('abricate_results', 'klebo_plasmidfinder_summary_report.html'),
            'klebo_ecoh_summary_report.html': ('abricate_results', 'klebo_ecoh_summary_report.html')
        }

        # Sample-centric reporter expects additional files (e.g., TSV summaries)
        self.sample_centric_summary_files = {
            'klebo_amrfinder_summary.tsv': ('klebo_amrfinder_results', 'klebo_amrfinder_summary.tsv'),
            'klebo_amrfinder_statistics_summary.tsv': ('klebo_amrfinder_results', 'klebo_amrfinder_statistics_summary.tsv'),
            'mutation_summary.tsv': ('klebo_amrfinder_results', 'mutation_summary.tsv'),
            'klebo_kaptive_summary.tsv': ('kaptive_results', 'klebo_kaptive_summary.tsv'),
            'mlst_summary.tsv': ('mlst_results', 'mlst_summary.tsv'),
            'klebo_fasta_qc_summary.tsv': ('fasta_qc_results', 'klebo_fasta_qc_summary.tsv'),
        }
        # Add ABRicate TSV files
        abricate_tsv_files = [
            'klebo_card_abricate_summary.tsv',
            'klebo_resfinder_abricate_summary.tsv',
            'klebo_vfdb_abricate_summary.tsv',
            'klebo_argannot_abricate_summary.tsv',
            'klebo_megares_abricate_summary.tsv',
            'klebo_ecoli_vf_abricate_summary.tsv',
            'klebo_bacmet2_abricate_summary.tsv',
            'klebo_plasmidfinder_abricate_summary.tsv',
            'klebo_ncbi_abricate_summary.tsv',
            'klebo_ecoh_abricate_summary.tsv'
        ]
        for fname in abricate_tsv_files:
            self.sample_centric_summary_files[fname] = ('abricate_results', fname)

        self.module_subtitles = {
            'QC': 'Sequence Quality Control & Statistics',
            'MLST': 'Multi-Locus Sequence Typing',
            'Kaptive': 'K/O Locus Typing',
            'ABRicate': 'Comprehensive Resistance & Virulence Screening',
            'AMR': 'Antimicrobial Resistance Gene Detection',
            'Ultimate Reporter': 'Gene‑centric Integrated Analysis',
            'Sample Centric Reporter': 'Isolate‑centric Interactive Analysis'
        }

    def setup_colors(self):
        self.color_info = Color.CYAN
        self.color_success = Color.BRIGHT_GREEN
        self.color_warning = Color.BRIGHT_YELLOW
        self.color_error = Color.BRIGHT_RED
        self.color_highlight = Color.BRIGHT_CYAN
        self.color_banner = Color.BRIGHT_MAGENTA
        self.color_module = Color.BRIGHT_BLUE
        self.color_sample = Color.GREEN
        self.color_file = Color.YELLOW
        self.color_reset = Color.RESET

    def print_color(self, message: str, color: str = Color.RESET, bold: bool = False):
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")

    def print_header(self, title: str, subtitle: str = ""):
        print()
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * 20}{title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}{' ' * 22}{subtitle}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print()

    def print_info(self, message: str):
        print(f"{self.color_info}[INFO]{Color.RESET} {message}")

    def print_success(self, message: str):
        print(f"{self.color_success}✓{Color.RESET} {message}")

    def print_warning(self, message: str):
        print(f"{self.color_warning}⚠️{Color.RESET} {message}")

    def print_error(self, message: str):
        print(f"{self.color_error}✗{Color.RESET} {message}")

    def print_command(self, command: str):
        print(f"{Color.DIM}{Color.WHITE}  $ {command}{Color.RESET}")

    def _cleanup_temp_dirs(self):
        """Remove all temporary directories on exit."""
        if self.temp_dirs:
            for td in self.temp_dirs:
                try:
                    shutil.rmtree(td, ignore_errors=True)
                    if self.logger:
                        self.logger.info(f"Cleaned up temp dir: {td}")
                except Exception as e:
                    if self.logger:
                        self.logger.warning(f"Could not remove temp dir {td}: {e}")

    def _get_scientific_quotes(self):
        return [
            {"quote": "Science is organised knowledge.", "author": "Herbert Spencer", "theme": "knowledge"},
            {"quote": "The science of today is the technology of tomorrow.", "author": "Edward Teller", "theme": "technology"},
            {"quote": "Nature is the source of all true knowledge.", "author": "Leonardo da Vinci", "theme": "nature"},
            {"quote": "Biology is the most powerful technology ever created.", "author": "Freeman Dyson", "theme": "biology"},
            {"quote": "Genomics is a lens on biology.", "author": "Eric Lander", "theme": "genomics"},
            {"quote": "Every microbe has its own story.", "author": "Anonymous", "theme": "microbiology"},
            {"quote": "Data beats emotions.", "author": "Sean Rad", "theme": "data"},
            {"quote": "Code is poetry.", "author": "WordPress", "theme": "programming"},
            {"quote": "Sequence today, understand tomorrow.", "author": "Anonymous", "theme": "sequencing"},
            {"quote": "Microbes rule the world.", "author": "Paul Stamets", "theme": "microbiology"},
            {"quote": "In every drop, a universe.", "author": "Antonie van Leeuwenhoek", "theme": "microscopy"},
            {"quote": "Genes are the language of life.", "author": "Francis Collins", "theme": "genetics"},
            {"quote": "Resistance is not futile.", "author": "Antibiotic Researcher", "theme": "resistance"},
            {"quote": "Evolution in a petri dish.", "author": "Richard Lenski", "theme": "evolution"},
            {"quote": "Small things, big impact.", "author": "Microbiologist", "theme": "microbes"},
            {"quote": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson", "theme": "science"},
            {"quote": "In science, there are no shortcuts to truth.", "author": "Karl Popper", "theme": "science"},
            {"quote": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie", "theme": "understanding"},
            {"quote": "The most exciting phrase to hear in science is not 'Eureka!' but 'That's funny...'", "author": "Isaac Asimov", "theme": "discovery"},
            {"quote": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur", "theme": "global"},
            {"quote": "DNA is like a computer program but far, far more advanced than any software ever created.", "author": "Bill Gates", "theme": "genomics"},
        ]

    def display_random_quote(self):
        if not self.quotes:
            return
        quote_data = random.choice(self.quotes)
        quote = quote_data["quote"]
        author = quote_data["author"]
        theme = quote_data.get("theme", "science")
        quote_color = random.choice(self.quote_colors)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        theme_icons = {"microbiology": "🦠", "discovery": "🔬", "knowledge": "📚", "medicine": "⚕️",
                       "science": "🧪", "research": "🔍", "exploration": "🚀", "curiosity": "🤔",
                       "practice": "🛠️", "motivation": "💪", "nature": "🌿", "inquiry": "❓",
                       "technology": "💻", "understanding": "🧠", "perspective": "👁️", "innovation": "💡",
                       "recognition": "🏆", "purpose": "🎯", "biology": "🧬", "genomics": "🧬",
                       "data": "📊", "programming": "💻", "sequencing": "🧬", "microscopy": "🔬",
                       "genetics": "🧬", "resistance": "🛡️", "evolution": "🔄", "microbes": "🦠",
                       "global": "🌍", "conservation": "🌱", "time": "⏳", "collaboration": "🤝"}
        icon = theme_icons.get(theme, "💭")
        print()
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}[{current_time}] {icon} SCIENTIFIC INSIGHT: {Color.RESET}")
        print()
        print(f"{quote_color}   \"{quote}\"{Color.RESET}")
        print(f"{Color.BOLD}{Color.WHITE}   — {author}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}   Theme: {theme.capitalize()}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print()

    def update_amr_database(self, force: bool = False) -> bool:
        amr_module_path = self.base_dir / "modules" / "kleb_amr_module"
        amr_script = amr_module_path / "klebo_amrfinder.py"
        if not amr_script.exists():
            self.print_error(f"AMR script not found at: {amr_script}")
            return False
        self.print_info("Updating AMRfinderPlus database...")
        flag = "--force-update" if force else "--update-db"
        cmd = [sys.executable, str(amr_script), flag]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0:
            self.print_success("AMR database updated successfully.")
            version_cmd = [sys.executable, str(amr_script), "--db-version"]
            version_result = subprocess.run(version_cmd, capture_output=True, text=True, cwd=amr_module_path)
            if version_result.returncode == 0:
                self.print_info(f"New database version: {version_result.stdout.strip()}")
            return True
        else:
            self.print_error("AMR database update failed.")
            if result.stderr:
                print(result.stderr)
            return False

    def find_fasta_files(self, input_path: str) -> List[Path]:
        self.print_info(f"Searching for files with pattern: {input_path}")
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and
                           f.lower().endswith(('.fna', '.fasta', '.fa', '.fn')) and
                           not Path(f).name.startswith('.')]
            self.print_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)

        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fn']:
            self.print_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]

        if input_path_obj.is_dir():
            patterns = [f"{input_path}/*.fna", f"{input_path}/*.fasta", f"{input_path}/*.fa", f"{input_path}/*.fn"]
            fasta_files = []
            for pattern in patterns:
                matched_files = glob.glob(pattern)
                for file_path in matched_files:
                    path = Path(file_path)
                    if path.is_file() and not path.name.startswith('.'):
                        fasta_files.append(path)
            fasta_files = sorted(list(set(fasta_files)))
            if fasta_files:
                self.print_success(f"Found {len(fasta_files)} FASTA files in directory")
            else:
                self.print_warning(f"No FASTA files found in directory: {input_path}")
            return fasta_files

        self.print_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return "*.fna"
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            return f"*{list(extensions)[0]}"
        return "*"

    def run_module_in_temp(self, module_name: str, fasta_files: List[Path],
                           cmd_str: str, result_subdir: str = None) -> Tuple[bool, str]:
        global _should_exit
        if _should_exit:
            return False, "Execution interrupted by user."

        module_orig = self.base_dir / "modules" / module_name
        if not module_orig.exists():
            return False, f"Module directory not found: {module_orig}"

        temp_dir = tempfile.mkdtemp(prefix=f"kleboscope_{module_name}_")
        self.temp_dirs.add(temp_dir)
        self.logger.info(f"Temporary directory for {module_name}: {temp_dir}")

        try:
            # Copy module
            shutil.copytree(module_orig, Path(temp_dir) / module_name, dirs_exist_ok=True)
            module_work_dir = Path(temp_dir) / module_name

            # Copy fasta files
            for f in fasta_files:
                shutil.copy2(f, module_work_dir / f.name)

            self.logger.info(f"Copied {len(fasta_files)} files to {module_name} module")
            pattern = self.get_file_pattern(fasta_files)
            self.logger.info(f"Running {module_name} analysis with pattern: {pattern}")

            # Run the command
            result = subprocess.run(cmd_str, shell=True, cwd=module_work_dir, capture_output=True, text=True)
            if result.stdout:
                self.logger.info(f"STDOUT:\n{result.stdout}")
            if result.stderr:
                self.logger.info(f"STDERR:\n{result.stderr}")

            if _should_exit:
                return False, "Execution interrupted after module completion."

            if result.returncode != 0:
                self.logger.error(f"{module_name} failed with return code {result.returncode}")
                return False, f"{module_name} failed with return code {result.returncode}"

            # Copy results
            if result_subdir:
                src = module_work_dir / result_subdir
                if src.exists():
                    dst = self.user_output_dir / result_subdir
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Results copied to {dst}")
                else:
                    self.logger.warning(f"Result directory not found: {src}")

            return True, f"{module_name} completed successfully"

        except Exception as e:
            self.logger.error(f"Exception in {module_name}: {e}\n{traceback.format_exc()}")
            return False, f"Exception in {module_name}: {e}"
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.temp_dirs.discard(temp_dir)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def run_qc(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"python klebo_fasta_qc.py '{pattern}'"
        return self.run_module_in_temp("kleb_qc_module", fasta_files, cmd, self.output_dirs['qc'])

    def run_mlst(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"python klebo_mlst.py -i '{pattern}' -o {self.output_dirs['mlst']} -db db -sc . --batch"
        return self.run_module_in_temp("kleb_mlst_module", fasta_files, cmd, self.output_dirs['mlst'])

    def run_kaptive(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"python klebo_kaptive.py -i '{pattern}' -o {self.output_dirs['kaptive']}"
        return self.run_module_in_temp("kleb_kaptive_module", fasta_files, cmd, self.output_dirs['kaptive'])

    def run_abricate(self, fasta_files: List[Path], output_dir: Path, threads: int,
                    min_id: Optional[float] = None, min_cov: Optional[float] = None) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"python kleb_abricate_module.py '{pattern}'"
        if min_id is not None:
            cmd += f" --min-id {min_id}"
        if min_cov is not None:
            cmd += f" --min-cov {min_cov}"
        return self.run_module_in_temp("kleb_abricate_module", fasta_files, cmd, "abricate_results")

    def run_amr(self, fasta_files: List[Path], output_dir: Path, threads: int,
                min_identity: Optional[float] = None, min_coverage: Optional[float] = None,
                skip_mutations: bool = False) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = f"python klebo_amrfinder.py '{pattern}'"
        if min_identity is not None:
            cmd += f" --min-identity {min_identity}"
        if min_coverage is not None:
            cmd += f" --min-coverage {min_coverage}"
        if skip_mutations:
            cmd += " --skip-mutations"
        return self.run_module_in_temp("kleb_amr_module", fasta_files, cmd, "klebo_amrfinder_results")

    def run_summary(self, output_dir: Path) -> Tuple[bool, str]:
        global _should_exit
        if _should_exit:
            return False, "Execution interrupted by user."

        summary_module_orig = self.base_dir / "modules" / "kleb_gene_centric_module"
        if not summary_module_orig.exists():
            return False, f"Summary module not found: {summary_module_orig}"

        # ✅ Create temporary directory (like other modules)
        temp_dir = tempfile.mkdtemp(prefix="kleboscope_summary_")
        self.temp_dirs.add(temp_dir)
        self.logger.info(f"Temporary directory for summary: {temp_dir}")

        try:
            # Copy module to temp
            shutil.copytree(summary_module_orig, Path(temp_dir) / "kleb_gene_centric_module", dirs_exist_ok=True)
            summary_work_dir = Path(temp_dir) / "kleb_gene_centric_module"

            # Copy summary files
            copied = 0
            missing = 0
            for target_name, (subdir, filename) in self.summary_files.items():
                if subdir:
                    source = output_dir / subdir / filename
                else:
                    source = output_dir / filename
                if source.exists():
                    shutil.copy2(source, summary_work_dir / target_name)
                    copied += 1
                else:
                    missing += 1

            if missing > 0:
                self.logger.warning(f"Copied {copied} files, {missing} missing.")

            # Run reporter in temp dir
            self.logger.info("Running gene‑centric ultimate reporter...")
            cmd = [sys.executable, str(summary_work_dir / "kleboscope_ultimate_gene_centric_reporter.py"), "-i", "."]
            result = subprocess.run(cmd, cwd=summary_work_dir, capture_output=True, text=True)

            # Copy results back
            summary_source = summary_work_dir / self.output_dirs['summary']
            summary_target = output_dir / self.output_dirs['summary']
            if summary_source.exists():
                if summary_target.exists():
                    shutil.rmtree(summary_target)
                shutil.copytree(summary_source, summary_target)
                self.logger.info(f"Gene‑centric reports copied to: {summary_target}")

            return result.returncode == 0, "Gene‑centric summary completed"

        except Exception as e:
            self.logger.error(f"Exception in summary: {e}\n{traceback.format_exc()}")
            return False, f"Exception in summary: {e}"
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.temp_dirs.discard(temp_dir)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def run_sample_centric_reporter(self, output_dir: Path) -> Tuple[bool, str]:
        """Copy summary files and run the sample‑centric reporter in a temporary directory."""
        global _should_exit
        if _should_exit:
            return False, "Execution interrupted by user."

        sample_module_orig = self.base_dir / "modules" / "kleb_sample_centric_module"
        if not sample_module_orig.exists():
            return False, f"Sample‑centric module not found: {sample_module_orig}"

        # ✅ Create temporary directory (like other modules)
        temp_dir = tempfile.mkdtemp(prefix="kleboscope_sample_centric_")
        self.temp_dirs.add(temp_dir)
        self.logger.info(f"Temporary directory for sample‑centric reporter: {temp_dir}")

        try:
            # Copy module to temp
            shutil.copytree(sample_module_orig, Path(temp_dir) / "kleb_sample_centric_module", dirs_exist_ok=True)
            sample_work_dir = Path(temp_dir) / "kleb_sample_centric_module"

            # Copy all required summary files (from summary_files and sample_centric_summary_files)
            all_files = {}
            all_files.update(self.summary_files)
            all_files.update(self.sample_centric_summary_files)

            copied = 0
            missing = 0
            for target_name, (subdir, filename) in all_files.items():
                if subdir:
                    source = output_dir / subdir / filename
                else:
                    source = output_dir / filename
                if source.exists():
                    shutil.copy2(source, sample_work_dir / target_name)
                    copied += 1
                else:
                    missing += 1

            if missing > 0:
                self.logger.warning(f"Copied {copied} files, {missing} missing for sample‑centric reporter.")

            # Run sample‑centric reporter in temp dir
            self.logger.info("Running sample‑centric ultimate reporter...")
            cmd = [sys.executable, str(sample_work_dir / "kleboscope_ultimate_sample_centric_report.py"), "-i", "."]
            result = subprocess.run(cmd, cwd=sample_work_dir, capture_output=True, text=True)
            if result.stdout:
                self.logger.info(result.stdout)
            if result.stderr:
                self.logger.info(result.stderr)

            if _should_exit:
                return False, "Execution interrupted after reporter completion."

            # Copy results back to output directory
            sample_source = sample_work_dir / self.output_dirs['sample_centric']
            sample_target = output_dir / self.output_dirs['sample_centric']
            if sample_source.exists():
                if sample_target.exists():
                    shutil.rmtree(sample_target)
                shutil.copytree(sample_source, sample_target)
                self.logger.info(f"Sample‑centric reports copied to: {sample_target}")
                files = list(sample_target.glob("*"))
                html_count = len([f for f in files if f.suffix == '.html'])
                json_count = len([f for f in files if f.suffix == '.json'])
                csv_count = len([f for f in files if f.suffix == '.csv'])
                self.logger.info(f"📊 {html_count} HTML, {json_count} JSON, {csv_count} CSV files")
            else:
                self.logger.warning(f"Sample‑centric reports directory not found: {sample_source}")

            return result.returncode == 0, "Sample‑centric summary completed"

        except Exception as e:
            self.logger.error(f"Exception in sample‑centric reporter: {e}\n{traceback.format_exc()}")
            return False, f"Exception in sample‑centric reporter: {e}"
        finally:
            if not self.keep_temp:
                shutil.rmtree(temp_dir, ignore_errors=True)
                self.temp_dirs.discard(temp_dir)
                self.logger.info(f"Removed temporary directory: {temp_dir}")

    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                              skip_modules: Dict[str, bool] = None,
                              skip_summary: bool = False,
                              skip_sample_centric: bool = False,
                              update_amr_db_only: bool = False,
                              force_update_amr: bool = False,
                              amr_min_identity: Optional[float] = None,
                              amr_min_coverage: Optional[float] = None,
                              amr_skip_mutations: bool = False,
                              abricate_min_id: Optional[float] = None,
                              abricate_min_cov: Optional[float] = None):
        if update_amr_db_only or force_update_amr:
            self.update_amr_database(force=force_update_amr)
            return

        if skip_modules is None:
            skip_modules = {}

        start_time = datetime.now()
        self.display_banner()
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        self.user_output_dir = output_path

        log_file = output_path / "kleboscope_run.log"
        self._setup_logging(log_file)

        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found! Analysis stopped.")
            return
        self.print_success(f"Starting analysis of {len(fasta_files)} K. pneumoniae samples")

        # Create output subdirectories
        for subdir in self.output_dirs.values():
            (output_path / subdir).mkdir(exist_ok=True)

        self.print_header("ANALYSIS PLAN", "Modules to be executed")
        plan = [
            ("QC", not skip_modules.get('qc', False)),
            ("MLST", not skip_modules.get('mlst', False)),
            ("Kaptive", not skip_modules.get('kaptive', False)),
            ("ABRicate", not skip_modules.get('abricate', False)),
            ("AMR", not skip_modules.get('amr', False)),
            ("Gene‑centric Reporter", not skip_summary),
            ("Sample‑centric Reporter", not skip_sample_centric),
        ]
        for analysis, enabled in plan:
            if enabled:
                print(f"   {Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET} - {analysis}")
            else:
                print(f"   {Color.YELLOW}⏸️  SKIPPED{Color.RESET} - {analysis}")
        if amr_min_identity is not None:
            print(f"   AMR min identity: {amr_min_identity}")
        if amr_min_coverage is not None:
            print(f"   AMR min coverage: {amr_min_coverage}")
        if amr_skip_mutations:
            print("   AMR mutation reporting: SKIPPED")
        if abricate_min_id is not None:
            print(f"   ABRicate min identity: {abricate_min_id}")
        if abricate_min_cov is not None:
            print(f"   ABRicate min coverage: {abricate_min_cov}")
        print()

        # First batch: QC, MLST, Kaptive
        tasks = []
        if not skip_modules.get('qc', False):
            tasks.append(("QC", self.run_qc, (fasta_files, output_path, threads)))
        if not skip_modules.get('mlst', False):
            tasks.append(("MLST", self.run_mlst, (fasta_files, output_path, threads)))
        if not skip_modules.get('kaptive', False):
            tasks.append(("Kaptive", self.run_kaptive, (fasta_files, output_path, threads)))

        if tasks:
            self.print_info(f"Running {len(tasks)} analyses in parallel...")
            results = {}
            with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
                future_to_name = {executor.submit(func, *args): name for name, func, args in tasks}
                for future in as_completed(future_to_name):
                    name = future_to_name[future]
                    try:
                        success, output = future.result()
                        results[name] = (success, output)
                    except Exception as e:
                        results[name] = (False, f"Exception: {e}")

            for name, _, _ in tasks:
                success, output = results.get(name, (False, "No result"))
                subtitle = self.module_subtitles.get(name, "")
                self.print_header(f"{name} Analysis", subtitle)
                if success:
                    self.print_success(f"✅ {name} completed")
                else:
                    self.print_error(f"❌ {name} failed")
                self.display_random_quote()
        else:
            self.print_info("No analyses in first batch (all skipped).")

        # ABRicate
        if not skip_modules.get('abricate', False):
            self.print_header("ABRICATE ANALYSIS", self.module_subtitles['ABRicate'])
            success, output = self.run_abricate(fasta_files, output_path, threads,
                                                min_id=abricate_min_id, min_cov=abricate_min_cov)
            if success:
                self.print_success("✅ ABRicate completed")
            else:
                self.print_error("❌ ABRicate failed")
            self.display_random_quote()
        else:
            self.print_info("Skipping ABRicate analysis.")

        # AMR
        if not skip_modules.get('amr', False):
            self.print_header("AMR ANALYSIS", self.module_subtitles['AMR'])
            success, output = self.run_amr(fasta_files, output_path, threads,
                                           min_identity=amr_min_identity,
                                           min_coverage=amr_min_coverage,
                                           skip_mutations=amr_skip_mutations)
            if success:
                self.print_success("✅ AMR completed")
            else:
                self.print_error("❌ AMR failed")
            self.display_random_quote()
        else:
            self.print_info("Skipping AMR analysis.")

        # Gene‑centric summary
        if not skip_summary:
            self.print_header("ULTIMATE REPORTER (Gene‑Centric)", self.module_subtitles['Ultimate Reporter'])
            success, output = self.run_summary(output_path)
            if success:
                self.print_success("✅ Gene‑centric reporter completed")
            else:
                self.print_warning("Gene‑centric reporter had issues")
            self.display_random_quote()

        # Sample‑centric summary
        if not skip_sample_centric:
            self.print_header("SAMPLE‑CENTRIC REPORTER", self.module_subtitles['Sample Centric Reporter'])
            success, output = self.run_sample_centric_reporter(output_path)
            if success:
                self.print_success("✅ Sample‑centric reporter completed")
            else:
                self.print_warning("Sample‑centric reporter had issues")
            self.display_random_quote()

        # Final report
        analysis_time = datetime.now() - start_time
        self.print_header("ANALYSIS COMPLETE", f"Time elapsed: {str(analysis_time).split('.')[0]}")
        self.print_success(f"🎉 Analysis complete! Results in: {output_path}")

        for subdir in sorted(output_path.iterdir()):
            if subdir.is_dir():
                file_count = len(list(subdir.glob("*")))
                self.print_info(f"  📁 {subdir.name} ({file_count} files)")

        print()
        print(f"{Color.BRIGHT_YELLOW}📚 Please cite Kleboscope:{Color.RESET}")
        print(f"   Beckley, B., et al. Kleboscope: a gene‑centric, species‑optimized computational pipeline for comprehensive")
        print(f"   Klebsiella pneumoniae genomic surveillance. (2026).")
        print(f"   {Color.BRIGHT_CYAN}https://github.com/bbeckley-hub/kleboscope{Color.RESET}")
        print()
        print(f"{Color.BRIGHT_YELLOW}💡 Support & Contributions:{Color.RESET}")
        print(f"   • Issues & feature requests: https://github.com/bbeckley-hub/kleboscope/issues")
        print(f"   • Email: brownbeckley94@gmail.com")
        print()

        self.display_random_quote()
        self.logger.info(f"Analysis completed. Results in: {output_path}")

    def _setup_logging(self, log_file: Path):
        self.logger = logging.getLogger("Kleboscope")
        self.logger.setLevel(logging.INFO)
        if self.logger.handlers:
            self.logger.handlers.clear()
        fh = logging.FileHandler(log_file, mode='w')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.logger.info(f"Kleboscope v{__version__} started")
        self.logger.info(f"Log file: {log_file}")

    def display_banner(self):
        banner = f"""{Color.BOLD}{Color.BRIGHT_MAGENTA}
{'='*80}
{' '*20}🦠 KLEBOSCOPE - K. pneumoniae Genomic Analysis Pipeline v{__version__}
{'='*80}
{Color.RESET}{Color.BRIGHT_CYAN}
Complete K. pneumoniae genomic analysis pipeline
MLST | K/O Locus | AMR | Virulence | Plasmid | Quality Control | Summary Reports

Critical Genes Tracked:
🔴 Carbapenemases (KPC, NDM, OXA-48)
🟠 Colistin (mcr)
🟡 Tigecycline (tetX)
🟢 ICEKp Markers (ybt, clb, iro, rmp)
🔵 Virulence Plasmid Markers (iro, iuc, rmp, rmpA2)
🟣 Biocides & Heavy Metals (qac, sil, mer, ars, pco)
⚪ Adhesins (fim, mrk, ecp)
🔶 Secretion Systems (tss)
💧 Siderophores
🧪 Toxins
{Color.RESET}{Color.DIM}
Author: Brown Beckley | Email: brownbeckley94@gmail.com
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Version: {__version__}
GitHub: https://github.com/bbeckley-hub/kleboscope
Citation: Beckley, B., et al. Kleboscope: a gene‑centric, species‑optimized computational pipeline for comprehensive 
          Klebsiella pneumoniae genomic surveillance. (2026).
{'='*80}{Color.RESET}
"""
        print(banner)

    def print_colored_help(self):
        self.display_banner()
        print(f"{Color.BRIGHT_YELLOW}USAGE:{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope{Color.RESET} {Color.CYAN}-i INPUT -o OUTPUT{Color.RESET} [OPTIONS]")
        print(f"  {Color.GREEN}kleboscope --update-amr-db{Color.RESET}  (Update AMR database incrementally)")
        print(f"  {Color.GREEN}kleboscope --force-update-amr-db{Color.RESET} (Force full AMR database update)")
        print(f"  {Color.GREEN}kleboscope --version{Color.RESET} Show version and exit")
        print()
        print(f"{Color.BRIGHT_YELLOW}REQUIRED ARGUMENTS (for analysis):{Color.RESET}")
        print(f"  {Color.GREEN}-i, --input{Color.RESET} INPUT    Input FASTA file(s) – supports glob patterns like \"*.fna\"")
        print(f"  {Color.GREEN}-o, --output{Color.RESET} OUTPUT  Output directory for results\n")
        print(f"{Color.BRIGHT_YELLOW}OPTIONAL ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-h, --help{Color.RESET}            Show this help message")
        print(f"  {Color.GREEN}-t, --threads{Color.RESET} THREADS  Number of threads (default: 2)")
        print(f"  {Color.GREEN}--version{Color.RESET}               Show version and exit")
        print(f"  {Color.GREEN}--update-amr-db{Color.RESET}         Update AMRfinderPlus database incrementally and exit")
        print(f"  {Color.GREEN}--force-update-amr-db{Color.RESET}  Force full AMR database update (overwrites old) and exit")
        print(f"  {Color.GREEN}--keep-temp{Color.RESET}             Do not delete temporary directories (for debugging)")
        print()
        print(f"{Color.BRIGHT_YELLOW}AMR FINDER PLUS FLAGS:{Color.RESET}")
        print(f"  {Color.GREEN}--amr-min-identity{Color.RESET}   Minimum identity for AMR hits (0..1)")
        print(f"  {Color.GREEN}--amr-min-coverage{Color.RESET}   Minimum coverage for AMR hits (0..1)")
        print(f"  {Color.GREEN}--skip-amr-mutations{Color.RESET}     Disable point mutation reporting in AMR (enabled by default)")
        print()
        print(f"{Color.BRIGHT_YELLOW}ABRICATE FLAGS:{Color.RESET}")
        print(f"  {Color.GREEN}--abricate-min-id{Color.RESET}   Minimum identity for ABRicate hits (default: 80)")
        print(f"  {Color.GREEN}--abricate-min-cov{Color.RESET}   Minimum coverage for ABRicate hits (default: 80)")
        print()
        print(f"{Color.BRIGHT_YELLOW}SKIP OPTIONS:{Color.RESET}")
        print(f"  {Color.GREEN}--skip-qc{Color.RESET}               Skip QC analysis")
        print(f"  {Color.GREEN}--skip-mlst{Color.RESET}             Skip MLST analysis")
        print(f"  {Color.GREEN}--skip-kaptive{Color.RESET}          Skip Kaptive analysis")
        print(f"  {Color.GREEN}--skip-abricate{Color.RESET}         Skip ABRicate analysis")
        print(f"  {Color.GREEN}--skip-amr{Color.RESET}              Skip AMR analysis")
        print(f"  {Color.GREEN}--skip-summary{Color.RESET}          Skip gene‑centric ultimate reporter")
        print(f"  {Color.GREEN}--skip-sample-centric{Color.RESET}   Skip sample‑centric ultimate reporter\n")
        print(f"{Color.BRIGHT_YELLOW}EXAMPLES:{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope -i \"*.fna\" -o results{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope -i \"*.fasta\" -o results --threads 4 --amr-min-identity 0.95 --amr-min-coverage 0.9{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope -i genome.fna -o results --skip-qc --skip-abricate{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope -i \"*.fna\" -o results --abricate-min-id 85 --abricate-min-cov 90{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope --update-amr-db{Color.RESET}   # Manually update AMR database\n")
        print(f"{Color.BRIGHT_YELLOW}Supported FASTA formats:{Color.RESET} {Color.CYAN}.fna, .fasta, .fa, .fn{Color.RESET}")
        print(f"{Color.BRIGHT_YELLOW}Note:{Color.RESET} Run kleboscope --update-amr-db or --force-update-amr-db (Mandatory).")
        print(f"{Color.BRIGHT_YELLOW}RECOMMENDED:{Color.RESET} run abricate --setupdb to ensure ABRicate databases are present.")


def main():
    if '--version' in sys.argv:
        print(f"kleboscope {__version__}")
        sys.exit(0)

    if '-h' in sys.argv or '--help' in sys.argv:
        temp = KleboscopeOrchestrator()
        temp.print_colored_help()
        sys.exit(0)

    parser = argparse.ArgumentParser(
        description="Kleboscope: Complete K. pneumoniae typing pipeline with parallel execution",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-i', '--input', help='Input FASTA file(s) – can use glob patterns like "*.fna" or "*.fasta"')
    parser.add_argument('-o', '--output', help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Number of threads (default: 2)')

    parser.add_argument('--amr-min-identity', type=float, help='Minimum identity for AMR hits (0..1)')
    parser.add_argument('--amr-min-coverage', type=float, help='Minimum coverage for AMR hits (0..1)')
    parser.add_argument('--skip-amr-mutations', action='store_true', help='Disable point mutation reporting in AMR (enabled by default)')

    parser.add_argument('--abricate-min-id', type=float, help='Minimum identity for ABRicate hits (default: 80)')
    parser.add_argument('--abricate-min-cov', type=float, help='Minimum coverage for ABRicate hits (default: 80)')

    parser.add_argument('--skip-qc', action='store_true', help='Skip QC analysis')
    parser.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    parser.add_argument('--skip-kaptive', action='store_true', help='Skip Kaptive analysis')
    parser.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    parser.add_argument('--skip-amr', action='store_true', help='Skip AMR analysis')
    parser.add_argument('--skip-summary', action='store_true', help='Skip gene‑centric ultimate reporter')
    parser.add_argument('--skip-sample-centric', action='store_true', help='Skip sample‑centric ultimate reporter')

    parser.add_argument('--update-amr-db', action='store_true', help='Manually update AMRfinderPlus database incrementally and exit')
    parser.add_argument('--force-update-amr-db', action='store_true', help='Force full AMRfinderPlus database update (overwrites old) and exit')

    parser.add_argument('--keep-temp', action='store_true', help='Do not delete temporary directories (for debugging)')

    args = parser.parse_args()

    if args.update_amr_db or args.force_update_amr_db:
        orchestrator = KleboscopeOrchestrator()
        orchestrator.run_complete_analysis("", "", update_amr_db_only=True, force_update_amr=args.force_update_amr_db)
        sys.exit(0)

    if not args.input or not args.output:
        parser.error("When not using --update-amr-db or --force-update-amr-db, both -i/--input and -o/--output are required.")

    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'kaptive': args.skip_kaptive,
        'abricate': args.skip_abricate,
        'amr': args.skip_amr,
    }

    orchestrator = KleboscopeOrchestrator()
    orchestrator.keep_temp = args.keep_temp
    orchestrator.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        skip_summary=args.skip_summary,
        skip_sample_centric=args.skip_sample_centric,
        amr_min_identity=args.amr_min_identity,
        amr_min_coverage=args.amr_min_coverage,
        amr_skip_mutations=args.skip_amr_mutations,
        abricate_min_id=args.abricate_min_id,
        abricate_min_cov=args.abricate_min_cov
    )


if __name__ == "__main__":
    main()