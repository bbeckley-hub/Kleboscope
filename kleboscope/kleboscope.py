#!/usr/bin/env python3
"""
Kleboscope Main Orchestrator – Parallel Execution with Scientific Quotes
Complete K. pneumoniae typing & resistance pipeline
Author: Brown Beckley <brownbeckley94@gmail.com>
Version: 1.0.1
Date: 2026-04-16
Affiliation: University of Ghana Medical School – Department of Medical Biochemistry
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import random
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

__version__ = "1.0.1"

# =============================================================================
# Color class
# =============================================================================
class Color:
    """ANSI color codes for coloured output"""
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


# =============================================================================
# Kleboscope Orchestrator
# =============================================================================
class KleboscopeOrchestrator:
    """Kleboscope orchestrator with parallel first batch and sequential second batch"""

    def __init__(self):
        self.base_dir = Path(__file__).parent
        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]

        # Output directory names (final names in the user's output folder)
        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'mlst': 'mlst_results',
            'kaptive': 'kaptive_results',
            'abricate': 'klebo_abricate_results',
            'amr': 'klebo_amrfinder_results',
            'summary': 'KLEBOSCOPE_ULTIMATE_REPORTS'
        }

        # HTML files required for the summary module (exact names)
        self.summary_html_files = {
            'mlst_summary.html': 'mlst_results/mlst_summary.html',
            'klebo_fasta_qc_summary.html': 'fasta_qc_results/klebo_fasta_qc_summary.html',
            'klebo_kaptive_summary.html': 'kaptive_results/klebo_kaptive_summary.html',
            'klebo_amrfinder_summary_report.html': 'klebo_amrfinder_results/klebo_amrfinder_summary_report.html',
            'klebo_card_summary_report.html': 'klebo_abricate_results/klebo_card_summary_report.html',
            'klebo_ncbi_summary_report.html': 'klebo_abricate_results/klebo_ncbi_summary_report.html',
            'klebo_resfinder_summary_report.html': 'klebo_abricate_results/klebo_resfinder_summary_report.html',
            'klebo_vfdb_summary_report.html': 'klebo_abricate_results/klebo_vfdb_summary_report.html',
            'klebo_argannot_summary_report.html': 'klebo_abricate_results/klebo_argannot_summary_report.html',
            'klebo_megares_summary_report.html': 'klebo_abricate_results/klebo_megares_summary_report.html',
            'klebo_ecoli_vf_summary_report.html': 'klebo_abricate_results/klebo_ecoli_vf_summary_report.html',
            'klebo_bacmet2_summary_report.html': 'klebo_abricate_results/klebo_bacmet2_summary_report.html',
            'klebo_plasmidfinder_summary_report.html': 'klebo_abricate_results/klebo_plasmidfinder_summary_report.html',
            'klebo_ecoh_summary_report.html': 'klebo_abricate_results/klebo_ecoh_summary_report.html'
        }

    # --------------------------------------------------------------------------
    # Colour setup
    # --------------------------------------------------------------------------
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

    # --------------------------------------------------------------------------
    # Quotes
    # --------------------------------------------------------------------------
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

    # --------------------------------------------------------------------------
    # AMR Database Management (NEW)
    # --------------------------------------------------------------------------
    def update_amr_database(self) -> bool:
        """Run the AMR module's database update and return success."""
        amr_module_path = self.base_dir / "modules" / "kleb_amr_module"
        amr_script = amr_module_path / "klebo_amrfinder.py"
        
        if not amr_script.exists():
            self.print_error(f"AMR script not found at: {amr_script}")
            return False
        
        self.print_info("Updating AMRfinderPlus database...")
        cmd = [sys.executable, str(amr_script), "--update-db"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        
        if result.returncode == 0:
            self.print_success("AMR database updated successfully.")
            # Also show the new version
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

    def ensure_amr_database(self) -> bool:
        """Check if AMR database exists; if not, attempt to update."""
        amr_module_path = self.base_dir / "modules" / "kleb_amr_module"
        amr_script = amr_module_path / "klebo_amrfinder.py"
        if not amr_script.exists():
            self.print_error("AMR script not found, cannot check database.")
            return False
        
        # Run --db-version to see if database is present
        cmd = [sys.executable, str(amr_script), "--db-version"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0 and "Unknown" not in result.stdout and "No database" not in result.stdout and "None" not in result.stdout:
            self.print_success(f"AMR database already present: {result.stdout.strip()}")
            return True
        else:
            self.print_warning("AMR database not found or outdated. Attempting automatic update...")
            return self.update_amr_database()

    # --------------------------------------------------------------------------
    # File discovery
    # --------------------------------------------------------------------------
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

    # --------------------------------------------------------------------------
    # Copy files to module directory
    # --------------------------------------------------------------------------
    def copy_fasta_to_module(self, fasta_files: List[Path], module_path: Path):
        for fasta_file in fasta_files:
            target = module_path / fasta_file.name
            if not target.exists():
                shutil.copy2(fasta_file, target)

    # --------------------------------------------------------------------------
    # Module runners – all use subprocess.run to capture output silently
    # --------------------------------------------------------------------------
    def run_qc(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        qc_module = self.base_dir / "modules" / "kleb_qc_module"
        script = qc_module / "klebo_fasta_qc.py"
        if not script.exists():
            return False, f"Error: QC script not found: {script}"
        self.copy_fasta_to_module(fasta_files, qc_module)
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, str(script), pattern]
        short = f"[INFO] Running QC analysis with pattern: {pattern}\n"
        short += f"  $ python {script.name} {pattern}\n"
        result = subprocess.run(cmd, cwd=qc_module, capture_output=True, text=True)
        if result.returncode != 0:
            short += "⚠️ QC analysis had warnings\n"
            if result.stderr:
                short += f"stderr (first 5 lines):\n{chr(10).join(result.stderr.strip().split(chr(10))[:5])}\n"
        else:
            short += "✓ QC analysis completed!\n"
        # Copy results
        qc_source = qc_module / self.output_dirs['qc']
        qc_target = output_dir / self.output_dirs['qc']
        if qc_source.exists():
            if qc_target.exists():
                shutil.rmtree(qc_target)
            shutil.copytree(qc_source, qc_target)
            short += f"✓ QC results copied to: {qc_target}\n"
        else:
            short += f"⚠️ QC results directory not found: {qc_source} – skipping copy\n"
        return result.returncode == 0, short

    def run_mlst(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        mlst_module = self.base_dir / "modules" / "kleb_mlst_module"
        script = mlst_module / "klebo_mlst.py"
        if not script.exists():
            return False, f"Error: MLST script not found: {script}"
        self.copy_fasta_to_module(fasta_files, mlst_module)
        pattern = self.get_file_pattern(fasta_files)
        output_subdir = self.output_dirs['mlst']
        cmd = [
            sys.executable, str(script),
            "-i", pattern,
            "-o", output_subdir,
            "-db", "db/pubmlst",
            "-sc", ".",
            "--batch"
        ]
        short = f"[INFO] Running MLST analysis\n"
        short += f"  $ python {script.name} -i {pattern} -o {output_subdir} -db db/pubmlst -sc . --batch\n"
        result = subprocess.run(cmd, cwd=mlst_module, capture_output=True, text=True)
        if result.returncode != 0:
            short += "⚠️ MLST analysis had warnings\n"
            if result.stderr:
                short += f"stderr (first 5 lines):\n{chr(10).join(result.stderr.strip().split(chr(10))[:5])}\n"
        else:
            short += "✓ MLST analysis completed!\n"
        mlst_source = mlst_module / output_subdir
        mlst_target = output_dir / output_subdir
        if mlst_source.exists():
            if mlst_target.exists():
                shutil.rmtree(mlst_target)
            shutil.copytree(mlst_source, mlst_target)
            short += f"✓ MLST results copied to: {mlst_target}\n"
        else:
            short += f"⚠️ MLST results directory not found: {mlst_source} – skipping copy\n"
        return result.returncode == 0, short

    def run_kaptive(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        kaptive_module = self.base_dir / "modules" / "kleb_serotype_module"
        script = kaptive_module / "klebo_kaptive.py"
        if not script.exists():
            return False, f"Error: Kaptive script not found: {script}"
        self.copy_fasta_to_module(fasta_files, kaptive_module)
        pattern = self.get_file_pattern(fasta_files)
        output_subdir = self.output_dirs['kaptive']
        cmd = [sys.executable, str(script), "-i", pattern, "-o", output_subdir]
        short = f"[INFO] Running Kaptive analysis\n"
        short += f"  $ python {script.name} -i {pattern} -o {output_subdir}\n"
        result = subprocess.run(cmd, cwd=kaptive_module, capture_output=True, text=True)
        if result.returncode != 0:
            short += "⚠️ Kaptive analysis had warnings\n"
            if result.stderr:
                short += f"stderr (first 5 lines):\n{chr(10).join(result.stderr.strip().split(chr(10))[:5])}\n"
        else:
            short += "✓ Kaptive analysis completed!\n"
        kaptive_source = kaptive_module / output_subdir
        kaptive_target = output_dir / output_subdir
        if kaptive_source.exists():
            if kaptive_target.exists():
                shutil.rmtree(kaptive_target)
            shutil.copytree(kaptive_source, kaptive_target)
            short += f"✓ Kaptive results copied to: {kaptive_target}\n"
        else:
            short += f"⚠️ Kaptive results directory not found: {kaptive_source} – skipping copy\n"
        return result.returncode == 0, short

    def run_abricate(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        abricate_module = self.base_dir / "modules" / "kleb_abricate_module"
        script = abricate_module / "kleb_abricate_module.py"
        if not script.exists():
            return False, f"Error: ABRicate script not found: {script}"
        self.copy_fasta_to_module(fasta_files, abricate_module)
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, str(script), pattern]
        short = f"[INFO] Running ABRicate analysis\n"
        short += f"  $ python {script.name} {pattern}\n"
        result = subprocess.run(cmd, cwd=abricate_module, capture_output=True, text=True)
        if result.returncode != 0:
            short += "⚠️ ABRicate analysis had warnings\n"
            if result.stderr:
                short += f"stderr (first 5 lines):\n{chr(10).join(result.stderr.strip().split(chr(10))[:5])}\n"
        else:
            short += "✓ ABRicate analysis completed!\n"
        abricate_source = abricate_module / "abricate_results"
        abricate_target = output_dir / self.output_dirs['abricate']
        if abricate_source.exists():
            if abricate_target.exists():
                shutil.rmtree(abricate_target)
            shutil.copytree(abricate_source, abricate_target)
            short += f"✓ ABRicate results copied to: {abricate_target}\n"
        else:
            short += f"⚠️ ABRicate results directory not found: {abricate_source} – skipping copy\n"
        return result.returncode == 0, short

    def run_amr(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        amr_module = self.base_dir / "modules" / "kleb_amr_module"
        script = amr_module / "klebo_amrfinder.py"
        if not script.exists():
            return False, f"Error: AMR script not found: {script}"
        
        # Ensure database is present before running
        if not self.ensure_amr_database():
            return False, "AMR database is missing and could not be updated automatically.\nPlease run: kleboscope --update-amr-db"
        
        self.copy_fasta_to_module(fasta_files, amr_module)
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, str(script), pattern]
        short = f"[INFO] Running AMRfinder analysis\n"
        short += f"  $ python {script.name} {pattern}\n"
        result = subprocess.run(cmd, cwd=amr_module, capture_output=True, text=True)
        if result.returncode != 0:
            short += "⚠️ AMR analysis had warnings\n"
            if result.stderr:
                short += f"stderr (first 5 lines):\n{chr(10).join(result.stderr.strip().split(chr(10))[:5])}\n"
        else:
            short += "✓ AMR analysis completed!\n"
        amr_source = amr_module / "klebo_amrfinder_results"
        amr_target = output_dir / self.output_dirs['amr']
        if amr_source.exists():
            if amr_target.exists():
                shutil.rmtree(amr_target)
            shutil.copytree(amr_source, amr_target)
            short += f"✓ AMR results copied to: {amr_target}\n"
        else:
            short += f"⚠️ AMR results directory not found: {amr_source} – skipping copy\n"
        return result.returncode == 0, short

    def run_summary(self, output_dir: Path) -> Tuple[bool, str]:
        summary_module = self.base_dir / "modules" / "kleb_summary_module"
        script = summary_module / "kleboscope_ultimate_reporter.py"
        if not script.exists():
            return False, f"Error: Summary script not found: {script}"
        output = "Copying required HTML files to summary module...\n"
        copied = 0
        missing = 0
        for target_name, source_rel in self.summary_html_files.items():
            source = output_dir / source_rel
            if source.exists():
                shutil.copy2(source, summary_module / target_name)
                copied += 1
                output += f"  ✓ {target_name}\n"
            else:
                output += f"  ✗ {target_name} (missing at {source_rel})\n"
                missing += 1
        if missing > 0:
            output += f"⚠️ Copied {copied} files, {missing} missing. Some analysis may be incomplete.\n"
        output += "[INFO] Running ultimate reporter...\n"
        cmd = [sys.executable, str(script), "-i", "."]
        output += f"  $ python {script.name} -i .\n"
        result = subprocess.run(cmd, cwd=summary_module, capture_output=True, text=True)
        output += result.stdout
        if result.stderr:
            output += "\n" + result.stderr
        if result.returncode == 0:
            output += "✓ Ultimate reporter completed!\n"
        else:
            output += "⚠️ Ultimate reporter had warnings\n"
        summary_source = summary_module / self.output_dirs['summary']
        summary_target = output_dir / self.output_dirs['summary']
        if summary_source.exists():
            if summary_target.exists():
                shutil.rmtree(summary_target)
            shutil.copytree(summary_source, summary_target)
            output += f"✓ Ultimate reports copied to: {summary_target}\n"
            files = list(summary_target.glob("*"))
            html_count = len([f for f in files if f.suffix == '.html'])
            json_count = len([f for f in files if f.suffix == '.json'])
            csv_count = len([f for f in files if f.suffix == '.csv'])
            output += f"  📊 {html_count} HTML, {json_count} JSON, {csv_count} CSV files\n"
        else:
            output += f"⚠️ Ultimate reports directory not found: {summary_source}\n"
        return result.returncode == 0, output

    # --------------------------------------------------------------------------
    # Cleanup
    # --------------------------------------------------------------------------
    def cleanup_module(self, module_path: Path, fasta_files: List[Path]):
        try:
            for fasta_file in fasta_files:
                temp = module_path / fasta_file.name
                if temp.exists():
                    temp.unlink()
            dirs_to_remove = [
                self.output_dirs['qc'],
                self.output_dirs['mlst'],
                self.output_dirs['kaptive'],
                "abricate_results",
                "klebo_amrfinder_results"
            ]
            for dir_name in dirs_to_remove:
                dir_path = module_path / dir_name
                if dir_path.exists():
                    shutil.rmtree(dir_path)
            for html in module_path.glob("*.html"):
                html.unlink()
        except Exception as e:
            self.print_warning(f"Partial cleanup issue in {module_path.name}: {e}")

    # --------------------------------------------------------------------------
    # Main execution
    # --------------------------------------------------------------------------
    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                              skip_modules: Dict[str, bool] = None, skip_summary: bool = False,
                              update_amr_db_only: bool = False):
        # If only updating AMR database, do that and exit
        if update_amr_db_only:
            self.update_amr_database()
            return
        
        if skip_modules is None:
            skip_modules = {}
        start_time = datetime.now()
        self.display_banner()
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found! Analysis stopped.")
            return
        self.print_success(f"Starting analysis of {len(fasta_files)} K. pneumoniae samples")

        # Create output subdirectories
        for subdir in self.output_dirs.values():
            (output_path / subdir).mkdir(exist_ok=True)

        # Display plan
        self.print_header("ANALYSIS PLAN", "Modules to be executed")
        plan = [
            ("QC Analysis", not skip_modules.get('qc', False)),
            ("MLST Analysis", not skip_modules.get('mlst', False)),
            ("Kaptive Analysis", not skip_modules.get('kaptive', False)),
            ("ABRicate Analysis", not skip_modules.get('abricate', False)),
            ("AMR Analysis", not skip_modules.get('amr', False)),
            ("Ultimate Reporter", not skip_summary),
        ]
        for analysis, enabled in plan:
            if enabled:
                print(f"   {Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET} - {analysis}")
            else:
                print(f"   {Color.YELLOW}⏸️  SKIPPED{Color.RESET} - {analysis}")
        print()

        # =====================================================================
        # PARALLEL FIRST BATCH: QC, MLST, Kaptive
        # =====================================================================
        tasks = []
        if not skip_modules.get('qc', False):
            tasks.append(("QC", self.run_qc))
        if not skip_modules.get('mlst', False):
            tasks.append(("MLST", self.run_mlst))
        if not skip_modules.get('kaptive', False):
            tasks.append(("Kaptive", self.run_kaptive))

        if tasks:
            self.print_info(f"Running {len(tasks)} analyses in parallel...")
            results = {}
            with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
                future_to_name = {executor.submit(run, fasta_files, output_path, threads): name for name, run in tasks}
                for future in as_completed(future_to_name):
                    name = future_to_name[future]
                    try:
                        success, output = future.result()
                        results[name] = (success, output)
                    except Exception as e:
                        results[name] = (False, f"Exception: {e}")

            # Print results in order
            for name, _ in tasks:
                success, output = results.get(name, (False, "No result"))
                self.print_header(f"{name} Analysis")
                print(output.strip())
                self.display_random_quote()
                if success:
                    self.print_success(f"✅ {name} completed")
                else:
                    self.print_error(f"❌ {name} failed")
        else:
            self.print_info("No analyses in first batch (all skipped).")

        # =====================================================================
        # SEQUENTIAL SECOND BATCH: ABRicate, then AMR (with short output)
        # =====================================================================
        if not skip_modules.get('abricate', False):
            self.print_header("ABRICATE ANALYSIS", "Comprehensive Resistance & Virulence Screening")
            success, output = self.run_abricate(fasta_files, output_path, threads)
            print(output.strip())
            self.display_random_quote()
            if success:
                self.print_success("✅ ABRicate completed")
            else:
                self.print_error("❌ ABRicate failed")
            self.print_info("ABRicate analysis completed.")
        else:
            self.print_info("Skipping ABRicate analysis.")

        if not skip_modules.get('amr', False):
            self.print_header("AMR ANALYSIS", "Antimicrobial Resistance Gene Detection")
            success, output = self.run_amr(fasta_files, output_path, threads)
            print(output.strip())
            self.display_random_quote()
            if success:
                self.print_success("✅ AMR completed")
            else:
                self.print_error("❌ AMR failed")
            self.print_info("AMR analysis completed.")
        else:
            self.print_info("Skipping AMR analysis.")

        # =====================================================================
        # FINAL SUMMARY
        # =====================================================================
        if not skip_summary:
            self.print_info("Copying files to summary module and running ultimate reporter...")
            summary_module = self.base_dir / "modules" / "kleb_summary_module"
            for html in summary_module.glob("*.html"):
                html.unlink()
            self.print_header("ULTIMATE REPORTER", "Gene‑centric Integrated Analysis")
            success, output = self.run_summary(output_path)
            print(output.strip())
            self.display_random_quote()
            if not success:
                self.print_warning("Ultimate reporter had issues")

        # Clean up modules
        for module_dir in ["kleb_qc_module", "kleb_mlst_module", "kleb_serotype_module", "kleb_abricate_module", "kleb_amr_module"]:
            module_path = self.base_dir / "modules" / module_dir
            if module_path.exists():
                self.cleanup_module(module_path, fasta_files)

        # Final summary
        analysis_time = datetime.now() - start_time
        self.print_header("ANALYSIS COMPLETE", f"Time elapsed: {str(analysis_time).split('.')[0]}")
        self.print_success(f"🎉 Analysis complete! Results in: {output_path}")

        for subdir in sorted(output_path.iterdir()):
            if subdir.is_dir():
                file_count = len(list(subdir.glob("*")))
                self.print_info(f"  📁 {subdir.name} ({file_count} files)")

        self.display_random_quote()

    # --------------------------------------------------------------------------
    # Banner and colored help
    # --------------------------------------------------------------------------
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
{'='*80}{Color.RESET}
"""
        print(banner)

    def print_colored_help(self):
        self.display_banner()
        print(f"{Color.BRIGHT_YELLOW}USAGE:{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope{Color.RESET} {Color.CYAN}-i INPUT -o OUTPUT{Color.RESET} [OPTIONS]")
        print(f"  {Color.GREEN}kleboscope --update-amr-db{Color.RESET}  (Update AMR database only)")
        print()
        print(f"{Color.BRIGHT_YELLOW}REQUIRED ARGUMENTS (for analysis):{Color.RESET}")
        print(f"  {Color.GREEN}-i, --input{Color.RESET} INPUT    Input FASTA file(s) – supports glob patterns like \"*.fna\"")
        print(f"  {Color.GREEN}-o, --output{Color.RESET} OUTPUT  Output directory for results\n")
        print(f"{Color.BRIGHT_YELLOW}OPTIONAL ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-h, --help{Color.RESET}            Show this help message")
        print(f"  {Color.GREEN}-t, --threads{Color.RESET} THREADS  Number of threads (default: 2)")
        print(f"  {Color.GREEN}--version{Color.RESET}               Show version and exit")
        print(f"  {Color.GREEN}--update-amr-db{Color.RESET}         Update AMRfinderPlus database and exit (first‑time users)\n")
        print(f"{Color.BRIGHT_YELLOW}SKIP OPTIONS:{Color.RESET}")
        print(f"  {Color.GREEN}--skip-qc{Color.RESET}               Skip QC analysis")
        print(f"  {Color.GREEN}--skip-mlst{Color.RESET}             Skip MLST analysis")
        print(f"  {Color.GREEN}--skip-kaptive{Color.RESET}          Skip Kaptive analysis")
        print(f"  {Color.GREEN}--skip-abricate{Color.RESET}         Skip ABRicate analysis")
        print(f"  {Color.GREEN}--skip-amr{Color.RESET}              Skip AMR analysis")
        print(f"  {Color.GREEN}--skip-summary{Color.RESET}          Skip ultimate reporter generation\n")
        print(f"{Color.BRIGHT_YELLOW}EXAMPLES:{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope -i \"*.fna\" -o results{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope -i \"*.fasta\" -o results --threads 4{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope -i genome.fna -o results --skip-qc --skip-abricate{Color.RESET}")
        print(f"  {Color.GREEN}kleboscope --update-amr-db{Color.RESET}   # Mandatory first step for AMR analysis\n")
        print(f"{Color.BRIGHT_YELLOW}Supported FASTA formats:{Color.RESET} {Color.CYAN}.fna, .fasta, .fa, .fn{Color.RESET}")


# =============================================================================
# Main entry point
# =============================================================================
def main():
    # Handle --version early
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
    parser.add_argument('--skip-qc', action='store_true', help='Skip QC analysis')
    parser.add_argument('--skip-mlst', action='store_true', help='Skip MLST analysis')
    parser.add_argument('--skip-kaptive', action='store_true', help='Skip Kaptive analysis')
    parser.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate analysis')
    parser.add_argument('--skip-amr', action='store_true', help='Skip AMR analysis')
    parser.add_argument('--skip-summary', action='store_true', help='Skip ultimate reporter generation')
    parser.add_argument('--update-amr-db', action='store_true', help='Update AMRfinderPlus database and exit (mandatory before first AMR run)')

    args = parser.parse_args()

    # If only updating the AMR database, no need for input/output
    if args.update_amr_db:
        orchestrator = KleboscopeOrchestrator()
        orchestrator.run_complete_analysis("", "", update_amr_db_only=True)
        sys.exit(0)

    # Otherwise, input and output are required
    if not args.input or not args.output:
        parser.error("When not using --update-amr-db, both -i/--input and -o/--output are required.")

    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'kaptive': args.skip_kaptive,
        'abricate': args.skip_abricate,
        'amr': args.skip_amr,
    }

    orchestrator = KleboscopeOrchestrator()
    try:
        orchestrator.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip_modules=skip_modules,
            skip_summary=args.skip_summary
        )
    except KeyboardInterrupt:
        print(f"\n{Color.BRIGHT_RED}❌ Analysis interrupted by user{Color.RESET}")
    except Exception as e:
        print(f"\n{Color.BRIGHT_RED}💥 Critical error: {e}{Color.RESET}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()