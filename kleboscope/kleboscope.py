#!/usr/bin/env python3
"""
Kleboscope Main Orchestrator – Parallel Execution with Scientific Quotes
Complete K. pneumoniae typing & resistance pipeline
Author: Brown Beckley <brownbeckley94@gmail.com>
Version: 2.0.0 
Affiliation: University of Ghana Medical School – Department of Medical Biochemistry
"""

import os
import sys
import glob
import argparse
import subprocess
#import shutil
import random
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

__version__ = "2.0.0"

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
        self.base_dir = Path(__file__).parent.resolve()
        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]

        # Output directory names 
        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'mlst': 'mlst_results',
            'kaptive': 'kaptive_results',
            'abricate': 'klebo_abricate_results',
            'amr': 'klebo_amrfinder_results',
            'summary': 'KLEBOSCOPE_ULTIMATE_REPORTS'
        }

        # Centralized configurations to reduce repetition and improve maintainability
        # {pattern} = Absolute path glob, {output} = Absolute output dir, {threads} = Allocated CPU cores
        self.module_configs = {
            'qc': {
                'dir': 'kleb_qc_module', 
                'script': 'klebo_fasta_qc.py', 
                'args': ['-i', '{pattern}', '-o', '{output}', '-t', '{threads}']
            },
            'mlst': {
                'dir': 'kleb_mlst_module', 
                'script': 'klebo_mlst.py', 
                'args': ['-i', '{pattern}', '-o', '{output}', '-db', 'db/pubmlst', '-sc', '.', '--batch']
            },
            'kaptive': {
                'dir': 'kleb_serotype_module', 
                'script': 'klebo_kaptive.py', 
                'args': ['-i', '{pattern}', '-o', '{output}']
            },
            'abricate': {
                'dir': 'kleb_abricate_module', 
                'script': 'kleb_abricate_module.py', 
                'args': ['-i', '{pattern}', '-o', '{output}', '-t', '{threads}']
            },
            'amr': {
                'dir': 'kleb_amr_module', 
                'script': 'klebo_amrfinder.py', 
                'args': ['{pattern}', '-o', '{output}', '--cpus', '{threads}']
            }
        }

    # --------------------------------------------------------------------------
    # Colour setup & Output Logging
    # --------------------------------------------------------------------------
    def setup_colors(self):
        self.color_info = Color.CYAN
        self.color_success = Color.BRIGHT_GREEN
        self.color_warning = Color.BRIGHT_YELLOW
        self.color_error = Color.BRIGHT_RED

    def print_header(self, title: str, subtitle: str = ""):
        print(f"\n{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * 20}{title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}{' ' * 22}{subtitle}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}\n")

    def print_info(self, message: str): print(f"{self.color_info}[INFO]{Color.RESET} {message}")
    def print_success(self, message: str): print(f"{self.color_success}✓{Color.RESET} {message}")
    def print_warning(self, message: str): print(f"{self.color_warning}⚠️{Color.RESET} {message}")
    def print_error(self, message: str): print(f"{self.color_error}✗{Color.RESET} {message}")

    # --------------------------------------------------------------------------
    # Quotes
    # --------------------------------------------------------------------------
    def _get_scientific_quotes(self):
        return [
            {"quote": "Science is organised knowledge.", "author": "Herbert Spencer", "theme": "knowledge"},
            {"quote": "The science of today is the technology of tomorrow.", "author": "Edward Teller", "theme": "technology"},
            {"quote": "Biology is the most powerful technology ever created.", "author": "Freeman Dyson", "theme": "biology"},
            {"quote": "Genomics is a lens on biology.", "author": "Eric Lander", "theme": "genomics"},
            {"quote": "Resistance is not futile.", "author": "Antibiotic Researcher", "theme": "resistance"},
            {"quote": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson", "theme": "science"}
        ]

    def display_random_quote(self):
        if not self.quotes: return
        quote_data = random.choice(self.quotes)
        quote_color = random.choice(self.quote_colors)
        print(f"\n{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] 🧪 SCIENTIFIC INSIGHT: {Color.RESET}\n")
        print(f"{quote_color}   \"{quote_data['quote']}\"{Color.RESET}")
        print(f"{Color.BOLD}{Color.WHITE}   — {quote_data['author']}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}\n")

    # --------------------------------------------------------------------------
    # Absolute File Discovery
    # --------------------------------------------------------------------------
    # help reducing  I/O errors and path issues by ensuring all file paths are absolute and resolved at the start of the analysis
    def find_fasta_files(self, input_path: str) -> List[Path]:
        self.print_info(f"Searching for files with pattern: {input_path}")
        
        # Resolve to absolute paths immediately
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f).resolve() for f in matched_files if Path(f).is_file() and
                           f.lower().endswith(('.fna', '.fasta', '.fa', '.fn'))]
        else:
            input_path_obj = Path(input_path).resolve()
            if input_path_obj.is_file():
                fasta_files = [input_path_obj]
            elif input_path_obj.is_dir():
                fasta_files = [Path(f).resolve() for f in glob.glob(f"{input_path_obj}/*") 
                               if Path(f).is_file() and str(f).lower().endswith(('.fna', '.fasta', '.fa', '.fn'))]
            else:
                fasta_files = []

        if fasta_files:
            self.print_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        
        self.print_error(f"Input path not found or empty: {input_path}")
        return []

    def get_absolute_file_pattern(self, fasta_files: List[Path]) -> str:
        """Returns a glob string based on the absolute path of the directory."""
        if not fasta_files:
            return "*.fna"
        # fix to handle one input
        if len(fasta_files) == 1:
            return str(fasta_files[0])
        parent_dir = fasta_files[0].parent
        extensions = set(f.suffix.lower() for f in fasta_files)
        ext = f"*{list(extensions)[0]}" if len(extensions) == 1 else "*"
        return str(parent_dir / ext)

    # --------------------------------------------------------------------------
    # Module Runner
    # --------------------------------------------------------------------------
    def _execute_module(self, module_key: str, display_name: str, fasta_files: List[Path], base_output_dir: Path, threads: int) -> Tuple[bool, str]:
        config = self.module_configs[module_key]
        module_path = self.base_dir / "modules" / config['dir']
        script = module_path / config['script']
        
        if not script.exists():
            return False, f"Error: Script not found: {script}"
            
        output_dir = base_output_dir / self.output_dirs[module_key]
        output_dir.mkdir(parents=True, exist_ok=True)
        
        abs_pattern = self.get_absolute_file_pattern(fasta_files)
        
        cmd = [sys.executable, str(script)]
        for arg in config['args']:
            cmd.append(arg.format(
                pattern=abs_pattern,
                output=str(output_dir),
                threads=threads
            ))
            
        short = f"[INFO] Running {display_name} analysis\n"
        short += f"  $ {' '.join(cmd)}\n"
        
        result = subprocess.run(cmd, cwd=module_path, capture_output=True, text=True)
        
        if result.returncode != 0:
            short += f"⚠️ {display_name} analysis had warnings or errors\n"
            if result.stderr:
                short += f"stderr (first 5 lines):\n{chr(10).join(result.stderr.strip().split(chr(10))[:5])}\n"
        else:
            short += f"✓ {display_name} analysis completed successfully in {output_dir}!\n"
            
        return result.returncode == 0, short

    def run_summary(self, output_dir: Path) -> Tuple[bool, str]:
        """Runs the ultimate reporter, assuming it pulls from the defined output_dir"""
        summary_module = self.base_dir / "modules" / "kleb_summary_module"
        script = summary_module / "kleboscope_ultimate_reporter.py"
        
        if not script.exists():
            return False, f"Error: Summary script not found: {script}"
            
        final_summary_dir = output_dir / self.output_dirs['summary']
        final_summary_dir.mkdir(parents=True, exist_ok=True)
        
        output = "[INFO] Running ultimate reporter...\n"
        cmd = [sys.executable, str(script), "-i", str(output_dir), "-o", str(final_summary_dir)]
        output += f"  $ {' '.join(cmd)}\n"
        
        result = subprocess.run(cmd, cwd=summary_module, capture_output=True, text=True)
        output += result.stdout
        
        if result.stderr: output += "\n" + result.stderr
        
        if result.returncode == 0:
            output += f"✓ Ultimate reports generated at: {final_summary_dir}\n"
        else:
            output += "⚠️ Ultimate reporter had warnings\n"
            
        return result.returncode == 0, output

    # --------------------------------------------------------------------------
    # Main execution
    # --------------------------------------------------------------------------
    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                              skip_modules: Dict[str, bool] = None, skip_summary: bool = False):
        if skip_modules is None: skip_modules = {}
        start_time = datetime.now()
        self.display_banner()
        
        output_path = Path(output_dir).resolve()
        output_path.mkdir(parents=True, exist_ok=True)

        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found! Analysis stopped.")
            return

        # Display plan
        self.print_header("ANALYSIS PLAN", "Modules to be executed")
        plan = [
            ("qc", "QC", not skip_modules.get('qc', False)),
            ("mlst", "MLST", not skip_modules.get('mlst', False)),
            ("kaptive", "Kaptive", not skip_modules.get('kaptive', False)),
            ("abricate", "ABRicate", not skip_modules.get('abricate', False)),
            ("amr", "AMR", not skip_modules.get('amr', False)),
            ("summary", "Ultimate Reporter", not skip_summary),
        ]
        
        for key, name, enabled in plan:
            status = f"{Color.BRIGHT_GREEN}✅ ENABLED" if enabled else f"{Color.YELLOW}⏸️  SKIPPED"
            print(f"   {status}{Color.RESET} - {name} Analysis")
        print()

        # =====================================================================
        # PARALLEL FIRST BATCH: QC, MLST, Kaptive
        # =====================================================================
        tasks = [(k, n) for k, n, e in plan[:3] if e]
        
        if tasks:
            # Prevent system thrashing by dividing threads among tasks
            threads_per_task = max(1, threads // len(tasks))
            self.print_info(f"Running {len(tasks)} analyses in parallel (allocating {threads_per_task} threads per task)...")
            
            results = {}
            with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
                future_to_name = {
                    executor.submit(self._execute_module, key, name, fasta_files, output_path, threads_per_task): name 
                    for key, name in tasks
                }
                for future in as_completed(future_to_name):
                    name = future_to_name[future]
                    try:
                        results[name] = future.result()
                    except Exception as e:
                        results[name] = (False, f"Exception: {e}")

            # Print results
            for _, name in tasks:
                success, output = results.get(name, (False, "No result"))
                self.print_header(f"{name} Analysis")
                print(output.strip())
                self.display_random_quote()
                if success: self.print_success(f"✅ {name} completed")
                else: self.print_error(f"❌ {name} failed")
        else:
            self.print_info("No analyses in first batch (all skipped).")

        # =====================================================================
        # SEQUENTIAL SECOND BATCH: ABRicate, then AMR
        # =====================================================================
        sequential_tasks = [(k, n) for k, n, e in plan[3:5] if e]
        
        for key, name in sequential_tasks:
            self.print_header(f"{name.upper()} ANALYSIS")
            # Give full thread pool to sequential tasks
            success, output = self._execute_module(key, name, fasta_files, output_path, threads)
            print(output.strip())
            self.display_random_quote()
            if success: self.print_success(f"✅ {name} completed")
            else: self.print_error(f"❌ {name} failed")

        # =====================================================================
        # FINAL SUMMARY
        # =====================================================================
        if not skip_summary:
            self.print_header("ULTIMATE REPORTER", "Gene‑centric Integrated Analysis")
            success, output = self.run_summary(output_path)
            print(output.strip())
            self.display_random_quote()

        # Final summary
        analysis_time = datetime.now() - start_time
        self.print_header("ANALYSIS COMPLETE", f"Time elapsed: {str(analysis_time).split('.')[0]}")
        self.print_success(f"🎉 Analysis complete! Results in: {output_path}")

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
        print(f"  {Color.GREEN}kleboscope{Color.RESET} {Color.CYAN}-i INPUT -o OUTPUT{Color.RESET} [OPTIONS]\n")
        print(f"{Color.BRIGHT_YELLOW}REQUIRED ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-i, --input{Color.RESET} INPUT    Input FASTA file(s) – supports glob patterns like \"*.fna\"")
        print(f"  {Color.GREEN}-o, --output{Color.RESET} OUTPUT  Output directory for results\n")
        print(f"{Color.BRIGHT_YELLOW}OPTIONAL ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-h, --help{Color.RESET}            Show this help message")
        print(f"  {Color.GREEN}-t, --threads{Color.RESET} THREADS  Number of threads (default: 2)")
        print(f"  {Color.GREEN}--version{Color.RESET}               Show version and exit\n")
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
        print(f"  {Color.GREEN}kleboscope -i genome.fna -o results --skip-qc --skip-abricate{Color.RESET}\n")
        print(f"{Color.BRIGHT_YELLOW}Supported FASTA formats:{Color.RESET} {Color.CYAN}.fna, .fasta, .fa, .fn{Color.RESET}")
# =============================================================================
# Main entry point
# =============================================================================
def main():
    if '--version' in sys.argv:
        print(f"kleboscope {__version__}")
        sys.exit(0)

    if '-h' in sys.argv or '--help' in sys.argv:
        KleboscopeOrchestrator().print_colored_help()
        sys.exit(0)

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-t', '--threads', type=int, default=2)
    parser.add_argument('--skip-qc', action='store_true')
    parser.add_argument('--skip-mlst', action='store_true')
    parser.add_argument('--skip-kaptive', action='store_true')
    parser.add_argument('--skip-abricate', action='store_true')
    parser.add_argument('--skip-amr', action='store_true')
    parser.add_argument('--skip-summary', action='store_true')

    args = parser.parse_args()

    skip_modules = {
        'qc': args.skip_qc, 'mlst': args.skip_mlst, 'kaptive': args.skip_kaptive,
        'abricate': args.skip_abricate, 'amr': args.skip_amr,
    }

    orchestrator = KleboscopeOrchestrator()
    try:
        orchestrator.run_complete_analysis(args.input, args.output, args.threads, skip_modules, args.skip_summary)
    except KeyboardInterrupt:
        print(f"\n{Color.BRIGHT_RED}❌ Analysis interrupted by user{Color.RESET}")
    except Exception as e:
        print(f"\n{Color.BRIGHT_RED}💥 Critical error: {e}{Color.RESET}")
        sys.exit(1)

if __name__ == "__main__":
    main()