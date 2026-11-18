#!/usr/bin/env python3
"""
Kleboscope Main Orchestrator - Complete K. pneumoniae Typing Pipeline
Version: 1.3.0

Runs all analysis and reporting modules in isolated temporary directories and
consolidates the final deliverables into Kleboscope_final_report/. Supports
HPC, Docker, Singularity, and conda installations via environment overrides
and user-local database caches.

Modules:
    analysis  : qc_module, mlst_module, kaptive_module, abricate_module,
                amr_module, mge_module
    reporting : gene_centric_module, sample_centric_module

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import tempfile
import logging
import traceback
import signal
import atexit
import random
import urllib.request
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Set, Tuple

__version__ = "1.3.0"

AUTHOR_NAME = "Brown Beckley"
AUTHOR_EMAIL = "brownbeckley94@gmail.com"
AUTHOR_AFFILIATION = (
    "University of Ghana Medical School - Department of Medical Biochemistry"
)


class Color:
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'


class ColoredHelpFormatter(argparse.HelpFormatter):
    """ANSI-aware help formatter for Kleboscope.

    The parser description and epilog are pre-formatted strings containing
    ANSI escape sequences, so they are emitted verbatim instead of being
    re-wrapped (textwrap counts invisible escape bytes toward line width,
    which scattered the text across wrong columns). Argument help strings
    contain no escape sequences, so they keep the standard argparse
    wrapping. Option names are colored only after padding is applied,
    keeping the option/description columns aligned.
    """

    def _fill_text(self, text, width, indent):
        return ''.join(indent + line + '\n' for line in text.splitlines())

    def _format_usage(self, usage, actions, groups, prefix):
        text = super()._format_usage(usage, actions, groups, prefix)
        if text.startswith('usage: '):
            text = text[len('usage: '):]
            return f"{Color.BOLD}{Color.GREEN}usage:{Color.RESET} {text}"
        return text

    def start_section(self, heading):
        super().start_section(
            f"{Color.BOLD}{Color.BRIGHT_GREEN}{heading}{Color.RESET}")

    def _format_action(self, action):
        help_position = min(self._action_max_length + 2,
                            self._max_help_position)
        help_width = max(self._width - help_position, 11)
        action_width = help_position - self._current_indent - 2
        header = self._format_action_invocation(action)
        colored_header = f"{Color.BRIGHT_CYAN}{header}{Color.RESET}"
        indent = ' ' * self._current_indent

        if not action.help:
            return f"{indent}{colored_header}\n"

        help_lines = self._split_lines(self._expand_help(action), help_width)
        if len(header) <= action_width:
            pad = action_width - len(header) + 2
            out = [f"{indent}{colored_header}{' ' * pad}{help_lines[0]}"]
        else:
            out = [f"{indent}{colored_header}",
                   f"{' ' * help_position}{help_lines[0]}"]
        out.extend(f"{' ' * help_position}{line}" for line in help_lines[1:])
        for subaction in self._iter_indented_subactions(action):
            out.append(self._format_action(subaction))
        return '\n'.join(out) + '\n'


DIR_QC = "fasta_qc_results"
DIR_MLST = "mlst_results"
DIR_KAPTIVE = "kaptive_results"
DIR_ABRICATE = "abricate_results"
DIR_AMR = "amrfinder_results"
DIR_MGE = "mge_results"
DIR_GENE_CENTRIC = "KLEBOSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS"
DIR_SAMPLE_CENTRIC = "KLEBOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS"
DIR_FINAL = "Kleboscope_final_report"

ALL_MODULE_DIRS = [DIR_QC, DIR_MLST, DIR_KAPTIVE, DIR_ABRICATE, DIR_AMR, DIR_MGE]

CITATION_LINES = [
    "Beckley B, et al. Kleboscope: a gene-centric, species-optimized",
    "computational pipeline for Klebsiella pneumoniae genomic surveillance.",
    "https://github.com/bbeckley-hub/kleboscope",
]


class KleboscopeOrchestrator:
    """Coordinates all Kleboscope analysis and reporting modules.

    Runs analysis modules in isolated temporary directories, collects
    batch summaries into a consolidated final report, and prunes
    intermediate files according to the selected report level.

    Full subprocess stdout/stderr is captured into the run log. On
    non-zero exit codes the complete output is also echoed to the
    console so failures are actionable without opening the log.
    """

    def __init__(self, quiet: bool = False, keep_temp: bool = False,
                 tmpdir: Optional[str] = None):
        self.base_dir = Path(__file__).parent
        self.quiet = quiet
        self.keep_temp = keep_temp
        self.logger: Optional[logging.Logger] = None
        self.user_output_dir: Optional[Path] = None
        self.temp_dirs: List[str] = []
        self.module_times: Dict[str, float] = {}
        self.running = True
        self.report_level = 'summary'
        self.warnings: List[str] = []
        self.errors: List[str] = []

        if tmpdir:
            self.tmp_root = Path(tmpdir)
            self.tmp_root.mkdir(parents=True, exist_ok=True)
        else:
            self.tmp_root = None

        self.quotes = self._get_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE,
        ]
        self.module_subtitles = {
            'qc': 'Sequence Quality Control and Species Confirmation',
            'mlst': 'Multi-Locus Sequence Typing',
            'kaptive': 'K/O Locus Typing',
            'abricate': 'Comprehensive Resistance and Virulence Screening',
            'amr': 'AMRFinderPlus Gene Detection with Point Mutations',
            'mge': 'Mobile Genetic Element Profiling',
            'gene_centric': 'Gene-Centric Integrated Report',
            'sample_centric': 'Isolate-Centric Interactive Report',
        }
        self.output_dirs = {
            'qc': DIR_QC,
            'mlst': DIR_MLST,
            'kaptive': DIR_KAPTIVE,
            'abricate': DIR_ABRICATE,
            'amr': DIR_AMR,
            'mge': DIR_MGE,
        }

        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
        atexit.register(self._cleanup_temp_dirs)

    # ------------------------------------------------------------------
    # Signal handling and temp management
    # ------------------------------------------------------------------
    def _signal_handler(self, signum, frame):
        self.print_warning(f"Received signal {signum}. Cleaning up and exiting...")
        self.running = False
        self._cleanup_temp_dirs()
        sys.exit(1)

    def _cleanup_temp_dirs(self):
        if self.keep_temp:
            return
        for d in list(self.temp_dirs):
            if d and Path(d).exists():
                shutil.rmtree(d, ignore_errors=True)
                if self.logger:
                    self.logger.info(f"Removed temporary directory: {d}")
        self.temp_dirs = []

    def _mkdtemp(self, prefix: str) -> str:
        """Create a temp dir honouring --tmpdir, KLEBOSCOPE_TMPDIR, TMPDIR."""
        if self.tmp_root:
            return tempfile.mkdtemp(prefix=prefix, dir=str(self.tmp_root))
        env_tmp = os.environ.get("KLEBOSCOPE_TMPDIR")
        if env_tmp:
            Path(env_tmp).mkdir(parents=True, exist_ok=True)
            return tempfile.mkdtemp(prefix=prefix, dir=env_tmp)
        return tempfile.mkdtemp(prefix=prefix)

    # ------------------------------------------------------------------
    # Console helpers
    # ------------------------------------------------------------------
    def setup_logging(self, output_dir: Path):
        log_file = output_dir / "kleboscope_run.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file, mode='w')],
        )
        self.logger = logging.getLogger("Kleboscope")
        self.logger.info(f"Kleboscope v{__version__} started")
        self.logger.info(f"Log file: {log_file}")
        self.user_output_dir = output_dir

    def print_info(self, m):
        if not self.quiet:
            print(f"{Color.CYAN}[INFO]{Color.RESET} {m}")
        if self.logger:
            self.logger.info(m)

    def print_success(self, m):
        if not self.quiet:
            print(f"{Color.BRIGHT_GREEN}✓{Color.RESET} {m}")
        if self.logger:
            self.logger.info(f"SUCCESS: {m}")

    def print_warning(self, m):
        if not self.quiet:
            print(f"{Color.BRIGHT_YELLOW}⚠{Color.RESET} {m}")
        if self.logger:
            self.logger.warning(m)
        self.warnings.append(m)

    def print_error(self, m):
        print(f"{Color.BRIGHT_RED}✗{Color.RESET} {m}")
        if self.logger:
            self.logger.error(m)
        self.errors.append(m)

    def print_header(self, title: str, subtitle: str = ""):
        if not self.quiet:
            width = 80
            print()
            print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'=' * width}{Color.RESET}")
            pad = max(0, (width - len(title)) // 2)
            print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * pad}{title}{Color.RESET}")
            if subtitle:
                sub_pad = max(0, (width - len(subtitle)) // 2)
                print(f"{Color.DIM}{Color.WHITE}{' ' * sub_pad}{subtitle}{Color.RESET}")
            print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'=' * width}{Color.RESET}")
            print()
        if self.logger:
            self.logger.info(f"=== {title} ===")
            if subtitle:
                self.logger.info(f"   {subtitle}")

    def display_module_header(self, title: str, subtitle: str = ""):
        if self.quiet:
            return
        width = 80
        line = f"{Color.BOLD}{Color.BRIGHT_MAGENTA}{'─' * width}{Color.RESET}"
        print(line)
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}  ▸ {title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}    {subtitle}{Color.RESET}")
        print(line)

    def display_banner(self):
        if self.quiet:
            return
        print(f"""{Color.BOLD}{Color.BRIGHT_MAGENTA}
{'=' * 80}
{' ' * 20}KLEBOSCOPE v{__version__}
{' ' * 12}K. pneumoniae Genomic Analysis Pipeline
{'=' * 80}{Color.RESET}{Color.BRIGHT_CYAN}
Complete K. pneumoniae genomic analysis:
  MLST | K/O Locus | ABRicate | AMRFinderPlus | MGE | QC | Reports

Author: {AUTHOR_NAME} <{AUTHOR_EMAIL}>
Affiliation: {AUTHOR_AFFILIATION}
{Color.RESET}{Color.DIM}{'=' * 80}{Color.RESET}
""")

    def _get_quotes(self) -> List[str]:
        return [
            "Science is organised knowledge. — Herbert Spencer",
            "The science of today is the technology of tomorrow. — Edward Teller",
            "Genomics is a lens on biology. — Eric Lander",
            "Every microbe has its own story.",
            "Sequence today, understand tomorrow.",
            "Genes are the language of life. — Francis Collins",
            "Small things, big impact.",
            "In science, there are no shortcuts to truth. — Karl Popper",
            "Nothing in life is to be feared, it is only to be understood. — Marie Curie",
            "Science knows no country, because knowledge belongs to humanity. — Louis Pasteur",
        ]

    def display_random_quote(self):
        if self.quiet or not self.quotes:
            return
        quote = random.choice(self.quotes)
        color = random.choice(self.quote_colors)
        print()
        print(f"{Color.DIM}{'─' * 80}{Color.RESET}")
        print(f"{color}   \"{quote}\"{Color.RESET}")
        print(f"{Color.DIM}{'─' * 80}{Color.RESET}")
        print()

    def _run_timed(self, label: str, fn, *args, **kwargs):
        start = datetime.now()
        result = fn(*args, **kwargs)
        elapsed = datetime.now() - start
        elapsed_secs = elapsed.total_seconds()
        elapsed_str = str(elapsed).split('.')[0]
        self.module_times[label] = elapsed_secs
        if not self.quiet:
            ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            print(f"{Color.DIM}[t]{Color.RESET} "
                  f"{Color.DIM}[{ts}]{Color.RESET} "
                  f"{label} finished in {Color.BRIGHT_WHITE}{elapsed_str}{Color.RESET}")
        if self.logger:
            self.logger.info(f"TIMING: {label} = {elapsed_str}")
        return result

    def _log_subprocess_result(self, name: str, result) -> bool:
        """Log full subprocess output; echo to console on failure."""
        if result.stdout:
            self.logger.info(f"=== {name} STDOUT ===\n{result.stdout}")
        if result.stderr:
            self.logger.warning(f"=== {name} STDERR ===\n{result.stderr}")
        if result.returncode != 0:
            print(f"\n{Color.BRIGHT_RED}--- {name} FAILED (rc={result.returncode}) ---{Color.RESET}")
            if result.stdout:
                print(f"{Color.DIM}STDOUT:{Color.RESET}\n{result.stdout}")
            if result.stderr:
                print(f"{Color.DIM}STDERR:{Color.RESET}\n{result.stderr}")
            print(f"{Color.BRIGHT_RED}--- END {name} ---{Color.RESET}\n")
            return False
        return True

    # ------------------------------------------------------------------
    # Module path resolution
    # ------------------------------------------------------------------
    def get_module_path(self, module_name: str) -> Path:
        """Resolve a module directory, honouring environment overrides."""
        env_override = os.environ.get(f'KLEBOSCOPE_{module_name.upper()}_MODULE_PATH')
        if env_override:
            return Path(env_override)
        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "kleboscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def _get_system_module_path(self, module_name: str) -> Path:
        if hasattr(sys, 'prefix'):
            share_path = Path(sys.prefix) / "share" / "kleboscope" / "modules" / module_name
            if share_path.exists():
                return share_path
        return self.base_dir / "modules" / module_name

    def get_shared_data_root(self) -> Path:
        """Writable user-local data root for read-only module installs."""
        root = Path(os.environ.get(
            "KLEBOSCOPE_SHARED_DATA",
            Path.home() / ".local" / "share" / "kleboscope",
        ))
        root.mkdir(parents=True, exist_ok=True)
        return root

    # ------------------------------------------------------------------
    # MLST database resolution
    # ------------------------------------------------------------------
    @staticmethod
    def _mlst_db_has_scheme(db_dir: Optional[Path]) -> bool:
        return db_dir is not None and (db_dir / "pubmlst" / "klebsiella").exists()

    def get_mlst_db_dir(self, for_update: bool = False) -> Path:
        """Resolve the MLST database directory.

        Resolution order:
          1. KLEBOSCOPE_MLST_DB env override
          2. User-local ~/.local/share/kleboscope/mlst_module/db
          3. Packaged module db/

        When for_update=True, always return a writable path.
        """
        env_db = os.environ.get("KLEBOSCOPE_MLST_DB")
        if env_db:
            p = Path(env_db)
            if for_update:
                p.mkdir(parents=True, exist_ok=True)
                return p
            if self._mlst_db_has_scheme(p):
                return p

        user_module = Path(os.environ.get(
            "KLEBOSCOPE_MLST_MODULE_PATH",
            self.get_shared_data_root() / "mlst_module",
        ))
        user_db = user_module / "db"

        if for_update:
            user_db.mkdir(parents=True, exist_ok=True)
            return user_db

        if self._mlst_db_has_scheme(user_db):
            self.print_info(f"Using user-local MLST database: {user_db}")
            return user_db

        module_db = self.get_module_path("mlst_module") / "db"
        if self._mlst_db_has_scheme(module_db):
            return module_db

        self.print_warning(
            "No populated MLST database found. Run `kleboscope --pull-mlst-db`."
        )
        user_db.mkdir(parents=True, exist_ok=True)
        return user_db

    # ------------------------------------------------------------------
    # MGE database resolution
    # ------------------------------------------------------------------
    @staticmethod
    def _mge_db_has_data(db_dir: Optional[Path]) -> bool:
        """Return True if the directory contains the beatrix DIAMOND index."""
        if db_dir is None or not db_dir.exists():
            return False
        candidates = [
            db_dir / "beatrix-1-6_v1_all" / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "mobileOG-db-beatrix-1.6.dmnd",
            db_dir / "beatrix-1-6_v1_all.dmnd",
        ]
        return any(c.exists() for c in candidates)

    def get_mge_data_dir(self, for_update: bool = False) -> Optional[Path]:
        """Return the mobileOG-db data directory.

        Resolution order:
          1. KLEBOSCOPE_MGE_DATA_DIR env override
          2. User-local ~/.local/share/kleboscope/mge_module/mobileOG-db/
          3. Packaged <module>/mobileOG-db/ (Docker / Singularity)
        """
        env_dir = os.environ.get("KLEBOSCOPE_MGE_DATA_DIR")
        if env_dir:
            env_path = Path(env_dir)
            if for_update:
                env_path.mkdir(parents=True, exist_ok=True)
                return env_path
            if self._mge_db_has_data(env_path):
                return env_path

        user_module = Path(os.environ.get(
            "KLEBOSCOPE_MGE_MODULE_PATH",
            self.get_shared_data_root() / "mge_module",
        ))
        user_dir = user_module / "mobileOG-db"

        if for_update:
            user_dir.mkdir(parents=True, exist_ok=True)
            return user_dir
        if self._mge_db_has_data(user_dir):
            if self.logger:
                self.logger.info(f"Using user-local mobileOG-db: {user_dir}")
            return user_dir

        module_db = self.get_module_path("mge_module") / "mobileOG-db"
        if self._mge_db_has_data(module_db):
            return module_db

        return None

    # ------------------------------------------------------------------
    # Exec-bit restoration
    # ------------------------------------------------------------------
    @staticmethod
    def _restore_exec_bits(root: Path, logger=None) -> int:
        """Restore +x on runnable files (zip extraction / NFS strip it)."""
        if not root.exists():
            return 0
        executable_exts = {".pl", ".sh", ".py", ".rb"}
        count = 0
        for f in root.rglob("*"):
            if not f.is_file():
                continue
            try:
                parent_is_bin = f.parent.name in ("bin", "scripts")
                is_script = f.suffix.lower() in executable_exts
                if not (parent_is_bin or is_script):
                    continue
                current = f.stat().st_mode
                new_mode = current | 0o111
                if new_mode != current:
                    f.chmod(new_mode)
                    count += 1
            except Exception:
                pass
        if logger and count > 0:
            logger.info(f"Restored execute bit on {count} file(s) under {root}")
        return count

    # ------------------------------------------------------------------
    # MLST database: pull
    # ------------------------------------------------------------------
    def _ensure_mlst_blast_index(self, target_path: Path) -> bool:
        """Build db/blast/mlst.fa if missing (required by bin/mlst)."""
        blast_fa = target_path / "db" / "blast" / "mlst.fa"
        if blast_fa.exists():
            return True

        self.print_info("BLST blast index missing; building db/blast/mlst.fa...")
        for candidate in ("scripts/mlst-make_blast_db.sh", "bin/mlst-make_blast_db.sh"):
            script = target_path / candidate
            if script.exists():
                result = subprocess.run(
                    ["bash", str(script)], cwd=target_path,
                    capture_output=True, text=True,
                )
                if result.stdout:
                    self.logger.info(result.stdout)
                if result.stderr:
                    self.logger.warning(result.stderr)
                if result.returncode == 0 and blast_fa.exists():
                    self.print_success("MLST blast index built.")
                    return True
                self.print_warning(
                    f"{candidate} did not produce {blast_fa} "
                    f"(makeblastdb must be on PATH)."
                )
                return False

        self.print_warning(
            f"No make_blast_db script found; run it manually in {target_path} "
            f"or MLST will fail with 'db/blast/mlst.fa: No such file'."
        )
        return False

    def pull_mlst_database(self) -> bool:
        """Install the Klebsiella MLST database.

        Writes into the packaged module directory if writable, otherwise
        ~/.local/share/kleboscope/mlst_module/. Prefers the module's
        database.sh installer; falls back to the prebuilt GitHub archive.
        """
        module_path = self.get_module_path("mlst_module")

        if module_path.exists() and os.access(module_path, os.W_OK):
            target_path = module_path
            self.print_info(f"Using writable module directory: {target_path}")
        else:
            target_path = Path(os.environ.get(
                "KLEBOSCOPE_MLST_MODULE_PATH",
                self.get_shared_data_root() / "mlst_module",
            ))
            target_path.mkdir(parents=True, exist_ok=True)
            self.print_info(
                f"System MLST module is read-only; using user-local path: {target_path}"
            )

        system_path = self._get_system_module_path("mlst_module")
        target_db = target_path / "db"

        if self._mlst_db_has_scheme(target_db):
            self.print_success(f"MLST database already present at {target_db}")
            self._ensure_mlst_blast_index(target_path)
            os.environ['KLEBOSCOPE_MLST_MODULE_PATH'] = str(target_path)
            os.environ['KLEBOSCOPE_MLST_DB'] = str(target_db)
            return True

        installer = module_path / "database.sh"
        if installer.exists():
            self.print_info(f"Running {installer.name} to install MLST database...")
            env = os.environ.copy()
            env['KLEBOSCOPE_MLST_DB'] = str(target_db)
            env['MLST_DATADIR'] = str(target_db)
            result = subprocess.run(
                ["bash", str(installer)], cwd=module_path,
                capture_output=True, text=True, env=env,
            )
            if not self._log_subprocess_result("MLST database install", result):
                self.print_error("MLST database installation failed.")
                return False
            if not self._mlst_db_has_scheme(target_db):
                self.print_error(
                    f"Installation completed but no scheme found at {target_db}"
                )
                return False
            self._restore_exec_bits(target_path, logger=self.logger)
            self._ensure_mlst_blast_index(target_path)
            self.print_success(f"MLST database installed to {target_db}")
            os.environ['KLEBOSCOPE_MLST_MODULE_PATH'] = str(target_path)
            os.environ['KLEBOSCOPE_MLST_DB'] = str(target_db)
            return True

        return self._download_mlst_from_github(target_path, system_path)

    def _download_mlst_from_github(self, target_path: Path,
                                   system_path: Path) -> bool:
        """Download a prebuilt MLST database archive from GitHub."""
        url = os.environ.get(
            "KLEBOSCOPE_MLST_PULL_URL",
            "https://github.com/bbeckley-hub/mlst/archive/refs/heads/master.zip",
        )
        extract_dir = None
        try:
            self.print_info(f"Downloading MLST database from {url}")
            with urllib.request.urlopen(url, timeout=60) as resp:
                data = resp.read()
            self.print_info(f"Downloaded {len(data) / 1024:.1f} KB")

            extract_dir = Path(self._mkdtemp(prefix="kleboscope_mlst_extract_"))
            zip_path = extract_dir / "mlst.zip"
            zip_path.write_bytes(data)
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(extract_dir)

            candidates = list(extract_dir.glob("mlst-*"))
            if not candidates:
                self.print_error("Extracted repository root not found in archive")
                return False
            src_root = candidates[0]

            for dirname in ('bin', 'db', 'perl5', 'scripts'):
                src = src_root / dirname
                if src.exists():
                    dst = target_path / dirname
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                    self.logger.info(f"Copied {dirname}/ to {dst}")
                else:
                    self.print_warning(f"'{dirname}' not found in repository, skipping.")

            self._restore_exec_bits(target_path, logger=self.logger)

            if target_path != system_path:
                for script in ('klebo_mlst.py', 'database.sh'):
                    src_script = system_path / script
                    if src_script.exists():
                        shutil.copy2(src_script, target_path / script)

            db_check = target_path / "db" / "pubmlst" / "klebsiella"
            if not db_check.exists():
                self.print_error(f"MLST scheme missing after pull: {db_check}")
                return False

            allele_files = sorted(db_check.glob("*.tfa"))
            if not allele_files:
                self.print_error(f"No allele files (*.tfa) found in {db_check}")
                return False
            self.print_info(f"Verified {len(allele_files)} allele loci in {db_check}")

            removed = 0
            pubmlst = target_path / "db" / "pubmlst"
            if pubmlst.exists():
                for item in pubmlst.iterdir():
                    if item.is_dir() and item.name != "klebsiella":
                        shutil.rmtree(item, ignore_errors=True)
                        removed += 1
            if removed:
                self.print_info(f"Pruned {removed} unused MLST scheme(s)")

            self._ensure_mlst_blast_index(target_path)

            os.environ['KLEBOSCOPE_MLST_MODULE_PATH'] = str(target_path)
            os.environ['KLEBOSCOPE_MLST_DB'] = str(target_path / "db")
            self.print_success(f"MLST database installed to {target_path / 'db'}")
            return True

        except Exception as exc:
            self.print_error(f"MLST download failed: {exc}")
            self.print_info(
                "If this URL is unavailable, download manually or set "
                "KLEBOSCOPE_MLST_PULL_URL to a valid archive URL."
            )
            return False
        finally:
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    def update_mlst_database(self) -> bool:
        """Refresh the Klebsiella MLST scheme via klebo_mlst.py --update-db."""
        module_path = self.get_module_path("mlst_module")
        if module_path.exists() and os.access(module_path, os.W_OK):
            target_module = module_path
        else:
            target_module = Path(os.environ.get(
                "KLEBOSCOPE_MLST_MODULE_PATH",
                self.get_shared_data_root() / "mlst_module",
            ))
            target_module.mkdir(parents=True, exist_ok=True)

        script = target_module / "klebo_mlst.py"
        if not script.exists():
            script = module_path / "klebo_mlst.py"
        if not script.exists():
            self.print_error(f"MLST script not found: {script}")
            return False

        db_dir = target_module / "db"
        db_dir.mkdir(parents=True, exist_ok=True)

        if not self._mlst_db_has_scheme(db_dir):
            self.print_error(
                f"MLST module tree is not populated at {target_module}. "
                "Run `kleboscope --pull-mlst-db` first."
            )
            return False

        os.environ['KLEBOSCOPE_MLST_MODULE_PATH'] = str(target_module)
        os.environ['KLEBOSCOPE_MLST_DB'] = str(db_dir)

        self.print_info(f"Updating MLST database at {db_dir}...")
        env = os.environ.copy()
        env['KLEBOSCOPE_MLST_DB'] = str(db_dir)
        result = subprocess.run(
            [sys.executable, str(script), "--update-db",
             "-db", str(db_dir), "-sc", str(target_module)],
            capture_output=True, text=True, cwd=target_module, env=env,
        )
        if not self._log_subprocess_result("MLST database update", result):
            self.print_error("MLST database update failed.")
            return False

        if not self._mlst_db_has_scheme(db_dir):
            self.print_error(f"Update finished but no scheme found at {db_dir}")
            return False

        self._ensure_mlst_blast_index(target_module)
        self.print_success(f"MLST database updated at {db_dir}")
        return True

    # ------------------------------------------------------------------
    # MGE database: pull
    # ------------------------------------------------------------------
    def pull_mge_database(self) -> bool:
        """Download mobileOG-db from the GitHub release.

        Docker and Singularity users do not need this — the data is baked
        into the image and get_mge_data_dir() finds it there automatically.
        """
        existing = self.get_mge_data_dir()
        if existing is not None:
            self.print_success(f"mobileOG-db already present at {existing}")
            os.environ["KLEBOSCOPE_MGE_DATA_DIR"] = str(existing)
            return True

        module_path = self.get_module_path("mge_module")
        if module_path.exists() and os.access(module_path, os.W_OK):
            target_root = module_path
            self.print_info(f"Using writable module directory: {target_root}")
        else:
            target_root = self.get_shared_data_root() / "mge_module"
            target_root.mkdir(parents=True, exist_ok=True)
            self.print_info(
                f"Module dir is read-only; using user-local path: {target_root}"
            )

        url = os.environ.get(
            "KLEBOSCOPE_MGE_PULL_URL",
            "https://github.com/bbeckley-hub/Kleboscope/releases/download/"
            "v1.3.0/kleboscope-mge-data-v1.3.0.zip",
        )
        extract_dir = None
        try:
            self.print_info(f"Downloading mobileOG-db from {url}")
            self.print_info(
                "Payload is approximately 300 MB, this may take a few minutes..."
            )
            temp_zip = Path(self._mkdtemp(prefix="kleboscope_mge_zip_")) / "mge.zip"
            urllib.request.urlretrieve(url, temp_zip)
            self.print_success("Download complete")

            extract_dir = Path(self._mkdtemp(prefix="kleboscope_mge_extract_"))
            with zipfile.ZipFile(temp_zip, "r") as zf:
                zf.extractall(extract_dir)
            self.logger.info(f"Extracted to {extract_dir}")

            beatrix_dir = None
            for root, dirs, _ in os.walk(extract_dir):
                if "beatrix-1-6_v1_all" in dirs:
                    beatrix_dir = Path(root) / "beatrix-1-6_v1_all"
                    break
            if beatrix_dir is None:
                self.print_error("Could not find beatrix-1-6_v1_all in archive")
                return False

            dmnd = beatrix_dir / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd.exists():
                self.print_error(f"DIAMOND index missing from archive: {dmnd.name}")
                return False

            target_beatrix = target_root / "mobileOG-db" / "beatrix-1-6_v1_all"
            if target_beatrix.exists():
                shutil.rmtree(target_beatrix)
            target_beatrix.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(beatrix_dir, target_beatrix)

            dmnd_target = target_beatrix / "mobileOG-db-beatrix-1.6.dmnd"
            if not dmnd_target.exists():
                self.print_error(
                    f"Install verification failed: {dmnd_target} not present"
                )
                return False

            size_mb = dmnd_target.stat().st_size / (1024 * 1024)
            self.print_success(
                f"mobileOG-db installed to {target_root / 'mobileOG-db'} "
                f"(DIAMOND index: {size_mb:.0f} MB)"
            )
            os.environ["KLEBOSCOPE_MGE_DATA_DIR"] = str(target_root / "mobileOG-db")
            os.environ["KLEBOSCOPE_MGE_MODULE_PATH"] = str(target_root)
            return True

        except Exception as exc:
            self.print_error(f"mobileOG-db download failed: {exc}")
            self.print_info(
                "If this URL is unavailable, place mobileOG-db manually at "
                "~/.local/share/kleboscope/mge_module/mobileOG-db/ or set "
                "KLEBOSCOPE_MGE_PULL_URL to a valid archive URL."
            )
            return False
        finally:
            if extract_dir and extract_dir.exists():
                shutil.rmtree(extract_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # AMR database management
    # ------------------------------------------------------------------
    def update_amr_database(self, force: bool = False) -> bool:
        """Update the AMRFinderPlus database."""
        module_path = self.get_module_path("amr_module")
        script = module_path / "klebo_amrfinder.py"
        if not script.exists():
            self.print_error(f"AMR script not found: {script}")
            return False

        if force:
            self.print_info("Force mode: module will overwrite existing database.")

        self.print_info("Updating AMRFinderPlus database...")
        flag = "--force-update" if force else "--update-db"
        result = subprocess.run(
            [sys.executable, str(script), flag],
            capture_output=True, text=True, cwd=module_path,
        )
        if not self._log_subprocess_result("AMR database update", result):
            self.print_error("AMR database update failed.")
            return False

        self.print_success("AMR database updated successfully.")
        ver_result = subprocess.run(
            [sys.executable, str(script), "--db-version"],
            capture_output=True, text=True, cwd=module_path,
        )
        if ver_result.returncode == 0 and ver_result.stdout:
            self.print_info(f"New database version: {ver_result.stdout.strip()}")
        return True

    def ensure_amr_database(self) -> bool:
        """Verify AMRFinderPlus database exists; auto-update if missing."""
        module_path = self.get_module_path("amr_module")
        script = module_path / "klebo_amrfinder.py"
        if not script.exists():
            self.print_error("AMR script not found, cannot check database.")
            return False

        result = subprocess.run(
            [sys.executable, str(script), "--db-version"],
            capture_output=True, text=True, cwd=module_path,
        )
        if (result.returncode == 0
                and "Unknown" not in result.stdout
                and "None" not in result.stdout):
            self.print_success(
                f"AMR database already present: {result.stdout.strip()}"
            )
            return True

        self.print_warning("AMR database missing; attempting automatic update...")
        return self.update_amr_database(force=False)

    # ------------------------------------------------------------------
    # FASTA discovery
    # ------------------------------------------------------------------
    def find_fasta_files(self, input_path: str) -> List[Path]:
        if '*' in input_path or '?' in input_path:
            matched = glob.glob(input_path)
            files = [Path(f) for f in matched
                     if Path(f).is_file()
                     and f.lower().endswith(('.fna', '.fasta', '.fa', '.fsa'))
                     and not Path(f).name.startswith('.')]
            self.print_info(f"Found {len(files)} FASTA files from pattern")
            return sorted(files)

        p = Path(input_path)
        if p.is_file() and p.suffix.lower() in ('.fna', '.fasta', '.fa', '.fsa'):
            return [p]
        if p.is_dir():
            files = []
            for ext in ('*.fna', '*.fasta', '*.fa', '*.fsa',
                        '*.FNA', '*.FASTA', '*.FA', '*.FSA'):
                files.extend(p.glob(ext))
            files = [f for f in files if f.is_file() and not f.name.startswith('.')]
            self.print_info(f"Found {len(files)} FASTA files in {p}")
            return sorted(set(files))
        self.print_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return "*.fna"
        exts = set(f.suffix.lower() for f in fasta_files)
        if len(exts) == 1:
            return f"*{list(exts)[0]}"
        return "*"

    # ------------------------------------------------------------------
    # Copy helpers
    # ------------------------------------------------------------------
    def _copy_module_results(self, src: Path, dst: Path,
                             skip_dirs: Optional[Set[str]] = None) -> None:
        """Copy module results to the user output dir, respecting report_level."""
        if not src.exists():
            return
        if dst.exists():
            shutil.rmtree(dst)
        dst.mkdir(parents=True, exist_ok=True)

        skip_dirs = skip_dirs or set()
        for item in src.iterdir():
            if item.name in skip_dirs:
                continue
            if item.is_file():
                shutil.copy2(item, dst / item.name)
            elif item.is_dir() and self.report_level == 'full':
                shutil.copytree(item, dst / item.name)

    # ------------------------------------------------------------------
    # Generic module runner
    # ------------------------------------------------------------------
    def run_module_in_temp(
        self,
        module_name: str,
        fasta_files: List[Path],
        cmd_parts: List[str],
        result_subdir: Optional[str] = None,
        extra_result_files: Optional[List[str]] = None,
        env_overrides: Optional[Dict[str, str]] = None,
    ) -> Tuple[bool, str]:
        """Copy module to a temp dir, run it there, copy results back.

        The module is copied to <temp>/<module_name>, FASTA files are
        copied next to it, and the command runs with cwd=<temp>/<module_name>
        so module-relative paths (db/, bin/, scripts/) resolve correctly.
        """
        module_orig = self.get_module_path(module_name)
        if not module_orig.exists():
            return False, f"Module directory not found: {module_orig}"

        temp_dir = Path(self._mkdtemp(prefix=f"kleboscope_{module_name}_"))
        self.temp_dirs.append(str(temp_dir))
        if self.logger:
            self.logger.info(f"Temp dir for {module_name}: {temp_dir}")

        try:
            work_dir = temp_dir / module_name
            shutil.copytree(module_orig, work_dir, dirs_exist_ok=True)
            self._restore_exec_bits(work_dir, logger=self.logger)

            for f in fasta_files:
                shutil.copy2(f, work_dir / f.name)

            script_name = Path(cmd_parts[1]).name
            script_path = work_dir / script_name
            if not script_path.exists():
                return False, f"Script not found in {module_name}: {script_name}"

            abs_cmd = [cmd_parts[0], str(script_path)] + cmd_parts[2:]

            env = os.environ.copy()
            if env_overrides:
                env.update(env_overrides)

            if self.logger:
                self.logger.info(f"Command: {' '.join(abs_cmd)}")
            result = subprocess.run(
                abs_cmd, cwd=work_dir, capture_output=True, text=True, env=env,
            )
            if not self._log_subprocess_result(module_name, result):
                return False, f"{module_name} failed (rc={result.returncode})"

            if result_subdir:
                src = work_dir / result_subdir
                if src.exists():
                    dst = self.user_output_dir / result_subdir
                    self._copy_module_results(src, dst)
                    if self.logger:
                        self.logger.info(f"Results copied to {dst}")
                else:
                    self.print_warning(
                        f"{module_name} ran but produced no {result_subdir}/ directory"
                    )

            if extra_result_files:
                for extra in extra_result_files:
                    src_extra = work_dir / extra
                    if src_extra.exists():
                        shutil.copy2(src_extra, self.user_output_dir / extra)
                        if self.logger:
                            self.logger.info(f"Copied {extra} to output")

            return True, f"{module_name} completed successfully"

        except Exception as exc:
            if self.logger:
                self.logger.error(
                    f"Exception in {module_name}: {exc}\n{traceback.format_exc()}"
                )
            return False, f"Exception: {exc}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)
                if self.logger:
                    self.logger.info(f"Removed temporary directory: {temp_dir}")

    # ------------------------------------------------------------------
    # Analysis modules
    # ------------------------------------------------------------------
    def run_qc_analysis(self, fasta_files: List[Path], output_dir: Path,
                        threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "klebo_fasta_qc.py", pattern,
               "-o", DIR_QC, "-c", str(threads)]
        return self.run_module_in_temp(
            "qc_module", fasta_files, cmd, result_subdir=DIR_QC,
        )

    def run_mlst_analysis(self, fasta_files: List[Path], output_dir: Path,
                          threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        db_dir = self.get_mlst_db_dir()
        cmd = [sys.executable, "klebo_mlst.py",
               "-i", pattern, "-o", DIR_MLST,
               "-db", str(db_dir), "-sc", "bin",
               "-c", str(threads), "--batch"]
        return self.run_module_in_temp(
            "mlst_module", fasta_files, cmd, result_subdir=DIR_MLST,
            env_overrides={
                'KLEBOSCOPE_MLST_DB': str(db_dir),
                'MLST_DATADIR': str(db_dir),
            },
        )

    def run_kaptive_analysis(self, fasta_files: List[Path], output_dir: Path,
                             threads: int) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "klebo_kaptive.py",
               "-i", pattern, "-o", DIR_KAPTIVE, "-c", str(threads)]
        return self.run_module_in_temp(
            "kaptive_module", fasta_files, cmd, result_subdir=DIR_KAPTIVE,
        )

    def run_abricate_analysis(self, fasta_files: List[Path], output_dir: Path,
                              threads: int, min_id: int = 80,
                              min_cov: int = 80) -> Tuple[bool, str]:
        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "kleb_abricate_module.py", pattern,
               "--min-id", str(min_id), "--min-cov", str(min_cov),
               "--cpus", str(threads)]
        return self.run_module_in_temp(
            "abricate_module", fasta_files, cmd, result_subdir=DIR_ABRICATE,
        )

    def run_amr_analysis(self, fasta_files: List[Path], output_dir: Path,
                         threads: int, min_identity: Optional[float] = None,
                         min_coverage: Optional[float] = None,
                         skip_mutations: bool = False,
                         force_update: bool = False) -> Tuple[bool, str]:
        if force_update:
            self.update_amr_database(force=True)
        elif not self.ensure_amr_database():
            return False, "AMR database missing and could not be updated."

        pattern = self.get_file_pattern(fasta_files)
        cmd = [sys.executable, "klebo_amrfinder.py", pattern,
               "--cpus", str(threads)]
        if min_identity is not None:
            cmd.extend(["--min-identity", str(min_identity)])
        if min_coverage is not None:
            cmd.extend(["--min-coverage", str(min_coverage)])
        if skip_mutations:
            cmd.append("--skip-mutations")
        return self.run_module_in_temp(
            "amr_module", fasta_files, cmd, result_subdir=DIR_AMR,
        )

    # ------------------------------------------------------------------
    # MGE summary collection
    # ------------------------------------------------------------------
    def _collect_mge_summaries(self, temp_dir: Path, dst: Path) -> List[str]:
        """Locate the MGEFinder batch summaries and stage them for reporters.

        MGEFinder.py may write the summary files at the top level of the
        run directory, inside the -d mge_results/ directory, or inside the
        per-sample annotations subtree (which the result copy deliberately
        skips). This searches the entire temporary run directory, copies
        every summary it finds into dst under its canonical name, and
        returns the canonical names that are still missing.
        """
        expected = [
            "kleboscope_mge_summary.tsv",
            "kleboscope_mge_all_hits.tsv",
            "kleboscope_mge_summary.json",
            "kleboscope_mge_summary.html",
        ]
        variants = {
            "kleboscope_mge_summary.tsv": "mge_summary.tsv",
            "kleboscope_mge_all_hits.tsv": "mge_all_hits.tsv",
            "kleboscope_mge_summary.json": "mge_summary.json",
            "kleboscope_mge_summary.html": "mge_summary.html",
        }
        excluded = {"mobileOG-db", "__pycache__"}

        def _search(name: str) -> List[Path]:
            hits = [
                p for p in temp_dir.rglob(name)
                if p.is_file() and not (excluded & set(p.parts))
            ]
            hits.sort(key=lambda p: len(p.parts))
            return hits

        missing = []
        for canonical in expected:
            hits = _search(canonical)
            if hits:
                if hits[0] != dst / canonical:
                    shutil.copy2(hits[0], dst / canonical)
                    self.logger.info(
                        f"Collected {canonical} from {hits[0]}"
                    )
                continue

            hits = _search(variants[canonical])
            if hits:
                shutil.copy2(hits[0], dst / canonical)
                self.print_warning(
                    f"MGE summary written as '{variants[canonical]}'; "
                    f"staged as '{canonical}' in {dst}"
                )
            else:
                missing.append(canonical)

        if "kleboscope_mge_summary.html" in missing:
            mge_src = temp_dir / "mge_results"
            if mge_src.exists():
                for html in sorted(mge_src.rglob("*.html")):
                    if html.is_file():
                        self.logger.info(
                            f"MGE html present under unexpected name: {html}"
                        )
        return missing

    def run_mge_analysis(self, fasta_files: List[Path], output_dir: Path,
                         threads: int) -> Tuple[bool, str]:
        """Run MGE profiling with mobileOG-db.

        The resolved mobileOG-db payload is copied into the temp working
        directory (no symlinks) so the pipeline works on filesystems that
        do not support symlinks. Actual filenames are passed to
        MGEFinder.py — it does not expand glob patterns. Batch summaries
        are collected from wherever MGEFinder writes them and staged in
        the consolidated mge_results/ directory for the reporters.
        """
        module_orig = self.get_module_path("mge_module")
        if not module_orig.exists():
            return False, f"mge_module not found: {module_orig}"

        mge_data = self.get_mge_data_dir()
        if mge_data is None:
            return False, (
                "mobileOG-db not found. Run `kleboscope --pull-mge-db` first, "
                "or use the Docker/Singularity image (which has it bundled)."
            )
        if self.logger:
            self.logger.info(f"Resolved mobileOG-db: {mge_data}")

        temp_dir = Path(self._mkdtemp(prefix="kleboscope_mge_"))
        self.temp_dirs.append(str(temp_dir))
        if self.logger:
            self.logger.info(f"Temp dir for mge: {temp_dir}")

        try:
            shutil.copytree(
                module_orig, temp_dir, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns(
                    "mge_results", "annotations",
                    "*.fna", "*.fasta", "*.fa", "*.fsa",
                    "mobileOG-db", "__pycache__", "*.pyc",
                ),
            )
            self._restore_exec_bits(temp_dir, logger=self.logger)

            beatrix_src = mge_data / "beatrix-1-6_v1_all"
            if not beatrix_src.exists() and mge_data.name == "beatrix-1-6_v1_all":
                beatrix_src = mge_data
            if not beatrix_src.exists():
                return False, f"beatrix-1-6_v1_all not found under {mge_data}"

            beatrix_dst = temp_dir / "mobileOG-db" / "beatrix-1-6_v1_all"
            beatrix_dst.parent.mkdir(parents=True, exist_ok=True)
            if beatrix_dst.exists():
                shutil.rmtree(beatrix_dst)
            shutil.copytree(beatrix_src, beatrix_dst)
            if self.logger:
                self.logger.info(f"Copied {beatrix_src} -> {beatrix_dst}")

            for f in fasta_files:
                shutil.copy2(f, temp_dir / f.name)

            script = temp_dir / "MGEFinder.py"
            if not script.exists():
                return False, "MGEFinder.py not found in mge_module"

            self.print_info(
                f"Running MGE analysis on {len(fasta_files)} genome(s)..."
            )
            mge_inputs = [f.name for f in fasta_files]
            cmd = ([sys.executable, str(script), "-i"] + mge_inputs +
                   ["-d", "mge_results",
                    "-db_dir", "mobileOG-db/beatrix-1-6_v1_all"])
            result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True)
            if not self._log_subprocess_result("MGE", result):
                return False, f"MGE analysis failed (rc={result.returncode})"

            src = temp_dir / "mge_results"
            if not src.exists():
                return False, f"MGE results directory not found: {src}"

            dst = self.user_output_dir / DIR_MGE
            self._copy_module_results(src, dst, skip_dirs={"annotations"})

            missing = self._collect_mge_summaries(temp_dir, dst)
            if missing:
                self.print_warning(
                    "MGE ran but these summary files are missing: "
                    f"{', '.join(missing)}. Reporters may not include MGE data."
                )

            self.print_success("MGE profiling completed successfully.")
            return True, "MGE completed successfully"

        except Exception as exc:
            if self.logger:
                self.logger.error(traceback.format_exc())
            return False, f"MGE exception: {exc}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Species gate
    # ------------------------------------------------------------------
    def _load_species_gate(self, output_dir: Path) -> Optional[Set[str]]:
        gate_file = output_dir / DIR_QC / "kleboscope_species_passed.txt"
        if not gate_file.exists():
            return None
        with open(gate_file, 'r', encoding='utf-8') as f:
            names = {line.strip() for line in f if line.strip()}
        self.logger.info(
            f"Species gate: {len(names)} sample(s) passed from {gate_file}"
        )
        return names

    def _filter_by_species_gate(self, fasta_files: List[Path],
                                gate: Optional[Set[str]]) -> List[Path]:
        if gate is None:
            self.print_warning(
                "Species gate unavailable. Running ALL samples through "
                "downstream modules."
            )
            return fasta_files
        passed = [f for f in fasta_files if f.name in gate]
        excluded = [f for f in fasta_files if f.name not in gate]
        if excluded:
            self.print_info(
                f"Species gate: {len(passed)} of {len(fasta_files)} samples passed"
            )
            for f in excluded:
                self.print_warning(f"  Excluded: {f.name}")
        else:
            self.print_success(
                f"Species gate: all {len(passed)} sample(s) passed"
            )
        return passed

    # ------------------------------------------------------------------
    # Reporter inputs
    # ------------------------------------------------------------------
    def _copy_required_files_to_temp(self, src_output_dir: Path,
                                     temp_work_dir: Path) -> None:
        """Copy every batch file the two reporters require."""
        required_files = [
            (DIR_QC, ["kleboscope_fasta_qc_summary.html",
                      "kleboscope_fasta_qc_summary.tsv",
                      "kleboscope_fasta_qc_summary.json"]),
            (DIR_MLST, ["kleboscope_mlst_summary.html",
                        "kleboscope_mlst_summary.tsv",
                        "kleboscope_mlst_summary.json"]),
            (DIR_KAPTIVE, ["kleboscope_kaptive_summary.html",
                           "kleboscope_kaptive_summary.tsv",
                           "kleboscope_kaptive_summary.json"]),
            (DIR_MGE, ["kleboscope_mge_summary.html",
                       "kleboscope_mge_all_hits.tsv",
                       "kleboscope_mge_summary.tsv",
                       "kleboscope_mge_summary.json"]),
            (DIR_AMR, ["klebo_amrfinder_summary.tsv",
                       "klebo_amrfinder_statistics_summary.tsv",
                       "klebo_amrfinder_summary_report.html",
                       "klebo_amrfinder_master_summary.json",
                       "mutation_summary.html",
                       "mutation_summary.tsv",
                       "mutation_master_summary.json"]),
        ]

        abricate_dbs = ['card', 'resfinder', 'argannot', 'vfdb',
                        'plasmidfinder', 'megares', 'ncbi', 'ecoh',
                        'ecoli_vf', 'bacmet2']
        abricate_files = []
        for db in abricate_dbs:
            abricate_files.append(f"klebo_{db}_abricate_summary.tsv")
            abricate_files.append(f"klebo_{db}_summary_report.html")
        abricate_files.append("klebo_abricate_master_summary.json")
        required_files.append((DIR_ABRICATE, abricate_files))

        for dirname, filenames in required_files:
            src_dir = src_output_dir / dirname
            if not src_dir.exists():
                self.logger.warning(f"Source directory not found: {src_dir}")
                continue
            for fname in filenames:
                src_file = src_dir / fname
                if src_file.exists():
                    shutil.copy2(src_file, temp_work_dir / fname)
                    self.logger.info(f"Copied {fname} to reporter work dir")
                else:
                    self.logger.warning(f"Required file not found: {src_file}")

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------
    def run_gene_centric_reporter(self, output_dir: Path) -> Tuple[bool, str]:
        module_orig = self.get_module_path("gene_centric_module")
        if not module_orig.exists():
            return False, f"Module not found: {module_orig}"

        temp_dir = Path(self._mkdtemp(prefix="kleboscope_gene_centric_"))
        self.temp_dirs.append(str(temp_dir))
        try:
            work_dir = temp_dir / "gene_centric_module"
            shutil.copytree(module_orig, work_dir, dirs_exist_ok=True)
            self._restore_exec_bits(work_dir, logger=self.logger)
            self._copy_required_files_to_temp(output_dir, work_dir)

            script = work_dir / "kleboscope_ultimate_gene_centric_reporter.py"
            if not script.exists():
                return False, (
                    "kleboscope_ultimate_gene_centric_reporter.py not found"
                )

            self.print_info("Running gene-centric reporter...")
            result = subprocess.run(
                [sys.executable, str(script), "-i", "."],
                cwd=work_dir, capture_output=True, text=True,
            )
            if not self._log_subprocess_result("Gene-centric reporter", result):
                return False, "Gene-centric reporter failed."

            results_dir = work_dir / DIR_GENE_CENTRIC
            if not results_dir.exists():
                return False, f"{DIR_GENE_CENTRIC} not produced"

            dst = output_dir / DIR_GENE_CENTRIC
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(results_dir, dst)
            self.print_success(f"Gene-centric reports copied to {dst}")
            return True, "Gene-centric reporter completed successfully"
        except Exception as exc:
            self.logger.error(traceback.format_exc())
            return False, f"Gene-centric reporter exception: {exc}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)

    def run_sample_centric_reporter(self, output_dir: Path) -> Tuple[bool, str]:
        module_orig = self.get_module_path("sample_centric_module")
        if not module_orig.exists():
            return False, f"Module not found: {module_orig}"

        temp_dir = Path(self._mkdtemp(prefix="kleboscope_sample_centric_"))
        self.temp_dirs.append(str(temp_dir))
        try:
            work_dir = temp_dir / "sample_centric_module"
            shutil.copytree(module_orig, work_dir, dirs_exist_ok=True)
            self._restore_exec_bits(work_dir, logger=self.logger)
            self._copy_required_files_to_temp(output_dir, work_dir)

            script = work_dir / "kleboscope_ultimate_sample_centric_report.py"
            if not script.exists():
                return False, (
                    "kleboscope_ultimate_sample_centric_report.py not found"
                )

            self.print_info("Running sample-centric reporter...")
            result = subprocess.run(
                [sys.executable, str(script), "-i", "."],
                cwd=work_dir, capture_output=True, text=True,
            )
            if not self._log_subprocess_result("Sample-centric reporter", result):
                return False, "Sample-centric reporter failed."

            results_dir = work_dir / DIR_SAMPLE_CENTRIC
            if not results_dir.exists():
                return False, f"{DIR_SAMPLE_CENTRIC} not produced"

            dst = output_dir / DIR_SAMPLE_CENTRIC
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(results_dir, dst)
            self.print_success(f"Sample-centric reports copied to {dst}")
            return True, "Sample-centric reporter completed successfully"
        except Exception as exc:
            self.logger.error(traceback.format_exc())
            return False, f"Sample-centric reporter exception: {exc}"
        finally:
            if not self.keep_temp and self.running:
                shutil.rmtree(temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Final consolidation + pruning
    # ------------------------------------------------------------------
    def copy_summary_results_to_final_directory(self, output_dir: Path) -> None:
        try:
            self.print_info("Consolidating results into final report directory...")
            final_dir = output_dir / DIR_FINAL
            final_dir.mkdir(parents=True, exist_ok=True)

            for reporter_dir in (DIR_GENE_CENTRIC, DIR_SAMPLE_CENTRIC):
                src = output_dir / reporter_dir
                if src.exists() and src.is_dir():
                    target = final_dir / reporter_dir
                    if target.exists():
                        shutil.rmtree(target)
                    shutil.copytree(src, target)
                    self.print_info(f"  ▸ Copied {reporter_dir}/")

            self.print_success(f"Final report assembled at {final_dir}")
        except Exception as exc:
            self.print_error(f"Error consolidating final report: {exc}")

    def _cleanup_top_level_duplicates(self, output_path: Path) -> None:
        removed = 0
        for dirname in (DIR_GENE_CENTRIC, DIR_SAMPLE_CENTRIC):
            d = output_path / dirname
            if d.exists():
                shutil.rmtree(d, ignore_errors=True)
                removed += 1
        if removed:
            self.logger.info(f"Removed {removed} duplicate top-level director(ies)")

    def _finalize_output_by_report_level(self, output_path: Path) -> None:
        if self.report_level == 'full':
            self.print_info("Report level: full (all files retained)")
            return

        final_dir = output_path / DIR_FINAL
        if not final_dir.exists() or not any(final_dir.iterdir()):
            self.print_warning(
                "Kleboscope_final_report/ is missing or empty; skipping prune."
            )
            return

        if self.report_level == 'final':
            removed = 0
            for name in ALL_MODULE_DIRS:
                d = output_path / name
                if d.exists():
                    shutil.rmtree(d, ignore_errors=True)
                    removed += 1
            self.print_success(
                f"Report level: final — removed {removed} module directories"
            )
            return

        if self.report_level == 'summary':
            removed = 0
            for name in ALL_MODULE_DIRS:
                d = output_path / name
                if not d.exists():
                    continue
                for item in list(d.iterdir()):
                    if item.is_dir():
                        shutil.rmtree(item, ignore_errors=True)
                        removed += 1
            if removed:
                self.print_success(
                    f"Report level: summary — removed {removed} per-sample folders"
                )
            else:
                self.print_success(
                    "Report level: summary — batch summaries retained"
                )

    # ------------------------------------------------------------------
    # Sequential driver
    # ------------------------------------------------------------------
    def run_sequential_analyses(
        self,
        fasta_files: List[Path],
        output_dir: Path,
        threads: int,
        skip_modules: Dict[str, bool],
        amr_min_identity: Optional[float],
        amr_min_coverage: Optional[float],
        amr_skip_mutations: bool,
        amr_force_update: bool,
        abricate_min_id: int,
        abricate_min_cov: int,
    ) -> Dict[str, bool]:
        results: Dict[str, bool] = {}
        downstream_files = fasta_files

        if not skip_modules.get('qc', False):
            self.display_module_header("FASTA QC", self.module_subtitles['qc'])
            ok, _ = self._run_timed(
                "QC", self.run_qc_analysis, fasta_files, output_dir, threads,
            )
            results['QC'] = ok
            if ok:
                gate = self._load_species_gate(output_dir)
                downstream_files = self._filter_by_species_gate(fasta_files, gate)
                if not downstream_files:
                    self.print_error(
                        "No samples passed species gate. Downstream modules skipped."
                    )
                    for name in ('MLST', 'Kaptive', 'ABRicate', 'AMRFinderPlus', 'MGE'):
                        results[name] = False
                    return results
            else:
                self.print_warning(
                    "FASTA QC failed; species gate unavailable. All samples flow through."
                )
            self.display_random_quote()

        if not skip_modules.get('mlst', False):
            self.display_module_header("MLST", self.module_subtitles['mlst'])
            results['MLST'], _ = self._run_timed(
                "MLST", self.run_mlst_analysis,
                downstream_files, output_dir, threads,
            )
            self.display_random_quote()

        if not skip_modules.get('kaptive', False):
            self.display_module_header("Kaptive", self.module_subtitles['kaptive'])
            results['Kaptive'], _ = self._run_timed(
                "Kaptive", self.run_kaptive_analysis,
                downstream_files, output_dir, threads,
            )
            self.display_random_quote()

        if not skip_modules.get('abricate', False):
            self.display_module_header("ABRicate", self.module_subtitles['abricate'])
            results['ABRicate'], _ = self._run_timed(
                "ABRicate", self.run_abricate_analysis,
                downstream_files, output_dir, threads,
                abricate_min_id, abricate_min_cov,
            )
            self.display_random_quote()

        if not skip_modules.get('amr', False):
            self.display_module_header(
                "AMRFinderPlus", self.module_subtitles['amr'])
            results['AMRFinderPlus'], _ = self._run_timed(
                "AMRFinderPlus", self.run_amr_analysis,
                downstream_files, output_dir, threads,
                amr_min_identity, amr_min_coverage,
                amr_skip_mutations, amr_force_update,
            )
            self.display_random_quote()

        if not skip_modules.get('mge', False):
            self.display_module_header("MGE", self.module_subtitles['mge'])
            results['MGE'], _ = self._run_timed(
                "MGE", self.run_mge_analysis,
                downstream_files, output_dir, threads,
            )
            self.display_random_quote()

        return results

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------
    def run_complete_analysis(
        self,
        input_path: str,
        output_dir: str,
        threads: int = 2,
        skip_modules: Optional[Dict[str, bool]] = None,
        skip_gene_centric: bool = False,
        skip_sample_centric: bool = False,
        update_amr_only: bool = False,
        force_update_amr_only: bool = False,
        pull_mlst_only: bool = False,
        update_mlst_only: bool = False,
        pull_mge_only: bool = False,
        amr_min_identity: Optional[float] = None,
        amr_min_coverage: Optional[float] = None,
        amr_skip_mutations: bool = False,
        amr_force_update: bool = False,
        abricate_min_id: int = 80,
        abricate_min_cov: int = 80,
        clean_output: bool = False,
        report_level: str = 'summary',
    ) -> bool:
        """Execute the complete Kleboscope pipeline (or a db-management mode)."""
        if skip_modules is None:
            skip_modules = {}
        self.report_level = report_level

        db_only_modes = [
            (update_amr_only, "update_amr", self.update_amr_database, {}),
            (force_update_amr_only, "update_amr", self.update_amr_database,
             {'force': True}),
            (pull_mlst_only, "mlst_pull", self.pull_mlst_database, {}),
            (update_mlst_only, "mlst_update", self.update_mlst_database, {}),
            (pull_mge_only, "mge_pull", self.pull_mge_database, {}),
        ]
        for enabled, label, fn, kwargs in db_only_modes:
            if enabled:
                temp_log = Path(self._mkdtemp(prefix=f"kleboscope_{label}_"))
                self.temp_dirs.append(str(temp_log))
                self.setup_logging(temp_log)
                ok = fn(**kwargs)
                if not self.keep_temp:
                    shutil.rmtree(temp_log, ignore_errors=True)
                return ok

        start_time = datetime.now()
        self.display_banner()

        output_path = Path(output_dir)
        if clean_output and output_path.exists():
            self.print_info(f"Cleaning output directory: {output_path}")
            shutil.rmtree(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        self.setup_logging(output_path)

        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found. Analysis stopped.")
            return False

        self.print_success(
            f"Starting analysis of {len(fasta_files)} K. pneumoniae genomes"
        )

        for sub in self.output_dirs.values():
            (output_path / sub).mkdir(exist_ok=True)

        if not skip_modules.get('mge', False) and self.get_mge_data_dir() is None:
            self.print_warning(
                "MGE is enabled but mobileOG-db is not installed. "
                "Run `kleboscope --pull-mge-db` first, or add --skip-mge."
            )

        plan = [
            ("FASTA QC",       not skip_modules.get('qc', False)),
            ("MLST",           not skip_modules.get('mlst', False)),
            ("Kaptive",        not skip_modules.get('kaptive', False)),
            ("ABRicate",       not skip_modules.get('abricate', False)),
            ("AMRFinderPlus",  not skip_modules.get('amr', False)),
            ("MGE",            not skip_modules.get('mge', False)),
            ("Gene-Centric",   not skip_gene_centric),
            ("Sample-Centric", not skip_sample_centric),
        ]
        self.print_header("Analysis Plan", "Modules to be executed")
        for analysis, enabled in plan:
            status = (f"{Color.BRIGHT_GREEN}✓ ENABLED{Color.RESET}" if enabled
                      else f"{Color.YELLOW}⏸ SKIPPED{Color.RESET}")
            print(f"   {status} - {analysis}")
        print(f"\n   Report level: {self.report_level}\n")

        analysis_results = self.run_sequential_analyses(
            fasta_files, output_path, threads, skip_modules,
            amr_min_identity, amr_min_coverage,
            amr_skip_mutations, amr_force_update,
            abricate_min_id, abricate_min_cov,
        )

        if not skip_gene_centric:
            self.display_module_header(
                "Gene-Centric Reporter", self.module_subtitles['gene_centric'])
            analysis_results['Gene-Centric'], _ = self._run_timed(
                "Gene-Centric Reporter",
                self.run_gene_centric_reporter, output_path,
            )
            self.display_random_quote()

        if not skip_sample_centric:
            self.display_module_header(
                "Sample-Centric Reporter", self.module_subtitles['sample_centric'])
            analysis_results['Sample-Centric'], _ = self._run_timed(
                "Sample-Centric Reporter",
                self.run_sample_centric_reporter, output_path,
            )
            self.display_random_quote()

        self.copy_summary_results_to_final_directory(output_path)
        self._cleanup_top_level_duplicates(output_path)
        self._finalize_output_by_report_level(output_path)

        elapsed = datetime.now() - start_time
        elapsed_str = str(elapsed).split('.')[0]

        self.print_header("Analysis Complete", f"Time elapsed: {elapsed_str}")

        if self.module_times:
            for label, secs in sorted(self.module_times.items(), key=lambda x: -x[1]):
                mins = int(secs // 60)
                rem = secs % 60
                time_str = f"{mins}m {rem:.1f}s" if mins else f"{rem:.1f}s"
                print(f"   {Color.BRIGHT_CYAN}{label:<28}{Color.RESET} {time_str}")

        successful = sum(analysis_results.values())
        total = len(analysis_results)
        if successful == total:
            self.print_success(f"All {total} analyses completed successfully.")
        else:
            self.print_warning(
                f"{successful}/{total} analyses completed successfully."
            )

        print("\n📁 FINAL OUTPUT:")
        for item in sorted(output_path.glob('*')):
            if item.is_dir():
                count = len(list(item.rglob('*')))
                print(f"   📂 {item.name}/ ({count} items)")
            else:
                print(f"   📄 {item.name}")

        print("\n📚 Citation:")
        for line in CITATION_LINES:
            print(f"   {line}")
        print()

        self.display_random_quote()

        if self.warnings:
            print(f"{Color.BRIGHT_YELLOW}Warnings encountered:{Color.RESET}")
            for w in self.warnings:
                self.print_info(f"  - {w}")
        if self.errors:
            print(f"{Color.BRIGHT_RED}Errors encountered:{Color.RESET}")
            for e in self.errors:
                self.print_info(f"  - {e}")

        return successful == total


def main():
    parser = argparse.ArgumentParser(
        description=(
            f"{Color.BRIGHT_YELLOW}Kleboscope: Complete K. pneumoniae "
            f"Typing Pipeline{Color.RESET}"
        ),
        formatter_class=ColoredHelpFormatter,
        epilog=f"""
{Color.BOLD}{Color.BRIGHT_GREEN}EXAMPLES:{Color.RESET}
  kleboscope -i genome.fna -o results/
  kleboscope -i "*.fna" -o batch_results --threads 8
  kleboscope -i "*.fna" -o results --report-level final
  kleboscope -i "*.fna" -o results --tmpdir /scratch/user/kleboscope_tmp

{Color.BOLD}{Color.BRIGHT_GREEN}DATABASE SETUP (First-time users):{Color.RESET}
  1. kleboscope --update-amr-db       (AMRFinderPlus database)
  2. kleboscope --pull-mlst-db        (MLST scheme, from GitHub)
  3. kleboscope --pull-mge-db         (mobileOG-db, ~300 MB — conda users only)
  4. abricate --setupdb               (ABRicate databases)

{Color.BRIGHT_YELLOW}MLST DATABASE MANAGEMENT:{Color.RESET}
  --pull-mlst-db     : Download pre-built P. aeruginosa MLST scheme from the
                       GitHub repo. No API key required.
  --update-mlst-db   : Refresh the MLST scheme directly from PubMLST using
                       mlstdb. Requires a one-time PubMLST API key setup:
                         1. Register at https://pubmlst.org/site-accounts
                         2. Generate an API key in your profile
                         3. Run: mlstdb connect --db pubmlst --api-key
                         4. Paste your key when prompted

{Color.BOLD}{Color.BRIGHT_GREEN}REPORT LEVELS:{Color.RESET}
  final   - Only Kleboscope_final_report/ is kept.
  summary - (default) Final report + batch summaries; per-sample folders pruned.
  full    - Everything is kept, including per-sample HTMLs.

{Color.BOLD}{Color.BRIGHT_GREEN}ENVIRONMENT OVERRIDES (HPC / Docker):{Color.RESET}
  KLEBOSCOPE_SHARED_DATA          Writable data root (default: ~/.local/share/kleboscope)
  KLEBOSCOPE_TMPDIR               Temp directory root (or use --tmpdir)
  KLEBOSCOPE_MLST_DB              MLST database directory
  KLEBOSCOPE_MLST_MODULE_PATH     MLST module root (bin/, db/, scripts/)
  KLEBOSCOPE_MLST_PULL_URL        Override MLST pull URL
  KLEBOSCOPE_MGE_DATA_DIR         mobileOG-db data directory
  KLEBOSCOPE_MGE_MODULE_PATH      MGE module root
  KLEBOSCOPE_MGE_PULL_URL         Override mobileOG-db pull URL

{Color.BOLD}{Color.BRIGHT_GREEN}SUPPORTED FASTA FORMATS:{Color.RESET} .fna, .fasta, .fa

{Color.BOLD}{Color.BRIGHT_GREEN}CONTACT:{Color.RESET}
  {AUTHOR_NAME} <{AUTHOR_EMAIL}>
  {AUTHOR_AFFILIATION}
"""
    )

    parser.add_argument('-i', '--input',
                        help='Input FASTA file(s) - can use glob patterns like "*.fna"')
    parser.add_argument('-o', '--output',
                        help='Output directory for all results')
    parser.add_argument('-t', '--threads', type=int, default=2,
                        help='Number of threads (default: 2)')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress all non-error console output')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Do not delete temporary directories (for debugging)')
    parser.add_argument('--tmpdir',
                        help='Override temporary directory root (HPC/Docker)')
    parser.add_argument('--clean-output', action='store_true',
                        help='Delete output directory before analysis')

    parser.add_argument('--report-level',
                        choices=['final', 'summary', 'full'],
                        default='summary',
                        help='Controls which files are kept after the run '
                             '(default: summary).')

    parser.add_argument('-v', '--version', action='version',
                        version=f'kleboscope {__version__}')

    parser.add_argument('--update-amr-db', action='store_true',
                        help='Update AMRFinderPlus database (incremental) and exit')
    parser.add_argument('--force-update-amr-db', action='store_true',
                        help='Force complete AMRFinderPlus database update and exit')
    parser.add_argument('--pull-mlst-db', action='store_true',
                        help='Install the Klebsiella MLST database '
                             '(runs database.sh or downloads prebuilt archive) and exit')
    parser.add_argument('--update-mlst-db', action='store_true',
                        help='Update the Klebsiella MLST database from PubMLST '
                             '(requires mlstdb) and exit')
    parser.add_argument('--pull-mge-db', action='store_true',
                        help='Download the mobileOG-db payload (~300 MB) from GitHub '
                             'and exit')

    parser.add_argument('--amr-min-identity', type=float,
                        help='Minimum identity for AMR hits (0..1)')
    parser.add_argument('--amr-min-coverage', type=float,
                        help='Minimum coverage for AMR hits (0..1)')
    parser.add_argument('--skip-amr-mutations', action='store_true',
                        help='Disable point mutation reporting in AMR')
    parser.add_argument('--amr-force-update', action='store_true',
                        help='Force AMR database update before analysis')

    parser.add_argument('--abricate-min-id', type=int, default=80,
                        help='Minimum identity for ABRicate hits (default: 80)')
    parser.add_argument('--abricate-min-cov', type=int, default=80,
                        help='Minimum coverage for ABRicate hits (default: 80)')

    skip = parser.add_argument_group('Skip options')
    skip.add_argument('--skip-qc', action='store_true', help='Skip FASTA QC')
    skip.add_argument('--skip-mlst', action='store_true', help='Skip MLST')
    skip.add_argument('--skip-kaptive', action='store_true', help='Skip Kaptive')
    skip.add_argument('--skip-abricate', action='store_true', help='Skip ABRicate')
    skip.add_argument('--skip-amr', action='store_true', help='Skip AMRFinderPlus')
    skip.add_argument('--skip-mge', action='store_true', help='Skip MGE profiling')
    skip.add_argument('--skip-gene-centric', action='store_true',
                      help='Skip gene-centric reporter')
    skip.add_argument('--skip-sample-centric', action='store_true',
                      help='Skip sample-centric reporter')

    args = parser.parse_args()

    if args.abricate_min_id < 0 or args.abricate_min_id > 100:
        parser.error("--abricate-min-id must be between 0 and 100")
    if args.abricate_min_cov < 0 or args.abricate_min_cov > 100:
        parser.error("--abricate-min-cov must be between 0 and 100")

    orch = KleboscopeOrchestrator(
        quiet=args.quiet,
        keep_temp=args.keep_temp,
        tmpdir=args.tmpdir,
    )

    if args.report_level == 'final' and args.skip_gene_centric and args.skip_sample_centric:
        parser.error(
            "--report-level final requires at least one reporter to run "
            "(--skip-gene-centric and --skip-sample-centric both set would "
            "produce an empty final report)."
        )

    try:
        ok = orch.run_complete_analysis(
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            skip_modules={
                'qc':       args.skip_qc,
                'mlst':     args.skip_mlst,
                'kaptive':  args.skip_kaptive,
                'abricate': args.skip_abricate,
                'amr':      args.skip_amr,
                'mge':      args.skip_mge,
            },
            skip_gene_centric=args.skip_gene_centric,
            skip_sample_centric=args.skip_sample_centric,
            update_amr_only=args.update_amr_db,
            force_update_amr_only=args.force_update_amr_db,
            pull_mlst_only=args.pull_mlst_db,
            update_mlst_only=args.update_mlst_db,
            pull_mge_only=args.pull_mge_db,
            amr_min_identity=args.amr_min_identity,
            amr_min_coverage=args.amr_min_coverage,
            amr_skip_mutations=args.skip_amr_mutations,
            amr_force_update=args.amr_force_update,
            abricate_min_id=args.abricate_min_id,
            abricate_min_cov=args.abricate_min_cov,
            clean_output=args.clean_output,
            report_level=args.report_level,
        )
        sys.exit(0 if ok else 1)
    except KeyboardInterrupt:
        orch.print_error("Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        orch.print_error(f"Critical error: {e}")
        if not args.quiet:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()