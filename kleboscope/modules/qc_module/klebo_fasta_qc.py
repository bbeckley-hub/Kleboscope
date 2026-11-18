#!/usr/bin/env python3
"""
Kleboscope FASTA QC - Quality Control for Klebsiella pneumoniae assemblies
with fastANI-based species confirmation.

Per-sample output (1 file each):
    <sample>/<sample>_fasta_qc_report.html   (passed or species-unknown only)

Batch-level outputs:
    kleboscope_fasta_qc_summary.tsv          (full audit — Yes / No / NA)
    kleboscope_fasta_qc_summary.html         (passing/unknown only + exclusion count)
    kleboscope_fasta_qc_summary.json         (full audit + species filter policy)
    kleboscope_species_passed.txt            (filenames to forward)
    kleboscope_species_excluded.txt          (filenames that failed species check)

Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Version: 2.0.0 (species gate, single per-sample HTML, three-file batch summary)
"""

import argparse
import glob
import json
import logging
import os
import re
import statistics
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


KLEBOSCOPE_BANNER = (
    "██╗  ██╗ ██╗      ███████╗ ██████╗  ██████╗ ███████╗  ██████╗ ██████╗  ██████╗ ███████╗\n"
    "██║ ██╔╝ ██║      ██╔════╝ ██╔══██╗██╔═══██╗██╔════╝ ██╔════╝██╔═══██╗ ██╔══██╗██╔════╝\n"
    "█████╔╝  ██║      █████╗   ██████╔╝██║   ██║███████╗ ██║     ██║   ██║ ██████╔╝█████╗\n"
    "██╔═██╗  ██║      ██╔══╝   ██╔══██╗██║   ██║╚════██║ ██║     ██║   ██║ ██╔═══╝ ██╔══╝\n"
    "██║  ██╗ ███████╗ ███████╗ ██████╔╝╚██████╔╝███████║ ╚██████╗╚██████╔╝ ██║     ███████╗\n"
    "╚═╝  ╚═╝ ╚══════╝ ╚══════╝ ╚═════╝  ╚═════╝ ╚══════╝  ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝"
)


def _render_footer(title: str, current_time: str) -> str:
    """Return the shared KLEBOSCOPE footer used across every module."""
    return f"""
<div class="footer">
    <p><strong>KLEBOSCOPE</strong> - {title}</p>
    <p class="timestamp">Generated: {current_time}</p>
    <div class="authorship">
        <p><strong>Technical Support & Inquiries:</strong></p>
        <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
        <p>Email: brownbeckley94@gmail.com</p>
        <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
    </div>
</div>"""


class KleboscopeFASTAQC:
    """FASTA QC for K. pneumoniae with fastANI species confirmation."""

    def __init__(self, cpus: Optional[int] = None, ref_dir: Optional[str] = None):
        """Initialise the analyzer.

        Args:
            cpus: CPU core budget. Defaults to all available cores.
            ref_dir: Directory containing reference genomes. Defaults to
                ``./ref_db`` next to this script.
        """
        self.logger = self._setup_logging()
        self.cpus = cpus or os.cpu_count() or 1

        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "Kleboscope turns genomic complexity into actionable insights for AMR surveillance.", "author": "Brown Beckley"},
        ]

        self.metadata = {
            "tool_name": "Kleboscope FASTA QC",
            "version": "1.3.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School - Department of Medical Biochemistry",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "biopython_version": self._get_biopython_version(),
        }

        self.thresholds = {
            "gc_normal": (55, 58),
            "short_seq": 100,
            "long_seq": 1_000_000,
            "max_n_run": 100,
            "max_homopolymer": 20,
            "ambiguous_critical": 5.0,
            "ambiguous_warning": 1.0,
        }

        if ref_dir is None:
            ref_dir = Path(__file__).resolve().parent / "ref_db"
        self.ref_dir = Path(ref_dir)
        self.species_refs = self._load_references()

        if self.species_refs and not any(
            self._is_kpneumoniae(sp) for sp in self.species_refs
        ):
            self.logger.warning(
                "No K. pneumoniae reference detected in ref_db/. "
                "Species check will flag ALL samples as non-K. pneumoniae. "
                "Add klebsiella_pneumoniae_*.fna to ref_db/ to correct this."
            )

    @staticmethod
    def _setup_logging() -> logging.Logger:
        """Return the module logger."""
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
        )
        return logging.getLogger(__name__)

    @staticmethod
    def _get_biopython_version() -> str:
        """Return the installed BioPython version string."""
        try:
            import Bio
            return Bio.__version__
        except Exception:
            return "Unknown"

    def _get_random_quote(self) -> Dict[str, str]:
        """Return a random science quote."""
        import random
        return random.choice(self.science_quotes)

    # ------------------------------------------------------------------
    # Species confirmation helpers
    # ------------------------------------------------------------------
    def _load_references(self) -> Dict[str, Path]:
        """Return a ``{species_name: path}`` mapping for every FASTA in ref_dir."""
        refs: Dict[str, Path] = {}
        if not self.ref_dir.exists():
            self.logger.warning(
                "Reference directory %s not found. Species check disabled.",
                self.ref_dir,
            )
            return refs

        for f in sorted(self.ref_dir.glob("*.fna")):
            stem = f.stem.replace("_", " ")
            stem = " ".join(stem.split())
            species = stem.title()
            refs[species] = f

        if refs:
            self.logger.info(
                "Loaded %d reference genome(s): %s",
                len(refs), ", ".join(refs.keys()),
            )
        return refs

    @staticmethod
    def _is_kpneumoniae(species_name: str) -> bool:
        """Return True if the species name refers to K. pneumoniae (not a close relative)."""
        s = species_name.lower()
        if "pneumoniae" not in s:
            return False
        return not any(
            x in s
            for x in ("quasipneumoniae", "variicola", "oxytoca", "michiganensis")
        )

    @staticmethod
    def _species_status(species_check: Optional[Dict]) -> str:
        """Return 'Yes', 'No', or 'NA' for the species confirmation status."""
        if not species_check or "error" in species_check:
            return "NA"
        return "Yes" if species_check.get("passed", False) else "No"

    def _run_species_ani(self, query_fasta: str, threads: int = 4) -> Optional[Dict[str, Any]]:
        """Run fastANI against every reference genome and return the best match.

        Pass criteria: best match is K. pneumoniae AND ANI >= 95%.
        Contamination is flagged only when the second-best match is
        non-K. pneumoniae, ANI > 90%, and within 5% of the best hit.
        """
        if not self.species_refs:
            return None

        results: List[tuple] = []
        for species, ref_path in self.species_refs.items():
            tmp_fd, tmp_path = tempfile.mkstemp(suffix=".txt")
            os.close(tmp_fd)
            try:
                cmd = [
                    "fastANI",
                    "-q", query_fasta,
                    "-r", str(ref_path),
                    "-o", tmp_path,
                    "-t", str(threads),
                ]
                subprocess.run(cmd, capture_output=True, text=True, check=True)
                with open(tmp_path, "r", encoding="utf-8") as handle:
                    line = handle.readline().strip()
                if line:
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        try:
                            ani = float(parts[2])
                            results.append((species, ani))
                        except ValueError:
                            self.logger.warning("Could not parse ANI from: %s", line)
            except subprocess.CalledProcessError:
                self.logger.warning(
                    "fastANI failed for %s vs %s",
                    Path(query_fasta).name, species,
                )
            except FileNotFoundError:
                self.logger.error(
                    "fastANI executable not found on PATH. Species check disabled."
                )
                return None
            finally:
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

        if not results:
            return None

        results.sort(key=lambda x: x[1], reverse=True)
        best_species, best_ani = results[0]
        second_species = results[1][0] if len(results) > 1 else None
        second_ani = results[1][1] if len(results) > 1 else None

        passed = best_ani >= 95.0 and self._is_kpneumoniae(best_species)

        contamination = False
        if (
            second_ani is not None
            and second_ani > 90.0
            and (best_ani - second_ani) < 5.0
            and not self._is_kpneumoniae(second_species)
        ):
            contamination = True

        return {
            "best_match": best_species,
            "ani_percent": round(best_ani, 2),
            "passed": passed,
            "contamination_suspected": contamination,
            "second_best_match": second_species,
            "second_best_ani": round(second_ani, 2) if second_ani is not None else None,
            "all_matches": {sp: round(ani, 2) for sp, ani in results},
        }

    # ------------------------------------------------------------------
    # Per-file analysis
    # ------------------------------------------------------------------
    def analyze_file(self, fasta_file: str) -> Dict[str, Any]:
        """Run every QC metric plus species confirmation on one FASTA file."""
        try:
            self.logger.info("Analyzing %s ...", os.path.basename(fasta_file))

            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            if not sequences:
                return {
                    "filename": os.path.basename(fasta_file),
                    "filepath": fasta_file,
                    "status": "error",
                    "error": "No sequences found",
                }

            seq_lengths = [len(seq) for seq in sequences]
            total_length = sum(seq_lengths)
            sorted_lengths = sorted(seq_lengths, reverse=True)

            n50 = self._calculate_nx(sorted_lengths, total_length, 50)
            n75 = self._calculate_nx(sorted_lengths, total_length, 75)
            n90 = self._calculate_nx(sorted_lengths, total_length, 90)
            l50 = self._calculate_lx(sorted_lengths, total_length, 50)
            l75 = self._calculate_lx(sorted_lengths, total_length, 75)
            l90 = self._calculate_lx(sorted_lengths, total_length, 90)

            base_counts: Counter = Counter()
            gc_contents: List[float] = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                base_counts.update(seq_str)
                gc_contents.append(gc_fraction(seq_str) * 100)

            total_bases = sum(base_counts.values())
            a_percent = (base_counts.get("A", 0) / total_bases) * 100 if total_bases else 0
            t_percent = (base_counts.get("T", 0) / total_bases) * 100 if total_bases else 0
            g_percent = (base_counts.get("G", 0) / total_bases) * 100 if total_bases else 0
            c_percent = (base_counts.get("C", 0) / total_bases) * 100 if total_bases else 0
            gc_percent = g_percent + c_percent
            at_percent = a_percent + t_percent

            ambiguous_bases = sum(
                base_counts.get(b, 0)
                for b in ("N", "Y", "R", "W", "S", "K", "M", "B", "D", "H", "V")
            )
            ambiguous_percent = (
                (ambiguous_bases / total_bases) * 100 if total_bases else 0
            )

            sequences_with_n = sum(
                1 for seq in sequences if "N" in str(seq.seq).upper()
            )
            total_n_bases = base_counts.get("N", 0)

            max_n_run = 0
            n_runs: List[int] = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                for match in re.finditer(r"N+", seq_str):
                    run_len = len(match.group())
                    n_runs.append(run_len)
                    if run_len > max_n_run:
                        max_n_run = run_len

            homopolymers: List[Dict[str, Any]] = []
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                for base in ("A", "T", "G", "C"):
                    for match in re.finditer(f"{base}+", seq_str):
                        if len(match.group()) > 4:
                            homopolymers.append({
                                "base": base,
                                "length": len(match.group()),
                                "position": match.start(),
                            })
            max_homopolymer = (
                max(h["length"] for h in homopolymers) if homopolymers else 0
            )

            seq_hashes = set()
            duplicate_sequences = 0
            for seq in sequences:
                h = hash(str(seq.seq))
                if h in seq_hashes:
                    duplicate_sequences += 1
                else:
                    seq_hashes.add(h)

            short_sequences = sum(
                1 for length in seq_lengths if length < self.thresholds["short_seq"]
            )
            long_sequences = sum(
                1 for length in seq_lengths if length > self.thresholds["long_seq"]
            )
            length_distribution = self._create_length_bins(seq_lengths)
            kp_status = self._check_kpneumoniae_specific(
                gc_percent, gc_contents, total_length, seq_lengths
            )

            species_result = self._run_species_ani(fasta_file, self.cpus)
            species_check = (
                species_result
                if species_result
                else {"error": "Species check not available"}
            )
            species_status = self._species_status(species_check)

            results: Dict[str, Any] = {
                "filename": os.path.basename(fasta_file),
                "filepath": fasta_file,
                "file_size_mb": os.path.getsize(fasta_file) / (1024 ** 2),
                "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "status": "success",
                "total_sequences": len(sequences),
                "total_length": total_length,
                "total_bases": total_bases,
                "longest_sequence": max(seq_lengths),
                "shortest_sequence": min(seq_lengths),
                "mean_length": statistics.mean(seq_lengths),
                "median_length": statistics.median(seq_lengths),
                "n50": n50,
                "n75": n75,
                "n90": n90,
                "l50": l50,
                "l75": l75,
                "l90": l90,
                "gc_percent": gc_percent,
                "at_percent": at_percent,
                "a_percent": a_percent,
                "t_percent": t_percent,
                "g_percent": g_percent,
                "c_percent": c_percent,
                "ambiguous_percent": ambiguous_percent,
                "sequences_with_n": sequences_with_n,
                "total_n_bases": total_n_bases,
                "max_n_run": max_n_run,
                "homopolymers_count": len(homopolymers),
                "max_homopolymer": max_homopolymer,
                "duplicate_sequences": duplicate_sequences,
                "short_sequences": short_sequences,
                "long_sequences": long_sequences,
                "gc_distribution": {
                    "mean": statistics.mean(gc_contents),
                    "median": statistics.median(gc_contents),
                    "min": min(gc_contents),
                    "max": max(gc_contents),
                },
                "length_distribution": length_distribution,
                "base_counts": dict(base_counts),
                "n_runs": n_runs,
                "kp_status": kp_status,
                "species_check": species_check,
                "species_status": species_status,
                "warnings": self._generate_warnings(
                    gc_percent, ambiguous_percent, max_n_run, max_homopolymer,
                    short_sequences, long_sequences, duplicate_sequences,
                    len(sequences), kp_status, species_result,
                ),
            }

            self.logger.info(
                "%s: %d sequences, %s bp, N50 %s, GC %.1f%%, species %s",
                os.path.basename(fasta_file), len(sequences),
                f"{total_length:,}", f"{n50:,}", gc_percent, species_status,
            )
            return results

        except Exception as exc:
            self.logger.error("Error analyzing %s: %s", fasta_file, exc)
            return {
                "filename": os.path.basename(fasta_file),
                "status": "error",
                "error": str(exc),
            }

    def _check_kpneumoniae_specific(
        self,
        gc_percent: float,
        gc_contents: List[float],
        total_length: int,
        seq_lengths: List[int],
    ) -> Dict[str, Any]:
        """Return K. pneumoniae expected-range assessment."""
        genome_size_ok = 5_000_000 <= total_length <= 5_600_000
        gc_ok = self.thresholds["gc_normal"][0] <= gc_percent <= self.thresholds["gc_normal"][1]

        if len(seq_lengths) <= 50:
            assembly_quality = "Good"
        elif len(seq_lengths) <= 100:
            assembly_quality = "Moderate"
        else:
            assembly_quality = "Fragmented"

        return {
            "genome_size_mbp": total_length / 1_000_000,
            "genome_size_status": "Normal" if genome_size_ok else "Atypical",
            "gc_status": "Normal" if gc_ok else "Atypical",
            "assembly_quality": assembly_quality,
            "contig_count": len(seq_lengths),
            "expected_gc_range": f"{self.thresholds['gc_normal'][0]}-{self.thresholds['gc_normal'][1]}%",
            "expected_genome_size": "5.0-5.6 Mbp",
        }

    @staticmethod
    def _calculate_nx(sorted_lengths: List[int], total_length: int, x: int) -> int:
        """Return Nx statistic for x in [0, 100]."""
        if not sorted_lengths:
            return 0
        target = total_length * (x / 100)
        cumulative = 0
        for length in sorted_lengths:
            cumulative += length
            if cumulative >= target:
                return length
        return sorted_lengths[-1]

    @staticmethod
    def _calculate_lx(sorted_lengths: List[int], total_length: int, x: int) -> int:
        """Return Lx statistic for x in [0, 100]."""
        if not sorted_lengths:
            return 0
        target = total_length * (x / 100)
        cumulative = 0
        for i, length in enumerate(sorted_lengths, 1):
            cumulative += length
            if cumulative >= target:
                return i
        return len(sorted_lengths)

    @staticmethod
    def _create_length_bins(lengths: List[int]) -> Dict[str, int]:
        """Return per-bin counts for the given sequence lengths."""
        bins = {
            "< 100 bp": 0,
            "100-500 bp": 0,
            "500-1k bp": 0,
            "1k-5k bp": 0,
            "5k-10k bp": 0,
            "10k-50k bp": 0,
            "50k-100k bp": 0,
            "100k-500k bp": 0,
            "500k-1M bp": 0,
            "> 1M bp": 0,
        }
        for length in lengths:
            if length < 100:
                bins["< 100 bp"] += 1
            elif length < 500:
                bins["100-500 bp"] += 1
            elif length < 1_000:
                bins["500-1k bp"] += 1
            elif length < 5_000:
                bins["1k-5k bp"] += 1
            elif length < 10_000:
                bins["5k-10k bp"] += 1
            elif length < 50_000:
                bins["10k-50k bp"] += 1
            elif length < 100_000:
                bins["50k-100k bp"] += 1
            elif length < 500_000:
                bins["100k-500k bp"] += 1
            elif length < 1_000_000:
                bins["500k-1M bp"] += 1
            else:
                bins["> 1M bp"] += 1
        return bins

    def _generate_warnings(
        self,
        gc_percent: float,
        ambiguous_percent: float,
        max_n_run: int,
        max_homopolymer: int,
        short_sequences: int,
        long_sequences: int,
        duplicate_sequences: int,
        total_sequences: int,
        kp_status: Dict[str, Any],
        species_result: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, str]]:
        """Return a list of warning dictionaries for the sample."""
        warnings: List[Dict[str, str]] = []

        low, high = self.thresholds["gc_normal"]
        if gc_percent < low:
            warnings.append({
                "level": "warning",
                "message": (
                    f"Low GC content ({gc_percent:.1f}%) - below "
                    f"K. pneumoniae range ({low}-{high}%)"
                ),
            })
        elif gc_percent > high:
            warnings.append({
                "level": "warning",
                "message": (
                    f"High GC content ({gc_percent:.1f}%) - above "
                    f"K. pneumoniae range ({low}-{high}%)"
                ),
            })

        if ambiguous_percent > self.thresholds["ambiguous_critical"]:
            warnings.append({
                "level": "danger",
                "message": (
                    f"High ambiguous bases ({ambiguous_percent:.2f}%) - "
                    "may indicate poor quality"
                ),
            })
        elif ambiguous_percent > self.thresholds["ambiguous_warning"]:
            warnings.append({
                "level": "warning",
                "message": f"Elevated ambiguous bases ({ambiguous_percent:.2f}%)",
            })

        if max_n_run > 100:
            warnings.append({
                "level": "danger",
                "message": (
                    f"Very long N-run detected ({max_n_run} bases) - "
                    "may indicate assembly gaps"
                ),
            })
        elif max_n_run > 10:
            warnings.append({
                "level": "warning",
                "message": f"Long N-run detected ({max_n_run} bases)",
            })

        if max_homopolymer > 20:
            warnings.append({
                "level": "danger",
                "message": (
                    f"Very long homopolymer ({max_homopolymer} bases) - "
                    "may cause sequencing errors"
                ),
            })
        elif max_homopolymer > 10:
            warnings.append({
                "level": "warning",
                "message": f"Long homopolymer ({max_homopolymer} bases)",
            })

        if short_sequences > total_sequences * 0.5:
            warnings.append({
                "level": "danger",
                "message": (
                    f"Many short sequences ({short_sequences}) - "
                    "may indicate poor assembly"
                ),
            })
        elif short_sequences > total_sequences * 0.1:
            warnings.append({
                "level": "warning",
                "message": f"Some short sequences ({short_sequences})",
            })

        if long_sequences > 0:
            warnings.append({
                "level": "warning",
                "message": (
                    f"Very long sequences detected ({long_sequences}) - "
                    "may indicate contamination"
                ),
            })

        dup_pct = (duplicate_sequences / total_sequences) * 100 if total_sequences else 0
        if dup_pct > 10:
            warnings.append({
                "level": "warning",
                "message": (
                    f"Duplicate sequences detected "
                    f"({duplicate_sequences}, {dup_pct:.1f}%)"
                ),
            })

        if kp_status.get("genome_size_status") == "Atypical":
            warnings.append({
                "level": "warning",
                "message": (
                    f"Atypical genome size "
                    f"({kp_status['genome_size_mbp']:.2f} Mbp) for "
                    f"K. pneumoniae (expected: "
                    f"{kp_status['expected_genome_size']})"
                ),
            })

        if kp_status.get("assembly_quality") == "Fragmented":
            warnings.append({
                "level": "warning",
                "message": (
                    f"High contig count ({kp_status['contig_count']}) - "
                    "consider additional assembly polishing"
                ),
            })

        if species_result and "error" not in species_result:
            if not species_result.get("passed", False):
                warnings.append({
                    "level": "danger",
                    "message": (
                        f"Species mismatch: best match is "
                        f"{species_result.get('best_match', 'ND')} "
                        f"(ANI {species_result.get('ani_percent', 'ND')}%). "
                        f"Expected K. pneumoniae (ANI >= 95%). "
                        f"Sample will be EXCLUDED from downstream modules."
                    ),
                })
            if species_result.get("contamination_suspected", False):
                warnings.append({
                    "level": "danger",
                    "message": (
                        f"Contamination suspected: second best match is "
                        f"{species_result.get('second_best_match', 'unknown')} "
                        f"(ANI {species_result.get('second_best_ani', 'ND')}%), "
                        f"which is not K. pneumoniae."
                    ),
                })

        return warnings

    # ------------------------------------------------------------------
    # Per-sample HTML report
    # ------------------------------------------------------------------
    def create_individual_html_report(self, results: Dict[str, Any], output_dir: str) -> str:
        """Write the per-sample FASTA QC HTML report and return its path."""
        filename_no_ext = Path(results["filename"]).stem
        sample_dir = os.path.join(output_dir, filename_no_ext)
        os.makedirs(sample_dir, exist_ok=True)

        html_file = os.path.join(
            sample_dir, f"{filename_no_ext}_fasta_qc_report.html"
        )
        random_quote = self._get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        species_status = results.get("species_status", "NA")
        species_check = results.get("species_check", {}) or {}
        kp_status = results.get("kp_status", {})

        warnings_html = ""
        for warning in results.get("warnings", []):
            color = "#dc2626" if warning["level"] == "danger" else "#f59e0b"
            warnings_html += (
                f'<div style="background: rgba(220,38,38,0.1); padding: 10px; '
                f'border-radius: 6px; margin: 5px 0; border-left: 4px solid {color};">'
                f'<div style="font-weight: bold; color: {color};">'
                f'{warning["level"].upper()}</div>'
                f'<div>{warning["message"]}</div></div>'
            )

        exclusion_banner = ""
        if species_status == "No":
            exclusion_banner = (
                '<div class="report-section" style="border-left: 4px solid #dc2626; '
                'background: #fef2f2;">'
                '<h2 style="color: #dc2626; border-bottom: 3px solid #dc2626;">'
                'Sample Excluded from Downstream Analysis</h2>'
                '<p style="color: #991b1b; font-weight: bold;">This sample failed '
                'the fastANI species confirmation (best match is not '
                'K. pneumoniae or ANI is below 95%). It will NOT be forwarded '
                'to MLST, Kaptive, AMR, ABRicate, or MGE modules.</p></div>'
            )

        composition_html = f"""
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px; margin-top: 15px;">
            <div style="background: linear-gradient(135deg, #10b981 0%, #059669 100%); padding: 15px; border-radius: 8px; text-align: center; color: white;">
                <div style="font-size: 20px; font-weight: bold;">{results['gc_percent']:.1f}%</div>
                <div style="font-size: 12px; opacity: 0.9;">GC Content</div>
            </div>
            <div style="background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); padding: 15px; border-radius: 8px; text-align: center; color: white;">
                <div style="font-size: 20px; font-weight: bold;">{results['at_percent']:.1f}%</div>
                <div style="font-size: 12px; opacity: 0.9;">AT Content</div>
            </div>
            <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                <div style="font-size: 20px; font-weight: bold;">{results['a_percent']:.1f}%</div>
                <div style="font-size: 12px;">Adenine (A)</div>
            </div>
            <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                <div style="font-size: 20px; font-weight: bold;">{results['t_percent']:.1f}%</div>
                <div style="font-size: 12px;">Thymine (T)</div>
            </div>
            <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                <div style="font-size: 20px; font-weight: bold;">{results['g_percent']:.1f}%</div>
                <div style="font-size: 12px;">Guanine (G)</div>
            </div>
            <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; text-align: center; color: #0369a1; border: 1px solid #bae6fd;">
                <div style="font-size: 20px; font-weight: bold;">{results['c_percent']:.1f}%</div>
                <div style="font-size: 12px;">Cytosine (C)</div>
            </div>
        </div>
        """

        length_dist_html = ""
        for length_range, count in results["length_distribution"].items():
            pct = (count / results["total_sequences"] * 100) if results["total_sequences"] else 0
            length_dist_html += (
                f"<tr><td>{length_range}</td>"
                f"<td>{count:,}</td>"
                f"<td>{pct:.1f}%</td></tr>"
            )

        kp_html = ""
        if kp_status:
            size_color = "#10b981" if kp_status["genome_size_status"] == "Normal" else "#f59e0b"
            gc_color = "#10b981" if kp_status["gc_status"] == "Normal" else "#f59e0b"
            assembly_color = (
                "#10b981" if kp_status["assembly_quality"] == "Good"
                else "#f59e0b" if kp_status["assembly_quality"] == "Moderate"
                else "#dc2626"
            )
            kp_html = f"""
            <div style="margin-top: 20px; padding: 20px; background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%); border-radius: 10px; color: white;">
                <h3 style="color: white; margin-bottom: 15px;">Klebsiella pneumoniae Specific Analysis</h3>
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px;">
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Genome Size</div>
                        <div style="font-size: 18px; font-weight: bold; color: {size_color};">{kp_status['genome_size_mbp']:.2f} Mbp</div>
                        <div style="font-size: 11px; opacity: 0.7;">Status: {kp_status['genome_size_status']}</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">GC Content</div>
                        <div style="font-size: 18px; font-weight: bold; color: {gc_color};">{results['gc_percent']:.1f}%</div>
                        <div style="font-size: 11px; opacity: 0.7;">Expected: {kp_status['expected_gc_range']}</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Assembly Quality</div>
                        <div style="font-size: 18px; font-weight: bold; color: {assembly_color};">{kp_status['assembly_quality']}</div>
                        <div style="font-size: 11px; opacity: 0.7;">Contigs: {kp_status['contig_count']}</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Expected Size</div>
                        <div style="font-size: 18px; font-weight: bold;">{kp_status['expected_genome_size']}</div>
                        <div style="font-size: 11px; opacity: 0.7;">Typical K. pneumoniae</div>
                    </div>
                </div>
            </div>
            """

        species_html = ""
        if species_check and "error" not in species_check:
            pass_color = "#10b981" if species_check.get("passed") else "#dc2626"
            contam_color = "#f59e0b" if species_check.get("contamination_suspected") else "#10b981"
            species_html = f"""
            <div style="margin-top: 20px; padding: 20px; background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%); border-radius: 10px; color: white;">
                <h3 style="color: white; margin-bottom: 15px;">Species Confirmation (fastANI)</h3>
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px;">
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Best Match</div>
                        <div style="font-size: 18px; font-weight: bold;">{species_check.get('best_match', 'ND')}</div>
                        <div style="font-size: 11px;">ANI: {species_check.get('ani_percent', 'ND')}%</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">K. pneumoniae Confirmed</div>
                        <div style="font-size: 18px; font-weight: bold; color: {pass_color};">{'Yes' if species_check.get('passed') else 'No'}</div>
                        <div style="font-size: 11px;">Threshold >= 95%</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Contamination Suspected</div>
                        <div style="font-size: 18px; font-weight: bold; color: {contam_color};">{'Yes' if species_check.get('contamination_suspected') else 'No'}</div>
                        <div style="font-size: 11px;">Second best &gt;90% and not K. pneumoniae</div>
                    </div>
                    <div>
                        <div style="font-size: 12px; opacity: 0.8;">Second Best</div>
                        <div style="font-size: 18px; font-weight: bold;">{species_check.get('second_best_match', 'ND') or 'ND'}</div>
                        <div style="font-size: 11px;">ANI: {species_check.get('second_best_ani', 'ND') if species_check.get('second_best_ani') is not None else 'ND'}%</div>
                    </div>
                </div>
                <details style="margin-top: 15px;">
                    <summary style="cursor: pointer; color: #3b82f6;">All ANI values</summary>
                    <ul>
            """
            for sp, ani in species_check.get("all_matches", {}).items():
                species_html += f"<li>{sp}: {ani}%</li>"
            species_html += "</ul></details></div>"
        elif species_status == "NA":
            species_html = (
                '<div style="margin-top: 20px; padding: 20px; background: '
                'linear-gradient(135deg, #1e293b 0%, #0f172a 100%); '
                'border-radius: 10px; color: white;">'
                '<h3 style="color: white; margin-bottom: 15px;">Species Confirmation (fastANI)</h3>'
                '<p style="color: #fbbf24; font-weight: bold;">'
                'Species check could not be performed (fastANI not installed or '
                'ref_db/ empty). Sample will be retained in downstream modules '
                'because species status is unknown.</p></div>'
            )

        footer = _render_footer("FASTA QC Analysis Module", current_time)

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KLEBOSCOPE - {results['filename']} FASTA QC Report</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color: #ffffff; padding: 20px; min-height: 100vh;
}}
.container {{ max-width: 1600px; margin: 0 auto; }}
.ascii-container {{
    background: rgba(0, 0, 0, 0.7); padding: 20px; border-radius: 15px;
    margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
    border: 2px solid rgba(0, 119, 190, 0.3); overflow-x: auto; text-align: center;
}}
.ascii-art {{
    font-family: 'Cascadia Mono', 'Consolas', 'Monaco', 'Menlo',
                 'DejaVu Sans Mono', 'Liberation Mono', 'Courier New', monospace;
    font-size: 12px; line-height: 1.05; white-space: pre;
    color: #00d4ff; text-shadow: 0 0 8px rgba(0, 212, 255, 0.4);
    display: inline-block; padding: 8px 0; text-align: left;
    font-variant-ligatures: none;
    font-feature-settings: "liga" 0, "calt" 0, "kern" 0;
    font-kerning: none; letter-spacing: 0;
    text-rendering: geometricPrecision;
}}
.quote-container {{
    background: rgba(255, 255, 255, 0.1); backdrop-filter: blur(10px);
    padding: 20px; border-radius: 10px; margin-bottom: 30px; text-align: center;
    min-height: 100px; display: flex; flex-direction: column;
    justify-content: center; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
    border: 1px solid rgba(255, 255, 255, 0.2); transition: opacity 0.5s ease-in-out;
}}
.quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
.quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
.report-section {{
    background: rgba(255, 255, 255, 0.95); color: #1f2937;
    padding: 25px; border-radius: 10px; margin-bottom: 20px;
    box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
}}
.report-section h2 {{
    color: #0b3d5f; border-bottom: 3px solid #0077be;
    padding-bottom: 10px; margin-bottom: 20px; font-size: 24px;
}}
.metrics-grid {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 20px; margin-top: 15px;
}}
.metric-card {{
    background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
    color: white; padding: 20px; border-radius: 8px;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
}}
.metric-label {{ font-size: 13px; opacity: 0.9; margin-bottom: 5px; text-transform: uppercase; letter-spacing: 0.5px; }}
.metric-value {{ font-size: 22px; font-weight: bold; word-break: break-word; }}
.stat-table {{
    width: 100%; border-collapse: collapse; margin: 20px 0;
    background: white; border-radius: 8px; overflow: hidden;
    box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
}}
.stat-table th {{
    background: #005a8c; color: white; padding: 12px;
    text-align: left; font-weight: bold;
}}
.stat-table td {{ padding: 10px; border-bottom: 1px solid #e5e7eb; }}
.stat-table tr:nth-child(even) {{ background-color: #f8fafc; }}
.stat-table tr:hover {{ background-color: #e0f2fe; }}
.footer {{
    text-align: center; margin-top: 30px; padding: 20px;
    background: rgba(0, 0, 0, 0.3); border-radius: 10px; font-size: 14px;
}}
.timestamp {{ color: #fbbf24; font-weight: bold; }}
.authorship {{
    margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1);
    border-radius: 8px; font-size: 12px;
}}
@media (max-width: 768px) {{
    .ascii-art {{ font-size: 6px; }}
    .metrics-grid {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>
<div class="container">
<div class="ascii-container">
<div class="ascii-art">{KLEBOSCOPE_BANNER}</div>
</div>

<div class="quote-container" id="quoteContainer">
    <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
    <div class="quote-author" id="quoteAuthor">- {random_quote['author']}</div>
</div>

{exclusion_banner}

<div class="report-section">
    <h2>Sample Information</h2>
    <div class="metrics-grid">
        <div class="metric-card"><div class="metric-label">Sample</div><div class="metric-value" style="font-size: 16px;">{results['filename']}</div></div>
        <div class="metric-card"><div class="metric-label">Analysis Date</div><div class="metric-value" style="font-size: 14px;">{results['analysis_date']}</div></div>
        <div class="metric-card"><div class="metric-label">File Size</div><div class="metric-value">{results['file_size_mb']:.2f} MB</div></div>
    </div>
</div>

<div class="report-section">
    <h2>FASTA QC Summary</h2>
    <div class="metrics-grid">
        <div class="metric-card"><div class="metric-label">Total Sequences</div><div class="metric-value">{results['total_sequences']:,}</div></div>
        <div class="metric-card"><div class="metric-label">Total Length</div><div class="metric-value">{results['total_length']:,} bp</div></div>
        <div class="metric-card"><div class="metric-label">N50</div><div class="metric-value">{results['n50']:,} bp</div></div>
        <div class="metric-card"><div class="metric-label">GC Content</div><div class="metric-value">{results['gc_percent']:.1f}%</div></div>
        <div class="metric-card"><div class="metric-label">AT Content</div><div class="metric-value">{results['at_percent']:.1f}%</div></div>
        <div class="metric-card"><div class="metric-label">Ambiguous Bases</div><div class="metric-value">{results['ambiguous_percent']:.2f}%</div></div>
    </div>
</div>

<div class="report-section">
    <h2>Basic Statistics</h2>
    <table class="stat-table">
        <thead><tr><th>Metric</th><th>Value</th><th>Description</th></tr></thead>
        <tbody>
            <tr><td><strong>Total Sequences</strong></td><td>{results['total_sequences']:,}</td><td>Number of sequences</td></tr>
            <tr><td><strong>Total Length</strong></td><td>{results['total_length']:,} bp</td><td>Sum of all sequence lengths</td></tr>
            <tr><td><strong>Total Bases</strong></td><td>{results['total_bases']:,} bp</td><td>Total bases counted</td></tr>
            <tr><td><strong>Longest Sequence</strong></td><td>{results['longest_sequence']:,} bp</td><td>Longest sequence</td></tr>
            <tr><td><strong>Shortest Sequence</strong></td><td>{results['shortest_sequence']:,} bp</td><td>Shortest sequence</td></tr>
            <tr><td><strong>Mean Length</strong></td><td>{results['mean_length']:,.0f} bp</td><td>Average length</td></tr>
            <tr><td><strong>Median Length</strong></td><td>{results['median_length']:,} bp</td><td>Median length</td></tr>
            <tr><td><strong>N50</strong></td><td>{results['n50']:,} bp</td><td>50% of bases in longer sequences</td></tr>
            <tr><td><strong>N75</strong></td><td>{results['n75']:,} bp</td><td>75% of bases in longer sequences</td></tr>
            <tr><td><strong>N90</strong></td><td>{results['n90']:,} bp</td><td>90% of bases in longer sequences</td></tr>
            <tr><td><strong>L50</strong></td><td>{results['l50']:,}</td><td>Count at N50</td></tr>
            <tr><td><strong>L75</strong></td><td>{results['l75']:,}</td><td>Count at N75</td></tr>
            <tr><td><strong>L90</strong></td><td>{results['l90']:,}</td><td>Count at N90</td></tr>
        </tbody>
    </table>
</div>

<div class="report-section">
    <h2>Nucleotide Composition</h2>
    {composition_html}
</div>

<div class="report-section">
    <h2>Length Distribution</h2>
    <table class="stat-table">
        <thead><tr><th>Length Range</th><th>Count</th><th>Percentage</th></tr></thead>
        <tbody>{length_dist_html}</tbody>
    </table>
</div>

{kp_html}
{species_html}

<div class="report-section">
    <h2>Quality Warnings</h2>
    {warnings_html if warnings_html else '<p style="color: #10b981; font-weight: bold;">No quality warnings detected</p>'}
</div>

<div class="report-section">
    <h2>Additional Statistics</h2>
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-top: 15px;">
        <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;"><div style="font-size: 16px; font-weight: bold;">{results['sequences_with_n']:,}</div><div style="font-size: 12px;">Sequences with Ns</div></div>
        <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;"><div style="font-size: 16px; font-weight: bold;">{results['total_n_bases']:,}</div><div style="font-size: 12px;">Total N Bases</div></div>
        <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;"><div style="font-size: 16px; font-weight: bold;">{results['max_n_run']:,}</div><div style="font-size: 12px;">Max N-run</div></div>
        <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;"><div style="font-size: 16px; font-weight: bold;">{results['homopolymers_count']:,}</div><div style="font-size: 12px;">Homopolymers</div></div>
        <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;"><div style="font-size: 16px; font-weight: bold;">{results['max_homopolymer']:,}</div><div style="font-size: 12px;">Max Homopolymer</div></div>
        <div style="background: #f0f9ff; padding: 15px; border-radius: 8px; color: #0369a1; border: 1px solid #bae6fd;"><div style="font-size: 16px; font-weight: bold;">{results['duplicate_sequences']:,}</div><div style="font-size: 12px;">Duplicate Sequences</div></div>
    </div>
</div>
{footer}
</div>

<script>
const quotes = {json.dumps(self.science_quotes)};
const quoteContainer = document.getElementById('quoteContainer');
const quoteText = document.getElementById('quoteText');
const quoteAuthor = document.getElementById('quoteAuthor');
function displayQuote() {{
    quoteContainer.style.opacity = '0';
    setTimeout(() => {{
        const q = quotes[Math.floor(Math.random() * quotes.length)];
        quoteText.textContent = '"' + q.text + '"';
        quoteAuthor.textContent = '- ' + q.author;
        quoteContainer.style.opacity = '1';
    }}, 500);
}}
displayQuote();
setInterval(displayQuote, 10000);
</script>
</body>
</html>"""

        with open(html_file, "w", encoding="utf-8") as handle:
            handle.write(html_content)

        self.logger.info("HTML report written: %s", html_file)
        return html_file

    # ------------------------------------------------------------------
    # Batch summary
    # ------------------------------------------------------------------
    def create_summary_report(
        self, all_results: List[Dict[str, Any]], output_dir: str
    ) -> None:
        """Write the batch TSV, JSON, and HTML summary plus species filter files."""
        successful = [r for r in all_results if r.get("status") == "success"]
        if not successful:
            self.logger.warning("No successful analyses; skipping summary outputs.")
            return

        tsv_file = os.path.join(output_dir, "kleboscope_fasta_qc_summary.tsv")
        self._write_summary_tsv(successful, tsv_file)

        json_summary = self._build_json_summary(successful)
        json_file = os.path.join(output_dir, "kleboscope_fasta_qc_summary.json")
        with open(json_file, "w", encoding="utf-8") as handle:
            json.dump(json_summary, handle, indent=2, default=str)

        display = [
            r for r in successful
            if r.get("species_status", self._species_status(r.get("species_check"))) != "No"
        ]
        excluded_count = len(successful) - len(display)

        if display:
            html_file = os.path.join(output_dir, "kleboscope_fasta_qc_summary.html")
            self._write_summary_html(display, html_file, excluded_count)
            self.logger.info("HTML summary written: %s", html_file)
        else:
            self.logger.warning(
                "All samples failed species check; HTML summary not generated. "
                "See kleboscope_species_excluded.txt."
            )

        self._write_species_filter_files(successful, output_dir)

        self.logger.info("TSV summary written: %s", tsv_file)
        self.logger.info("JSON summary written: %s", json_file)

    @staticmethod
    def _write_summary_tsv(results: List[Dict[str, Any]], output_path: str) -> None:
        """Write the batch TSV, one row per sample including species audit columns."""
        headers = [
            "Filename", "Total Sequences", "Total Length", "Total Bases",
            "GC Content (%)", "AT Content (%)", "N50", "N75", "N90",
            "Median Length", "Mean Length", "Longest Sequence", "Shortest Sequence",
            "Ambiguous Bases (%)", "Sequences with Ns", "Max N-run",
            "Homopolymers", "Max Homopolymer", "Duplicate Sequences",
            "Short Sequences (<100 bp)", "Long Sequences (>1M bp)",
            "File Size (MB)", "Warnings", "K. pneumoniae Status",
            "Best Species", "ANI (%)", "K. pneumoniae Confirmed",
            "Contamination Suspected",
        ]
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write("\t".join(headers) + "\n")
            for result in results:
                kp_status = result.get("kp_status", {})
                status_str = (
                    f"GC:{kp_status.get('gc_status', 'Unknown')},"
                    f"Size:{kp_status.get('genome_size_status', 'Unknown')}"
                )
                sc = result.get("species_check", {}) or {}
                species_confirmed = result.get(
                    "species_status",
                    KleboscopeFASTAQC._species_status(sc),
                )
                contam = "Yes" if sc.get("contamination_suspected") else "No"
                row = [
                    result["filename"],
                    str(result["total_sequences"]),
                    str(result["total_length"]),
                    str(result["total_bases"]),
                    f"{result['gc_percent']:.2f}",
                    f"{result['at_percent']:.2f}",
                    str(result["n50"]),
                    str(result["n75"]),
                    str(result["n90"]),
                    str(result["median_length"]),
                    f"{result['mean_length']:.0f}",
                    str(result["longest_sequence"]),
                    str(result["shortest_sequence"]),
                    f"{result['ambiguous_percent']:.2f}",
                    str(result["sequences_with_n"]),
                    str(result["max_n_run"]),
                    str(result["homopolymers_count"]),
                    str(result["max_homopolymer"]),
                    str(result["duplicate_sequences"]),
                    str(result["short_sequences"]),
                    str(result["long_sequences"]),
                    f"{result['file_size_mb']:.2f}",
                    str(len(result.get("warnings", []))),
                    status_str,
                    sc.get("best_match", "ND"),
                    str(sc.get("ani_percent", "ND")),
                    species_confirmed,
                    contam,
                ]
                handle.write("\t".join(row) + "\n")

    def _build_json_summary(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Build the single consolidated batch JSON."""
        passed_names = [
            r["filename"] for r in results
            if r.get("species_status", self._species_status(r.get("species_check")))
            in ("Yes", "NA")
        ]
        excluded_names = [
            r["filename"] for r in results
            if r.get("species_status", self._species_status(r.get("species_check")))
            == "No"
        ]

        return {
            "metadata": {
                "tool": self.metadata["tool_name"],
                "version": self.metadata["version"],
                "biopython_version": self.metadata["biopython_version"],
                "analysis_date": self.metadata["analysis_date"],
                "total_files": len(results),
                "total_sequences": sum(r["total_sequences"] for r in results),
                "total_length": sum(r["total_length"] for r in results),
            },
            "statistics": {
                "gc_content_range": {
                    "min": min(r["gc_percent"] for r in results),
                    "max": max(r["gc_percent"] for r in results),
                    "mean": statistics.mean(r["gc_percent"] for r in results),
                    "median": statistics.median(r["gc_percent"] for r in results),
                },
                "n50_range": {
                    "min": min(r["n50"] for r in results),
                    "max": max(r["n50"] for r in results),
                    "mean": statistics.mean(r["n50"] for r in results),
                    "median": statistics.median(r["n50"] for r in results),
                },
                "files_with_warnings": sum(1 for r in results if r.get("warnings")),
                "total_warnings": sum(len(r.get("warnings", [])) for r in results),
                "files_with_high_gc": sum(
                    1 for r in results
                    if r["gc_percent"] > self.thresholds["gc_normal"][1]
                ),
                "files_with_low_gc": sum(
                    1 for r in results
                    if r["gc_percent"] < self.thresholds["gc_normal"][0]
                ),
                "files_confirmed_kpneumoniae": sum(
                    1 for r in results if r.get("species_status") == "Yes"
                ),
                "files_species_check_na": sum(
                    1 for r in results if r.get("species_status") == "NA"
                ),
                "files_failed_species_check": sum(
                    1 for r in results if r.get("species_status") == "No"
                ),
                "files_with_contamination_suspected": sum(
                    1 for r in results
                    if (r.get("species_check") or {}).get("contamination_suspected", False)
                ),
            },
            "species_filter": {
                "passed_files": passed_names,
                "excluded_files": excluded_names,
                "policy": (
                    "Yes -> forward, NA -> forward (check unavailable), "
                    "No -> exclude"
                ),
            },
            "files": results,
        }

    def _write_species_filter_files(
        self, results: List[Dict[str, Any]], output_dir: str
    ) -> None:
        """Write passed and excluded filename lists for downstream orchestration."""
        passed_file = os.path.join(output_dir, "kleboscope_species_passed.txt")
        excluded_file = os.path.join(output_dir, "kleboscope_species_excluded.txt")

        passed: List[str] = []
        excluded: List[str] = []
        for r in results:
            status = r.get(
                "species_status", self._species_status(r.get("species_check"))
            )
            if status == "No":
                excluded.append(r["filename"])
            else:
                passed.append(r["filename"])

        with open(passed_file, "w", encoding="utf-8") as handle:
            handle.write("\n".join(passed))
            if passed:
                handle.write("\n")

        with open(excluded_file, "w", encoding="utf-8") as handle:
            handle.write("\n".join(excluded))
            if excluded:
                handle.write("\n")

        self.logger.info(
            "Species filter: %d passed, %d excluded (%s)",
            len(passed), len(excluded), passed_file,
        )
        if excluded:
            self.logger.warning(
                "Samples excluded from downstream modules: %s",
                ", ".join(excluded),
            )

    def _write_summary_html(
        self,
        results: List[Dict[str, Any]],
        html_file: str,
        excluded_count: int,
    ) -> None:
        """Write the batch HTML summary for the given passing/unknown samples."""
        random_quote = self._get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        rows = []
        for result in results:
            warning_count = len(result.get("warnings", []))
            warning_class = (
                "none" if warning_count == 0
                else "low" if warning_count < 3
                else "high"
            )
            sc = result.get("species_check", {}) or {}
            species_status = result.get(
                "species_status", self._species_status(sc)
            )
            best_species = sc.get("best_match", "ND")
            ani = sc.get("ani_percent", "ND")
            confirmed = "Yes" if species_status == "Yes" else "NA"
            contam = "Yes" if sc.get("contamination_suspected") else "No"
            rows.append(
                f"<tr>"
                f"<td><strong>{result['filename']}</strong></td>"
                f"<td>{result['total_sequences']:,}</td>"
                f"<td>{result['total_length']:,}</td>"
                f"<td>{result['total_bases']:,}</td>"
                f"<td>{result['gc_percent']:.1f}</td>"
                f"<td>{result['at_percent']:.1f}</td>"
                f"<td>{result['n50']:,}</td>"
                f"<td>{result['n75']:,}</td>"
                f"<td>{result['n90']:,}</td>"
                f"<td>{result['median_length']:,}</td>"
                f"<td>{result['mean_length']:.0f}</td>"
                f"<td>{result['longest_sequence']:,}</td>"
                f"<td>{result['shortest_sequence']:,}</td>"
                f"<td>{result['ambiguous_percent']:.2f}</td>"
                f"<td>{result['sequences_with_n']:,}</td>"
                f"<td>{result['max_n_run']:,}</td>"
                f"<td>{result['homopolymers_count']:,}</td>"
                f"<td>{result['max_homopolymer']:,}</td>"
                f"<td>{result['duplicate_sequences']:,}</td>"
                f"<td>{result['short_sequences']:,}</td>"
                f"<td>{result['long_sequences']:,}</td>"
                f"<td>{result['file_size_mb']:.2f}</td>"
                f"<td>{warning_count}</td>"
                f"<td>{best_species}</td>"
                f"<td>{ani}%</td>"
                f"<td>{confirmed}</td>"
                f"<td>{contam}</td>"
                f"</tr>"
            )
        table_rows = "\n".join(rows)

        total_sequences = sum(r["total_sequences"] for r in results)
        total_length = sum(r["total_length"] for r in results)
        total_warnings = sum(len(r.get("warnings", [])) for r in results)

        excluded_note = ""
        if excluded_count > 0:
            excluded_note = (
                '<div class="report-section" style="border-left: 4px solid #94a3b8; '
                'background: #f1f5f9;">'
                '<p style="color: #475569; font-size: 13px; margin: 0;">'
                f'<strong>Note:</strong> {excluded_count} additional sample(s) '
                'failed the fastANI species check and are not shown here. See '
                '<code>kleboscope_species_excluded.txt</code> for the list.</p></div>'
            )

        footer = _render_footer("FASTA QC Batch Summary", current_time)

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KLEBOSCOPE - FASTA QC Batch Summary</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%);
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    color: #ffffff; padding: 20px; min-height: 100vh;
}}
.container {{ max-width: 1800px; margin: 0 auto; }}
.ascii-container {{
    background: rgba(0, 0, 0, 0.7); padding: 20px; border-radius: 15px;
    margin-bottom: 20px; box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
    border: 2px solid rgba(0, 119, 190, 0.3); overflow-x: auto; text-align: center;
}}
.ascii-art {{
    font-family: 'Cascadia Mono', 'Consolas', 'Monaco', 'Menlo',
                 'DejaVu Sans Mono', 'Liberation Mono', 'Courier New', monospace;
    font-size: 12px; line-height: 1.05; white-space: pre;
    color: #00d4ff; text-shadow: 0 0 8px rgba(0, 212, 255, 0.4);
    display: inline-block; padding: 8px 0; text-align: left;
    font-variant-ligatures: none;
    font-feature-settings: "liga" 0, "calt" 0, "kern" 0;
    font-kerning: none; letter-spacing: 0;
    text-rendering: geometricPrecision;
}}
.quote-container {{
    background: rgba(255, 255, 255, 0.1); backdrop-filter: blur(10px);
    padding: 20px; border-radius: 10px; margin-bottom: 30px; text-align: center;
    min-height: 100px; display: flex; flex-direction: column;
    justify-content: center; box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
    border: 1px solid rgba(255, 255, 255, 0.2); transition: opacity 0.5s ease-in-out;
}}
.quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
.quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
.report-section {{
    background: rgba(255, 255, 255, 0.95); color: #1f2937;
    padding: 25px; border-radius: 10px; margin-bottom: 20px;
    box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
}}
.report-section h2 {{
    color: #0b3d5f; border-bottom: 3px solid #0077be;
    padding-bottom: 10px; margin-bottom: 20px; font-size: 24px;
}}
.summary-table {{
    width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 12px;
}}
.summary-table th {{
    background: #005a8c; color: white; padding: 12px; text-align: left;
    font-weight: bold; position: sticky; top: 0;
}}
.summary-table td {{ padding: 10px; border-bottom: 1px solid #e5e7eb; }}
.summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
.summary-table tr:hover {{ background-color: #e0f2fe; }}
.stats-grid {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 15px; margin-bottom: 20px;
}}
.stat-card {{
    background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
    color: white; padding: 15px; border-radius: 8px; text-align: center;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
}}
.stat-value {{ font-size: 24px; font-weight: bold; margin-bottom: 5px; }}
.stat-label {{ font-size: 12px; opacity: 0.9; }}
.footer {{
    text-align: center; margin-top: 30px; padding: 20px;
    background: rgba(0, 0, 0, 0.3); border-radius: 10px; font-size: 14px;
}}
.timestamp {{ color: #fbbf24; font-weight: bold; }}
.authorship {{
    margin-top: 15px; padding: 15px; background: rgba(255, 255, 255, 0.1);
    border-radius: 8px; font-size: 12px;
}}
.table-container {{
    overflow-x: auto; margin: 20px 0; max-height: 600px; overflow-y: auto;
}}
@media (max-width: 768px) {{
    .ascii-art {{ font-size: 6px; }}
    .summary-table {{ font-size: 11px; }}
}}
</style>
</head>
<body>
<div class="container">
<div class="ascii-container">
<div class="ascii-art">{KLEBOSCOPE_BANNER}</div>
</div>

<div class="quote-container" id="quoteContainer">
    <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
    <div class="quote-author" id="quoteAuthor">- {random_quote['author']}</div>
</div>

<div class="report-section">
    <h2>FASTA QC Batch Summary</h2>
    <div class="stats-grid">
        <div class="stat-card"><div class="stat-value">{len(results)}</div><div class="stat-label">SAMPLES SHOWN</div></div>
        <div class="stat-card"><div class="stat-value">{total_sequences:,}</div><div class="stat-label">TOTAL SEQUENCES</div></div>
        <div class="stat-card"><div class="stat-value">{total_length:,}</div><div class="stat-label">TOTAL BASES</div></div>
        <div class="stat-card"><div class="stat-value">{total_warnings}</div><div class="stat-label">TOTAL WARNINGS</div></div>
    </div>
    <p><strong>Analysis Date:</strong> {self.metadata['analysis_date']}</p>
    <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
    <p><strong>BioPython Version:</strong> {self.metadata['biopython_version']}</p>
</div>

{excluded_note}

<div class="report-section">
    <h2>Detailed FASTA QC Results</h2>
    <p style="color: #666; margin-bottom: 15px;">
        Species columns show the fastANI best match and confirmation status.
        Samples that explicitly failed the species check are not listed here.
    </p>
    <div class="table-container">
        <table class="summary-table">
            <thead>
                <tr>
                    <th>Filename</th><th>Sequences</th><th>Total Length</th><th>Total Bases</th>
                    <th>GC %</th><th>AT %</th><th>N50</th><th>N75</th><th>N90</th>
                    <th>Median</th><th>Mean</th><th>Longest</th><th>Shortest</th>
                    <th>Ambiguous %</th><th>Seqs with Ns</th><th>Max N-run</th>
                    <th>Homopolymers</th><th>Max Homopolymer</th><th>Duplicates</th>
                    <th>Short (&lt;100 bp)</th><th>Long (&gt;1M bp)</th>
                    <th>File Size (MB)</th><th>Warnings</th>
                    <th>Best Species</th><th>ANI %</th>
                    <th>K. pneumoniae</th><th>Contamination</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>
    </div>
</div>
{footer}
</div>

<script>
const quotes = {json.dumps(self.science_quotes)};
const quoteContainer = document.getElementById('quoteContainer');
const quoteText = document.getElementById('quoteText');
const quoteAuthor = document.getElementById('quoteAuthor');
function displayQuote() {{
    quoteContainer.style.opacity = '0';
    setTimeout(() => {{
        const q = quotes[Math.floor(Math.random() * quotes.length)];
        quoteText.textContent = '"' + q.text + '"';
        quoteAuthor.textContent = '- ' + q.author;
        quoteContainer.style.opacity = '1';
    }}, 500);
}}
displayQuote();
setInterval(displayQuote, 10000);
</script>
</body>
</html>"""

        with open(html_file, "w", encoding="utf-8") as handle:
            handle.write(html_content)

    # ------------------------------------------------------------------
    # Batch driver
    # ------------------------------------------------------------------
    def process_files(self, pattern: str, output_dir: str = "fasta_qc_results") -> List[Dict[str, Any]]:
        """Run FASTA QC on every matching file and write all outputs."""
        print("\n" + "=" * 80)
        print("KLEBOSCOPE FASTA QC - K. pneumoniae Quality Control")
        print("=" * 80)
        print(f"Author: Brown Beckley | Email: brownbeckley94@gmail.com")
        print(f"Affiliation: University of Ghana Medical School - Department of Medical Biochemistry")
        print("=" * 80)
        print(f"Output directory: {output_dir}")
        print(f"Using {self.cpus} CPU core(s)")
        if self.species_refs:
            print(f"Species references: {', '.join(self.species_refs.keys())}")
        else:
            print("Species check disabled (no references found in ref_db/)")
        print("=" * 80)

        fasta_exts = [
            ".fna", ".fasta", ".fa", ".faa", ".fn",
            ".fna.gz", ".fasta.gz", ".fa.gz",
        ]
        files: List[str] = []
        files.extend(glob.glob(pattern))
        if not files:
            for ext in fasta_exts:
                files.extend(glob.glob(f"{pattern}{ext}"))
        files = list(set(f for f in files if os.path.exists(f)))

        if not files:
            self.logger.error("No FASTA files found matching pattern: %s", pattern)
            self.logger.info("Supported extensions: %s", ", ".join(fasta_exts))
            return []

        self.logger.info("Found %d FASTA file(s)", len(files))
        os.makedirs(output_dir, exist_ok=True)

        all_results: List[Dict[str, Any]] = []
        successful = 0

        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            futures = {executor.submit(self.analyze_file, f): f for f in files}
            for future in as_completed(futures):
                file = futures[future]
                try:
                    result = future.result()
                    all_results.append(result)
                    if result.get("status") == "success":
                        successful += 1
                        if result.get("species_status") != "No":
                            self.create_individual_html_report(result, output_dir)
                        else:
                            sc = result.get("species_check") or {}
                            self.logger.warning(
                                "%s: excluded from HTML (best match: %s, ANI %s%%)",
                                result["filename"],
                                sc.get("best_match", "ND"),
                                sc.get("ani_percent", "ND"),
                            )
                    else:
                        self.logger.error(
                            "%s: %s",
                            result["filename"],
                            result.get("error", "unknown error"),
                        )
                except Exception as exc:
                    self.logger.error("%s: %s", os.path.basename(file), exc)
                    all_results.append({
                        "filename": os.path.basename(file),
                        "status": "error",
                        "error": str(exc),
                    })

        if successful > 0:
            self.create_summary_report(all_results, output_dir)

        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80)
        print(f"Total files: {len(files)}")
        print(f"Successful analyses: {successful}")
        print(f"Failed analyses: {len(files) - successful}")

        successful_results = [r for r in all_results if r.get("status") == "success"]
        excluded = [
            r["filename"] for r in successful_results
            if r.get("species_status") == "No"
        ]
        kept = len(successful_results) - len(excluded)
        print(f"Species check passed/unknown: {kept}")
        print(f"Species check failed (excluded): {len(excluded)}")
        if excluded:
            print(f"Excluded samples: {', '.join(excluded)}")

        print("\nOutput files:")
        print(f"  Per-sample HTML: {output_dir}/*/*_fasta_qc_report.html")
        if successful > 1:
            print(f"  Summary TSV: {output_dir}/kleboscope_fasta_qc_summary.tsv")
            print(f"  Summary JSON: {output_dir}/kleboscope_fasta_qc_summary.json")
            if kept > 0:
                print(f"  Summary HTML: {output_dir}/kleboscope_fasta_qc_summary.html")
        print(f"  Species passed: {output_dir}/kleboscope_species_passed.txt")
        print(f"  Species excluded: {output_dir}/kleboscope_species_excluded.txt")
        print("=" * 80)
        return all_results


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description=(
            "Kleboscope FASTA QC - K. pneumoniae quality control with "
            "HTML reports and fastANI species confirmation"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python klebo_fasta_qc.py "*.fna"
  python klebo_fasta_qc.py "genomes/*.fasta" --output qc_results
  python klebo_fasta_qc.py "*.fna" --ref_dir /path/to/ref_db

Supported FASTA extensions: .fna, .fasta, .fa, .faa, .fn

Per-sample output (1 file each):
  <sample>/<sample>_fasta_qc_report.html  (passed/unknown only)

Batch-level outputs:
  kleboscope_fasta_qc_summary.tsv         (full audit)
  kleboscope_fasta_qc_summary.html        (passed/unknown only + exclusion count)
  kleboscope_fasta_qc_summary.json        (full audit)
  kleboscope_species_passed.txt           (forward to downstream modules)
  kleboscope_species_excluded.txt         (drop)
        """,
    )
    parser.add_argument("pattern", help='File pattern (e.g., "*.fna", "genomes/*.fasta")')
    parser.add_argument("--output", "-o", default="fasta_qc_results", help="Output directory")
    parser.add_argument("--cpus", "-c", type=int, default=None, help="CPU core budget")
    parser.add_argument("--ref_dir", default=None, help="Reference directory (default: ./ref_db)")
    args = parser.parse_args()

    try:
        qc = KleboscopeFASTAQC(cpus=args.cpus, ref_dir=args.ref_dir)
        qc.process_files(args.pattern, args.output)
    except Exception as exc:
        print(f"FASTA QC analysis failed: {exc}")
        sys.exit(1)


if __name__ == "__main__":
    main()