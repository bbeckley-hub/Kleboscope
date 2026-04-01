#!/usr/bin/env python3
"""
Kleboscope FASTA QC - Comprehensive Quality Control for K. pneumoniae
Author: Brown Beckley
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Version: 2.0.0
"""

import os
import sys
import glob
import json
import statistics
import re
import logging
import argparse
import random
from pathlib import Path
from datetime import datetime
from collections import Counter
from typing import List, Dict, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

# BioPython imports
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

class KleboscopeFASTAQC:
    """Comprehensive FASTA QC for K. pneumoniae"""
    
    def __init__(self, cpus: int = 1):
        self.logger = self._setup_logging()
        self.cpus = cpus
        self.logger.info(f"Initialized with {self.cpus}") # value passed by kleboscope.py
        
        self.metadata = {
            "tool_name": "Kleboscope FASTA QC",
            "version": "2.0.0", 
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "biopython_version": self._get_biopython_version()
        }
        
        # K. pneumoniae specific thresholds
        self.thresholds = {
            'gc_normal': (55, 58),
            'short_seq': 100,
            'long_seq': 1000000,
            'max_n_run': 100,
            'max_homopolymer': 20,
            'ambiguous_critical': 5.0,
            'ambiguous_warning': 1.0,
        }
        
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "Kleboscope turns genomic complexity into actionable insights.", "author": "Brown Beckley"}
        ]
        
        self.ascii_art = """
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
"""

    def _setup_logging(self):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        return logging.getLogger(__name__)
    
    def _get_biopython_version(self) -> str:
        try:
            import Bio
            return Bio.__version__
        except ImportError:
            return "Unknown"

    # =========================================================================
    # CORE LOGIC
    # =========================================================================
    def analyze_file(self, fasta_file: str) -> Dict[str, Any]:
        """Comprehensive analysis of a FASTA file"""
        try:
            filename = Path(fasta_file).name
            self.logger.info(f"🔬 Analyzing {filename}...")
            
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            if not sequences:
                return {'filename': filename, 'status': 'error', 'error': 'No sequences found'}
            
            seq_lengths = sorted([len(seq) for seq in sequences], reverse=True)
            total_length = sum(seq_lengths)
            
            # N and L stats
            n50, l50 = self._calc_nx_lx(seq_lengths, total_length, 50)
            n75, l75 = self._calc_nx_lx(seq_lengths, total_length, 75)
            n90, l90 = self._calc_nx_lx(seq_lengths, total_length, 90)
            
            # Base Composition & Ambiguous Runs
            base_counts = Counter()
            gc_contents = []
            max_n_run = 0
            homopolymers = []
            
            for seq in sequences:
                seq_str = str(seq.seq).upper()
                base_counts.update(seq_str)
                gc_contents.append(gc_fraction(seq_str) * 100)
                
                # N-runs
                runs = [len(m.group()) for m in re.finditer(r'N+', seq_str)]
                if runs: max_n_run = max(max_n_run, max(runs))
                
                # Homopolymers
                for base in ['A', 'T', 'G', 'C']:
                    homopolymers.extend([len(m.group()) for m in re.finditer(f'{base}+', seq_str) if len(m.group()) > 4])
            
            total_bases = sum(base_counts.values())
            pct = lambda k: (base_counts.get(k, 0) / total_bases * 100) if total_bases else 0
            
            gc_percent = pct('G') + pct('C')
            ambig_bases = sum(base_counts.get(b, 0) for b in ['N', 'Y', 'R', 'W', 'S', 'K', 'M', 'B', 'D', 'H', 'V'])
            ambig_percent = (ambig_bases / total_bases * 100) if total_bases else 0
            
            kp_status = self._check_kpneumoniae_specific(gc_percent, total_length, len(seq_lengths))
            
            results = {
                'filename': filename,
                'file_size_mb': os.path.getsize(fasta_file) / (1024 ** 2),
                'analysis_date': self.metadata['analysis_date'],
                'status': 'success',
                'total_sequences': len(sequences),
                'total_length': total_length,
                'total_bases': total_bases,
                'longest_sequence': seq_lengths[0],
                'shortest_sequence': seq_lengths[-1],
                'mean_length': statistics.mean(seq_lengths),
                'median_length': statistics.median(seq_lengths),
                'n50': n50, 'l50': l50, 'n75': n75, 'l75': l75, 'n90': n90, 'l90': l90,
                'gc_percent': gc_percent, 'at_percent': pct('A') + pct('T'),
                'a_percent': pct('A'), 't_percent': pct('T'), 'g_percent': pct('G'), 'c_percent': pct('C'),
                'ambiguous_percent': ambig_percent,
                'sequences_with_n': sum(1 for seq in sequences if 'N' in str(seq.seq).upper()),
                'total_n_bases': base_counts.get('N', 0),
                'max_n_run': max_n_run,
                'homopolymers_count': len(homopolymers),
                'max_homopolymer': max(homopolymers) if homopolymers else 0,
                'duplicate_sequences': len(sequences) - len(set(hash(str(s.seq)) for s in sequences)),
                'short_sequences': sum(1 for l in seq_lengths if l < self.thresholds['short_seq']),
                'long_sequences': sum(1 for l in seq_lengths if l > self.thresholds['long_seq']),
                'length_distribution': self._create_length_bins(seq_lengths),
                'kp_status': kp_status
            }
            
            results['warnings'] = self._generate_warnings(results)
            self.logger.info(f"✅ {filename}: {len(sequences)} seqs, N50: {n50:,}, GC: {gc_percent:.1f}%")
            return results
            
        except Exception as e:
            self.logger.error(f"❌ Error analyzing {fasta_file}: {e}")
            return {'filename': Path(fasta_file).name, 'status': 'error', 'error': str(e)}

    def _calc_nx_lx(self, sorted_lengths: List[int], total_length: int, x: int) -> tuple:
        target = total_length * (x / 100)
        cumulative = 0
        for i, length in enumerate(sorted_lengths, 1):
            cumulative += length
            if cumulative >= target:
                return length, i
        return 0, 0
        
    def _check_kpneumoniae_specific(self, gc_percent: float, total_length: int, contigs: int) -> Dict:
        size_ok = 5000000 <= total_length <= 5600000
        gc_ok = self.thresholds['gc_normal'][0] <= gc_percent <= self.thresholds['gc_normal'][1]
        return {
            'genome_size_mbp': total_length / 1000000,
            'genome_size_status': 'Normal' if size_ok else 'Atypical',
            'gc_status': 'Normal' if gc_ok else 'Atypical',
            'assembly_quality': "Good" if contigs <= 50 else "Moderate" if contigs <= 100 else "Fragmented",
            'contig_count': contigs,
            'expected_gc_range': f"{self.thresholds['gc_normal'][0]}-{self.thresholds['gc_normal'][1]}%",
            'expected_genome_size': "5.0-5.6 Mbp"
        }

    def _create_length_bins(self, lengths: List[int]) -> Dict[str, int]:
        bins = {'< 100 bp': 0, '100-500 bp': 0, '500-1k bp': 0, '1k-5k bp': 0, '5k-10k bp': 0, '10k-50k bp': 0, '50k-100k bp': 0, '100k-500k bp': 0, '500k-1M bp': 0, '> 1M bp': 0}
        for l in lengths:
            if l < 100: bins['< 100 bp'] += 1
            elif l < 500: bins['100-500 bp'] += 1
            elif l < 1000: bins['500-1k bp'] += 1
            elif l < 5000: bins['1k-5k bp'] += 1
            elif l < 10000: bins['5k-10k bp'] += 1
            elif l < 50000: bins['10k-50k bp'] += 1
            elif l < 100000: bins['50k-100k bp'] += 1
            elif l < 500000: bins['100k-500k bp'] += 1
            elif l < 1000000: bins['500k-1M bp'] += 1
            else: bins['> 1M bp'] += 1
        return bins

    def _generate_warnings(self, r: Dict) -> List[Dict]:
        w = []
        low, high = self.thresholds['gc_normal']
        if r['gc_percent'] < low or r['gc_percent'] > high: w.append({'level': 'warning', 'message': f"GC content ({r['gc_percent']:.1f}%) outside Kp range ({low}-{high}%)"})
        if r['ambiguous_percent'] > self.thresholds['ambiguous_critical']: w.append({'level': 'danger', 'message': f"High ambiguous bases ({r['ambiguous_percent']:.2f}%)"})
        if r['max_n_run'] > 100: w.append({'level': 'danger', 'message': f"Long N-run ({r['max_n_run']} bases)"})
        if r['max_homopolymer'] > 20: w.append({'level': 'danger', 'message': f"Long homopolymer ({r['max_homopolymer']} bases)"})
        if r['short_sequences'] > r['total_sequences'] * 0.5: w.append({'level': 'danger', 'message': f"Many short sequences ({r['short_sequences']})"})
        if r['kp_status']['genome_size_status'] == 'Atypical': w.append({'level': 'warning', 'message': f"Atypical size ({r['kp_status']['genome_size_mbp']:.2f} Mbp)"})
        if r['kp_status']['assembly_quality'] == 'Fragmented': w.append({'level': 'warning', 'message': f"Fragmented assembly ({r['kp_status']['contig_count']} contigs)"})
        return w

    # =========================================================================
    # HTML TEMPLATING
    # =========================================================================
    def _get_base_html(self, title: str, body_content: str) -> str:
        quote = random.choice(self.science_quotes)
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ background: linear-gradient(135deg, #0b3d5f 0%, #1b6b8f 50%, #2a9d8f 100%); font-family: 'Segoe UI', Tahoma, Geneva, sans-serif; color: #fff; padding: 20px; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{ background: rgba(0,0,0,0.7); padding: 20px; border-radius: 15px; margin-bottom: 20px; border: 2px solid rgba(0,119,190,0.3); }}
        .ascii-art {{ font-family: monospace; font-size: 10px; color: #0077be; white-space: pre; overflow-x: auto; }}
        .quote-container {{ background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px; text-align: center; border: 1px solid rgba(255,255,255,0.2); margin-bottom: 20px; }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{ background: rgba(255,255,255,0.95); color: #1f2937; padding: 25px; border-radius: 10px; margin-bottom: 20px; }}
        .report-section h2 {{ color: #1e3a8a; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; margin-bottom: 20px; }}
        .metrics-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 20px; }}
        .metric-card {{ background: linear-gradient(135deg, #0077be 0%, #005a8c 100%); color: white; padding: 20px; border-radius: 8px; text-align: center; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .summary-table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px; background: white; }}
        .summary-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; }}
        .summary-table td {{ padding: 12px; border-bottom: 1px solid #e5e7eb; }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .footer {{ text-align: center; padding: 20px; background: rgba(0,0,0,0.8); border-radius: 10px; margin-top: 30px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container"><div class="ascii-art">{self.ascii_art}</div></div>
            <div class="quote-container"><div class="quote-text">"{quote['text']}"</div><div class="quote-author">— {quote['author']}</div></div>
        </div>
        {body_content}
        <div class="footer"><p>KLEBOSCOPE FASTA QC | Generated: {self.metadata['analysis_date']}</p></div>
    </div>
</body>
</html>"""

    def create_individual_html_report(self, r: Dict, output_dir: str):
        sample_dir = Path(output_dir) / Path(r['filename']).stem
        sample_dir.mkdir(parents=True, exist_ok=True)
        
        warnings_html = ''.join([f'<div style="background: {"#fee2e2" if w["level"]=="danger" else "#fef3c7"}; padding: 10px; margin: 5px 0; border-left: 4px solid {"#dc2626" if w["level"]=="danger" else "#f59e0b"};"><strong>{w["level"].upper()}</strong>: {w["message"]}</div>' for w in r.get('warnings', [])])
        
        content = f"""
        <div class="report-section">
            <h2>📊 Sample Information</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{r['filename']}</div><div>Filename</div></div>
                <div class="metric-card"><div class="metric-value">{r['file_size_mb']:.2f} MB</div><div>File Size</div></div>
            </div>
        </div>
        <div class="report-section">
            <h2>🎯 FASTA QC Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{r['total_sequences']:,}</div><div>Sequences</div></div>
                <div class="metric-card"><div class="metric-value">{r['total_length']:,}</div><div>Total Length</div></div>
                <div class="metric-card"><div class="metric-value">{r['n50']:,}</div><div>N50</div></div>
                <div class="metric-card"><div class="metric-value">{r['gc_percent']:.1f}%</div><div>GC Content</div></div>
            </div>
        </div>
        <div class="report-section">
            <h2>⚠️ Quality Warnings</h2>
            {warnings_html if warnings_html else '<p style="color: #10b981; font-weight: bold;">✅ No quality warnings detected</p>'}
        </div>
        """
        
        html_file = sample_dir / f"{Path(r['filename']).stem}_fasta_qc_report.html"
        with open(html_file, 'w', encoding='utf-8') as f: f.write(self._get_base_html(f"QC Report - {r['filename']}", content))
        
        with open(sample_dir / f"{Path(r['filename']).stem}_fasta_qc_report.json", 'w') as f: json.dump(r, f, indent=2, default=str)

    def create_summary_report(self, all_results: List[Dict], output_dir: str):
        successful = [r for r in all_results if r.get('status') == 'success']
        if not successful: return
        
        # TSV
        with open(Path(output_dir) / "klebo_fasta_qc_summary.tsv", 'w') as f:
            headers = ['Filename', 'Sequences', 'Length', 'GC%', 'N50', 'Warnings', 'Assembly']
            f.write('\t'.join(headers) + '\n')
            for r in successful: f.write(f"{r['filename']}\t{r['total_sequences']}\t{r['total_length']}\t{r['gc_percent']:.2f}\t{r['n50']}\t{len(r.get('warnings',[]))}\t{r['kp_status']['assembly_quality']}\n")
            
        # JSON
        with open(Path(output_dir) / "klebo_fasta_qc_summary.json", 'w') as f: json.dump({'files': successful}, f, indent=2, default=str)
        
        # HTML
        rows = ''.join([f'<tr><td>{r["filename"]}</td><td>{r["total_sequences"]}</td><td>{r["total_length"]:,}</td><td>{r["gc_percent"]:.1f}%</td><td>{r["n50"]:,}</td><td>{len(r.get("warnings",[]))}</td><td>{r["kp_status"]["assembly_quality"]}</td></tr>' for r in successful])
        content = f"""
        <div class="report-section">
            <h2>📊 Batch QC Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{len(successful)}</div><div>Files Processed</div></div>
                <div class="metric-card"><div class="metric-value">{sum(r['total_sequences'] for r in successful):,}</div><div>Total Sequences</div></div>
            </div>
            <div style="overflow-x:auto;">
                <table class="summary-table">
                    <thead><tr><th>Filename</th><th>Sequences</th><th>Length</th><th>GC%</th><th>N50</th><th>Warnings</th><th>Assembly Quality</th></tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>
        </div>
        """
        with open(Path(output_dir) / "klebo_fasta_qc_summary.html", 'w', encoding='utf-8') as f: f.write(self._get_base_html("Batch QC Summary", content))

    # =========================================================================
    # BATCH PROCESSOR
    # =========================================================================
    def process_files(self, pattern: str, output_dir: str) -> List[Dict[str, Any]]:
        self.logger.info(f"Using {self.cpus} CPU cores")
        
        files = list(set(glob.glob(pattern) + glob.glob(f"{pattern}.fasta") + glob.glob(f"{pattern}.fna")))
        if not files:
            self.logger.error(f"❌ No FASTA files found matching pattern: {pattern}")
            return []
            
        os.makedirs(output_dir, exist_ok=True)
        all_results = []
        
        with ThreadPoolExecutor(max_workers=self.cpus) as executor:
            futures = {executor.submit(self.analyze_file, f): f for f in files}
            for future in as_completed(futures):
                res = future.result()
                all_results.append(res)
                if res['status'] == 'success': self.create_individual_html_report(res, output_dir)
                
        self.create_summary_report(all_results, output_dir)
        return all_results

def main():
    parser = argparse.ArgumentParser(description='Kleboscope FASTA QC - Orchestrator Compliant')
    parser.add_argument('-i', '--input', required=True, help='File pattern for FASTA files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of CPU cores')
    args = parser.parse_args()
    
    try:
        qc = KleboscopeFASTAQC(cpus=args.threads)
        qc.process_files(args.input, args.output)
    except Exception as e:
        print(f"❌ FASTA QC analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()