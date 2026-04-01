#!/usr/bin/env python3
"""
Kleboscope MLST Module - Orchestrator Compliant
Author: Brown Beckley
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Version: 2.0.0 
"""

#import os
#import sys
import json
import glob
import argparse
import subprocess
import random
from pathlib import Path
from typing import Dict, List
from datetime import datetime

class KleboscopeMLSTAnalyzer:
    def __init__(self, database_dir: Path, script_dir: Path):
        self.database_dir = database_dir
        self.script_dir = script_dir
        self.mlst_bin = script_dir / "bin" / "mlst"
        
        self._load_lineage_db()
        
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"text": "Kleboscope turns genomic complexity into actionable insights for AMR surveillance.", "author": "Brown Beckley"}
        ]
        
        self.ascii_art = """
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
"""

    def _load_lineage_db(self):
        db_path = Path(__file__).parent / "klebo_mlst_lineages.json"
        try:
            with open(db_path, 'r') as f:
                self.lineage_db = json.load(f)
        except FileNotFoundError:
            print(f"⚠️ Warning: {db_path.name} not found. Lineage data will be unavailable.")
            self.lineage_db = {}

    def find_fasta_files(self, input_path: str) -> List[Path]:
        fasta_patterns = [input_path, f"{input_path}.fasta", f"{input_path}.fa", f"{input_path}.fna"]
        genome_files = []
        for p in fasta_patterns:
            genome_files.extend(glob.glob(p))
        return [Path(f) for f in sorted(list(set(genome_files))) if Path(f).is_file()]

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "klebsiella") -> Dict:
        print(f"🔬 Processing: {input_file.name}")
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        raw_output_file = sample_output_dir / "mlst_raw_output.txt"
        
        if not self.mlst_bin.exists():
            print(f"❌ MLST binary not found at: {self.mlst_bin}")
            return self._handle_error(input_file.name, sample_output_dir)
            
        mlst_cmd = ["perl", str(self.mlst_bin), str(input_file), "--scheme", scheme, "--csv", "--nopath"]
        
        try:
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)
            with open(raw_output_file, 'w') as f:
                f.write(f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
                
            mlst_results = self.parse_mlst_csv(result.stdout, input_file.name)
            st = mlst_results.get('st', 'ND')
            mlst_results.update(self.get_lineage_info(st))
            mlst_results.update(self.get_identity_coverage(st))
            
            self.generate_output_files(mlst_results, sample_output_dir)
            print(f"✅ Completed: {input_file.name} -> ST{st}")
            return mlst_results
        except subprocess.CalledProcessError:
            print(f"❌ MLST failed for {input_file.name}")
            return self._handle_error(input_file.name, sample_output_dir)

    def _handle_error(self, sample_name: str, output_dir: Path) -> Dict:
        res = self.get_fallback_results(sample_name)
        self.generate_output_files(res, output_dir)
        return res

    def parse_mlst_csv(self, stdout: str, sample_name: str) -> Dict:
        lines = stdout.strip().split('\n')
        if not lines: return self.get_empty_results(sample_name)
        
        result_line = next((line.strip() for line in reversed(lines) if line.strip() and ',' in line and not line.startswith('[')), None)
        if not result_line: return self.get_empty_results(sample_name)
        
        parts = result_line.split(',')
        if len(parts) < 3: return self.get_empty_results(sample_name)
        
        st = parts[2]
        alleles = {}
        allele_parts = []
        for i in range(3, len(parts)):
            if '(' in parts[i] and ')' in parts[i]:
                gene, allele = parts[i].split('(')
                allele = allele.rstrip(')')
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
                
        return {
            "sample": sample_name, "st": st, "scheme": "klebsiella", "alleles": alleles,
            "allele_profile": '-'.join(allele_parts),
            "confidence": "HIGH" if st and st not in ['-', 'ND'] else "LOW",
            "mlst_assigned": bool(st and st not in ['-', 'ND'])
        }

    def get_lineage_info(self, st: str) -> Dict:
        if st in self.lineage_db:
            return self.lineage_db[st]
        if st.isdigit():
            return {"clonal_complex": f"Unknown (ST{st})", "classification": "Not in database", "geographic_distribution": "Unknown", "clinical_significance": f"ST{st} is not currently in the Kleboscope database", "common_virulence": ["Unknown"], "outbreak_potential": "UNKNOWN", "typical_capsule": ["Unknown"], "resistance_profile": ["Unknown"]}
        return {"clonal_complex": "Not Assigned", "classification": "MLST typing failed", "geographic_distribution": "N/A", "clinical_significance": "Could not determine sequence type", "common_virulence": ["Cannot determine"], "outbreak_potential": "UNKNOWN", "typical_capsule": ["Unknown"], "resistance_profile": ["Unknown"]}

    def get_identity_coverage(self, st: str) -> Dict:
        if st and st not in ['-', 'ND', 'UNKNOWN']:
            return {"identity": "100%", "coverage": "100%", "mlst_status": "Assigned"}
        return {"identity": "Not Assigned", "coverage": "Not Assigned", "mlst_status": "Not Assigned"}

    def get_empty_results(self, sample_name: str) -> Dict:
        return {"sample": sample_name, "st": "ND", "scheme": "klebsiella", "alleles": {}, "allele_profile": "", "confidence": "LOW", "mlst_assigned": False}

    def get_fallback_results(self, sample_name: str) -> Dict:
        return {"sample": sample_name, "st": "UNKNOWN", "scheme": "klebsiella", "alleles": {}, "allele_profile": "", "confidence": "LOW", "mlst_assigned": False, "error": "MLST analysis failed"}

    # =========================================================================
    # HTML TEMPLATING & REPORT GENERATION
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
        .allele-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); gap: 15px; }}
        .allele-card {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 15px; border-radius: 8px; text-align: center; }}
        .virulence-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 10px; }}
        .virulence-card {{ background: #e0f2fe; color: #0369a1; padding: 12px; border-radius: 6px; text-align: center; font-weight: bold; border-left: 4px solid #0ea5e9; }}
        .profile-box {{ background: #f8fafc; padding: 15px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #3b82f6; }}
        .summary-table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px; background: white; }}
        .summary-table th {{ background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%); color: white; padding: 12px; text-align: left; }}
        .summary-table td {{ padding: 12px; border-bottom: 1px solid #e5e7eb; }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .st-cell {{ font-weight: bold; color: #1e40af; }}
        .allele-cell {{ font-family: monospace; background-color: #f0f9ff; color: #0369a1; font-weight: bold; }}
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
        <div class="footer"><p>KLEBOSCOPE MLST Analyzer | Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p></div>
    </div>
</body>
</html>"""

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        if 'identity' not in mlst_results: mlst_results.update(self.get_identity_coverage(mlst_results.get('st', 'ND')))
        
        # 1. HTML Report
        alleles_html = ''.join([f'<div class="allele-card"><div style="font-size: 12px; opacity: 0.9;">{g}</div><div style="font-size: 18px;">{a}</div></div>' for g, a in mlst_results.get('alleles', {}).items()])
        virulence_html = ''.join([f'<div class="virulence-card">{v}</div>' for v in mlst_results.get('common_virulence', [])[:6]])
        
        content = f"""
        <div class="report-section">
            <h2>📊 Sample Information</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{mlst_results['sample']}</div><div>Sample Name</div></div>
                <div class="metric-card"><div class="metric-value">Klebsiella</div><div>Scheme</div></div>
            </div>
        </div>
        <div class="report-section">
            <h2>🎯 MLST Results</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">ST{mlst_results['st']}</div><div>Sequence Type</div></div>
                <div class="metric-card"><div class="metric-value">{mlst_results.get('identity', 'N/A')}</div><div>Identity</div></div>
            </div>
            <h3>Allele Profile</h3><div class="profile-box"><code style="font-size: 16px; color: #1e40af; font-weight: bold;">{mlst_results['allele_profile']}</code></div>
            <h3>Individual Alleles</h3><div class="allele-grid">{alleles_html}</div>
        </div>
        <div class="report-section">
            <h2>🌍 Lineage Information</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{mlst_results.get('clonal_complex', 'Unknown')}</div><div>Clonal Complex</div></div>
                <div class="metric-card"><div class="metric-value">{mlst_results.get('outbreak_potential', 'UNKNOWN')}</div><div>Outbreak Potential</div></div>
            </div>
            <h3>Clinical Significance</h3><div class="profile-box"><p>{mlst_results.get('clinical_significance', 'Further analysis required.')}</p></div>
            <div style="display: flex; gap: 20px; flex-wrap: wrap;">
                <div style="flex: 1; min-width: 300px;"><h3>Characteristics</h3><p><strong>Capsule:</strong> {', '.join(mlst_results.get('typical_capsule', ['Unknown']))}</p><p><strong>Resistance:</strong> {', '.join(mlst_results.get('resistance_profile', ['Unknown']))}</p></div>
                <div style="flex: 1; min-width: 300px;"><h3>Common Virulence</h3><div class="virulence-grid">{virulence_html}</div></div>
            </div>
        </div>
        """
        with open(output_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(self._get_base_html(f"MLST Report - {mlst_results['sample']}", content))
            
        # 2. TSV Report
        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write("Sample\tST\tScheme\tMLST_Status\tIdentity\tCoverage\tClonal_Complex\tClassification\n")
            f.write(f"{mlst_results['sample']}\t{mlst_results['st']}\t{mlst_results['scheme']}\t{mlst_results.get('mlst_status', 'N/A')}\t{mlst_results.get('identity', 'N/A')}\t{mlst_results.get('coverage', 'N/A')}\t{mlst_results.get('clonal_complex', 'Unknown')}\t{mlst_results.get('classification', 'Unknown')}\n")

    def create_mlst_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        print("📊 Creating MLST batch summary files...")
        
        # HTML Summary
        assigned = sum(1 for r in all_results.values() if r.get('mlst_status') == 'Assigned')
        rows = ''
        for name, r in all_results.items():
            color = '#10b981' if r.get('mlst_status') == 'Assigned' else '#dc2626'
            rows += f'<tr><td><strong>{name}</strong></td><td class="st-cell">ST{r.get("st", "ND")}</td><td style="color: {color}">{r.get("mlst_status", "N/A")}</td><td>{r.get("identity", "N/A")}</td><td class="allele-cell">{r.get("allele_profile", "")}</td><td>{r.get("clonal_complex", "Unknown")}</td><td>{r.get("outbreak_potential", "UNKNOWN")}</td></tr>'
            
        content = f"""
        <div class="report-section">
            <h2>📊 MLST Summary - All Samples</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-value">{len(all_results)}</div><div>TOTAL SAMPLES</div></div>
                <div class="metric-card"><div class="metric-value">{assigned}</div><div>ASSIGNED</div></div>
                <div class="metric-card"><div class="metric-value">{len(all_results)-assigned}</div><div>NOT ASSIGNED</div></div>
            </div>
            <div style="overflow-x: auto;">
                <table class="summary-table">
                    <thead><tr><th>Sample</th><th>ST</th><th>MLST Status</th><th>Identity</th><th>Allele Profile</th><th>Clonal Complex</th><th>Outbreak Potential</th></tr></thead>
                    <tbody>{rows}</tbody>
                </table>
            </div>
        </div>
        """
        with open(output_dir / "mlst_summary.html", 'w', encoding='utf-8') as f:
            f.write(self._get_base_html("MLST Batch Summary", content))
            
        # JSON Summary
        json_data = {
            "metadata": {"analysis_date": datetime.now().isoformat(), "total_samples": len(all_results), "scheme": "klebsiella"},
            "samples": {name: {"sequence_type": res.get('st', 'ND'), "mlst_status": res.get('mlst_status'), "allele_profile": res.get('allele_profile'), "clonal_complex": res.get('clonal_complex')} for name, res in all_results.items()}
        }
        with open(output_dir / "mlst_summary.json", 'w') as f:
            json.dump(json_data, f, indent=2)

    def run_mlst_batch(self, input_path: str, output_dir: Path, scheme: str = "klebsiella") -> Dict[str, Dict]:
        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}
            
        print(f"📁 Found {len(fasta_files)} FASTA files")
        results = {f.name: self.run_mlst_single(f, output_dir, scheme) for f in fasta_files}
        self.create_mlst_summary(results, output_dir)
        return results

def main():
    parser = argparse.ArgumentParser(description='Kleboscope MLST Analyzer - Orchestrator Compliant')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA pattern')
    parser.add_argument('-o', '--output-dir', required=True, dest='output_dir', help='Output directory')
    parser.add_argument('-db', '--database-dir', required=True, help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True, help='Script directory (contains bin/mlst)')
    parser.add_argument('-s', '--scheme', default='klebsiella', help='MLST scheme')
    parser.add_argument('--batch', action='store_true', help='Process multiple files')
    
    args = parser.parse_args()
    
    analyzer = KleboscopeMLSTAnalyzer(Path(args.database_dir), Path(args.script_dir))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.batch:
        results = analyzer.run_mlst_batch(args.input, output_dir, args.scheme)
        print(f"🎉 Batch MLST completed! Processed {len(results)} samples")
    else:
        for fasta_file in analyzer.find_fasta_files(args.input):
            analyzer.run_mlst_single(fasta_file, output_dir, args.scheme)

if __name__ == "__main__":
    main()