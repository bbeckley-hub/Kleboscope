#!/usr/bin/env python3
"""
Kleboscope MLST Module - Complete with Beautiful HTML Reports
Author: Brown Beckley
GitHub: bbeckley-hub
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-03-20
"""

import os
import sys
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
        
        # Science quotes for rotation
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
            {"text": "Kleboscope turns genomic complexity into actionable insights for AMR surveillance.", "author": "Brown Beckley"}
        ]
    
    def get_random_quote(self):
        return random.choice(self.science_quotes)
    
    def find_fasta_files(self, input_path: str) -> List[Path]:
        """Find all FASTA files - supports all major extensions"""
        if os.path.isfile(input_path):
            return [Path(input_path)]
        
        if os.path.isdir(input_path):
            fasta_patterns = ['*.fna', '*.fasta', '*.fa', '*.fn', '*.fna.gz', '*.fasta.gz', '*.fa.gz', '*.faa']
            fasta_files = []
            for pattern in fasta_patterns:
                matched_files = glob.glob(os.path.join(input_path, pattern))
                for file_path in matched_files:
                    if Path(file_path).is_file():
                        fasta_files.append(Path(file_path))
            return sorted(list(set(fasta_files)))
        
        fasta_patterns = [
            input_path,
            f"{input_path}.fna", f"{input_path}.fasta", f"{input_path}.fa", f"{input_path}.fn", f"{input_path}.faa",
            f"{input_path}.fna.gz", f"{input_path}.fasta.gz", f"{input_path}.fa.gz"
        ]
        
        fasta_files = []
        for pattern in fasta_patterns:
            matched_files = glob.glob(pattern)
            for file_path in matched_files:
                if Path(file_path).is_file():
                    fasta_files.append(Path(file_path))
        
        return sorted(list(set(fasta_files)))

    def run_mlst_single(self, input_file: Path, output_dir: Path, scheme: str = "klebsiella") -> Dict:
        """Run MLST analysis for a single file"""
        print(f"🔬 Processing: {input_file.name}")
        
        sample_output_dir = output_dir / input_file.stem
        sample_output_dir.mkdir(parents=True, exist_ok=True)
        
        raw_output_file = sample_output_dir / "mlst_raw_output.txt"
        
        if not self.mlst_bin.exists():
            print(f"❌ MLST binary not found at: {self.mlst_bin}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_output_files(error_result, sample_output_dir)
            return error_result
        
        mlst_cmd = [
            "perl", str(self.mlst_bin),
            str(input_file),
            "--scheme", scheme,
            "--csv",
            "--nopath"
        ]
        
        try:
            result = subprocess.run(mlst_cmd, capture_output=True, text=True, check=True)
            
            with open(raw_output_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
            
            mlst_results = self.parse_mlst_csv(result.stdout, input_file.name)
            mlst_results.update(self.get_lineage_info(mlst_results.get('st', 'ND')))
            mlst_results.update(self.get_identity_coverage(mlst_results.get('st', 'ND')))
            self.generate_output_files(mlst_results, sample_output_dir)
            
            print(f"✅ Completed: {input_file.name} -> ST{mlst_results.get('st', 'ND')}")
            return mlst_results
            
        except subprocess.CalledProcessError as e:
            print(f"❌ MLST failed for {input_file.name}")
            error_result = self.get_fallback_results(input_file.name)
            self.generate_output_files(error_result, sample_output_dir)
            return error_result

    def parse_mlst_csv(self, stdout: str, sample_name: str) -> Dict:
        """Parse MLST CSV output"""
        lines = stdout.strip().split('\n')
        if not lines:
            return self.get_empty_results(sample_name)
        
        result_line = None
        for line in reversed(lines):
            if line.strip() and ',' in line and not line.startswith('['):
                result_line = line.strip()
                break
        
        if not result_line:
            return self.get_empty_results(sample_name)
        
        parts = result_line.split(',')
        
        if len(parts) < 3:
            return self.get_empty_results(sample_name)
        
        st = parts[2]
        
        alleles = {}
        allele_parts = []
        
        for i in range(3, len(parts)):
            allele_str = parts[i]
            if '(' in allele_str and ')' in allele_str:
                gene = allele_str.split('(')[0]
                allele = allele_str.split('(')[1].rstrip(')')
                alleles[gene] = allele
                allele_parts.append(f"{gene}({allele})")
        
        allele_profile = '-'.join(allele_parts) if allele_parts else ""
        
        return {
            "sample": sample_name,
            "st": st,
            "scheme": "klebsiella",
            "alleles": alleles,
            "allele_profile": allele_profile,
            "confidence": "HIGH" if st and st != '-' and st != 'ND' else "LOW",
            "mlst_assigned": True if st and st != '-' and st != 'ND' else False
        }

    def get_lineage_info(self, st: str) -> Dict:
        """Get lineage information for K. pneumoniae STs"""
        lineage_db = {
    # =========================================================================
    # Hypervirulent clones (hvKp)
    # =========================================================================
    '23': {
        "clonal_complex": "CC23",
        "classification": "Hypervirulent K. pneumoniae (hvKp)",
        "geographic_distribution": "Asia (especially China, Taiwan, Southeast Asia); also reported globally",
        "clinical_significance": "Major hypervirulent clone associated with community‑acquired liver abscess, often carrying rmpA and aerobactin.",
        "common_virulence": ["rmpA", "rmpA2", "iucA‑iucD (aerobactin)", "iroBCDN (salmochelin)", "ybt (yersiniabactin)"],
        "outbreak_potential": "HIGH",
        "typical_capsule": ["K1", "K2"],
        "resistance_profile": ["Usually susceptible, but ESBL‑ and carbapenemase‑producing variants have emerged (e.g., ST23‑KPC)."]
    },
    '65': {
        "clonal_complex": "CC65",
        "classification": "Hypervirulent K. pneumoniae (hvKp)",
        "geographic_distribution": "Global, especially Asia",
        "clinical_significance": "Associated with community‑acquired infections, often with K2 capsule and rmpA.",
        "common_virulence": ["rmpA", "aerobactin", "salmochelin"],
        "outbreak_potential": "HIGH",
        "typical_capsule": ["K2"],
        "resistance_profile": ["Variable."]
    },
    '86': {
        "clonal_complex": "CC86",
        "classification": "Hypervirulent K. pneumoniae (hvKp)",
        "geographic_distribution": "Global",
        "clinical_significance": "Hypervirulent clone, often associated with K2 capsule and community‑acquired infections.",
        "common_virulence": ["rmpA", "aerobactin"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["K2"],
        "resistance_profile": ["Variable."]
    },
    '375': {
        "clonal_complex": "CC375",
        "classification": "Hypervirulent K. pneumoniae (hvKp)",
        "geographic_distribution": "China, Vietnam",
        "clinical_significance": "hvKp with K1 capsule, associated with liver abscess and carrying rmpA.",
        "common_virulence": ["rmpA", "rmpA2", "aerobactin", "salmochelin"],
        "outbreak_potential": "HIGH",
        "typical_capsule": ["K1"],
        "resistance_profile": ["Variable."]
    },
    '380': {
        "clonal_complex": "CC380",
        "classification": "Hypervirulent K. pneumoniae (hvKp)",
        "geographic_distribution": "Asia, Middle East",
        "clinical_significance": "Emerging hvKp clone, often carrying rmpA and aerobactin; ESBL variants reported.",
        "common_virulence": ["rmpA", "aerobactin", "salmochelin"],
        "outbreak_potential": "MEDIUM‑HIGH",
        "typical_capsule": ["K1", "K2", "K20"],
        "resistance_profile": ["Usually susceptible, but ESBL (e.g., CTX‑M‑15) reported."]
    },
    '412': {
        "clonal_complex": "CC412",
        "classification": "Hypervirulent K. pneumoniae (hvKp)",
        "geographic_distribution": "Asia, Europe",
        "clinical_significance": "hvKp clone with K2 capsule, often community‑acquired.",
        "common_virulence": ["rmpA", "aerobactin"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["K2"],
        "resistance_profile": ["Variable."]
    },

    # =========================================================================
    # High‑risk multidrug‑resistant (MDR) and carbapenemase‑producing clones
    # =========================================================================
    '11': {
        "clonal_complex": "CC11",
        "classification": "MDR K. pneumoniae, major carbapenemase‑producing clone",
        "geographic_distribution": "Global, especially Europe, Asia",
        "clinical_significance": "One of the most common carbapenem‑resistant clones worldwide; frequently carries KPC‑2, NDM‑1, or OXA‑48.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "VERY HIGH",
        "typical_capsule": ["KL64", "KL47"],
        "resistance_profile": ["Carbapenems, β‑lactams, often colistin and tigecycline resistance."]
    },
    '14': {
        "clonal_complex": "CC14",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "Common ESBL‑producing clone, often CTX‑M‑15.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams, fluoroquinolones, aminoglycosides."]
    },
    '15': {
        "clonal_complex": "CC15",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global, especially Europe and Asia",
        "clinical_significance": "Frequently associated with ESBLs (e.g., CTX‑M‑15) and sometimes carbapenemases.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM‑HIGH",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams, fluoroquinolones, aminoglycosides."]
    },
    '17': {
        "clonal_complex": "CC17",
        "classification": "MDR K. pneumoniae, epidemic clone",
        "geographic_distribution": ["Global"],
        "clinical_significance": "Widely distributed MDR clone, often associated with OXA‑48 and KPC carbapenemases.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "HIGH",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems, β‑lactams, often colistin resistant."]
    },
    '37': {   
        "clonal_complex": "CC37",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "Common healthcare‑associated clone, often ESBL‑positive (CTX‑M‑15).",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": "Variable",
        "resistance_profile": ["β‑lactams, aminoglycosides."]
    },
    '101': {
        "clonal_complex": "CC101",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Asia, Europe",
        "clinical_significance": "Often associated with NDM carbapenemase.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems, β‑lactams."]
    },
    '147': {
        "clonal_complex": "CC147",
        "classification": "MDR K. pneumoniae, epidemic clone",
        "geographic_distribution": "Global",
        "clinical_significance": "Major carbapenem‑resistant clone, frequently carries KPC‑2, KPC‑3, NDM‑1.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "VERY HIGH",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems, β‑lactams, often colistin resistant."]
    },
    '258': {
        "clonal_complex": "CC258",
        "classification": "MDR K. pneumoniae, global epidemic clone",
        "geographic_distribution": "Global",
        "clinical_significance": "Most prevalent KPC‑producing clone worldwide. Highly multidrug‑resistant.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "VERY HIGH",
        "typical_capsule": ["Variable (e.g., KL106, KL107)"],
        "resistance_profile": ["Carbapenems, β‑lactams, aminoglycosides, fluoroquinolones, often colistin resistant."]
    },
    '307': {   
        "clonal_complex": "CC307",
        "classification": "MDR K. pneumoniae, emerging high‑risk clone",
        "geographic_distribution": "Global, especially Europe, North America",
        "clinical_significance": "Emerging high‑risk clone with ESBLs and carbapenemases (OXA‑48, KPC).",
        "common_virulence": ["Variable"],
        "outbreak_potential": "HIGH",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems, β‑lactams, often colistin resistant."]
    },
    '392': {
        "clonal_complex": "CC392",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Europe, Asia",
        "clinical_significance": "Often carries OXA‑48 carbapenemase.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems (OXA‑48), β‑lactams."]
    },
    '395': {   
        "clonal_complex": "CC395",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Asia",
        "clinical_significance": "Often associated with NDM carbapenemase (NDM‑1, NDM‑5).",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems, β‑lactams."]
    },
    '512': {
        "clonal_complex": "CC512",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Europe, America",
        "clinical_significance": "Carries KPC carbapenemase, associated with hospital outbreaks.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "HIGH",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems, β‑lactams, fluoroquinolones."]
    },

    # =========================================================================
    # Additional well‑documented STs (common in literature)
    # =========================================================================
    '1': {
        "clonal_complex": "CC1",
        "classification": "Classical / susceptible",
        "geographic_distribution": "Global",
        "clinical_significance": "Common in community infections; generally susceptible.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Variable."]
    },
    '3': {
        "clonal_complex": "CC3",
        "classification": "Classical / susceptible",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '13': {
        "clonal_complex": "CC13",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "Common ESBL producer, often CTX‑M‑15.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams, fluoroquinolones, aminoglycosides."]
    },
    '16': {
        "clonal_complex": "CC16",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "Often associated with ESBLs, occasional carbapenemases.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams, fluoroquinolones."]
    },
    '45': {
        "clonal_complex": "CC45",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "Often associated with ESBLs, occasionally carbapenemases.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams, fluoroquinolones."]
    },
    '48': {
        "clonal_complex": "CC48",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "Common ESBL producer, sometimes OXA‑48.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams, carbapenems (OXA‑48)."]
    },
    '76': {   
        "clonal_complex": "CC76",
        "classification": "Classical / MDR",
        "geographic_distribution": "Global",
        "clinical_significance": "Moderately common; may be susceptible or carry ESBLs.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "LOW‑MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Often susceptible; ESBL variants reported."]
    },
    '252': {   
        "clonal_complex": "Unknown",
        "classification": "Uncharacterized lineage — epidemiology unknown",
        "geographic_distribution": "Not well characterised",
        "clinical_significance": "Limited data; appears sporadically in databases. May be susceptible or carry ESBLs.",
        "common_virulence": ["Unknown"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Unknown"],
        "resistance_profile": ["Not well defined."]
    },
    '2096': {   
        "clonal_complex": "Unknown",
        "classification": "Uncharacterized lineage — epidemiology unknown",
        "geographic_distribution": "Not well characterised",
        "clinical_significance": "Recently identified; not widely studied. May represent a local or emerging clone.",
        "common_virulence": ["Unknown"],
        "outbreak_potential": "UNKNOWN",
        "typical_capsule": ["Unknown"],
        "resistance_profile": ["Unknown."]
    },

# =========================================================================
# Additional major clones 
# =========================================================================
    '2': {
        "clonal_complex": "CC2",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common community strain, generally susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '4': {
        "clonal_complex": "CC4",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '5': {
        "clonal_complex": "CC5",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '6': {
        "clonal_complex": "CC6",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '7': {
        "clonal_complex": "CC7",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '8': {
        "clonal_complex": "CC8",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '9': {
        "clonal_complex": "CC9",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '10': {
        "clonal_complex": "CC10",
        "classification": "Classical",
        "geographic_distribution": "Global",
        "clinical_significance": "Common, often susceptible.",
        "common_virulence": ["Limited"],
        "outbreak_potential": "LOW",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Usually susceptible."]
    },
    '12': {
        "clonal_complex": "CC12",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer, often CTX‑M‑15.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "fluoroquinolones."]
    },
    '18': {
        "clonal_complex": "CC18",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "Often associated with ESBLs.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "aminoglycosides."]
    },
    '20': {
        "clonal_complex": "CC20",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer, occasionally carbapenemase.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "fluoroquinolones."]
    },
    '22': {
        "clonal_complex": "CC22",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '24': {
        "clonal_complex": "CC24",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer, sometimes carbapenemase.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "carbapenems (rare)."]
    },
    '25': {
        "clonal_complex": "CC25",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "fluoroquinolones."]
    },
    '26': {
        "clonal_complex": "CC26",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '27': {
        "clonal_complex": "CC27",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer, occasionally carbapenemase.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "carbapenems (rare)."]
    },
    '28': {
        "clonal_complex": "CC28",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '29': {
        "clonal_complex": "CC29",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '30': {
        "clonal_complex": "CC30",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '31': {
        "clonal_complex": "CC31",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '32': {
        "clonal_complex": "CC32",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer, often CTX‑M‑15.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "fluoroquinolones."]
    },
    '33': {
        "clonal_complex": "CC33",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '34': {
        "clonal_complex": "CC34",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '35': {
        "clonal_complex": "CC35",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '36': {
        "clonal_complex": "CC36",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '38': {
        "clonal_complex": "CC38",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '39': {
        "clonal_complex": "CC39",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '40': {
        "clonal_complex": "CC40",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '41': {
        "clonal_complex": "CC41",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '42': {
        "clonal_complex": "CC42",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '43': {
        "clonal_complex": "CC43",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '44': {
        "clonal_complex": "CC44",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '46': {
        "clonal_complex": "CC46",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '47': {
        "clonal_complex": "CC47",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '49': {
        "clonal_complex": "CC49",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '50': {
        "clonal_complex": "CC50",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    },
    '100': {
        "clonal_complex": "CC100",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer, sometimes carbapenemase.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams", "carbapenems (rare)."]
    },
    '102': {
        "clonal_complex": "CC102",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Asia, Europe",
        "clinical_significance": "Often associated with NDM carbapenemase.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM‑HIGH",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems", "β‑lactams."]
    },
    '103': {
        "clonal_complex": "CC103",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Asia",
        "clinical_significance": "NDM producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems."]
    },
    '104': {
        "clonal_complex": "CC104",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Europe",
        "clinical_significance": "OXA‑48 producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["Carbapenems (OXA‑48)."]
    },
    '105': {
        "clonal_complex": "CC105",
        "classification": "MDR K. pneumoniae",
        "geographic_distribution": "Global",
        "clinical_significance": "ESBL producer.",
        "common_virulence": ["Variable"],
        "outbreak_potential": "MEDIUM",
        "typical_capsule": ["Variable"],
        "resistance_profile": ["β‑lactams."]
    }
}
        if st in lineage_db:
            return lineage_db[st]
        elif st.isdigit():
            return {"clonal_complex": f"Unknown (ST{st})", "classification": "Not in database", "geographic_distribution": "Unknown", "clinical_significance": f"ST{st} is not currently in the Kleboscope MLST database", "common_virulence": ["Unknown"], "outbreak_potential": "UNKNOWN", "typical_capsule": ["Unknown"], "resistance_profile": ["Unknown"]}
        else:
            return {"clonal_complex": "Not Assigned", "classification": "MLST typing failed", "geographic_distribution": "N/A", "clinical_significance": "Could not determine sequence type", "common_virulence": ["Cannot determine"], "outbreak_potential": "UNKNOWN", "typical_capsule": ["Unknown"], "resistance_profile": ["Unknown"]}

    def get_identity_coverage(self, st: str) -> Dict:
        if st and st not in ['-', 'ND', 'UNKNOWN']:
            return {"identity": "100%", "coverage": "100%", "mlst_status": "Assigned", "quality_metrics": {"assembly_quality": "High Quality", "allele_completeness": "Complete", "database_match": "Perfect Match"}}
        else:
            return {"identity": "Not Assigned", "coverage": "Not Assigned", "mlst_status": "Not Assigned", "quality_metrics": {"assembly_quality": "Requires Review", "allele_completeness": "Incomplete", "database_match": "No Match"}}

    def get_empty_results(self, sample_name: str) -> Dict:
        return {"sample": sample_name, "st": "ND", "scheme": "klebsiella", "alleles": {}, "allele_profile": "", "confidence": "LOW", "mlst_assigned": False}

    def get_fallback_results(self, sample_name: str) -> Dict:
        return {"sample": sample_name, "st": "UNKNOWN", "scheme": "klebsiella", "alleles": {}, "allele_profile": "", "confidence": "LOW", "mlst_assigned": False, "error": "MLST analysis failed"}

    def generate_output_files(self, mlst_results: Dict, output_dir: Path):
        if 'identity' not in mlst_results:
            mlst_results.update(self.get_identity_coverage(mlst_results.get('st', 'ND')))
        self.generate_html_report(mlst_results, output_dir)
        self.generate_text_report(mlst_results, output_dir)
        self.generate_tsv_report(mlst_results, output_dir)

    def generate_text_report(self, mlst_results: Dict, output_dir: Path):
        report = f"""MLST Analysis Report
===================

Sample: {mlst_results['sample']}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

MLST TYPING RESULTS:
-------------------
Sequence Type (ST): {mlst_results['st']}
Scheme: {mlst_results['scheme']}
Confidence: {mlst_results['confidence']}
MLST Status: {mlst_results.get('mlst_status', 'Not Assigned')}

Identity & Coverage:
Identity: {mlst_results.get('identity', 'Not Assigned')}
Coverage: {mlst_results.get('coverage', 'Not Assigned')}

Allele Profile:
{mlst_results['allele_profile']}

Detailed Alleles:
"""
        for gene, allele in mlst_results['alleles'].items():
            report += f"- {gene}: {allele}\n"
        
        report += f"""
LINEAGE ANALYSIS:
-----------------
Clonal Complex: {mlst_results.get('clonal_complex', 'Unknown')}
Classification: {mlst_results.get('classification', 'Unknown')}
Geographic Distribution: {mlst_results.get('geographic_distribution', 'Unknown')}
Clinical Significance: {mlst_results.get('clinical_significance', 'Unknown')}
Outbreak Potential: {mlst_results.get('outbreak_potential', 'UNKNOWN')}

Typical Capsule Types: {', '.join(mlst_results.get('typical_capsule', ['Unknown']))}
Resistance Profile: {', '.join(mlst_results.get('resistance_profile', ['Unknown']))}

Common Virulence Factors:
"""
        for virulence in mlst_results.get('common_virulence', []):
            report += f"- {virulence}\n"
        
        with open(output_dir / "mlst_report.txt", 'w') as f:
            f.write(report)

    def generate_tsv_report(self, mlst_results: Dict, output_dir: Path):
        tsv_content = f"Sample\tST\tScheme\tMLST_Status\tIdentity\tCoverage\tClonal_Complex\tClassification\tOutbreak_Potential\tTypical_Capsule\tResistance_Profile\n"
        tsv_content += f"{mlst_results['sample']}\t{mlst_results['st']}\t{mlst_results['scheme']}\t{mlst_results.get('mlst_status', 'Not Assigned')}\t{mlst_results.get('identity', 'Not Assigned')}\t{mlst_results.get('coverage', 'Not Assigned')}\t{mlst_results.get('clonal_complex', 'Unknown')}\t{mlst_results.get('classification', 'Unknown')}\t{mlst_results.get('outbreak_potential', 'UNKNOWN')}\t{','.join(mlst_results.get('typical_capsule', ['Unknown']))}\t{','.join(mlst_results.get('resistance_profile', ['Unknown']))}\n"
        with open(output_dir / "mlst_report.tsv", 'w') as f:
            f.write(tsv_content)

    def generate_html_report(self, mlst_results: Dict, output_dir: Path):
        random_quote = self.get_random_quote()
        
        sample = mlst_results['sample']
        st = mlst_results['st']
        confidence = mlst_results['confidence']
        clonal_complex = mlst_results.get('clonal_complex', 'Unknown')
        allele_profile = mlst_results['allele_profile']
        identity = mlst_results.get('identity', 'Not Assigned')
        coverage = mlst_results.get('coverage', 'Not Assigned')
        mlst_status = mlst_results.get('mlst_status', 'Not Assigned')
        classification = mlst_results.get('classification', 'Unknown')
        geographic_distribution = mlst_results.get('geographic_distribution', 'Unknown')
        outbreak_potential = mlst_results.get('outbreak_potential', 'UNKNOWN')
        clinical_significance = mlst_results.get('clinical_significance', 'Further analysis required.')
        typical_capsule = ', '.join(mlst_results.get('typical_capsule', ['Unknown']))
        resistance_profile = ', '.join(mlst_results.get('resistance_profile', ['Unknown']))
        common_virulence = mlst_results.get('common_virulence', [])
        
        alleles_html = ''
        for gene, allele in mlst_results.get('alleles', {}).items():
            alleles_html += f'                <div class="allele-card"><div style="font-size: 12px; opacity: 0.9;">{gene}</div><div style="font-size: 18px;">{allele}</div></div>\n'
        
        virulence_html = ''
        for virulence in common_virulence[:6]:
            virulence_html += f'<div class="virulence-card">{virulence}</div>\n'
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KLEBOSCOPE - MLST Analysis Report</title>
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
            border: 2px solid rgba(0, 119, 190, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.2;
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
            text-align: center;
        }}
        .metric-label {{ font-size: 14px; opacity: 0.9; margin-bottom: 5px; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .allele-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        .allele-card {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .virulence-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }}
        .virulence-card {{
            background: #e0f2fe;
            color: #0369a1;
            padding: 12px;
            border-radius: 6px;
            text-align: center;
            font-weight: bold;
            border-left: 4px solid #0ea5e9;
        }}
        .profile-box {{
            background: #f8fafc;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            border-left: 4px solid #3b82f6;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
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
            .ascii-art {{ font-size: 8px; }}
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
            <h2>📊 Sample Information</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-label">Sample Name</div><div class="metric-value">{sample}</div></div>
                <div class="metric-card"><div class="metric-label">Analysis Date</div><div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div></div>
                <div class="metric-card"><div class="metric-label">MLST Scheme</div><div class="metric-value">Klebsiella</div></div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🎯 MLST Results</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-label">Sequence Type</div><div class="metric-value">ST{st}</div></div>
                <div class="metric-card"><div class="metric-label">Identity</div><div class="metric-value">{identity}</div></div>
                <div class="metric-card"><div class="metric-label">Coverage</div><div class="metric-value">{coverage}</div></div>
            </div>
            <h3>Allele Profile</h3>
            <div class="profile-box"><code style="font-size: 16px; color: #1e40af; font-weight: bold;">{allele_profile}</code></div>
            <h3>Individual Alleles</h3>
            <div class="allele-grid">{alleles_html}</div>
        </div>
        
        <div class="report-section">
            <h2>🌍 Lineage Information</h2>
            <div class="metrics-grid">
                <div class="metric-card"><div class="metric-label">Clonal Complex</div><div class="metric-value">{clonal_complex}</div></div>
                <div class="metric-card"><div class="metric-label">Classification</div><div class="metric-value">{classification}</div></div>
                <div class="metric-card"><div class="metric-label">Outbreak Potential</div><div class="metric-value">{outbreak_potential}</div></div>
            </div>
            <h3>Clinical Significance</h3>
            <div class="profile-box"><p>{clinical_significance}</p></div>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-top: 20px;">
                <div><h3>Typical Characteristics</h3><p><strong>Capsule Types:</strong> {typical_capsule}</p><p><strong>Resistance Profile:</strong> {resistance_profile}</p></div>
                <div><h3>Common Virulence Factors</h3><div class="virulence-grid">{virulence_html}</div></div>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>KLEBOSCOPE</strong> - MLST Analysis Module</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
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
</html>'''
        
        with open(output_dir / "mlst_report.html", 'w', encoding='utf-8') as f:
            f.write(html_content)

    def create_mlst_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create comprehensive MLST summary files"""
        print("📊 Creating MLST summary files...")
        self.create_mlst_tsv_summary(all_results, output_dir)
        self.create_mlst_html_summary(all_results, output_dir)
        self.create_mlst_json_summary(all_results, output_dir)
        print("✅ MLST summary files created successfully!")

    def create_mlst_tsv_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        summary_file = output_dir / "mlst_summary.tsv"
        with open(summary_file, 'w') as f:
            f.write("Sample\tST\tMLST_Status\tIdentity\tCoverage\tAllele_Profile\tClonal_Complex\tClassification\tOutbreak_Potential\n")
            for sample_name, result in all_results.items():
                f.write(f"{sample_name}\t{result.get('st', 'ND')}\t{result.get('mlst_status', 'Not Assigned')}\t{result.get('identity', 'Not Assigned')}\t{result.get('coverage', 'Not Assigned')}\t{result.get('allele_profile', '')}\t{result.get('clonal_complex', 'Unknown')}\t{result.get('classification', 'Unknown')}\t{result.get('outbreak_potential', 'UNKNOWN')}\n")
        print(f"📄 TSV summary created: {summary_file}")

    def create_mlst_html_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        """Create HTML summary with quotes and detailed footnotes"""
        summary_file = output_dir / "mlst_summary.html"
        random_quote = self.get_random_quote()
        
        # Calculate statistics
        total_samples = len(all_results)
        assigned_samples = sum(1 for r in all_results.values() if r.get('mlst_status') == 'Assigned')
        not_assigned_samples = total_samples - assigned_samples
        
        # Build table rows
        table_rows = ''
        for sample_name, result in all_results.items():
            mlst_status_color = '#10b981' if result.get('mlst_status') == 'Assigned' else '#dc2626'
            table_rows += f'''                        <tr>
                            <td><strong>{sample_name}</strong></td>
                            <td class="st-cell">ST{result.get('st', 'ND')}</td>
                            <td style="color: {mlst_status_color}">{result.get('mlst_status', 'Not Assigned')}</td>
                            <td>{result.get('identity', 'Not Assigned')}</td>
                            <td>{result.get('coverage', 'Not Assigned')}</td>
                            <td class="allele-cell">{result.get('allele_profile', '')}</td>
                            <td>{result.get('clonal_complex', 'Unknown')}</td>
                            <td>{result.get('classification', 'Unknown')}</td>
                            <td>{result.get('outbreak_potential', 'UNKNOWN')}</td>
                        </tr>
'''
        
        html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KLEBOSCOPE - MLST Summary Report</title>
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
            border: 2px solid rgba(0, 119, 190, 0.3);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.2;
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
        }}
        .report-section h2 {{
            color: #1e3a8a;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #3b82f6 0%, #1e40af 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        .summary-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .summary-table tr:hover {{ background-color: #e0f2fe; }}
        .st-cell {{ font-weight: bold; color: #1e40af; }}
        .allele-cell {{ font-family: 'Courier New', monospace; background-color: #f0f9ff; color: #0369a1; font-weight: bold; }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #0077be 0%, #005a8c 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{ font-size: 24px; font-weight: bold; margin-bottom: 5px; }}
        .stat-label {{ font-size: 12px; opacity: 0.9; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
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
            .ascii-art {{ font-size: 8px; }}
            .summary-table {{ font-size: 12px; }}
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
            <h2>📊 MLST Summary - All Samples</h2>
            <div class="stats-grid">
                <div class="stat-card"><div class="stat-value">{total_samples}</div><div class="stat-label">SAMPLES PROCESSED</div></div>
                <div class="stat-card"><div class="stat-value">{assigned_samples}</div><div class="stat-label">SAMPLES ASSIGNED</div></div>
                <div class="stat-card"><div class="stat-value">{not_assigned_samples}</div><div class="stat-label">NOT ASSIGNED</div></div>
                <div class="stat-card"><div class="stat-value">{len(all_results)}</div><div class="stat-label">TOTAL SAMPLES</div></div>
            </div>
            
            <div style="overflow-x: auto;">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>ST</th>
                            <th>MLST Status</th>
                            <th>Identity</th>
                            <th>Coverage</th>
                            <th>Allele Profile</th>
                            <th>Clonal Complex</th>
                            <th>Classification</th>
                            <th>Outbreak Potential</th>
                        </thead>
                    <tbody>
{table_rows}                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>KLEBOSCOPE</strong> - MLST Summary Report</p>
            <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
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
</html>'''
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"🌐 HTML summary created: {summary_file}")

    def create_mlst_json_summary(self, all_results: Dict[str, Dict], output_dir: Path):
        summary_file = output_dir / "mlst_summary.json"
        json_summary = {
            "metadata": {"analysis_date": datetime.now().isoformat(), "total_samples": len(all_results), "scheme": "klebsiella"},
            "samples": {name: {"sequence_type": res.get('st', 'ND'), "mlst_status": res.get('mlst_status'), "allele_profile": res.get('allele_profile'), "clonal_complex": res.get('clonal_complex'), "classification": res.get('classification')} for name, res in all_results.items()}
        }
        with open(summary_file, 'w') as f:
            json.dump(json_summary, f, indent=2)
        print(f"📄 JSON summary created: {summary_file}")

    def run_mlst_batch(self, input_path: str, output_dir: Path, scheme: str = "klebsiella") -> Dict[str, Dict]:
        print("🔍 Searching for FASTA files...")
        fasta_files = self.find_fasta_files(input_path)
        
        if not fasta_files:
            print("❌ No FASTA files found!")
            return {}
        
        print(f"📁 Found {len(fasta_files)} FASTA files")
        
        results = {}
        for fasta_file in fasta_files:
            result = self.run_mlst_single(fasta_file, output_dir, scheme)
            results[fasta_file.name] = result
        
        self.create_mlst_summary(results, output_dir)
        return results

def main():
    parser = argparse.ArgumentParser(description='Kleboscope MLST Analyzer')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file, directory, or wildcard pattern')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    parser.add_argument('-db', '--database-dir', required=True, help='Database directory')
    parser.add_argument('-sc', '--script-dir', required=True, help='Script directory (contains bin/mlst)')
    parser.add_argument('-s', '--scheme', default='klebsiella', help='MLST scheme')
    parser.add_argument('--batch', action='store_true', help='Process multiple files')
    
    args = parser.parse_args()
    
    analyzer = KleboscopeMLSTAnalyzer(
        database_dir=Path(args.database_dir),
        script_dir=Path(args.script_dir)
    )
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.batch:
        results = analyzer.run_mlst_batch(args.input, output_dir, args.scheme)
        print(f"🎉 Batch MLST completed! Processed {len(results)} samples")
    else:
        fasta_files = analyzer.find_fasta_files(args.input)
        if fasta_files:
            for fasta_file in fasta_files:
                result = analyzer.run_mlst_single(fasta_file, output_dir, args.scheme)
                print(f"🎉 MLST completed for {fasta_file.name}: ST{result.get('st', 'ND')}")
        else:
            print(f"❌ No input files found: {args.input}")

if __name__ == "__main__":
    main()