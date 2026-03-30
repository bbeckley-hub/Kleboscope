# 🧬 Kleboscope

### **A gene‑centric, species‑optimized computational pipeline for comprehensive *Klebsiella pneumoniae* genomic surveillance**

#### **Complete K. pneumoniae typing, resistance, virulence, plasmid, and environmental marker analysis — from FASTA to actionable insights**

[![Version](https://img.shields.io/badge/version-1.0.0-blue)](https://github.com/bbeckley-hub/Kleboscope)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/kleboscope)](https://github.com/bbeckley-hub/Kleboscope/issues)

---

```text
██╗  ██╗██╗     ███████╗██████╗  ██████╗ ███████╗ ██████╗ ██████╗  ██████╗ ███████╗
██║ ██╔╝██║     ██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔═══██╗ ██╔══██╗██╔════╝
█████╔╝ ██║     █████╗  ██████╔╝██║   ██║███████╗██║     ██║   ██║ ██████╔╝█████╗  
██╔═██╗ ██║     ██╔══╝  ██╔══██╗██║   ██║╚════██║██║     ██║   ██║ ██╔═══╝ ██╔══╝  
██║  ██╗███████╗███████╗██████╔╝╚██████╔╝███████║╚██████╗╚██████╔╝ ██║     ███████╗
╚═╝  ╚═╝╚══════╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝  ╚═╝     ╚══════╝
```

---

## 📋 **Table of Contents**

- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [⚡ Quick Start](#-quick-start)
- [🔧 Installation](#-installation)
- [🚀 Usage Guide](#-usage-guide)
- [📁 Output Structure](#-output-structure)
- [🔍 Analytical Modules](#-analytical-modules)
- [📈 Performance](#-performance)
- [🔬 Validation](#-validation)
- [🔄 Alternative Tools](#-alternative-tools)
- [🤖 AI Integration Guide](#-ai-integration-guide)
- [🔮 Future Directions](#-future-directions)
- [❓ Frequently Asked Questions](#-frequently-asked-questions)
- [📚 Citation & Acknowledgements](#-citation--acknowledgements)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License & Third‑Party Notices](#-license--third-party-notices)

---

## 🎯 **Overview**

**Kleboscope** is an automated, locally‑executable computational pipeline designed specifically for comprehensive *Klebsiella pneumoniae* genomic surveillance. It addresses the growing threat of multidrug‑resistant and hypervirulent *K. pneumoniae* by integrating **seven essential analysis modules** into a single, cohesive workflow. Instead of listing genes per sample, Kleboscope presents a **gene‑centric** view: each gene is shown with **all genomes** that contain it, together with its frequency (`count (percentage%)`), enabling rapid pattern discovery.

### 🌍 **The Problem**
- **Fragmented Analysis**: *K. pneumoniae* typing requires separate tools for MLST, capsule typing, AMR detection, virulence screening, plasmid profiling, and quality control.
- **Time‑Consuming Workflows**: Combining outputs from multiple tools is error‑prone and delays outbreak responses.
- **Interpretation Gap**: Raw data lacks clinical context; identifying high‑risk clones (e.g., ST11‑KL64 with *blaKPC‑2* and *iuc*) requires manual correlation.
- **Missing Environmental Markers**: Existing tools ignore biocide and heavy metal resistance genes that co‑select for antibiotic resistance.

### 💡 **Our Solution**
Kleboscope delivers:
- **✅ Single‑command installation** via Conda 
- **✅ Parallel execution** of QC, MLST, and Kaptive for maximum speed
- **✅ Comprehensive gene‑centric HTML report** with interactive tables, search, export, and scrollable genome lists (no truncation)
- **✅ Tracking of critical resistance genes**: carbapenemases (*blaKPC*, *blaNDM*, *blaOXA‑48*), colistin (*mcr*), tigecycline (*tetX*), 16S rRNA methyltransferases
- **✅ Tracking of hypervirulence markers**: aerobactin (*iuc*), salmochelin (*iro*), yersiniabactin (*ybt*), colibactin (*clb*), regulators of hypermucoidy (*rmpA*, *rmpA2*)
- **✅ Environmental co‑selection markers**: biocide resistance (*qac*), heavy metal resistance (*sil*, *mer*, *ars*, *pco*, *czc*), mobile genetic elements, stress response genes
- **✅ Built‑in pattern discovery**: ST‑capsule associations, high‑risk combinations, ICEKp and virulence plasmid tracking
- **✅ AI integration guide** – prompts to help you mine data with large language models

**Perfect for**: Clinical laboratories, outbreak investigations, research studies, and public health surveillance.

---

## ✨ **Key Features**

### 🔬 **Core Analytical Modules**

| Module | Purpose | Key Outputs | Speed* |
|--------|---------|-------------|--------|
| **FASTA QC** | Assembly statistics, N50, GC%, contig quality | HTML, TSV, JSON | <30 sec |
| **MLST Typing** | Sequence type (Pasteur scheme) | ST, allele profile, clonal complex | <1 min |
| **Kaptive Capsule Typing** | K (capsule) and O (lipopolysaccharide) loci | K/O types, identity, coverage | 1-2 min |
| **ABRicate Screening** | 11 databases (CARD, ResFinder, VFDB, PlasmidFinder, BacMet2, etc.) | Gene‑centric tables per database | 2-3 min per DB |
| **AMRfinderPlus** | Acquired resistance genes & point mutations | Gene frequency, risk levels | 3-4 min |
| **Ultimate Reporter** | Integrated gene‑centric HTML report | Interactive tables, patterns, AI guide | <1 sec |

*Timings for a single genome on a laptop; parallel execution for batch processing*

### 🛡️ **Species‑Specific Innovations for *K. pneumoniae***

- **Comprehensive resistance tracking**:
  - Carbapenemases: *blaKPC*, *blaNDM*, *blaIMP*, *blaVIM*, *blaOXA‑48*, *blaGES*, *blaIMI*, *blaSME*, *blaFRI*
  - Colistin resistance: *mcr‑1* to *mcr‑10*, *pmrAB*, *lpx*, *arn*
  - Tigecycline resistance: *tet(X)* variants, efflux pumps (*adeABC*, *acrAB*)
  - 16S rRNA methyltransferases: *armA*, *rmtA*‑*rmtH*, *npmA*
  - β‑lactamases: all major ESBLs, AmpC, and narrow‑spectrum

- **Hypervirulence marker tracking**:
  - Aerobactin (*iucA‑D*, *iutA*)
  - Salmochelin (*iroB‑E*, *iroN*)
  - Yersiniabactin (*ybtA‑X*, *irp1‑2*)
  - Colibactin (*clbA‑S*)
  - Hypermucoidy regulators (*rmpA*, *rmpA2*, *rmpC*, *rmpD*)
  - Additional virulence factors: *peg‑344*, *allS*, *kfu*, *mrkD*

- **Environmental co‑selection markers** (BACMET2 database):
  - Biocide resistance: *qac*, *cep*, *form*, *oqx*
  - Heavy metal resistance: *sil*, *mer*, *ars*, *pco*, *czc*, *znt*, *cop*, *cad*, *chr*
  - Mobile genetic elements: *tra*, *mob*, *rep*, *int*, *tnp*
  - Stress response: *sox*, *mar*, *rob*, *rpo*

- **Gene‑centric ultimate report**:
  - Each gene displayed with **all genomes** that contain it
  - Frequency as `count (percentage%)`
  - Scrollable genome lists (no truncation)
  - Search, sort, export, and print functionality
  - Pattern discovery tabs: ST‑K/O associations, high‑risk combinations, ICEKp markers, virulence plasmid markers

### 🚀 **Performance Advantages**
- **Parallel first batch**: QC, MLST, and Kaptive run concurrently → ~4 minutes for 30 genomes
- **Optimal resource usage**: Auto‑detects CPU cores and RAM using `psutil`
- **Low memory footprint**: Runs comfortably on 4 GB RAM
- **Scales linearly**: 30 genomes in 5.5 hours; 100 genomes estimated in ~1 day (single thread)

---

## ⚡ **Quick Start**

### **Install in 60 seconds**

**Conda (recommended):**
```bash
conda create -n kleboscope -c conda-forge -c bioconda -c bbeckley-hub kleboscope
conda activate kleboscope
```

**From source (for developers):**
```bash
git clone https://github.com/bbeckley-hub/Kleboscope.git
cd kleboscope
pip install -e .
```

**Verify installation:**
```bash
kleboscope --version
```

### **Run your first analysis**

```bash
# Single genome
kleboscope -i genome.fna -o results

# Batch processing (30 genomes)
kleboscope -i "*.fna" -o results --threads 4
# Complete in ~2-4 hours on a laptop 🎉
```

---

## 🔧 **Installation**

### **System Requirements**
| Resource | Minimum | Recommended |
|----------|---------|-------------|
| **CPU Cores** | 2 | 4+ |
| **RAM** | 4 GB | 8 GB |
| **Storage** | 2 GB (plus genomes) | 10 GB+ |
| **OS** | Linux, macOS, WSL2 | Linux |

### **Dependencies**
All external tools and databases are **bundled** with Kleboscope as build dependencies, but you can inspect the required Python packages in the [`environment.yml`](environment.yml) file.

**Key dependencies:**  
- Python ≥3.9
- biopython, pandas, numpy, beautifulsoup4, psutil
- Bundled binaries: MLST, Kaptive, ABRicate, AMRfinderPlus

No separate installation of external tools is required.

---

## 🚀 **Usage Guide**

### **Basic Commands**

```bash
kleboscope -i "*.fna" -o results --threads 4
kleboscope -i genome.fna -o results --skip-qc --skip-amr
```

### **Skip Options**

| Option | Effect |
|--------|--------|
| `--skip-qc` | Skip FASTA QC |
| `--skip-mlst` | Skip MLST typing |
| `--skip-kaptive` | Skip Kaptive capsule typing |
| `--skip-abricate` | Skip ABRicate screening |
| `--skip-amr` | Skip AMRfinder analysis |
| `--skip-summary` | Skip ultimate reporter generation |

### **Real‑World Examples**

**Clinical laboratory: routine surveillance of 50 isolates**
```bash
kleboscope -i "*.fna" -o weekly_surveillance --threads 4
# Results in ~5.5 hours (depends on hardware); interactive HTML report ready
```

**Outbreak investigation: hypervirulence screening**
```bash
kleboscope -i "outbreak/*.fasta" -o urgent --skip-abricate --skip-amr
# Focus on MLST + Kaptive + virulence markers; results in ~5 minutes
```

---

## 📁 **Output Structure**

```
results/
├── fasta_qc_results/               # Individual QC reports + summary
├── mlst_results/                   # MLST reports + summary
├── kaptive_results/                # Kaptive reports + summary
├── klebo_abricate_results/         # Per-database summary HTML/TSV/JSON
├── klebo_amrfinder_results/        # AMRfinder reports + summary
└── KLEBOSCOPE_ULTIMATE_REPORTS/    # Gene‑centric ultimate report
    ├── kleboscope_ultimate_report.html   # Interactive report (open in browser!)
    ├── kleboscope_ultimate_report.json   # Complete data
    ├── kleboscope_samples.csv            # Sample overview
    ├── kleboscope_amr_genes.csv          # All AMR genes with genomes
    ├── kleboscope_virulence_genes.csv    # All virulence genes with genomes
    ├── kleboscope_environmental_markers.csv  # Biocide/heavy metal genes
    ├── kleboscope_patterns.csv           # High‑risk combos, ST‑K/O associations
    └── kleboscope_database_coverage.csv  # Database performance stats
```

---

## 🔍 **Analytical Modules**

### **1. FASTA QC**
- **Metrics**: total length, contig count, N50/N75/N90, L50/L75/L90, mean/median length, GC%, AT%, ambiguous bases, N‑runs, homopolymers, duplicate sequences
- **Output**: interactive HTML with warnings, TSV/JSON for downstream analysis

### **2. MLST Typing**
- **Database**: [PubMLST](https://pubmlst.org/kpneumoniae/) (Pasteur scheme)
- **Tool**: [MLST](https://github.com/tseemann/mlst) (v2.23)
- **Output**: ST, allele profile, clonal complex, classification, outbreak potential

### **3. Kaptive Capsule Typing**
- **Databases**: `kp_k` and `kp_o` (curated by [Kaptive team](https://github.com/katholt/Kaptive))
- **Tool**: [Kaptive](https://github.com/katholt/Kaptive) (v3.1.0)
- **Output**: K locus, O locus, identity, coverage, confidence

### **4. ABRicate Screening**
- **Tool**: [ABRicate](https://github.com/tseemann/abricate) (v1.2.0)
- **Databases (11)**:  

  | Database | Description | Link |
  |----------|-------------|------|
  | CARD | Comprehensive Antibiotic Resistance Database | [https://card.mcmaster.ca/](https://card.mcmaster.ca/) |
  | ResFinder | Acquired resistance genes | [https://cge.cbs.dtu.dk/services/ResFinder/](https://cge.cbs.dtu.dk/services/ResFinder/) |
  | ARG‑ANNOT | Antibiotic resistance gene database | [https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/) |
  | VFDB | Virulence Factor Database | [http://www.mgc.ac.cn/VFs/](http://www.mgc.ac.cn/VFs/) |
  | PlasmidFinder | Plasmid replicons | [https://cge.cbs.dtu.dk/services/PlasmidFinder/](https://cge.cbs.dtu.dk/services/PlasmidFinder/) |
  | NCBI | NCBI AMR reference | [https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/) |
  | MegaRes | Megaresistance database | [https://megares.meglab.org/](https://megares.meglab.org/) |
  | EcoH | E. coli O/H typing | [https://github.com/ssi-dk/ecoh](https://github.com/ssi-dk/ecoh) |
  | EcoLI_VF | E. coli virulence factors | [https://github.com/ssi-dk/ecoli_vf](https://github.com/ssi-dk/ecoli_vf) |
  | BacMet2 | Biocide and Metal Resistance Database | [http://bacmet.biomedicine.gu.se/](http://bacmet.biomedicine.gu.se/) |
  | NCBIfam | Additional resistance families from NCBI | Included in NCBI database |

- **Thresholds**: 80% identity and coverage
- **Output per database**:
  - “Genes by Genome” table
  - “Gene Frequency” table (gene, count, percentage, genomes)

### **5. AMRfinderPlus**
- **Tool**: [AMRfinderPlus](https://github.com/ncbi/amr) (v4.2.7) 
- **Database**: 2026‑01‑21.1 (included)
- **Organism**: *Klebsiella pneumoniae* (enables point‑mutation detection)
- **Output**: per‑sample reports + summary table with risk levels (Critical, High, Standard)

### **6. Ultimate Reporter**
- **Input**: All module summary HTML files
- **Processing**:
  - Parses all tables using BeautifulSoup
  - Normalises sample names (removes extensions, standardises GCF/GCA)
  - Builds unified gene‑centric structure
- **Output**: Interactive HTML report with:
  - Dashboard cards (samples, STs, capsule types, unique genes, high‑risk combos)
  - Sample overview (ST, K/O, gene counts)
  - MLST distribution table
  - Kaptive tables (K, O, K:O combinations)
  - AMR gene table (gene, category, database, risk, frequency, scrollable genome list)
  - Virulence gene table
  - Plasmid replicon table
  - Pattern discovery (ST‑K, ST‑O, K:O associations; high‑risk combinations; ICEKp and virulence plasmid markers)
  - Environmental markers (BACMET2 categories)
  - AI integration guide
  - Export buttons (CSV, JSON)

---

## 📈 **Performance**

| System | Genomes | Time |
|--------|---------|------|
| Laptop (2 cores, 4 GB) | 30 | 5h 16m |
| Workstation (16 cores, 16 GB) | 30 | ~3h (estimated) |

- **QC + MLST + Kaptive** (parallel): ~4 min
- **ABRicate** (11 DBs): ~1h
- **AMRfinder**: ~3h
- **Ultimate reporter**: <1 sec

Memory usage never exceeded 3 GB.

---

## 🔬 **Validation**

Kleboscope was validated on 30 publicly available *K. pneumoniae* genomes. **All results were in complete concordance with PubMLST, CGE, and Kaptive reference databases.**

### **Validation Highlights**
- **MLST**: 10 distinct STs, dominated by ST11 (50.0%)
- **Capsule types**: 11 K loci, 4 O loci; KL64:OL2α.1 most common (33.3%)
- **Critical resistance genes**: *blaKPC‑2* (56.7%), *blaNDM‑1* (10.0%), *blaOXA‑232* (10.0%)
- **Hypervirulence markers**: *iuc* (70%), *rmpA2* (70%), *ybt* (80%), *clb* (13.3% – ST23/ST2096)
- **Environmental markers**: *qacEdelta1* (56.7%), *sil* (50%), *mer* (70%)

For full validation data, see the interactive report (`kleboscope_ultimate_report.html`).

---

## 🔄 **Alternative Tools**

If Kleboscope does not fit your workflow, you may consider:

- **[Kleborate](https://github.com/katholt/Kleborate)** – the established tool for *K. pneumoniae* typing, resistance and virulence scoring (tabular output).  
- **[Bactopia](https://bactopia.github.io/)** – a flexible multi‑species pipeline for raw reads (requires Nextflow).  
- **[Nullarbor](https://github.com/tseemann/nullarbor)** – another multi‑species pipeline for raw reads.

Kleboscope complements these tools by offering a **gene‑centric interactive HTML report** with environmental markers, plasmid typing, and pattern discovery, all in a single, easy‑to‑use command.

---

## 🤖 **AI Integration Guide**

Kleboscope’s ultimate HTML report is designed to be AI‑friendly. Use any large language model (ChatGPT, Claude, Gemini) to gain deeper insights.

### **Quick Start**
1. **Open** `kleboscope_ultimate_report.html` in your browser
2. **Copy** any table or section
3. **Paste** into your AI chat
4. **Ask** questions like:

**For MLST and capsule types:**
- “Which sequence types are most common in this dataset?”
- “Show me all ST11 isolates and their capsule types.”
- “Are there any ST‑K associations that suggest high‑risk clones?”

**For resistance genes:**
- “Which samples carry carbapenemases? List them.”
- “What is the prevalence of *blaKPC‑2* and which STs does it associate with?”
- “Are there any isolates with both carbapenemase and colistin resistance genes?”

**For virulence markers:**
- “Which samples carry the aerobactin (*iuc*) system?”
- “Show me all hypervirulent candidates (ST23, *rmpA*, *iuc*).”
- “What are the clinical implications of a strain carrying both *rmpA* and *blaKPC‑2*?”

**For environmental markers:**
- “Which samples have *qacEΔ1*? What does this mean for disinfectant tolerance?”
- “List all isolates with silver (*sil*) or mercury (*mer*) resistance genes.”

**For pattern discovery:**
- “Identify high‑risk combinations (critical resistance + hypervirulence).”
- “What capsule types are associated with ST11? With ST23?”

### **Pro Tips**
- **Provide context**: “I’m analysing 30 *K. pneumoniae* genomes. Here is the gene frequency table…”
- **Ask for summaries**: “Summarise the resistance profile of this outbreak.”
- **Combine tables**: Copy the AMR table and the virulence table together to find correlations.

> *“AI accelerates pattern discovery, but always verify critical findings with domain experts.”*

---

## 🔮 **Future Directions**

We are actively developing Kleboscope and welcome community contributions. Planned enhancements include:

- **Raw read support** (FASTQ) with integrated assembly (Shovill) to eliminate the need for pre‑assembled genomes.
- **Web interface** similar to StaphScope Web, allowing non‑bioinformaticians to run the pipeline in a browser.
- **Real‑time database updates** for ABRicate and AMRfinder.
- **Machine learning module** for phenotype prediction and outbreak risk scoring.
- **Expanded environmental marker database** (e.g., disinfectant residues, heavy metal contamination).
- **Plugin system** for community‑contributed analysis modules.
- **Integration with public health databases** (e.g., NCBI Pathogen Detection) for large‑scale surveillance.

We invite collaboration on these fronts – see [Authors & Contact](#-authors--contact).

---

## ❓ **Frequently Asked Questions**

### **General Questions**

**Q: Is Kleboscope free?**  
A: Yes! Kleboscope is open‑source under the MIT license. Free for academic, clinical, and commercial use.

**Q: How is Kleboscope different from Kleborate?**  
A: Kleboscope offers:
- Gene‑centric, interactive HTML report (not just tabular)
- Environmental co‑selection markers (BACMET2)
- Plasmid replicon typing (PlasmidFinder)
- Pattern discovery (ST‑capsule associations, high‑risk combos, ICEKp/virulence plasmid tracking)
- AI integration guide
- Scrollable genome lists with no truncation

**Q: Can I use Kleboscope for clinical diagnostics?**  
A: Kleboscope is a research tool. While highly accurate, results should be validated with orthogonal methods for clinical decision‑making.

### **Technical Questions**

**Q: Why only assembled genomes?**  
A: The pipeline is optimised for assembled genomes, which are commonly available from public databases and clinical labs. Raw read support is planned.

**Q: How do I update databases?**  
A: Run `abricate --setupdb` for ABRicate databases. For AMRfinder, the bundled database is updated with each release. We will provide updates regularly.

**Q: Can I run Kleboscope on Windows?**  
A: Yes, via WSL2 (Windows Subsystem for Linux). Native Windows support is planned.

**Q: How do I handle very large datasets (1000+ genomes)?**  
A: Use the CLI with glob patterns; the pipeline scales linearly. On a cluster, you can increase `--threads` (AMRfinder uses multiple cores). Consider running modules separately with skip flags if needed.

### **Analysis Questions**

**Q: What does “Capsule null” mean in Kaptive results?**  
A: This indicates the assembly does not contain a complete K locus, often because of fragmentation or because the strain lacks the typical capsule. The O locus may still be present.

**Q: How is risk level assigned?**  
A: Risk levels are based on clinical relevance:
- **Critical**: carbapenemases, mcr, tetX, 16S rRNA methyltransferases, etc.
- **High**: ESBLs, colistin (point mutations), aminoglycoside resistance
- **Standard**: other resistance genes

**Q: Are virulence factors from other species filtered out?**  
A: Yes. The ultimate reporter uses a curated list of *K. pneumoniae*‑relevant genes. ABRicate databases are generic, but the report categories focus on known markers.

---

## 📚 **Citation & Acknowledgements**

### **Citing Kleboscope**
If you use Kleboscope in your research, please cite:

> Beckley,B. _et. al_, (2026). Kleboscope: A gene‑centric, species‑optimized computational pipeline for comprehensive *Klebsiella pneumoniae* genomic surveillance. *Nature Com.* (under preparation).

**Software citation:**
```bibtex
@software{kleboscope2026,
  author = {Brown Beckley and Vincent Amarh},
  title = {Kleboscope: A gene‑centric, species‑optimized computational pipeline for comprehensive Klebsiella pneumoniae genomic surveillance},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/Kleboscope}
}
```

### **Citing the Integrated Tools & Databases**

Kleboscope stands on the shoulders of many outstanding open‑source projects. When using Kleboscope, please also cite the tools and databases that make it possible:

#### **Core Tools**

| Tool | Citation |
|------|----------|
| **MLST** | Seemann T. (2018). mlst. GitHub. https://github.com/tseemann/mlst |
| **Kaptive** | Wick RR, Heinz E, Holt KE, Wyres KL. (2018). Kaptive web: user‑friendly capsule and lipopolysaccharide serotype prediction for *Klebsiella* genomes. *J Clin Microbiol*, 56(6):e00197-18. |
| **ABRicate** | Seemann T. (2018). ABRicate. GitHub. https://github.com/tseemann/abricate |
| **AMRfinderPlus** | Feldgarden M, et al. (2019). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep*, 11:12728. |

#### **Databases**

| Database | Citation |
|----------|----------|
| **PubMLST** | Jolley KA, Bray JE, Maiden MCJ. (2018). Open‑access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. *Wellcome Open Res*, 3:124. |
| **CARD** | Alcock BP, et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res*, 51(D1):D690–D699. |
| **ResFinder** | Bortolaia V, et al. (2020). ResFinder 4.0 for predictions of phenotypes from genotypes. *J Antimicrob Chemother*, 75(12):3491–3500. |
| **ARG‑ANNOT** | Gupta SK, et al. (2014). ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. *Antimicrob Agents Chemother*, 58(1):212–220. |
| **VFDB** | Chen L, et al. (2016). VFDB 2016: hierarchical and refined dataset for big data analysis—10 years on. *Nucleic Acids Res*, 44(D1):D694–D697. |
| **PlasmidFinder** | Carattoli A, et al. (2014). *In silico* detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother*, 58(7):3895–3903. |
| **MegaRes** | Bonin N, et al. (2023). MEGARes and AMR++, v3.0: an updated comprehensive database of antimicrobial resistance determinants and an improved software pipeline for classification using high‑throughput sequencing. *Nucleic Acids Res*, 51(D1):D744–D752. |
| **BacMet2** | Pal C, et al. (2014). BacMet: antibacterial biocide and metal resistance genes database. *Nucleic Acids Res*, 42(D1):D737–D743. |
| **EcoH / EcoLI_VF** | Joensen KG, et al. (2015). Rapid and easy *in silico* serotyping of *Escherichia coli* isolates by use of whole‑genome sequencing data. *J Clin Microbiol*, 53(8):2410–2426. |

### **Acknowledgements**

We thank the developers of all the tools and databases that Kleboscope integrates, and the open‑source community for their invaluable contributions. Special thanks to Torsten Seemann (MLST, ABRicate), the NCBI AMR team, the CGE group (PlasmidFinder, ResFinder), and the Kaptive team for their foundational work.

---

## 👥 **Authors & Contact**

**Brown Beckley** (Primary Developer)  
- University of Ghana Medical School  
- 📧 brownbeckley94@gmail.com  
- 🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)  
- LinkedIn: [@brownbeckley](https://www.linkedin.com/in/brown-beckley-190315319/)

**Vincent Amarh** (Co‑Author)  
- University of Ghana Medical School

### **Collaboration Opportunities**
We welcome collaborations on:
- *K. pneumoniae* epidemiology and outbreak studies
- Clinical validation of resistance/virulence markers
- Expanding the environmental marker database
- Development of a web interface
- Integration with public health surveillance systems
- Machine learning applications for phenotype prediction

---

## 📄 **License & Third‑Party Notices**

### **Core Kleboscope Code**
The Kleboscope pipeline code (workflow engine, report generation, HTML templates, Python modules) is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

### **Third‑Party Tools & Databases**
Kleboscope bundles several external tools and databases under their own licenses. By using Kleboscope, you agree to comply with the respective licenses of these components.

| Component | License |
|-----------|---------|
| MLST (tseemann) | GPL v2 |
| Kaptive | GPL v3 |
| ABRicate (tseemann) | GPL v2 |
| AMRfinderPlus (NCBI) | Public Domain |
| PlasmidFinder (CGE) | Free for academic use |
| ResFinder (CGE) | Free for academic use |
| CARD | CC BY 4.0 |
| VFDB | Open access |
| BacMet2 | Open access |
| PubMLST | Open access for research |

Full license texts and attribution details are available in the respective tool repositories and websites linked above.

---

<div align="center">

## **🚀 Ready to transform your *K. pneumoniae* surveillance?**

**From FASTA to actionable insights in one command.**

[![Get Started](https://img.shields.io/badge/GET_STARTED-Now-green?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/kleboscope#-quick-start)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/kleboscope/issues)

*Join the Fight Against Antimicrobial Resistance*

Antimicrobial resistance threatens modern medicine. We invite researchers, clinicians, and public health professionals to collaborate with us in expanding and validating our database, sharing regional epidemiological data, and advancing AMR surveillance.

**Together, we can enhance global AMR monitoring and develop more effective treatment strategies.**

</div>
