<p align="center">
  <img src="https://raw.githubusercontent.com/bbeckley-hub/Kleboscope/main/kleboscope.png" alt="Kleboscope Banner" width="100%">
</p>

<div align="center">

# 🧬 Kleboscope

### **A gene-centric, species-optimized computational pipeline for comprehensive *Klebsiella pneumoniae* genomic surveillance**

#### **Complete K. pneumoniae typing, resistance, virulence, and mobile genetic element analysis in minutes, not hours**

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
![Version](https://img.shields.io/badge/version-1.3.0-blue)
![Platforms](https://anaconda.org/bbeckley-hub/kleboscope/badges/platforms.svg)
![License](https://anaconda.org/bbeckley-hub/kleboscope/badges/license.svg)

[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/kleboscope)](https://hub.docker.com/r/bbeckleyhub/kleboscope)
[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/kleboscope/latest)](https://hub.docker.com/r/bbeckleyhub/kleboscope)
[![Docker Version](https://img.shields.io/docker/v/bbeckleyhub/kleboscope?sort=semver)](https://hub.docker.com/r/bbeckleyhub/kleboscope)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/kleboscope?label=Conda%20Downloads)](https://anaconda.org/bbeckley-hub/kleboscope)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/brown-beckley-190315319)
[![Stage](https://img.shields.io/badge/status-active-brightgreen)](#)

[![Powered by 🧠](https://img.shields.io/badge/powered%20by-science%20🔬-purple)](https://github.com/bbeckley-hub/Kleboscope)
[![Coffee](https://img.shields.io/badge/built%20with-%E2%98%95%20coffee-orange)](https://github.com/bbeckley-hub/Kleboscope)
[![Made with ❤️](https://img.shields.io/badge/made%20with-%E2%9D%A4%EF%B8%8F-red)](https://github.com/bbeckley-hub/Kleboscope)
[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)
[![Made for Research](https://img.shields.io/badge/made%20for-Research-0066cc.svg)](https://github.com/bbeckley-hub/Kleboscope)

[![Last Commit](https://img.shields.io/github/last-commit/bbeckley-hub/Kleboscope)](https://github.com/bbeckley-hub/Kleboscope/commits)
[![Contributors](https://img.shields.io/github/contributors/bbeckley-hub/Kleboscope)](https://github.com/bbeckley-hub/Kleboscope/graphs/contributors)
[![Security: bandit](https://img.shields.io/badge/security-bandit-yellow.svg)](https://github.com/PyCQA/bandit)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![CI](https://img.shields.io/github/actions/workflow/status/bbeckley-hub/Kleboscope/ci.yml?branch=main&label=CI)](https://github.com/bbeckley-hub/Kleboscope/actions)
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/bbeckley-hub/Kleboscope)

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/badge/issues-Open-orange)](https://github.com/bbeckley-hub/Kleboscope/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/Kleboscope)](https://github.com/bbeckley-hub/Kleboscope/stargazers)
![Profile Views](https://komarev.com/ghpvc/?username=bbeckley-hub&label=Profile%20Views&color=0e75b6&style=flat)
[![Google Scholar](https://img.shields.io/badge/Google%20Scholar-Profile-4285F4?style=flat&logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=CYNOsqIAAAAJ&hl=en)

![GitHub stats](https://github-readme-stats.vercel.app/api?username=bbeckley-hub&show_icons=true&theme=radical)
![Top Langs](https://github-readme-stats.vercel.app/api/top-langs/?username=bbeckley-hub&layout=compact&theme=radical)
[![GitHub Streak](https://streak-stats.demolab.com?user=bbeckley-hub&theme=radical&date_format=j%20M%5B%20Y%5D)](https://git.io/streak-stats)

**One command in:** 🖥️ **Command-line tool** for high-throughput, local, privacy-preserving analysis of assembl  **Two actionable reports out:** 🧬 gene-centric for epidemiology, 📋 sample-centric for the clinic

</div>

---

## 🎉 **What's New in v1.3.0 — MGE Profiling, Compare Tab, and a Hardened Orchestrator**

Version 1.3.0 is the biggest structural release of Kleboscope. The pipeline has been re-architected end to end, a new analysis dimension has been added, and the reports now answer a question the old reports could not: **how did these genes get here?**

### 🦠 **1. New MGE Module: Mobile Genetic Element Profiling**

By default, Docker runs containers as `root`, so any files written to bind-mounted directories will be owned by `root:root`, producing padlock icons and requiring `sudo chown`. **The fix is simple:** add `-u $(id -u):$(id -g)` to run the container with your host user's UID/GID. But first, the headline:

Kleboscope now integrates **mobileOG-db** to profile the mobile genetic elements that shuttle resistance and virulence between genomes.

- **Per-sample counts** across ten functional categories: Integrase, Transfer, Stability, Phage, Replication, IS-associated, ICE-associated, Plasmid-associated, Phage-associated, and Key MGE signatures.
- **Aggregate profiles** by sequence type and by K and O locus. A heavy MGE burden in one lineage but not another is a real signal.
- **Batch summaries** in HTML, TSV, and JSON, fully wired into both ultimate reporters.
- A new **`--pull-mge-db`** installer for conda users, and the Docker image ships with the complete mobileOG payload baked in, so container users never download a thing.

> **Why it matters:** you have AMR genes and virulence genes. MGE profiling tells you the postal system that delivered them. A plasmid-associated gene cluster appearing in two unrelated STs is a smoking gun for horizontal transfer.

### ⚖️ **2. The Compare Tab Comes to Kleboscope**

The feature StaphScope users asked for first has landed in the **gene-centric reporter**. Pick any two isolates and get:

1. **A verdict banner**: 🔴 near-identical (≥95% similarity), 🟠 high (≥85%), 🟡 moderate (≥60%), 🔵 distinct
2. **A similarity gauge**, an animated colored ring showing the overall match percentage
3. **A typing table** with every field side by side, green pill for match, red pill for difference
4. **Gene content split three ways**: 🟢 shared, 🔵 only in sample A, 🟣 only in sample B, across AMR, virulence, BACMET, plasmid, and mutation categories

Two isolates from the same ward with identical MLST, K-locus, and O-locus are very likely a transmission pair. The Compare tab turns "I think they're the same bug" into a structured, defensible observation. Click any two isolates in the gene-centric report and judge for yourself.

### ⚙️ **3. The Orchestrator, Rebuilt From the Ground Up**

- **Fully sequential, isolated execution**: every module is copied into its own temporary directory and runs there, so module-relative paths (db/, bin/, scripts/) always resolve correctly, in conda, on HPC, and in containers.
- **Robust database resolution chains** with environment overrides (`KLEBOSCOPE_MLST_DB`, `KLEBOSCOPE_MGE_DATA_DIR`, and friends) so the same code runs on a laptop and on a locked-down cluster node.
- **`--report-level` flag** (`final`, `summary`, `full`) to control which files are kept after a run and prevent disk exhaustion on 300+ sample cohorts.
- **`--tmpdir` and `KLEBOSCOPE_TMPDIR`** for HPC scratch redirection, with signal-handled and `atexit` cleanup.
- **`--quiet` mode**, per-module timing, a warnings and errors recap at the end of every run, and proper exit codes for pipeline integration.
- **Streamlined MLST reports**: sequence type, allele profile, identity, and coverage. The clonal complex annotation columns have been retired; the ST tells the story.

### 🔬 **4. MLST Database Management: PubMLST, Your Way**

```bash
kleboscope --pull-mlst-db      # install the Klebsiella scheme, no API key needed
kleboscope --update-mlst-db    # refresh directly from PubMLST, API key authenticated
```

The direct-from-source refresh uses **mlstdb** and the official PubMLST API, so your scheme is always current:

1. Register at [https://pubmlst.org/site-accounts](https://pubmlst.org/site-accounts)
2. Generate an API key in your profile
3. Run: `mlstdb connect --db pubmlst --api-key`
4. Paste your key when prompted

### 🎨 **5. Reports With One Visual Language**

All Kleboscope HTML reports now share a **consistent blue theme**, so a cohort report, an MLST summary, and an MGE profile look and feel like one product instead of five tools taped together. Interactive tables, search, export buttons, and click-to-expand gene lists behave identically everywhere.

### 🐳 **6. A Friendlier Docker Image**

The image has been rebuilt around the same philosophy as the rest of the ESKAPE suite: **pull and run, no padlock**. Every database (MLST scheme with BLAST index, mobileOG-db, AMRFinderPlus, ABRicate libraries) and every binary is baked in at build time, the container runs as a non-root user, and output files on mounted volumes are owned by you. Details in the [Docker & Singularity section](#-docker--singularity-usage--avoid-the-padlock-).

### 🔮 **Coming Soon: KLEOS, the Klebsiella Retrieval Database**

We are building **KLEOS**, a free, interactive, browser-based resource for typing-based genome retrieval, in the spirit of [SAURUS](https://bbeckley-hub.github.io/saurus/). The goal: put the public *K. pneumoniae* repository behind real-time filters for MLST, K and O locus, AMR genes, virulence markers, geography, and isolation source, then export a plain-text accession list ready for `datasets download genome accession --inputfile list.txt`. **Click. Filter. Export.** Watch the releases page.

---

## 🎯 **Overview**

**Kleboscope** is an automated, locally-executable computational pipeline designed specifically for comprehensive *Klebsiella pneumoniae* genomic surveillance. It addresses the growing threat of multidrug-resistant and hypervirulent *K. pneumoniae* by integrating **eight essential analysis modules** into a single, cohesive workflow, and consolidates the results into **two complementary report views**:

- **Gene-centric** – each gene is shown with **all genomes** that contain it, together with its frequency (`count (percentage%)`), enabling rapid cross-genome pattern discovery.
- **Sample-centric** – each isolate gets its own **interactive box** with typing badges (MLST, K-locus, O-locus, hypervirulence), per-database tables (AMR, Virulence, BACMET, Plasmids), and full mutation details. Perfect for clinical reports and patient-level investigations.

### 🌍 **The Problem**
- **Fragmented bioinformatics**: *K. pneumoniae* typing requires separate tools for MLST, capsule typing, AMR detection, virulence screening, plasmid profiling, MGE profiling, and quality control.
- **Time-consuming workflows**: combining outputs from multiple tools is error-prone and delays outbreak responses.
- **Interpretation gap**: raw data lacks clinical context; identifying high-risk clones (e.g., ST11-KL64 with *blaKPC-2* and *iuc*) requires manual correlation.
- **Missing environmental markers**: existing tools ignore biocide and heavy metal resistance genes that co-select for antibiotic resistance.

### 💡 **Our Solution**
Kleboscope delivers:
- ✅ **Single-command installation** via Conda or Docker.
- ✅ **MGE profiling** with mobileOG-db, the missing dimension of resistance dissemination.
- ✅ **Compare tab** for pairwise isolate comparison and cluster-level similarity heatmaps.
- ✅ **Two report types**: gene-centric for epidemiology, sample-centric for clinical decision-making.
- ✅ **100% local execution** with data privacy.
- ✅ **Tracking of critical resistance genes**: carbapenemases (*blaKPC*, *blaNDM*, *blaOXA-48*), colistin (*mcr*), tigecycline (*tetX*), 16S rRNA methyltransferases.
- ✅ **Tracking of hypervirulence markers**: aerobactin (*iuc*), salmochelin (*iro*), yersiniabactin (*ybt*), colibactin (*clb*), regulators of hypermucoidy (*rmpA*, *rmpA2*).
- ✅ **Environmental co-selection markers**: biocide resistance (*qac*), heavy metal resistance (*sil*, *mer*, *ars*, *pco*, *czc*), mobile genetic elements, stress response genes.
- ✅ **Built-in pattern discovery**: ST-capsule associations, high-risk combinations, ICEKp and virulence plasmid tracking.
- ✅ **AI integration guide** with prompts to help you mine your data with large language models.
- ✅ **HPC-friendly orchestrator**: isolated temp directories, comprehensive logging, modular skipping, report-level pruning, signal handling, and automatic cleanup.

**Perfect for**: clinical laboratories, outbreak investigations, research studies, and public health surveillance.

---

## ✨ **Key Features**

### 🔬 **Core Analytical Modules**

| Module | Purpose | Key Outputs | Speed* |
|--------|---------|-------------|--------|
| **FASTA QC** | Assembly statistics, N50, GC%, contig quality, species confirmation | HTML, TSV, JSON | <30 sec |
| **MLST Typing** | Sequence type (Pasteur scheme) | ST, allele profile, identity, coverage | <1 min |
| **Kaptive Capsule Typing** | K (capsule) and O (lipopolysaccharide) loci | K/O types, identity, coverage, confidence | 1-2 min |
| **ABRicate Screening** | 11 databases (CARD, ResFinder, VFDB, PlasmidFinder, BacMet2, and more) | Gene-centric tables per database | 2-3 min per DB |
| **AMRFinderPlus** | Acquired resistance genes and point mutations | Gene frequency, risk levels | 3-4 min |
| **MGE Profiling** | mobileOG-db signатure counts across functional categories | Per-sample MGE profiles | 2-3 min |
| **Gene-centric Reporter** | Integrated cross-genome HTML report with Compare tab | Interactive tables, patterns, similarities | <1 sec |
| **Sample-centric Reporter** | Per-isolate interactive HTML report | Isolate boxes, typing badges, mutations | <1 sec |

*Indicative timings on a laptop; parallel batch execution for cohorts.

### 🛡️ **Species-Specific Innovations for *K. pneumoniae***

- **Comprehensive resistance tracking**:
  - Carbapenemases: *blaKPC*, *blaNDM*, *blaIMP*, *blaVIM*, *blaOXA-48*, *blaGES*, *blaIMI*, *blaSME*
  - Colistin resistance: *mcr-1* to *mcr-10*, *pmrAB*, *lpx*, *arn*
  - Tigecycline resistance: *tet(X)* variants, efflux pumps (*acrAB*)
  - 16S rRNA methyltransferases: *armA*, *rmtA* to *rmtH*, *npmA*
  - β-lactamases: all major ESBLs, AmpC, and narrow-spectrum enzymes

- **Hypervirulence marker tracking**:
  - Aerobactin (*iucA-D*, *iutA*)
  - Salmochelin (*iroB-E*, *iroN*)
  - Yersiniabactin (*ybtA-X*, *irp1-2*)
  - Colibactin (*clbA-S*)
  - Hypermucoidy regulators (*rmpA*, *rmpA2*, *rmpC*, *rmpD*)
  - Additional markers: *peg-344*, *allS*, *kfu*, *mrkD*

- **Environmental co-selection markers** (BACMET2):
  - Biocide resistance: *qac*, *cep*, *oqx*
  - Heavy metal resistance: *sil*, *mer*, *ars*, *pco*, *czc*, *cop*, *cad*, *chr*
  - Mobile genetic elements: *tra*, *mob*, *rep*, *int*, *tnp*

- **Two complementary report views**:
  - **Gene-centric**: every gene with all genomes that carry it, frequency as `count (percentage%)`, scrollable genome lists, and the Compare tab for pairwise analysis.
  - **Sample-centric**: every isolate in its own box with typing badges, per-database tables, and mutation details.

---

## 🔥 **Feature Deep-Dives**

<details>
<summary><b>⚖️ The Compare Tab — Because "I Think They're the Same Bug" Isn't Good Enough (Click to expand)</b></summary>
<br>

Let's be honest about how outbreak investigations often go. Two genomes, two reports, two browser windows, and a mental diff of MLST, K-locus, O-locus, AMR genes, virulence genes, plasmids, and mutations, ending in *"yeah, I think they're the same."*

**The Compare tab replaces the eyeball.** Pick any two isolates in the gene-centric report, click compare, and get a verdict:

1. **A verdict banner**: 🔴 near-identical (≥95% similarity), 🟠 high (≥85%), 🟡 moderate (≥60%), 🔵 distinct
2. **A similarity gauge**: an animated colored ring showing the overall match percentage
3. **A typing table**: every field side by side, green pill for match, red pill for difference
4. **Gene content split three ways**: 🟢 shared, 🔵 only in A, 🟣 only in B, per category (AMR, Virulence, BACMET, Plasmids, Mutations)
5. **Full profile expanders**: click any category to see both samples' complete gene lists side by side

**Why this matters:**
- **Outbreak confirmation**: identical MLST + K-locus + O-locus across isolates from the same ward is a strong transmission signal your infection control team can act on.
- **Discordance analysis**: two isolates identical in every field except one, and that one field tells you which mobile element moved.
- **QC checks**: two genomes from the same sequencing run that come back identical are probably duplicates. Catch that now, not at revision time.
- **Teaching**: explaining "why isn't this hypervirulent?" to a clinician takes 30 seconds with a Compare tab.

**A caveat we state openly:** identical typing suggests a transmission pair, it does not prove one. These markers evolve slowly, and the gold standard remains core-genome SNP distance, which Kleboscope does not compute at this scope. We say *"possible transmission pair"* deliberately. If we ever lose that hedge, shout at us.
</details>

<details>
<summary><b>🦠 The MGE Tab — Where Resistance Gets a Postal Address (Click to expand)</b></summary>
<br>

You have AMR genes. You have virulence genes. You know *what* a genome carries. But *how did it get there?*

**Mobile genetic elements are the answer.** Plasmids, prophages, transposons, integrative conjugative elements: these are the physical vehicles that shuttle resistance and virulence genes between unrelated bacteria. If you do not understand the MGE landscape, you are reading the mail without ever seeing the postal system.

**What the MGE analysis gives you:**
- **Per-sample counts** across ten functional categories (Integrase, Transfer, Stability, Phage, Replication, IS-associated, ICE-associated, Plasmid-associated, Phage-associated, Key MGE signatures).
- **Aggregate profiles by typing** — mean MGE load per ST, per K-locus, per O-locus. High MGE burden in one lineage but not another is a real signal.
- **HTML, TSV, and JSON summaries** that flow straight into the gene-centric and sample-centric reports.

**What it helps with:**
- **Resistance dissemination**: a plasmid-associated gene cluster in two unrelated STs is a smoking gun for horizontal transfer.
- **Virulence acquisition**: hypervirulence payloads often ride on virulence plasmids and ICEKp.
- **Outlier detection**: an isolate with twice the cohort's MGE load is either something new or something broken. Either way, look at it.

**A word of honesty:** mobileOG-db counts protein *signatures*, not reconstructed elements. Multiple proteins from the same element produce multiple hits, and some categories also fire on chromosomal genes. That is why we use the word *"signatures"* everywhere. For element-level claims you need co-localization analysis, and that is a study, not a dashboard column.
</details>

<details>
<summary><b>📊 The Report-Level System — Because 300 Sample Cohorts Should Not Fill Your Disk (Click to expand)</b></summary>
<br>

Every Kleboscope module writes per-sample outputs. On a single genome that is a perk. On a 300-genome surveillance cohort it is a disk-filling hazard, especially on a shared cluster where you get a polite but firm email from the administrator.

**`--report-level` gives you control over what survives the run:**

| Level | What is kept | Use case |
|---|---|---|
| `final` | Only `Kleboscope_final_report/` | Clinical report delivery, the output goes straight to a collaborator |
| `summary` (default) | Final report plus batch summary files; per-sample folders are pruned | Surveillance cohorts; the summaries contain everything the reports need |
| `full` | Everything, including every per-sample HTML | Re-analysis, audits, debugging |

**Why the default is `summary`:** the two ultimate reports are generated from the batch summary TSV, JSON, and HTML files. The per-sample folders are inputs to those summaries, not a separate deliverable. Pruning them cuts output size by an order of magnitude on large cohorts without losing a single result that the reports use.

Combined with `--tmpdir` (or `KLEBOSCOPE_TMPDIR`) for scratch redirection, Kleboscope can run a 300-genome cohort on a node with strict home quotas and leave behind only what you actually want.
</details>

<details>
<summary><b>🔬 MLST Database Management — No API Key Required to Start (Click to expand)</b></summary>
<br>

Kleboscope ships a straightforward path to the current PubMLST Pasteur scheme for *K. pneumoniae*:

```bash
kleboscope --pull-mlst-db      # install the Klebsiella scheme, no API key needed
kleboscope --update-mlst-db    # refresh directly from PubMLST, API key authenticated
```

`--pull-mlst-db` installs the scheme into your data cache and builds the BLAST allele index automatically. It works with the module's own installer when available and falls back to a prebuilt archive, so there is nothing to configure on a fresh install.

`--update-mlst-db` uses **mlstdb** and the official PubMLST API:

1. Register at [https://pubmlst.org/site-accounts](https://pubmlst.org/site-accounts)
2. Generate an API key in your profile
3. Run: `mlstdb connect --db pubmlst --api-key`
4. Paste your key when prompted

**Where the database lives:**
- Conda and source installs: `~/.local/share/kleboscope/mlst_module/`
- Docker and Singularity: baked into the image, resolution is automatic
- Any environment override: `KLEBOSCOPE_MLST_DB` and `KLEBOSCOPE_MLST_MODULE_PATH`

The orchestrator checks the resolution chain on every MLST run and tells you clearly if anything is missing, with the exact command to fix it.
</details>

---

## 🌐 **ESKAPE Ecosystem**

Kleboscope is a core member of a growing suite of species-optimized pipelines built on the same architecture:

| Pipeline | Species | Repository |
|----------|---------|------------|
| **StaphScope** | *Staphylococcus aureus* | [bbeckley-hub/staphscope-typing-tool](https://github.com/bbeckley-hub/staphscope-typing-tool) |
| **Kleboscope** | *Klebsiella pneumoniae* | [bbeckley-hub/Kleboscope](https://github.com/bbeckley-hub/Kleboscope) |
| **PseudoScope** | *Pseudomonas aeruginosa* | [bbeckley-hub/pseudoscope](https://github.com/bbeckley-hub/pseudoscope) |
| **AcinetoScope** | *Acinetobacter baumannii* | [bbeckley-hub/acinetoscope](https://github.com/bbeckley-hub/acinetoscope) |

Shared architecture means shared reliability: identical orchestrator design, identical database-management philosophy, identical report styling, and the same Compare tab across species. The web platform for the suite is the [ESKAPE AMR Platform](https://eskape.bio), and the *S. aureus* retrieval database [SAURUS](https://bbeckley-hub.github.io/saurus/) previews what KLEOS will offer for *K. pneumoniae*.

---

## ⚡ **Quick Start (CLI)**

### Install in 60 seconds
```bash
# Method 1: Conda (Recommended)
conda create -n kleboscope -c conda-forge -c bioconda kleboscope -y
conda activate kleboscope

# Method 2: From source
git clone https://github.com/bbeckley-hub/Kleboscope.git
cd Kleboscope
conda env create -f environment.yml
conda activate kleboscope
pip install -e .

# Method 3: Docker
docker pull bbeckleyhub/kleboscope:1.3.0
```

### Run your first analysis
```bash
# Single genome, both reports
kleboscope -i genome.fna -o results/

# Batch processing
kleboscope -i "*.fna" -o batch_results --threads 8

# Just the final deliverable reports
kleboscope -i "*.fna" -o results --report-level final
```

### Update Databases (Conda and source users)
```bash
# ABRicate databases
abricate --setupdb

# AMR database
kleboscope --update-amr-db

# MLST database
kleboscope --pull-mlst-db

# mobileOG-db for MGE profiling (~300 MB)
kleboscope --pull-mge-db
```

> **Docker and Singularity users:** skip all of the above. Every database is baked into the image. Pull and run.

---

## 🐳 **Docker & Singularity Usage – avoid the padlock 🔓**

Docker and Singularity users get **everything pre-installed**. No first-run database setup, no `--pull-mlst-db`, no `--pull-mge-db`, no `--update-amr-db`, no `abricate --setupdb`. Every database (MLST scheme with its BLAST index, mobileOG-db, AMRFinderPlus, ABRicate DBs) and every tool (BLAST, DIAMOND, fastANI, ABRicate, AMRFinderPlus) is baked into the image. **Pull and run.**

### 🐳 **Docker (laptops, workstations, cloud VMs)**

Docker runs containers as `root` by default, so files written to bind-mounted directories come out owned by `root:root`. **Fix: add `-u $(id -u):$(id -g)`.**

```bash
# Pull the image
docker pull bbeckleyhub/kleboscope:1.3.0

# Run (no padlock, no sudo chown)
docker run --rm \
  -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  bbeckleyhub/kleboscope:1.3.0 \
  -i "/data/*.fna" -o /data/output -t 4
```

- `-u $(id -u):$(id -g)` → files owned by **you**
- `-v $(pwd):/data` → mounts your current dir as `/data` inside the container
- `-i "/data/*.fna"` → uses the **container path**, not your host path
- `-o /data/output` → lands at `./output/` on your host

<details>
<summary><b>Verify the image and open a debugging shell (Click to expand)</b></summary>

```bash
# Verify the container runs
docker run --rm bbeckleyhub/kleboscope:1.3.0 --help

# Confirm the databases are baked in (no download happens)
docker run --rm bbeckleyhub/kleboscope:1.3.0 --pull-mlst-db
docker run --rm bbeckleyhub/kleboscope:1.3.0 --pull-mge-db

# Interactive shell (for debugging)
docker run --rm -it bbeckleyhub/kleboscope:1.3.0 bash
```
</details>

---

### 🔬 **Singularity / Apptainer (HPC clusters) — recommended for HPC**

**HPC users should use Singularity, not Docker.** HPC login and compute nodes do not run the Docker daemon; Singularity (or its successor Apptainer) is the standard. Same image format, rootless by default, and no UID juggling: it runs as your real user automatically.

<details>
<summary><b>1. Load the runtime (Click to expand)</b></summary>

Module names vary by cluster. Check with `module avail`:

```bash
# Newer clusters — Apptainer
module load apptainer

# Older clusters — Singularity
module load containers/singularity/3.6.4

# Pulling the image also needs mksquashfs from squashfs-tools
module load squashfs-tools/4.4

singularity --version   # or: apptainer --version
```
</details>

<details>
<summary><b>2. Pull the image (one-time, needs internet on the login node) (Click to expand)</b></summary>

```bash
singularity pull kleboscope_1.3.0.sif docker://bbeckleyhub/kleboscope:1.3.0
```

Creates `kleboscope_1.3.0.sif` in the current directory. **`squashfs-tools` is only needed for the pull, not for running.**
</details>

<details>
<summary><b>3. Run the pipeline (Click to expand)</b></summary>

```bash
singularity run --bind $(pwd):/data \
  kleboscope_1.3.0.sif \
  -i "/data/*.fna" -o /data/output -t 16
```

- `--bind $(pwd):/data` — the HPC equivalent of Docker's `-v`
- **No `-u` needed** — Singularity runs as your real UID
- **No `-e HOME=` needed** — Singularity uses your real `$HOME`, so MLST and AMR caches persist between runs

> ⚠️ **Use `singularity run`, not `singularity exec`.** The image's entrypoint sets up `PATH` and `LD_LIBRARY_PATH` for the conda environment. `run` invokes it, `exec` skips it, which will break AMRFinderPlus.
</details>

<details>
<summary><b>4. SLURM job script (Click to expand)</b></summary>

```bash
#!/bin/bash
#SBATCH --job-name=kleboscope
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=kleboscope_%j.log

module load containers/singularity/3.6.4

cd $SLURM_SUBMIT_DIR

singularity run --bind $(pwd):/data \
  /path/to/kleboscope_1.3.0.sif \
  -i "/data/*.fna" \
  -o /data/output \
  -t ${SLURM_CPUS_PER_TASK}
```
</details>

<details>
<summary><b>5. HPC notes (Click to expand)</b></summary>

- **No `module load` needed at runtime** — only the login-node pull requires `squashfs-tools`.
- **Pull once, run everywhere** — the SIF is a single self-contained file. Share it on your scratch filesystem and every compute node can use it.
- **Small home quota?** Redirect the Singularity build cache to scratch:
  ```bash
  export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
  ```
- **Conda-install MLST / MGE databases** live under `$HOME/.local/share/kleboscope/` by default. On tight home quotas:
  ```bash
  export KLEBOSCOPE_SHARED_DATA=/scratch/$USER/kleboscope
  ```
- **Scratch-bound temporary directories** for large cohorts:
  ```bash
  export KLEBOSCOPE_TMPDIR=/scratch/$USER/kleboscope_tmp
  ```
</details>

<details>
<summary><b>6. Troubleshooting (Click to expand)</b></summary>

| Error | Fix |
|---|---|
| `mksquashfs: executable file not found` | `module load squashfs-tools/4.4` before `singularity pull` |
| `manifest unknown` | The Docker tag has not been pushed, or the name is wrong. Check https://hub.docker.com/r/bbeckleyhub/kleboscope/tags |
| `singularity: command not found` | `module avail 2>&1 \| grep -i singularity` to find the module name |
| `BUNDLED AMRFinderPlus not properly installed` | You used `singularity exec`. Use `singularity run` instead |
| `mobileOG-db not found` | Conda users: run `kleboscope --pull-mge-db`. Container users: the DB is baked in, so this means a broken image |
| `No populated MLST database found` | Conda users: run `kleboscope --pull-mlst-db` |
</details>

---

### 🤔 **Which should you use?**

| Platform | Recommendation | Command pattern | Output ownership |
|---|---|---|---|
| **Laptop / workstation** | Docker | `docker run --rm -u $(id -u):$(id -g) -v "$PWD:/data" …` | You |
| **HPC login + compute nodes** | **Singularity (or Apptainer)** | `singularity run --bind "$PWD:/data" kleboscope.sif …` | You (automatic) |
| **Cloud VM** | Either. Docker is simpler | see above | You |
| **Shared cluster without Docker** | Singularity | see above | You (automatic) |

**The rule:** Docker on machines where you control the daemon; Singularity on shared HPC systems where you do not.

---

## 🚀 **Usage Guide**

### **Basic Commands**

```bash
# Full analysis, both reports
kleboscope -i "*.fna" -o results --threads 8

# Single genome
kleboscope -i genome.fna -o results

# Custom scratch directory (HPC, Docker)
kleboscope -i "*.fna" -o results --tmpdir /scratch/user/kleboscope_tmp

# Custom ABRicate thresholds
kleboscope -i "*.fna" -o results --abricate-min-id 85 --abricate-min-cov 90

# Quiet mode for pipeline integration
kleboscope -i "*.fna" -o results --quiet
```

### **Skip Options**

| Option | Effect |
|--------|--------|
| `--skip-qc` | Skip FASTA QC |
| `--skip-mlst` | Skip MLST typing |
| `--skip-kaptive` | Skip Kaptive capsule typing |
| `--skip-abricate` | Skip ABRicate screening |
| `--skip-amr` | Skip AMRFinderPlus analysis |
| `--skip-mge` | Skip MGE profiling |
| `--skip-gene-centric` | Skip the gene-centric reporter |
| `--skip-sample-centric` | Skip the sample-centric reporter |

### **Report Level**

| Option | Effect |
|--------|--------|
| `--report-level final` | Keep only `Kleboscope_final_report/` |
| `--report-level summary` | Final report plus batch summaries; per-sample folders pruned (default) |
| `--report-level full` | Keep every file, including per-sample HTMLs |

### **AMRFinderPlus Flags**

| Option | Effect |
|--------|--------|
| `--amr-min-identity` | Minimum identity for AMR hits (0..1) |
| `--amr-min-coverage` | Minimum coverage for AMR hits (0..1) |
| `--skip-amr-mutations` | Disable point mutation reporting in AMR |
| `--amr-force-update` | Force AMR database update before analysis |

### **ABRicate Flags**

| Option | Effect |
|--------|--------|
| `--abricate-min-id` | Minimum identity for ABRicate hits (default: 80) |
| `--abricate-min-cov` | Minimum coverage for ABRicate hits (default: 80) |

### **Database Management Commands**

| Command | Effect |
|---------|--------|
| `--update-amr-db` | Update the AMRFinderPlus database and exit |
| `--force-update-amr-db` | Force a complete AMRFinderPlus database update and exit |
| `--pull-mlst-db` | Install the Klebsiella MLST scheme and build its BLAST index, no API key needed |
| `--update-mlst-db` | Refresh the MLST scheme from PubMLST via mlstdb (API key required) |
| `--pull-mge-db` | Download mobileOG-db (~300 MB) for conda and source installs |

### **Real-World Examples**

**Clinical laboratory: routine surveillance of 50 isolates**
```bash
kleboscope -i "*.fna" -o weekly_surveillance --threads 8
# Both reports ready, batch summaries only
```

**Outbreak investigation: rapid typing, no screening panels**
```bash
kleboscope -i "outbreak/*.fasta" -o urgent --skip-abricate --skip-amr --skip-mge --report-level final
# MLST + Kaptive + QC against the outbreak isolates, final report only
```

**Clinical report: patient-level analysis**
```bash
kleboscope -i "clinical/*.fna" -o patient_reports --skip-gene-centric
# Sample-centric only, ideal for clinical rounds
```

---

## 📁 **Output Structure**

With the default `--report-level summary`:

```
results/
├── fasta_qc_results/                            # QC reports + batch summary
├── mlst_results/                                # MLST reports + batch summary
├── kaptive_results/                             # Kaptive reports + batch summary
├── abricate_results/                            # Per-database ABRicate summaries
├── amrfinder_results/                           # AMRFinderPlus reports + summaries
├── mge_results/                                 # mobileOG-db MGE profiles + summaries ✨ NEW
└── Kleboscope_final_report/
    ├── KLEBOSCOPE_ULTIMATE_GENE_CENTRIC_REPORTS/
    │   ├── kleboscope_ultimate_gene_centric_report.html
    │   ├── kleboscope_ultimate_gene_centric_report.json
    │   ├── sample_overview.csv
    │   ├── amr_genes.csv
    │   ├── virulence_genes.csv
    │   ├── bacmet_genes.csv
    │   ├── plasmid_replicons.csv
    │   ├── mutations.csv
    │   └── pattern_discovery.csv
    └── KLEBOSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS/
        ├── kleboscope_ultimate_sample_centric_report.html
        ├── kleboscope_ultimate_sample_centric_report.json
        ├── sample_overview.csv
        ├── amr_genes.csv
        ├── virulence_genes.csv
        ├── bacmet_genes.csv
        ├── plasmid_replicons.csv
        └── mutations.csv
```

- `--report-level final` leaves only `Kleboscope_final_report/`.
- `--report-level full` keeps all per-sample folders inside the module directories as well.

**Gene-centric report**: epidemiology, cross-genome comparisons, pattern discovery, and the Compare tab.
**Sample-centric report**: clinical reports, patient-level investigations, and infection control rounds.

---

## 🔍 **Analytical Modules**

### **1. FASTA QC**
- **Metrics**: total length, contig count, N50/N75/N90, L50/L75/L90, mean/median length, GC%, ambiguous bases, and species confirmation via fastANI
- **Output**: interactive HTML with warnings, TSV/JSON for downstream analysis
- **Gate**: samples failing the species check are excluded downstream and reported clearly, so a contaminant never pollutes your cohort statistics

### **2. MLST Typing**
- **Database**: [PubMLST](https://pubmlst.org/organisms/klebsiella-pneumoniae) (Pasteur scheme)
- **Tool**: [MLST](https://github.com/tseemann/mlst)
- **Output**: ST, allele profile, identity, coverage
- **Management**: `--pull-mlst-db` to install, `--update-mlst-db` for API-key-authenticated refresh from PubMLST
- **Surveillance value**: identify epidemic clones (ST11, ST258, ST23, ST307) and track their spread

### **3. Kaptive Capsule Typing**
- **Databases**: `kp_k` and `kp_o` (curated by the Kaptive team)
- **Tool**: [Kaptive](https://github.com/klebgenomics/Kaptive)
- **Output**: K locus, O locus, identity, coverage, confidence
- **Clinical relevance**: K1 and K2 are associated with hypervirulence; O1/O2 influence serum resistance

### **4. ABRicate Screening**
- **Tool**: [ABRicate](https://github.com/tseemann/abricate)
- **Thresholds**: customisable with `--abricate-min-id` and `--abricate-min-cov` (default 80/80)
- **Databases (11)**:

  | Database | Description |
  |----------|-------------|
  | CARD | Comprehensive Antibiotic Resistance Database |
  | ResFinder | Acquired resistance genes |
  | ARG-ANNOT | Antibiotic resistance gene annotations |
  | VFDB | Virulence Factor Database |
  | PlasmidFinder | Plasmid replicons |
  | NCBI | NCBI AMR reference gene catalog |
  | MegaRes | Antimicrobial resistance determinants |
  | BacMet2 | Biocide and heavy metal resistance |
  | EcoH / EcoLI_VF | Cross-species screening panels |

### **5. AMRFinderPlus**
- **Tool**: [AMRFinderPlus](https://github.com/ncbi/amr)
- **Organism**: *Klebsiella pneumoniae* (enables point-mutation detection)
- **Output**: per-sample reports plus a summary with risk levels (Critical, High, Standard)
- **Mutations**: point mutations in *gyrA*, *parC*, *rpoB*, *mgrB*, *pmrAB*, *phoPQ*, and more

### **6. MGE Profiling** ✨ NEW
- **Database**: [mobileOG-db](https://github.com/clb21565/mobileOG-db)
- **Output**: per-sample counts across ten functional MGE categories, aggregate profiles by ST and K/O locus, batch HTML, TSV, and JSON summaries
- **Purpose**: connect resistance and virulence content to the mobile elements that carry it

### **7. Gene-centric Ultimate Reporter**
- **Dashboard cards**: samples, STs, capsule types, unique genes, high-risk combinations
- **Sample overview** with typing badges and gene counts
- **MLST, Kaptive, AMR, virulence, plasmid, and BACMET tables** with scrollable genome lists and frequencies
- **Compare tab** for pairwise isolate similarity with verdict banners and shared/unique gene splits
- **Pattern discovery**: ST-K associations, high-risk combinations, ICEKp and virulence plasmid markers
- **Exports**: CSV and JSON from every table

### **8. Sample-centric Ultimate Reporter**
- **Isolate boxes**: one interactive box per sample with typing badges (MLST, K-locus, O-locus, hypervirulence)
- **Per-database tables**: AMRFinderPlus, ABRicate panels, MGE, mutations
- **Mutations table**: gene, mutation, class, contig, coverage, identity, accession
- **Filters and search**: by sample name or database, with one-click CSV export
- **Perfect for**: clinical reports, patient-level investigations, infection control rounds

---

## 🤖 **AI Integration Guide**

Kleboscope reports are structured with clear tables and organized data that AI can easily understand.

1. **Install any AI browser extension** (ChatGPT, Claude, Gemini).
2. **Open your report**: `kleboscope_ultimate_gene_centric_report.html` or the sample-centric report.
3. **Select text** in any section (AMR Genes, MLST, MGE, Compare).
4. **Right-click → Ask AI** with your question.

**Example Questions:**
- *"Which STs carry carbapenemases in this cohort, and what are their K loci?"*
- *"List all isolates with both *blaKPC-2* and colistin resistance genes."*
- *"Which sample pairs in the Compare tab came back near-identical? Are any from the same ward?"*
- *"Which sequence type carries the highest mean MGE load, and what does that suggest?"*

---

## 🔄 **Alternative Tools**

If Kleboscope does not fit your workflow, consider these alternatives:

- **[Kleborate](https://github.com/katholt/Kleborate)**: well-established *K. pneumoniae* typing, resistance, and virulence scoring with tabular output.
- **[Bactopia](https://bactopia.github.io/)**: flexible multi-species pipeline for raw reads (FASTQ) built on Nextflow.
- **[Nullarbor](https://github.com/tseemann/nullarbor)**: multi-species pipeline for raw reads, popular in public health microbiology.
- **[Pathogenwatch](https://pathogenwatch.com/)**: web-based genomic surveillance including *K. pneumoniae*.
- **[AMRFinderPlus standalone](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/)**: if you only need AMR gene detection.

Each has its own strengths and may suit your requirements better depending on raw-read support, web interface preference, or output format.

---

## 📚 **Citation**

### **Citing Kleboscope**
If you use Kleboscope in your research, please cite:

> Beckley, B. (2026). Kleboscope: A gene-centric, species-optimized computational pipeline for comprehensive *Klebsiella pneumoniae* genomic surveillance. *Under review*.

```bibtex
@software{kleboscope2026,
  author = {Beckley, Brown},
  title = {Kleboscope: A gene-centric, species-optimized computational pipeline for comprehensive Klebsiella pneumoniae genomic surveillance},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/Kleboscope}
}
```

### **Acknowledgements**

Kleboscope stands on the shoulders of giants. We are deeply grateful to:

- **Torsten Seemann** for MLST and ABRicate, and countless foundational tools.
- **NCBI team** for AMRFinderPlus.
- **Kaptive team** (Wick, Wyres, Holt) for K and O locus databases.
- **PubMLST, CARD, ResFinder, VFDB, PlasmidFinder, ARG-ANNOT, MegaRes, BacMet** for essential databases.
- **mobileOG-db team** (Brown et al.) for the mobile genetic element resource.
- **ParBLiSS** team for fastANI.
- The broader **ESKAPE ecosystem**, including StaphScope, PseudoScope, AcinetoScope, and the [ESKAPE AMR Platform](https://eskape.bio), for sharing architecture, ideas, and a common mission.

> *"If we ever meet in person, the drinks are on me!" – Brown Beckley*

---

## 👥 **Authors & Contact**

**Brown Beckley** (Primary Developer)
- University of Ghana Medical School – Department of Medical Biochemistry
- 📧 brownbeckley94@gmail.com
- 🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)
- LinkedIn: [@brownbeckley](https://www.linkedin.com/in/brown-beckley-190315319/)

We welcome collaborations on:
- *K. pneumoniae* epidemiology and outbreak studies
- Clinical validation of resistance and virulence markers
- The KLEOS retrieval database
- Development of the ESKAPE AMR Platform web interface
- Integration with public health surveillance systems
- Machine learning applications for phenotype prediction

---

## 📄 **License**

### Core Kleboscope Code
The Kleboscope pipeline code (workflow engine, report generation, HTML templates, and Python modules) is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

### Third-Party Tools
Kleboscope executes several external bioinformatics tools installed as dependencies. Each tool is the property of its respective developers and is used under its own license.

<details>
<summary><b>📚 Third-Party Tool and Database Citations (Click to expand)</b></summary>

#### **MLST**
```bibtex
@software{seemann_mlst_2018,
  author = {Seemann, T.},
  title = {MLST: Scan contig files against traditional PubMLST typing schemes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/mlst}
}
```

#### **Kaptive**
```bibtex
@article{wick_kaptive_2018,
  author = {Wick, R. R. and Heinz, E. and Holt, K. E. and Wyres, K. L.},
  title = {Kaptive Web: User-Friendly Capsule and Lipopolysaccharide Serotype Prediction for Klebsiella Genomes},
  journal = {Journal of Clinical Microbiology},
  volume = {56},
  number = {6},
  pages = {e00197-18},
  year = {2018},
  doi = {10.1128/JCM.00197-18}
}
```

#### **PubMLST**
```bibtex
@article{jolley_pubmlst_2018,
  author = {Jolley, K. A. and Bray, J. E. and Maiden, M. C. J.},
  title = {Open-access bacterial population genomics: {BIGSdb} software, the {PubMLST.org} website and their applications},
  journal = {Wellcome Open Research},
  volume = {3},
  pages = {124},
  year = {2018},
  doi = {10.12688/wellcomeopenres.14826.1}
}
```

#### **ABRicate**
```bibtex
@software{seemann_abricate_2018,
  author = {Seemann, T.},
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/abricate}
}
```

#### **AMRFinderPlus**
```bibtex
@article{feldgarden_amrfinderplus_2021,
  author = {Feldgarden, M. et al.},
  title = {AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
  journal = {Scientific Reports},
  volume = {11},
  pages = {12728},
  year = {2021},
  doi = {10.1038/s41598-021-91456-0}
}
```

#### **mobileOG-db**
```bibtex
@article{brown_mobileogdb_2022,
  author = {Brown, C. L. and Mullet, J. and Hindi, F. and Stoll, J. E. and Gupta, S. and Choi, M. and Keenum, I. and Vikesland, P. and Pruden, A. and Zhang, L.},
  title = {mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements},
  journal = {Applied and Environmental Microbiology},
  volume = {88},
  number = {18},
  pages = {e00991-22},
  year = {2022},
  doi = {10.1128/aem.00991-22}
}
```

#### **fastANI**
```bibtex
@article{jain_fastani_2018,
  author = {Jain, C. and Rodriguez-R, L. M. and Phillippy, A. M. and Konstantinidis, K. T. and Aluru, S.},
  title = {High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries},
  journal = {Nature Communications},
  volume = {9},
  number = {1},
  pages = {5114},
  year = {2018},
  doi = {10.1038/s41467-018-07641-9}
}
```

#### **CARD**
```bibtex
@article{alcock_card_2023,
  author = {Alcock, B. P. et al.},
  title = {CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database},
  journal = {Nucleic Acids Research},
  volume = {51},
  number = {D1},
  pages = {D690-D699},
  year = {2023},
  doi = {10.1093/nar/gkac920}
}
```

#### **ResFinder**
```bibtex
@article{bortolaia_resfinder_2020,
  author = {Bortolaia, V. et al.},
  title = {ResFinder 4.0 for predictions of phenotypes from genotypes},
  journal = {Journal of Antimicrobial Chemotherapy},
  volume = {75},
  number = {12},
  pages = {3491-3500},
  year = {2020},
  doi = {10.1093/jac/dkaa345}
}
```

#### **PlasmidFinder**
```bibtex
@article{carattoli_plasmidfinder_2014,
  author = {Carattoli, A. et al.},
  title = {In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing},
  journal = {Antimicrobial Agents and Chemotherapy},
  volume = {58},
  number = {7},
  pages = {3895-3903},
  year = {2014},
  doi = {10.1128/AAC.02412-14}
}
```

#### **VFDB**
```bibtex
@article{chen_vfdb_2016,
  author = {Chen, L. et al.},
  title = {VFDB 2016: hierarchical and refined dataset for big data analysis--10 years on},
  journal = {Nucleic Acids Research},
  volume = {44},
  number = {D1},
  pages = {D694-D697},
  year = {2016},
  doi = {10.1093/nar/gkv1239}
}
```

#### **BacMet**
```bibtex
@article{pal_bacmet_2014,
  author = {Pal, C. et al.},
  title = {BacMet: antibacterial biocide and metal resistance genes database},
  journal = {Nucleic Acids Research},
  volume = {42},
  number = {D1},
  pages = {D737-D743},
  year = {2014},
  doi = {10.1093/nar/gkt1252}
}
```

</details>

<details>
<summary><b>🤖 AI Disclosure (Click to expand)</b></summary>

Kleboscope is developed with the assistance of AI tools, and we believe in being transparent about it.

**How we use AI:**
- **Code drafting and refactoring.** AI assistants help generate initial implementations, suggest refactors, and catch patterns we might miss.
- **Documentation.** Some prose in this README, module docstrings, and the HTML report text has been drafted with AI assistance and subsequently edited by human authors.
- **Testing scaffolds.** AI helps write boilerplate tests and mock fixtures.

**What we do NOT do:**
- **We do not ship code we don't understand.** Every line that lands in the repository has been reviewed by at least one human author.
- **We do not ship code we can't debug.** If a bug appears in code you can't reason about, that's a bug in the process, not a feature.
- **We do not use AI for scientific claims.** Biological interpretations, thresholds, and validation calls are made by the human authors based on peer-reviewed literature and reference genomes.
- **We do not train models on your data.** Kleboscope is a local tool. Nothing runs through a third-party service. Your genomes stay on your machine.

**Responsibility:**
The authors of Kleboscope take **full responsibility** for the correctness, behaviour, and safety of the released code. AI is a productivity tool in our hands, not an author on our paper. If something in Kleboscope is wrong, the fault is ours, not the tools we used to write it.

</details>

---

<div align="center">

## **🚀 Ready to transform your *K. pneumoniae* surveillance?**

| **Choose Your Tool** | |
|----------------------|-|
| 🖥️ **Command Line** | For high-throughput, local analysis |
| 🧬 **The ESKAPE Suite** | StaphScope · Kleboscope · PseudoScope · AcinetoScope |
| 🔮 **KLEOS Database** | Typing-based genome retrieval, coming soon |

[![Get Started CLI](https://img.shields.io/badge/GET_STARTED_CLI-Now-green?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/Kleboscope#quick-start-cli)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/Kleboscope/issues)

**From FASTA to insight. From fragments to patterns. From data to public health action.**

*Kleboscope: Precision surveillance for the antibiotic resistance era.*

⭐ **If you find this tool useful, please star the repository!** ⭐

*Join the Fight Against Antimicrobial Resistance*

Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in validating markers, sharing regional epidemiological data, and advancing AMR surveillance.

**Together, we can enhance global AMR monitoring and develop more effective treatment strategies.**

</div>
