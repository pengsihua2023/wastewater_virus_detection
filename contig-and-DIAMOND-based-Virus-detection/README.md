# Metagenome Assembly and Diamond Taxonomic Classification Workflow

**Version:** 3.0.0  
**Author:** Bioinformatics Pipeline  
**Last Updated:** October 2025

---

## Table of Contents

1. [Overview](#overview)
2. [Workflow Architecture](#workflow-architecture)
3. [Pipeline Steps](#pipeline-steps)
4. [Input Requirements](#input-requirements)
5. [Output Files](#output-files)
6. [Installation & Setup](#installation--setup)
7. [Usage](#usage)
8. [Result Interpretation](#result-interpretation)
9. [Troubleshooting](#troubleshooting)
10. [References](#references)

---

## Overview

This Nextflow pipeline performs **metagenome assembly and viral taxonomic classification** from paired-end sequencing reads. It combines two assembly strategies (MEGAHIT and metaSPAdes), performs gene prediction with Prodigal, and classifies sequences using Diamond BLASTP against the RVDB viral protein database with full NCBI taxonomy resolution.

### Key Features

- ✅ **Dual Assembly Strategy**: MEGAHIT (fast) + metaSPAdes (sensitive)
- ✅ **Quality Control**: Optional fastp preprocessing
- ✅ **Gene Prediction**: Prodigal in metagenome mode
- ✅ **Viral Classification**: Diamond BLASTP against RVDB database
- ✅ **Full Taxonomy Resolution**: Automatic TaxID to lineage conversion (Superkingdom → Species)
- ✅ **Comprehensive Reports**: Comparative analysis with phylum and family statistics
- ✅ **Reproducible**: Containerized tools (Apptainer/Singularity) + Conda environments
- ✅ **HPC-Ready**: SLURM workload manager integration

---

## Workflow Architecture

```
                        ┌─────────────────────────────────────────┐
                        │   Paired-End Sequencing Reads          │
                        │   (FASTQ files: R1 + R2)               │
                        └──────────────┬──────────────────────────┘
                                      │
                        ┌─────────────▼──────────────────────────┐
                        │   Step 1: Quality Control (Optional)   │
                        │   Tool: fastp                          │
                        │   • Adapter trimming                   │
                        │   • Quality filtering (Q20)            │
                        │   • Length filtering (≥50bp)           │
                        └──────────────┬──────────────────────────┘
                                      │
              ┌───────────────────────┴───────────────────────┐
              │                                               │
┌─────────────▼──────────────┐                  ┌─────────────▼──────────────┐
│   Step 2a: MEGAHIT Assembly│                  │   Step 2b: SPAdes Assembly │
│   Container: megahit:1.2.9 │                  │   Container: spades:3.15.5 │
│   • Fast k-mer based       │                  │   • de Bruijn graph based  │
│   • Memory efficient       │                  │   • High sensitivity       │
│   • Min contig: 1000 bp    │                  │   • Metagenomic mode       │
└─────────────┬──────────────┘                  └─────────────┬──────────────┘
              │                                               │
              │   Contigs (DNA)                              │   Contigs (DNA)
              │                                               │
┌─────────────▼──────────────┐                  ┌─────────────▼──────────────┐
│   Step 3a: Prodigal        │                  │   Step 3b: Prodigal        │
│   (MEGAHIT contigs)        │                  │   (SPAdes contigs)         │
│   Conda: prodigal=2.6.3    │                  │   Conda: prodigal=2.6.3    │
│   • Gene prediction        │                  │   • Gene prediction        │
│   • Metagenome mode (-p)   │                  │   • Metagenome mode (-p)   │
│   • Output: Proteins (FAA) │                  │   • Output: Proteins (FAA) │
└─────────────┬──────────────┘                  └─────────────┬──────────────┘
              │                                               │
              │   Predicted Proteins                         │   Predicted Proteins
              │                                               │
┌─────────────▼──────────────┐                  ┌─────────────▼──────────────┐
│   Step 4a: Diamond BLASTP  │                  │   Step 4b: Diamond BLASTP  │
│   (MEGAHIT proteins)       │                  │   (SPAdes proteins)        │
│   Conda: diamond=2.1.8     │                  │   Conda: diamond=2.1.8     │
│   • Query: Proteins        │                  │   • Query: Proteins        │
│   • Database: RVDB         │                  │   • Database: RVDB         │
│   • E-value: 1e-5          │                  │   • E-value: 1e-5          │
│   • Sensitive mode         │                  │   • Sensitive mode         │
│   • Output: TaxIDs         │                  │   • Output: TaxIDs         │
└─────────────┬──────────────┘                  └─────────────┬──────────────┘
              │                                               │
              │   Diamond Results (TaxIDs)                   │   Diamond Results (TaxIDs)
              │                                               │
              └───────────────────────┬───────────────────────┘
                                      │
                        ┌─────────────▼──────────────────────────┐
                        │   Step 5: Taxonomy Resolution &        │
                        │           Report Generation            │
                        │   Python Script (Pandas)               │
                        │                                        │
                        │   • Load NCBI Taxonomy (names.dmp +    │
                        │     nodes.dmp)                         │
                        │   • Map TaxIDs → Full Lineage          │
                        │     (Superkingdom → Species)           │
                        │   • Generate enhanced files with       │
                        │     organism names and taxonomy        │
                        │   • Create comparative statistics      │
                        │     (Phylum, Family levels)            │
                        │   • Generate comprehensive reports     │
                        └──────────────┬──────────────────────────┘
                                      │
                        ┌─────────────▼──────────────────────────┐
                        │   Final Outputs                        │
                        │                                        │
                        │   • MEGAHIT + Taxonomy (22 columns)    │
                        │   • SPAdes + Taxonomy (22 columns)     │
                        │   • Merged Report (TXT)                │
                        │   • Comparison Data (CSV)              │
                        └────────────────────────────────────────┘
```

---

## Pipeline Steps

### Step 1: Quality Control (Optional)

**Tool:** [fastp](https://github.com/OpenGene/fastp) v0.23.4

**Purpose:** Remove low-quality reads and adapters to improve assembly quality.

**Parameters:**
- Minimum quality score: Q20
- Maximum unqualified bases: 40%
- Minimum read length: 50 bp
- Automatic adapter detection and removal

**Input:** Paired-end FASTQ files (R1 + R2)  
**Output:** Clean reads (R1 + R2), HTML report, JSON statistics

**Skip this step:** Set `--skip_fastp` parameter

---

### Step 2: Metagenome Assembly

#### **2a. MEGAHIT Assembly**

**Tool:** [MEGAHIT](https://github.com/voutcn/megahit) v1.2.9  
**Container:** `docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1`

**Strategy:**
- k-mer based assembly (adaptive k-mer selection)
- Memory efficient (uses ~80% of available RAM)
- Fast (typically 1-2 hours for bacterial metagenomes)

**Parameters:**
- Memory usage: 80% of allocated RAM
- Minimum contig length: 1000 bp
- Threads: 16 CPUs
- Allocated memory: 64 GB
- Time limit: 12 hours

**Output:** Assembled contigs (FASTA format)

---

#### **2b. metaSPAdes Assembly**

**Tool:** [SPAdes](http://cab.spbu.ru/software/spades/) v3.15.5  
**Container:** `docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1`

**Strategy:**
- de Bruijn graph based assembly
- Metagenome-specific mode (`--meta`)
- Higher sensitivity for complex communities
- Error correction disabled (to avoid memory issues)

**Parameters:**
- Metagenome mode: Enabled
- Minimum contig length: 1000 bp
- Threads: 32 CPUs
- Allocated memory: 512 GB (large memory partition)
- Time limit: 48 hours

**Output:** Assembled contigs (FASTA format)

**Note:** SPAdes is more sensitive but requires more computational resources.

---

### Step 3: Gene Prediction

**Tool:** [Prodigal](https://github.com/hyattpd/Prodigal) v2.6.3  
**Environment:** Conda (`bioconda::prodigal=2.6.3`)

**Purpose:** Predict protein-coding genes from assembled DNA contigs.

**Mode:** Metagenome mode (`-p meta`)
- Optimized for short contigs
- Anonymous mode (no training)
- Suitable for diverse microbial communities

**Process:**
- Input: DNA contigs from MEGAHIT or SPAdes
- Output: 
  - Protein sequences (`.faa` - amino acids)
  - Gene sequences (`.fna` - nucleotides)

**Parameters:**
- Prediction mode: Metagenome (`-p meta`)
- Quiet mode: Enabled (`-q`)
- Threads: 4 CPUs
- Memory: 8 GB
- Time: 2 hours

---

### Step 4: Taxonomic Classification with Diamond

**Tool:** [Diamond](https://github.com/bbuchfink/diamond) v2.1.8  
**Environment:** Conda (`bioconda::diamond=2.1.8`)

**Database:** RVDB (Reference Viral DataBase) v30.0
- **Source:** [RVDB](https://hive.biochemistry.gwu.edu/rvdb)
- **Content:** Curated viral protein sequences
- **Size:** ~785,890 sequences, ~292 million amino acids
- **Taxonomy:** Enriched with NCBI TaxIDs

**Why RVDB?**
- Specialized for viral metagenomics
- High-quality curated sequences
- Regular updates
- Much smaller than NCBI nr (faster searches)

**Diamond BLASTP Algorithm:**
- Protein-protein alignment (BLASTP-like)
- ~20,000x faster than BLAST
- Sensitive mode for better sensitivity
- E-value threshold: 1e-5

**Parameters:**
- Search mode: Sensitive (`--sensitive`)
- E-value cutoff: 1e-5
- Max target sequences: 1 (best hit only)
- Threads: 32 CPUs
- Output format: Tabular (13 columns including staxids)

**Output Columns:**
```
1.  qseqid      - Query sequence ID (protein)
2.  sseqid      - Subject sequence ID (database hit)
3.  pident      - Percent identity
4.  length      - Alignment length
5.  mismatch    - Number of mismatches
6.  gapopen     - Number of gap openings
7.  qstart      - Query alignment start
8.  qend        - Query alignment end
9.  sstart      - Subject alignment start
10. send        - Subject alignment end
11. evalue      - E-value
12. bitscore    - Bit score
13. staxids     - Subject taxonomic ID
```

**Resource Allocation:**
- CPUs: 32
- Memory: 64 GB
- Time: 12 hours

---

### Step 5: Taxonomy Resolution and Report Generation

**Language:** Python 3 (embedded script)  
**Dependencies:** pandas v2.0.3

**Process Flow:**

#### **5.1 Load NCBI Taxonomy Database**

**Files Required:**
- `names.dmp`: Taxonomic names (~2.7M entries)
  - Maps TaxID → Scientific name
  - Example: `68887 → "Torque teno virus 1"`

- `nodes.dmp`: Taxonomic hierarchy (~3.5M nodes)
  - Defines parent-child relationships
  - Specifies taxonomic rank (superkingdom, phylum, class, etc.)

**Loading Time:** ~2-5 seconds

---

#### **5.2 TaxID to Lineage Resolution**

For each hit in Diamond results:

1. **Extract TaxID** from column 13 (staxids)
2. **Lookup organism name** in names.dmp
3. **Traverse taxonomy tree** from TaxID to root:
   - Start at current TaxID
   - Move up parent nodes
   - Collect names at each major rank
4. **Build lineage** with 9 taxonomic levels:
   - `organism_name`: Full scientific name
   - `superkingdom`: e.g., Viruses, Bacteria
   - `kingdom`
   - `phylum`: e.g., Pisuviricota
   - `class`
   - `order`
   - `family`: e.g., Anelloviridae
   - `genus`
   - `species`

**Handling Missing Ranks:** If a rank doesn't exist in the lineage, it's marked as "N/A"

**Example Resolution:**
```
TaxID 68887 →
  organism_name: Torque teno virus 1
  superkingdom:  Viruses
  kingdom:       N/A
  phylum:        N/A
  class:         N/A
  order:         N/A
  family:        Anelloviridae
  genus:         Alphatorquevirus
  species:       Torque teno virus 1
```

---

#### **5.3 Generate Enhanced Output Files**

**Format:** Tab-separated values (TSV)  
**Columns:** 22 (13 original + 9 taxonomy)

**File Naming:**
- `{sample}_megahit_with_taxonomy.txt`
- `{sample}_spades_with_taxonomy.txt`

---

#### **5.4 Generate Comparative Statistics**

The script calculates:

**Overall Statistics:**
- Total hits (alignment count)
- Unique query proteins
- Average identity (%)
- Average alignment length

**Phylum-Level Comparison:**
- Top 15 most abundant phyla
- Counts for MEGAHIT vs SPAdes
- Total counts across both assemblies

**Family-Level Comparison:**
- Top 15 most abundant families
- Counts for MEGAHIT vs SPAdes
- Useful for identifying dominant viral families

**Taxonomic ID Comparison:**
- Top 20 most abundant TaxIDs
- Detailed breakdown by assembler

**Unique Findings:**
- TaxIDs found only in SPAdes assembly
- TaxIDs found only in MEGAHIT assembly
- Useful for assembly method evaluation

---

#### **5.5 Generate Reports**

**Text Report** (`*_merged_report.txt`):
- Human-readable summary
- Formatted tables
- Statistical comparisons

**CSV Report** (`*_merged_report.csv`):
- Machine-readable format
- Full taxonomic comparison data
- Suitable for further analysis in R/Python/Excel

---

## Input Requirements

### Samplesheet Format

**File:** `samplesheet.csv`

**Required Columns:**
```csv
sample,fastq_1,fastq_2
```

**Example:**
```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
sample3,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz
```

**Notes:**
- Paths can be absolute or relative
- FASTQ files must be gzipped (`.fastq.gz`)
- Sample names must be unique
- No spaces in sample names (use underscores)

---

### Database Requirements

#### **1. RVDB Diamond Database**

**Required File:** `RVDB_prot_ref.dmnd`

**Location:** Specify in `run_metagenome_assembly_classification_en.sh`:
```bash
DIAMOND_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB/RVDB_prot_ref.dmnd"
```

**Building the Database:** See `RVDB Database Setup` section below.

---

#### **2. NCBI Taxonomy Files**

**Required Files:**
- `names.dmp`: Taxonomic names
- `nodes.dmp`: Taxonomic hierarchy

**Location:** Specify in `metagenome_assembly_classification_en.config`:
```groovy
taxonomy_names = '/path/to/names.dmp'
taxonomy_nodes = '/path/to/nodes.dmp'
```

**Download Instructions:**
```bash
cd /path/to/databases/RVDB

# Download NCBI taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extract
tar -xzvf taxdump.tar.gz

# This creates names.dmp and nodes.dmp
```

---

## Output Files

### Directory Structure

```
results/
├── fastp/                              # Quality control reports
│   ├── {sample}_fastp.html            # Visual QC report
│   └── {sample}_fastp.json            # QC statistics (JSON)
│
├── prodigal_megahit/                   # MEGAHIT gene predictions
│   └── {sample}_megahit_proteins.faa  # Predicted proteins
│
├── prodigal_spades/                    # SPAdes gene predictions
│   └── {sample}_spades_proteins.faa   # Predicted proteins
│
├── diamond_megahit/                    # MEGAHIT Diamond results
│   └── {sample}_megahit_diamond.txt   # Raw Diamond output (13 columns)
│
├── diamond_spades/                     # SPAdes Diamond results
│   └── {sample}_spades_diamond.txt    # Raw Diamond output (13 columns)
│
└── merged_reports/                     # Comprehensive analysis
    ├── {sample}_megahit_with_taxonomy.txt  # Enhanced MEGAHIT (22 columns)
    ├── {sample}_spades_with_taxonomy.txt   # Enhanced SPAdes (22 columns)
    ├── {sample}_merged_report.txt          # Text summary
    └── {sample}_merged_report.csv          # CSV comparison data
```

---

### Key Output Files

#### **1. Enhanced Diamond Results (with Taxonomy)**

**Files:**
- `{sample}_megahit_with_taxonomy.txt`
- `{sample}_spades_with_taxonomy.txt`

**Format:** Tab-separated, 22 columns

**Columns:**
```
Original 13 columns from Diamond + 9 taxonomy columns:

Diamond Columns (1-13):
qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, 
sstart, send, evalue, bitscore, staxids

Taxonomy Columns (14-22):
organism_name, superkingdom, kingdom, phylum, class, order, 
family, genus, species
```

**Example Row:**
```
k141_1_1  UGV43004.1  98.5  250  3  0  1  250  1  250  1.2e-50  450  68887  Torque teno virus 1  Viruses  N/A  N/A  N/A  N/A  Anelloviridae  Alphatorquevirus  Torque teno virus 1
```

---

#### **2. Merged Text Report**

**File:** `{sample}_merged_report.txt`

**Content:**

```
================================================================================
Diamond Comprehensive Analysis Report - MEGAHIT vs SPAdes Assembly Comparison
================================================================================

[Overall Statistics]
--------------------------------------------------------------------------------
SPAdes total hits:          12,345
MEGAHIT total hits:         15,234

SPAdes unique queries:      5,678
MEGAHIT unique queries:     6,789

SPAdes average identity:    87.45%
MEGAHIT average identity:   85.32%

SPAdes average length:      245.6 aa
MEGAHIT average length:     238.2 aa


[Phylum Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Phylum                         SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
N/A                           5,680           7,234           12,914
Pisuviricota                  2,340           3,120           5,460
Uroviricota                   1,890           2,456           4,346
...


[Family Level Comparison (Top 15)]
--------------------------------------------------------------------------------
Family                         SPAdes Count    MEGAHIT Count   Total          
--------------------------------------------------------------------------------
Anelloviridae                 3,450           4,890           8,340
Picornaviridae                1,230           1,560           2,790
Hepadnaviridae                890             1,120           2,010
...


[Taxonomic ID Comparison (Top 20)]
--------------------------------------------------------------------------------
Taxonomic ID         SPAdes Count    MEGAHIT Count   
--------------------------------------------------------------------------------
68887               1,250           1,680           
687329              980             1,340           
...


[Unique Findings]
--------------------------------------------------------------------------------

Taxonomic IDs found only in SPAdes: 234
  - 12345: 45 matches
  - 67890: 32 matches
  ...

Taxonomic IDs found only in MEGAHIT: 189
  - 23456: 67 matches
  - 78901: 54 matches
  ...

================================================================================
Analysis Complete
================================================================================
```

---

#### **3. CSV Comparison Data**

**File:** `{sample}_merged_report.csv`

**Columns:**
```
taxonomic_id, spades_count, megahit_count, total_count, 
spades_only, megahit_only, both
```

**Purpose:** Import into R/Python/Excel for custom analysis and visualization.

---

## Installation & Setup

### Prerequisites

**System Requirements:**
- Linux operating system
- SLURM workload manager
- Apptainer/Singularity (for containers)
- Conda/Mamba (for environment management)
- Nextflow ≥ 25.04

**Minimum Hardware:**
- CPUs: 32 cores
- RAM: 512 GB (for SPAdes large assemblies)
- Storage: 500 GB free space

---

### 1. Install Nextflow

```bash
# Download Nextflow
cd ~/bin
wget -qO- https://get.nextflow.io | bash

# Make executable
chmod +x nextflow

# Add to PATH
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Verify installation
nextflow -version
```

---

### 2. Setup Conda Environment

```bash
# Load Miniforge module (if using module system)
module load Miniforge3/24.11.3-0

# Create Nextflow environment
conda create -n nextflow_env \
    nextflow \
    fastp=0.23.4 \
    prodigal=2.6.3 \
    diamond=2.1.8 \
    pandas=2.0.3 \
    -c bioconda -c conda-forge

# Activate environment
conda activate nextflow_env
```

---

### 3. Prepare RVDB Database

See detailed instructions in `RVDB_DATABASE_SETUP.md` (provided separately).

**Quick Setup:**

```bash
# Navigate to database directory
cd /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB

# Download RVDB FASTA
wget https://hive.biochemistry.gwu.edu/cuts/rvdb/U-RVDBv30.0-prot_clustered.fasta.xz
xz -d U-RVDBv30.0-prot_clustered.fasta.xz

# Download NCBI taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz

# Build Diamond database with taxonomy
# (Use provided scripts: rebuild_rvdb_simple.sh or similar)
```

---

### 4. Clone Workflow Repository

```bash
cd /scratch/sp96859/Meta-genome-data-analysis/Apptainer
git clone <repository_url> Contig-Diamond-based
cd Contig-Diamond-based
```

---

## Usage

### Basic Usage

```bash
# Activate conda environment
conda activate nextflow_env

# Submit SLURM job
sbatch run_metagenome_assembly_classification_en.sh
```

---

### Advanced Usage

**Run with custom parameters:**

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet.csv \
    --outdir results \
    --diamond_db /path/to/RVDB_prot_ref.dmnd \
    --skip_fastp \
    --diamond_evalue 1e-10 \
    -resume
```

**Parameters:**
- `--input`: Path to samplesheet (required)
- `--outdir`: Output directory (default: `results`)
- `--diamond_db`: Path to Diamond database (required)
- `--skip_fastp`: Skip quality control
- `--skip_merge_reports`: Skip report generation
- `--diamond_evalue`: E-value threshold (default: 1e-5)
- `--diamond_max_target_seqs`: Max targets per query (default: 1)
- `-resume`: Resume from last checkpoint

---

### Monitoring Progress

```bash
# Check SLURM job status
squeue -u $USER

# View real-time log
tail -f slurm-<JOB_ID>.out

# Check Nextflow log
tail -f .nextflow.log
```

---

## Result Interpretation

### Understanding Taxonomy Levels

**Viral Classification Hierarchy:**
```
Superkingdom: Viruses
    ↓
Kingdom: [Often N/A for viruses]
    ↓
Phylum: e.g., Pisuviricota, Uroviricota
    ↓
Class: [Often N/A]
    ↓
Order: e.g., Picornavirales
    ↓
Family: e.g., Anelloviridae, Picornaviridae
    ↓
Genus: e.g., Alphatorquevirus
    ↓
Species: e.g., Torque teno virus 1
```

**Note:** Many viral lineages have "N/A" for intermediate ranks due to incomplete taxonomy.

---

### Comparing MEGAHIT vs SPAdes

**MEGAHIT:**
- ✅ Faster (1-2 hours)
- ✅ Lower memory usage
- ✅ Good for abundant organisms
- ⚠️ May miss low-abundance organisms

**SPAdes:**
- ✅ More sensitive
- ✅ Better for complex communities
- ✅ Longer contigs
- ⚠️ Slower (12-24 hours)
- ⚠️ Higher memory requirement

**Strategy:** Compare both results to get:
1. Shared hits: High-confidence identifications
2. MEGAHIT-only: Abundant, well-assembled
3. SPAdes-only: Low-abundance or complex regions

---

### Statistical Analysis Examples

**In R:**
```r
library(tidyverse)

# Load enhanced results
megahit <- read_tsv("sample_megahit_with_taxonomy.txt")
spades <- read_tsv("sample_spades_with_taxonomy.txt")

# Top families
megahit %>% 
  count(family, sort = TRUE) %>% 
  filter(family != "N/A") %>%
  head(10)

# Abundance by phylum
megahit %>%
  count(phylum, name = "count") %>%
  mutate(percentage = count / sum(count) * 100) %>%
  arrange(desc(count))

# Venn diagram of TaxIDs
library(VennDiagram)
megahit_taxids <- unique(megahit$staxids)
spades_taxids <- unique(spades$staxids)

venn.diagram(
  list(MEGAHIT = megahit_taxids, SPAdes = spades_taxids),
  filename = "taxid_venn.png"
)
```

---

## Troubleshooting

### Common Issues

#### **1. Diamond Error: "Options require taxonomy information"**

**Cause:** Database built without taxonomy  
**Solution:** Rebuild database with `--taxonmap`, `--taxonnodes`, `--taxonnames`  
See: RVDB database setup guide

---

#### **2. SPAdes Out of Memory**

**Cause:** Insufficient RAM for complex metagenomes  
**Solution:**
- Increase memory allocation in config (e.g., 768 GB)
- Use `--meta` mode only
- Reduce k-mer sizes
- Subsample reads

---

#### **3. Nextflow "Cannot find file" Error**

**Cause:** Incorrect paths in samplesheet  
**Solution:**
- Use absolute paths
- Verify files exist: `ls -lh /path/to/file.fastq.gz`
- Check for typos

---

#### **4. No Hits in Diamond Results**

**Possible Causes:**
- E-value too stringent → Try `1e-3` instead of `1e-5`
- Wrong database (bacterial db for viral samples)
- Poor quality assemblies → Check contig N50, total length
- Sample doesn't contain target organisms

---

## References

### Tools

1. **fastp:** Chen et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884-i890.

2. **MEGAHIT:** Li et al. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 31(10):1674-1676.

3. **SPAdes:** Bankevich et al. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. *J Comput Biol* 19(5):455-477.

4. **Prodigal:** Hyatt et al. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11:119.

5. **Diamond:** Buchfink et al. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods* 18:366-368.

6. **Nextflow:** Di Tommaso et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35:316-319.

### Databases

7. **RVDB:** Goodacre et al. (2018). A Reference Viral Database (RVDB) To Enhance Bioinformatics Analysis of High-Throughput Sequencing for Novel Virus Detection. *mSphere* 3(2):e00069-18.

8. **NCBI Taxonomy:** Schoch et al. (2020). NCBI Taxonomy: a comprehensive update on curation, resources and tools. *Database* 2020:baaa062.

---

## Support

For questions, issues, or contributions:
- Create an issue in the GitHub repository
- Contact: bioinformatics-support@institution.edu
- Documentation: See `/docs` directory

---

## License

This pipeline is released under the MIT License.

---

## Citation

If you use this pipeline in your research, please cite:

```
[Your Name] et al. (2025). Metagenome Assembly and Diamond Taxonomic 
Classification Workflow. GitHub repository: [URL]
```

And cite all individual tools used (see References section).

---

**Last Updated:** October 2025  
**Version:** 3.0.0  
**Status:** Production-ready ✅
