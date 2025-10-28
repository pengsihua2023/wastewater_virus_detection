# Viral Metagenome Analysis Workflow

**Version**: 5.1.0  
**Last Updated**: October 28, 2025  
**Pipeline**: Nextflow DSL2  

A comprehensive bioinformatics workflow for viral sequence identification from metagenomic data using dual viral classification tools (VirSorter2 + DeepVirFinder) with assembler comparison.

---

## 🌟 Key Features

- **Dual Assemblers**: Parallel assembly using MEGAHIT and metaSPAdes for comprehensive coverage
- **Dual Viral Classifiers**: VirSorter2 (rule-based + ML) and DeepVirFinder (deep learning)
- **Assembler Comparison**: Identify high-confidence consensus viral sequences across assemblers
- **Three-Level Validation**: Single tool → Tool consensus → Full consensus
- **Quality Control**: Automated read quality assessment with fastp
- **Reproducible**: Conda/Apptainer containers ensure reproducibility
- **Resume Support**: Continue from checkpoints with `-resume`
- **SLURM Compatible**: Optimized for HPC cluster execution

---

## 📋 Table of Contents

- [Workflow Overview](#workflow-overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Requirements](#input-requirements)
- [Output Structure](#output-structure)
- [Parameters](#parameters)
- [Understanding Results](#understanding-results)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)

---

## 🔬 Workflow Overview

### Pipeline Stages

```
Raw Reads (FASTQ)
    ↓
┌─────────────────┐
│  1. fastp QC    │ → Clean reads saved
└─────────────────┘
    ↓
┌─────────────────────────────────────┐
│  2. Parallel Assembly               │
│     • MEGAHIT (fast, complex)       │
│     • metaSPAdes (sensitive)        │
└─────────────────────────────────────┘
    ↓                    ↓
MEGAHIT Contigs    SPAdes Contigs (both saved)
    ↓                    ↓
┌─────────────────────────────────────┐
│  3. Dual Viral Classification       │
│     • VirSorter2 (ML + rules)       │
│     • DeepVirFinder (deep learning) │
└─────────────────────────────────────┘
    ↓                    ↓
┌─────────────────────────────────────┐
│  4. Merge Results per Assembler    │
│     • Tool consensus identification │
│     • Confidence scoring            │
└─────────────────────────────────────┘
    ↓                    ↓
MEGAHIT Viral List  SPAdes Viral List
    ↓                    ↓
┌─────────────────────────────────────┐
│  5. Assembler Comparison           │
│     • Cross-assembler validation   │
│     • ⭐ Final consensus list       │
└─────────────────────────────────────┘
    ↓
HIGH-CONFIDENCE VIRAL SEQUENCES ✨
```

### Validation Levels

1. **Level 1: Single Tool Detection**
   - Identified by VirSorter2 OR DeepVirFinder
   - Confidence: Medium (⭐)

2. **Level 2: Tool Consensus (Recommended)**
   - Identified by BOTH VirSorter2 AND DeepVirFinder
   - Confidence: High (⭐⭐)

3. **Level 3: Full Consensus (Highest)**
   - Tool consensus in BOTH MEGAHIT AND SPAdes assemblies
   - Confidence: Very High (⭐⭐⭐)
   - **→ Use these for downstream analysis!**

---

## 🛠️ Installation

### Prerequisites

- **Nextflow** ≥ 25.04
- **Conda/Mamba** (for package management)
- **Apptainer/Singularity** (for containers)
- **SLURM** (for HPC job submission)

### Step 1: Clone Repository

```bash
cd /your/project/directory
git clone <repository_url> VirSorter2-DeepVirFinder
cd VirSorter2-DeepVirFinder
```

### Step 2: Setup Conda Environment

```bash
# Load conda module (adjust for your system)
module load Miniforge3/24.11.3-0

# Create main workflow environment
conda create -n nextflow_env -c bioconda -c conda-forge \
    nextflow=25.04 \
    fastp \
    virsorter=2.2.4 \
    screed \
    pandas \
    numpy

# Activate environment
conda activate nextflow_env
```

### Step 3: Setup DeepVirFinder Environment

DeepVirFinder requires **Python 3.6** with specific older package versions:

```bash
# Create dedicated environment for DeepVirFinder
conda create -n dvf python=3.6 -y
conda activate dvf

# Install dependencies in specific order
conda install -c conda-forge -c bioconda h5py=2.10.0 -y
conda install -c conda-forge numpy=1.16.0 -y
pip install theano==1.0.5
pip install keras==2.2.4
conda install -c bioconda scikit-learn=0.22.1 biopython pandas -y

# Configure Keras to use Theano backend
mkdir -p ~/.keras
cat > ~/.keras/keras.json << EOF
{
    "image_data_format": "channels_last",
    "epsilon": 1e-07,
    "floatx": "float32",
    "backend": "theano"
}
EOF

# Clone DeepVirFinder
cd /your/tools/directory
git clone https://github.com/jessieren/DeepVirFinder
cd DeepVirFinder
```

### Step 4: Setup VirSorter2 Database

```bash
# Activate main environment
conda activate nextflow_env

# Setup database (~13GB download)
DB_PATH="/scratch/<your_id>/databases/virsorter2/db"
virsorter setup -d $DB_PATH -j 4

# Verify
ls -lh $DB_PATH
```

### Step 5: Update Configuration

Edit `run_metagenome_assembly_classification_en.sh`:

```bash
# Set your paths
VIRSORTER2_DB="/scratch/<your_id>/databases/virsorter2/db"
DEEPVIRFINDER_DIR="/path/to/DeepVirFinder"

# Adjust SLURM parameters
#SBATCH --partition=your_partition
#SBATCH --mem=256G
#SBATCH --time=72:00:00
```

---

## 🚀 Quick Start

### 1. Prepare Sample Sheet

Create `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

### 2. Run Workflow

```bash
# Activate environment
conda activate nextflow_env

# Submit to SLURM
sbatch run_metagenome_assembly_classification_en.sh

# Monitor progress
tail -f Viral_Classification_*.out

# Check job status
squeue -u $USER
```

### 3. Resume Failed/Incomplete Runs

```bash
# Workflow automatically resumes from last successful checkpoint
sbatch run_metagenome_assembly_classification_en.sh

# Nextflow detects cached results and only re-runs failed steps
```

---

## 📥 Input Requirements

### Required Files

1. **Sample Sheet** (`samplesheet.csv`)
   - Columns: `sample`, `fastq_1`, `fastq_2`
   - Paired-end FASTQ files (gzipped or uncompressed)

2. **VirSorter2 Database**
   - Location: Specified by `--virsorter2_db`
   - Size: ~13GB

3. **DeepVirFinder Installation**
   - Location: Specified by `--deepvirfinder_dir`
   - With pre-trained models

### Data Quality Recommendations

- **Read length**: ≥100 bp (150 bp recommended)
- **Coverage**: ≥1 million paired-end reads per sample
- **Quality**: Phred score ≥20 (fastp will filter)

---

## 📂 Output Structure

```
results/
├── clean_reads/                          # QC-filtered reads
│   ├── sample1_clean_R1.fastq.gz
│   └── sample1_clean_R2.fastq.gz
│
├── fastp/                                # Quality control reports
│   ├── sample1_fastp.html               ⭐ View QC metrics
│   └── sample1_fastp.json
│
├── assembly_megahit/                     # MEGAHIT assembled contigs
│   └── sample1_megahit_contigs.fa       ⭐ Contigs ≥1kb
│
├── assembly_spades/                      # SPAdes assembled contigs
│   └── sample1_spades_contigs.fa        ⭐ Contigs ≥1kb
│
├── virsorter2_megahit/                   # VirSorter2 results (MEGAHIT)
│   ├── sample1_megahit_vs2_final-viral-score.tsv
│   ├── sample1_megahit_vs2_final-viral-combined.fa
│   └── sample1_megahit_vs2_final-viral-boundary.tsv
│
├── virsorter2_spades/                    # VirSorter2 results (SPAdes)
│   └── sample1_spades_vs2_final-viral-score.tsv
│
├── deepvirfinder_megahit/                # DeepVirFinder results (MEGAHIT)
│   └── sample1_megahit_dvf_output.txt
│
├── deepvirfinder_spades/                 # DeepVirFinder results (SPAdes)
│   └── sample1_spades_dvf_output.txt
│
├── merged_viral_reports_megahit/         # Merged results per assembler
│   ├── sample1_megahit_viral_merged_report.txt
│   ├── sample1_megahit_viral_merged_report.csv   ⭐ Detailed data
│   └── sample1_megahit_viral_consensus.txt       ⭐ Tool consensus
│
├── merged_viral_reports_spades/
│   ├── sample1_spades_viral_merged_report.txt
│   ├── sample1_spades_viral_merged_report.csv
│   └── sample1_spades_viral_consensus.txt
│
└── assembler_comparison/                 # Final comparison
    ├── sample1_assembler_comparison.txt          ⭐ Summary report
    ├── sample1_assembler_comparison.csv
    └── sample1_consensus_viral_sequences.txt     ⭐⭐⭐ USE THIS!
```

### Most Important Outputs

1. **`assembler_comparison/sample_consensus_viral_sequences.txt`** ⭐⭐⭐
   - **Highest confidence viral sequences**
   - Validated by: Both tools + Both assemblers
   - **→ Use these for downstream analysis**

2. **`merged_viral_reports_*/sample_*_viral_merged_report.csv`** ⭐⭐
   - Detailed per-contig information
   - Scores from both tools
   - Tool consensus status

3. **`assembly_*/sample_*_contigs.fa`** ⭐
   - Assembled contigs for further analysis

---

## ⚙️ Parameters

### Essential Parameters

```bash
--input                 # Sample sheet CSV path
--outdir               # Output directory (default: results)
--virsorter2_db        # VirSorter2 database path
--deepvirfinder_dir    # DeepVirFinder installation directory
```

### Quality Control

```bash
--save_clean_reads     # Save fastp cleaned reads (default: true)
--fastp_qualified_quality_phred  # Min quality score (default: 20)
```

### Assembly Parameters

```bash
--megahit_min_contig_len      # MEGAHIT min contig length (default: 1000)
--spades_min_contig_len       # SPAdes min contig length (default: 1000)
```

### Viral Classification Thresholds

```bash
# VirSorter2
--virsorter2_min_length       # Min contig length (default: 1000 bp)
--virsorter2_min_score        # Min viral score (default: 0.5)

# DeepVirFinder
--deepvirfinder_min_length    # Min contig length (default: 1000 bp)
--deepvirfinder_pvalue        # P-value threshold (default: 0.05)
```

### Example Custom Run

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet.csv \
    --outdir results \
    --virsorter2_db /path/to/db \
    --deepvirfinder_dir /path/to/DeepVirFinder \
    --virsorter2_min_score 0.7 \
    --deepvirfinder_pvalue 0.01 \
    -resume
```

---

## 📊 Understanding Results

### Merged Report CSV Columns

```csv
sequence_name,identified_by,vs2_score,vs2_group,dvf_score,dvf_pvalue,consensus
k141_12345,Both,0.953,dsDNAphage,0.876,0.001,True
k141_67890,VirSorter2_only,0.721,ssDNA,,,False
NODE_123,DeepVirFinder_only,,,0.234,0.032,False
```

**Columns**:
- `sequence_name`: Contig identifier
- `identified_by`: Detection method (Both/VirSorter2_only/DeepVirFinder_only)
- `vs2_score`: VirSorter2 score (0-1, higher = more viral)
- `vs2_group`: Viral group (dsDNAphage, ssDNA, RNA, etc.)
- `dvf_score`: DeepVirFinder score (0-1)
- `dvf_pvalue`: Statistical significance (lower = more confident)
- `consensus`: TRUE if both tools agree

### Interpreting Scores

#### VirSorter2 Score
- **≥0.9**: High confidence viral
- **0.7-0.9**: Likely viral
- **0.5-0.7**: Possible viral (review manually)

#### DeepVirFinder P-value
- **<0.01**: High confidence viral
- **0.01-0.05**: Likely viral
- **>0.05**: Low confidence

### Typical Results Distribution

For a metagenomic sample with ~50,000 contigs:

```
Total Contigs: 50,000

VirSorter2:
├── Identified: 1,500 (3%)
└── High score (≥0.9): 600

DeepVirFinder:
├── Analyzed: 50,000
└── Significant (p<0.05): 2,000

Merged (Tool Consensus):
├── Both tools: 150-300 ⭐⭐
├── VirSorter2 only: 1,200-1,400
└── DeepVirFinder only: 1,700-1,800

Final (Assembler Consensus):
└── MEGAHIT + SPAdes: 50-150 ⭐⭐⭐
```

**Expected consensus rate**: 5-20% of single-tool identifications

---

## 🔧 Troubleshooting

### Common Issues

#### 1. ImportError: No module named 'numpy'

**Problem**: Python packages not in conda environment

**Solution**:
```bash
conda activate nextflow_env
conda install -c conda-forge pandas numpy -y
```

#### 2. ModuleNotFoundError: No module named 'h5py' (DeepVirFinder)

**Problem**: DeepVirFinder dependencies missing or wrong version

**Solution**:
```bash
conda activate dvf
conda install h5py=2.10.0 -y
pip install theano==1.0.5 keras==2.2.4
```

#### 3. Keras backend error

**Problem**: Keras using wrong backend (TensorFlow instead of Theano)

**Solution**:
```bash
# Set in shell
export KERAS_BACKEND=theano

# Or configure in ~/.keras/keras.json
echo '{"backend": "theano"}' > ~/.keras/keras.json
```

#### 4. VirSorter2 database not found

**Solution**:
```bash
virsorter setup -d /path/to/db -j 4
# Then update --virsorter2_db parameter
```

#### 5. Assembly failed: mount /lscratch error

**Problem**: Apptainer trying to mount non-existent paths

**Solution**: Already fixed in `metagenome_assembly_classification_en.config`:
```groovy
apptainer {
    autoMounts = true
    runOptions = '--bind /scratch:/scratch --bind /home:/home --containall'
}
```

#### 6. Zero consensus sequences (Sequence name mismatch bug)

**Problem**: Bug in v5.0.0 where VirSorter2 (`k141_12345||full`) and DeepVirFinder (`k141_12345 flag=...`) sequence names didn't match.

**Status**: ✅ **FIXED in v5.1.0**
- Sequence names now normalized before comparison
- See `BUG_FIX_SEQUENCE_NAME_MATCHING.md` for details

**To apply fix**:
```bash
# Option 1: Rerun merge/compare only (fast)
sbatch rerun_merge_and_compare.sh

# Option 2: Full rerun with resume
sbatch run_metagenome_assembly_classification_en.sh
```

### Performance Optimization

#### Speed up assembly
```bash
# Use MEGAHIT only (faster, less sensitive)
# Comment out SPAdes processes in workflow

# Or adjust resources
--megahit_memory 0.9  # Use more memory
```

#### Reduce computational load
```bash
# Increase minimum contig length
--virsorter2_min_length 2000
--deepvirfinder_min_length 2000

# More stringent thresholds
--virsorter2_min_score 0.7
--deepvirfinder_pvalue 0.01
```

### Getting Help

1. **Check logs**:
   ```bash
   cat Viral_Classification_*.out
   cat Viral_Classification_*.err
   ```

2. **Check Nextflow logs**:
   ```bash
   cat .nextflow.log
   ```

3. **Check individual process logs**:
   ```bash
   cd work/<hash>/<hash>/
   cat .command.out
   cat .command.err
   ```

---

## 📚 Downstream Analysis Examples

### Extract High-Confidence Viral Sequences

```bash
# Get sequence IDs
CONSENSUS_LIST="results/assembler_comparison/sample1_consensus_viral_sequences.txt"

# Extract from MEGAHIT assembly
grep -A 1 -Ff $CONSENSUS_LIST \
    results/assembly_megahit/sample1_megahit_contigs.fa \
    > high_confidence_viral_sequences.fa

# Or from merged VirSorter2 output
grep -A 1 -Ff $CONSENSUS_LIST \
    results/virsorter2_megahit/sample1_megahit_vs2_final-viral-combined.fa \
    > high_confidence_viral_sequences.fa
```

### Annotate Viral Sequences

```bash
# 1. Gene prediction
prodigal -i high_confidence_viral_sequences.fa \
    -a viral_proteins.faa \
    -d viral_genes.fna \
    -p meta

# 2. Functional annotation
diamond blastp -q viral_proteins.faa \
    -d /path/to/refseq_viral \
    -o viral_annotations.txt \
    --outfmt 6 --max-target-seqs 5

# 3. Taxonomic classification
# Using genomad or other viral taxonomy tools
```

### Statistical Analysis in R

```r
library(tidyverse)

# Load merged report
data <- read_csv("results/merged_viral_reports_megahit/sample1_megahit_viral_merged_report.csv")

# Visualize score distribution
ggplot(data, aes(x=vs2_score, fill=identified_by)) +
  geom_histogram(bins=50) +
  facet_wrap(~identified_by) +
  theme_minimal() +
  labs(title="VirSorter2 Score Distribution by Detection Method")

# Compare p-values
ggplot(data %>% filter(!is.na(dvf_pvalue)), 
       aes(x=-log10(dvf_pvalue), fill=consensus)) +
  geom_histogram(bins=50) +
  theme_minimal() +
  labs(title="DeepVirFinder Significance")
```

---

## 🏆 Best Practices

### Sample Preparation
1. **Viral enrichment**: Use virus-like particle (VLP) enrichment if possible
2. **Quality**: Ensure high-quality DNA extraction (avoid degradation)
3. **Sequencing depth**: ≥10M paired-end reads recommended

### Data Analysis
1. **Always check QC reports** (`fastp/*.html`)
2. **Prioritize consensus sequences** (assembler_comparison/*.txt)
3. **Validate unexpected results** manually (BLAST against NCBI)
4. **Keep assembly contigs** for re-analysis with different parameters

### Resource Management
1. **Use `-resume`** for interrupted runs
2. **Monitor memory usage** (adjust config if needed)
3. **Clean work directory** periodically:
   ```bash
   # After successful completion
   nextflow clean -f
   ```

### Version Control
1. **Record pipeline version** in your methods
2. **Save configuration files** with results
3. **Document parameter changes** from defaults

---

## 📖 Citations

### Tools Used in This Workflow

**VirSorter2**
```
Guo, J., Bolduc, B., Zayed, A. A., Varsani, A., Dominguez-Huerta, G., Delmont, T. O., ... & Sullivan, M. B. (2021).
VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses.
Microbiome, 9(1), 1-13.
```

**DeepVirFinder**
```
Ren, J., Ahlgren, N. A., Lu, Y. Y., Fuhrman, J. A., & Sun, F. (2017).
VirFinder: a novel k-mer based tool for identifying viral sequences from assembled metagenomic data.
Microbiome, 5(1), 1-20.

Ren, J., Song, K., Deng, C., Ahlgren, N. A., Fuhrman, J. A., Li, Y., ... & Sun, F. (2020).
Identifying viruses from metagenomic data using deep learning.
Quantitative Biology, 8(1), 64-77.
```

**MEGAHIT**
```
Li, D., Liu, C. M., Luo, R., Sadakane, K., & Lam, T. W. (2015).
MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.
Bioinformatics, 31(10), 1674-1676.
```

**SPAdes**
```
Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017).
metaSPAdes: a new versatile metagenomic assembler.
Genome research, 27(5), 824-834.
```

**fastp**
```
Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018).
fastp: an ultra-fast all-in-one FASTQ preprocessor.
Bioinformatics, 34(17), i884-i890.
```

**Nextflow**
```
Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017).
Nextflow enables reproducible computational workflows.
Nature biotechnology, 35(4), 316-319.
```

---

## 📝 Version History

### v5.1.0 (2025-10-28) - Current
- ✅ **Fixed**: Sequence name matching bug between VirSorter2 and DeepVirFinder
- ✅ **Added**: Sequence name normalization in merge scripts
- ✅ **Added**: Improved Python environment handling in rerun scripts
- ✅ **Updated**: Documentation with troubleshooting guide

### v5.0.0 (2025-10-27)
- Added assembler comparison step
- Added final consensus viral sequence list
- Added save_clean_reads parameter
- Removed Prodigal and Diamond steps (streamlined for viral focus)
- Enhanced result merging with pandas
- Improved error handling and logging

### v4.0.0 (2025-10-24)
- Initial release with dual assemblers
- Integrated VirSorter2 and DeepVirFinder
- Added comprehensive result merging
- Apptainer/Singularity support
- SLURM cluster compatibility

---

## 📧 Support

For issues, questions, or suggestions:
1. Check the [Troubleshooting](#troubleshooting) section
2. Review `BUG_FIX_SEQUENCE_NAME_MATCHING.md` for known issues
3. Check Nextflow logs (`.nextflow.log`)
4. Contact your system administrator for cluster-specific issues

---

## 📄 License

This workflow integrates multiple open-source tools. Please cite all tools used (see [Citations](#citations)).

---

## 🎯 Quick Reference Card

```bash
# Setup (one-time)
conda create -n nextflow_env -c bioconda nextflow fastp virsorter=2.2.4
conda create -n dvf python=3.6
virsorter setup -d /path/to/db

# Run
sbatch run_metagenome_assembly_classification_en.sh

# Resume
sbatch run_metagenome_assembly_classification_en.sh  # auto-resumes

# Monitor
tail -f Viral_Classification_*.out

# Results
cat results/assembler_comparison/*_assembler_comparison.txt
cat results/assembler_comparison/*_consensus_viral_sequences.txt
```

**Most Important Output**: 
`results/assembler_comparison/{sample}_consensus_viral_sequences.txt` ⭐⭐⭐

---

**Happy Viral Hunting! 🦠🔬**

