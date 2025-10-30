# Viral Metagenome Analysis Workflow Documentation

**Version**: 5.1.0  
**Last Updated**: 2025-10-30

---

## 📋 Workflow Overview

This workflow integrates two advanced viral identification tools, **VirSorter2** and **DeepVirFinder**, and provides the highest confidence viral sequence identification through **assembler comparison**.

### ✨ Key Features

- 🔬 **Dual-Tool Validation**: VirSorter2 (machine learning + rules) + DeepVirFinder (deep learning CNN)
- 🔄 **Dual Assemblers**: MEGAHIT (fast & efficient) + metaSPAdes (high quality) parallel assembly
- 🎯 **Three-Level Validation**: Single tool → Tool consensus → Assembler consensus
- 📊 **Automated Comparison**: Intelligent comparison of viral identification results from two assembly methods
- ⭐ **Final Consensus**: Generate high-confidence viral sequence list validated by four methods
- 💾 **Complete Preservation**: Save all intermediate results (clean reads, contigs, all analysis levels)

### Workflow Diagram

```
Raw Sequencing Data (FASTQ)
    ↓
1️⃣ Quality Control (fastp)
    ├─→ QC Reports
    └─→ Clean Reads (optional save)
    ↓
2️⃣ Parallel Metagenome Assembly
    ├── MEGAHIT (fast) → MEGAHIT contigs
    └── metaSPAdes (high quality) → SPAdes contigs
    ↓
3️⃣ Viral Sequence Identification (parallel)
    ├── VirSorter2 (MEGAHIT)
    ├── VirSorter2 (SPAdes)
    ├── DeepVirFinder (MEGAHIT)
    └── DeepVirFinder (SPAdes)
    ↓
4️⃣ Tool Result Integration (per assembler)
    ├── MEGAHIT: VirSorter2 ∩ DeepVirFinder
    └── SPAdes:  VirSorter2 ∩ DeepVirFinder
    ↓
5️⃣ Assembler Comparison ⭐ NEW!
    └── Consensus by existence (both assemblers detected viruses)
    ↓
🎯 Final High-Confidence Viral Sequence List
   (Four-way validation: 2 tools × 2 assemblers)
   
**Note**: Sequence IDs differ between assemblers (k141_XXX vs NODE_XXX), 
so "consensus" means both assemblers detected viruses, not identical sequence IDs.
```

---

## 🦠 Tool Description

### VirSorter2

| Item | Description |
|------|-------------|
| **Input** | Assembled contigs (FASTA format) |
| **Method** | Combines machine learning and rule-based approaches |
| **Database** | ~13GB comprehensive viral reference database |
| **Output** | `*_vs2_final-viral-score.tsv` (viral scores)<br>`*_vs2_final-viral-combined.fa` (viral contigs) |
| **Advantages** | • Identifies multiple viral types (dsDNA phage, NCLDV, RNA, ssDNA, etc.)<br>• Provides detailed confidence scores and viral type classification<br>• Based on large-scale reference database, excellent for known viruses |
| **Min Length** | Default 1000 bp |
| **Threshold** | max_score > 0.5 (adjustable) |

### DeepVirFinder

| Item | Description |
|------|-------------|
| **Input** | Assembled contigs (FASTA format) |
| **Method** | Convolutional Neural Network (CNN) for sequence feature learning |
| **Model** | Built-in pre-trained model (no additional download needed) |
| **Output** | `*_dvf_output.txt` (contains score and p-value) |
| **Advantages** | • No reference database required, suitable for novel/unknown viruses<br>• Provides statistical significance p-value<br>• Based on sequence patterns, independent of known viral features |
| **Min Length** | Default 1000 bp |
| **Threshold** | p-value < 0.05 (adjustable) |

### Three-Level Validation System 🎯

This workflow provides three levels of viral sequence validation:

#### Level 1: Single Tool Identification
- Location: `virsorter2_*/`, `deepvirfinder_*/`
- Confidence: ⭐⭐☆☆☆
- Description: All viral sequences identified by a single tool

#### Level 2: Tool Consensus (per assembler)
- Location: `merged_viral_reports_megahit/`, `merged_viral_reports_spades/`
- Output: `*_viral_consensus.txt`
- Confidence: ⭐⭐⭐⭐☆
- Description: Sequences identified by BOTH VirSorter2 AND DeepVirFinder

#### Level 3: Full Consensus (Recommended) ⭐
- Location: `assembler_comparison/`
- Output: `*_consensus_viral_sequences.txt`
- Confidence: ⭐⭐⭐⭐⭐
- Description: **Four-way validation by existence**
  - ✅ Both MEGAHIT and SPAdes assembled sequences from the sample
  - ✅ Both VirSorter2 and DeepVirFinder identified viruses from each assembler
  - ✅ **Consensus logic**: If both assemblers detected viruses, all high-confidence sequences are included
  - ⚠️ **Note**: Sequence IDs differ between assemblers (k141_XXX vs NODE_XXX), so consensus is based on both assemblers detecting viruses, not exact sequence ID matching

**💡 Recommendation**: Use Level 3 consensus sequences for downstream analysis to obtain the most reliable results!

---

## 🚀 Quick Start

### Prerequisites

#### 1. Conda Environment Setup

Two conda environments are required:

```bash
# Environment 1: nextflow_env (for running Nextflow and VirSorter2)
conda create -n nextflow_env python=3.9
conda activate nextflow_env
conda install -c bioconda nextflow virsorter=2.2.4 fastp screed -y

# Environment 2: dvf (for DeepVirFinder, Python 3.6)
# Note: See next step for detailed dependency installation
conda create -n dvf python=3.6
```

#### 2. DeepVirFinder Installation

```bash
# Activate dvf environment
conda activate dvf

# Install required dependencies (note version compatibility)
conda install -c conda-forge h5py=2.10.0 numpy -y
pip install theano==1.0.5
pip install keras==2.2.4

# Install other dependencies
conda install -c bioconda scikit-learn biopython pandas -y

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

# Download DeepVirFinder (if not already done)
cd /your/install/path/
git clone https://github.com/jessieren/DeepVirFinder.git
cd DeepVirFinder

# Test installation
python dvf.py -h
```

**⚠️ Important - DeepVirFinder Dependency Version Requirements**:

DeepVirFinder has **strict requirements** for dependency versions; using incorrect versions will cause failures:

| Dependency | Required Version | Reason |
|------------|-----------------|--------|
| **Python** | 3.6 | DeepVirFinder only supports Python 3.6 |
| **h5py** | 2.10.0 | h5py 3.x returns str instead of bytes, incompatible with Keras 2.2.4 |
| **keras** | 2.2.4 | keras 2.6+ uses TensorFlow 2.x format, incompatible with model files |
| **theano** | 1.0.5 | As Keras backend, cannot use TensorFlow |

**Common Failing Version Combinations**:
- ❌ h5py >= 3.0 + keras 2.2.4 → `AttributeError: 'str' object has no attribute 'decode'`
- ❌ keras >= 2.6 + TensorFlow 2.x → `TypeError: Parameter to MergeFrom()`
- ❌ Python 3.7+ → Various compatibility issues

**Recommended Installation Order**:
1. Create Python 3.6 environment
2. Install h5py=2.10.0
3. Install theano=1.0.5
4. Install keras=2.2.4
5. Configure Keras backend to theano
6. Verify all dependencies

#### 3. Database Setup

##### VirSorter2 Database (required, ~13GB)

```bash
# Activate nextflow_env environment
conda activate nextflow_env

# Download and setup database
virsorter setup -d /path/to/virsorter2/db -j 4

# Wait for download to complete (~13GB)
```

Or use the provided script:
```bash
bash setup_virsorter2_db.sh
```

---

### Configuration and Execution

#### Step 1: Prepare Sample Sheet

Create `samplesheet.csv` file:

```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
sample3,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz
```

#### Step 2: Configure Paths

Edit `run_metagenome_assembly_classification_en.sh`:

```bash
# VirSorter2 database path (required)
VIRSORTER2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/virsorter2/db"

# DeepVirFinder installation path (required)
DEEPVIRFINDER_DIR="/path/to/DeepVirFinder"
```

#### Step 3: Submit Job

```bash
# Submit to SLURM (default: full run from scratch, no -resume)
sbatch run_metagenome_assembly_classification_en.sh

# Or run locally (for testing)
bash run_metagenome_assembly_classification_en.sh
```

**About -resume parameter (optional):**
- The script does not enable `-resume` by default and performs a full run from scratch.
- To skip completed steps, manually add the `-resume` flag on the command line.
- With `-resume`, only failed or incomplete steps are re-run, saving time.

#### Step 4: Monitor Execution

```bash
# Check job status
squeue -u $USER

# View real-time log
tail -f Viral_Classification_*.out

# View Nextflow log
tail -f .nextflow.log
```

---

## 📊 Output Results

### Complete Directory Structure

```
results/
├── fastp/                          # Quality control reports
│   ├── *_fastp.html                # HTML quality reports
│   └── *_fastp.json                # JSON quality data
│
├── clean_reads/                    # Quality-controlled sequences (optional)
│   ├── *_clean_R1.fastq.gz         # Forward clean reads
│   └── *_clean_R2.fastq.gz         # Reverse clean reads
│
├── assembly_megahit/               # MEGAHIT assembly results
│   └── *_megahit_contigs.fa        # Assembled contigs
│
├── assembly_spades/                # SPAdes assembly results
│   └── *_spades_contigs.fa         # Assembled contigs
│
├── virsorter2_megahit/             # VirSorter2 results (MEGAHIT)
│   ├── *_megahit_vs2_final-viral-score.tsv
│   └── *_megahit_vs2_final-viral-combined.fa
│
├── virsorter2_spades/              # VirSorter2 results (SPAdes)
│   ├── *_spades_vs2_final-viral-score.tsv
│   └── *_spades_vs2_final-viral-combined.fa
│
├── deepvirfinder_megahit/          # DeepVirFinder results (MEGAHIT)
│   └── *_megahit_dvf_output.txt
│
├── deepvirfinder_spades/           # DeepVirFinder results (SPAdes)
│   └── *_spades_dvf_output.txt
│
├── merged_viral_reports_megahit/   # Tool integration (MEGAHIT) ⭐⭐⭐⭐
│   ├── *_megahit_viral_merged_report.txt
│   ├── *_megahit_viral_merged_report.csv
│   └── *_megahit_viral_consensus.txt
│
├── merged_viral_reports_spades/    # Tool integration (SPAdes) ⭐⭐⭐⭐
│   ├── *_spades_viral_merged_report.txt
│   ├── *_spades_viral_merged_report.csv
│   └── *_spades_viral_consensus.txt
│
└── assembler_comparison/           # Assembler comparison ⭐⭐⭐⭐⭐
    ├── *_assembler_comparison.txt      # Comprehensive comparison report
    ├── *_assembler_comparison.csv      # Detailed comparison data
    └── *_consensus_viral_sequences.txt # Final recommended list 🏆
```

### Key Output Files Explained

#### 1. VirSorter2 Output (`*_vs2_final-viral-score.tsv`)

| Column | Description |
|--------|-------------|
| seqname | Contig name |
| max_score | Maximum viral score (0-1) |
| max_score_group | Viral type (dsDNAphage, ssDNA, etc.) |
| length | Sequence length |
| hallmark | Number of hallmark genes |
| viral_gene | Number of viral genes |
| cellular_gene | Number of host genes |

**Interpretation**:
- `max_score > 0.9`: High confidence viral
- `max_score > 0.7`: Medium confidence
- `max_score > 0.5`: Low confidence

#### 2. DeepVirFinder Output (`*_dvf_output.txt`)

| Column | Description |
|--------|-------------|
| name | Contig name |
| len | Sequence length |
| score | Viral probability score (0-1) |
| pvalue | Statistical significance |

**Interpretation**:
- `pvalue < 0.01` and `score > 0.9`: High confidence
- `pvalue < 0.05` and `score > 0.7`: Medium confidence

#### 3. Tool Integration Report (`*_viral_merged_report.txt`) ⭐⭐⭐⭐

Contains:
- **Overall Statistics**: Number of viruses identified by each tool
- **Consensus Sequences**: Sequences identified by both tools (recommended)
- **Single Tool**: Sequences identified by only one tool
- **Detailed Comparison**: Dual-tool score comparison for each sequence

#### 4. Assembler Comparison Report (`*_assembler_comparison.txt`) ⭐⭐⭐⭐⭐

**This is the most important output!**

Contains:
- **Overall Statistics**:
  - Number of viruses identified by MEGAHIT
  - Number of viruses identified by SPAdes
  - Number of consensus viruses (highest confidence)
  - Assembler consistency percentage
  
- **Detailed Analysis**:
  - Sequences identified only by MEGAHIT
  - Sequences identified only by SPAdes
  - Consensus sequences from both
  
- **Recommendations**:
  - Prioritize consensus sequences for downstream analysis

Example output:
```
====================================================================================================
Assembler Comparison Report - Viral Identification Results
MEGAHIT vs metaSPAdes
Sample: sample1
====================================================================================================

[Overall Statistics]
----------------------------------------------------------------------------------------------------
MEGAHIT identified viral sequences:    150
SPAdes identified viral sequences:     180
Total unique viral sequences:          330

Consensus viral sequences (both assemblies detected): 330
MEGAHIT viral sequences:                150
SPAdes viral sequences:                180

**Note**: Consensus here means both assemblers detected viruses (existence consensus).
Sequence IDs differ between assemblers (k141_XXX vs NODE_XXX), so all high-confidence 
sequences from both assemblers are included in the consensus.

====================================================================================================
[Recommendation]
----------------------------------------------------------------------------------------------------
High-confidence viral sequences (identified by both assemblers): 330
Recommend prioritizing these consensus sequences for downstream analysis.
```

#### 5. Final Consensus Sequence List (`*_consensus_viral_sequences.txt`) 🏆

**This is the most important file!**

Format:
```
# High-confidence viral sequences from both MEGAHIT and SPAdes
# Both assemblers detected viruses (consensus by existence)
# Sample: sample1
# Total sequences: 330
# Note: Sequence IDs differ between assemblers, consensus means both detected viruses
#
k141_123456
NODE_1003_length_1759_cov_6.676643
k141_234567
NODE_1029_length_1742_cov_1.900415
...
```

**Important Note**: 
- Sequence IDs from MEGAHIT (k141_XXX) and SPAdes (NODE_XXX_length_XXX_cov_XXX) are different
- "Consensus" means both assemblers detected viruses, not that sequence IDs match
- All sequences from both assemblers are included if both detected viruses

**Uses**:
- Downstream viral genome analysis
- Extract viral sequences for annotation
- Construct viral phylogenetic trees
- Perform viral function prediction

---

## 🔬 Downstream Analysis Examples

### 1. Extract Final Consensus Viral Sequences

```bash
# Extract sequences using final consensus list (recommended)
seqkit grep -f results/assembler_comparison/sample1_consensus_viral_sequences.txt \
    results/assembly_megahit/sample1_megahit_contigs.fa \
    > sample1_high_confidence_viruses.fa

# Count sequences
echo "Total high-confidence viral sequences:"
grep -c ">" sample1_high_confidence_viruses.fa
```

### 2. Viral Sequence Annotation

```bash
# Use Prokka for gene prediction and annotation
prokka --outdir sample1_annotation \
       --prefix sample1_viruses \
       --metagenome \
       sample1_high_confidence_viruses.fa

# Or use DRAM-v for viral functional annotation
DRAM-v.py annotate -i sample1_high_confidence_viruses.fa \
                    -o sample1_dram_annotation
```

### 3. Viral Classification (optional)

```bash
# Use vConTACT2 for viral classification
vcontact2 --raw-proteins proteins.faa \
          --rel-mode 'Diamond' \
          --proteins-fp gene_to_genome.csv \
          --db 'ProkaryoticViralRefSeq94-Merged' \
          --output-dir vcontact2_results
```

### 4. Viral Abundance Analysis

```bash
# Map clean reads to viral contigs
bowtie2-build sample1_high_confidence_viruses.fa viral_index

bowtie2 -x viral_index \
        -1 results/clean_reads/sample1_clean_R1.fastq.gz \
        -2 results/clean_reads/sample1_clean_R2.fastq.gz \
        -S sample1_viral_mapping.sam \
        -p 8

# Calculate abundance
samtools view -bS sample1_viral_mapping.sam | samtools sort -o sample1_viral.bam
samtools index sample1_viral.bam
samtools idxstats sample1_viral.bam > sample1_viral_abundance.txt
```

---

## ⚙️ Parameter Adjustment

### Common Parameter Modifications

In `run_metagenome_assembly_classification_en.sh` or directly in command line:

```bash
nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet.csv \
    --outdir results \
    --virsorter2_db /path/to/virsorter2/db \
    --deepvirfinder_dir /path/to/DeepVirFinder \
    \
    # Quality control parameters
    --save_clean_reads true \              # Whether to save clean reads
    \
    # Assembly parameters
    --megahit_min_contig_len 1000 \        # MEGAHIT minimum contig length
    --megahit_memory 0.8 \                 # MEGAHIT memory usage ratio
    --spades_memory 512 \                  # SPAdes memory (GB)
    \
    # VirSorter2 parameters
    --virsorter2_min_length 1000 \         # Minimum sequence length
    --virsorter2_min_score 0.5 \           # Minimum viral score
    \
    # DeepVirFinder parameters
    --deepvirfinder_min_length 1000 \      # Minimum sequence length
    --deepvirfinder_pvalue 0.05            # p-value threshold
```

### Skip Specific Steps

```bash
# Skip fastp QC (if data already QC'd)
--skip_fastp

# Skip VirSorter2
--skip_virsorter2

# Skip DeepVirFinder
--skip_deepvirfinder

# Skip result merging
--skip_merge_reports
```

### Resource Configuration

Edit `metagenome_assembly_classification_en.config`:

```groovy
process {
    // Adjust MEGAHIT resources
    withName: 'MEGAHIT_ASSEMBLY' {
        cpus = 16
        memory = '128 GB'
        time = '24h'
    }
    
    // Adjust SPAdes resources
    withName: 'SPADES_ASSEMBLY' {
        cpus = 32
        memory = '512 GB'
        time = '48h'
    }
    
    // Adjust DeepVirFinder resources
    withName: 'DEEPVIRFINDER_*' {
        cpus = 8
        memory = '32 GB'
        time = '12h'
    }
}
```

---

## 📈 Performance Optimization Recommendations

### 1. Data Preprocessing

```bash
# If data is large, subsample for testing first
seqtk sample -s100 raw_R1.fastq.gz 1000000 > test_R1.fastq.gz
seqtk sample -s100 raw_R2.fastq.gz 1000000 > test_R2.fastq.gz
```

### 2. Parallelization

```bash
# Fully utilize multi-core CPU
# Modify cpus parameter in config file
# MEGAHIT: 8-16 cores
# SPAdes: 16-32 cores
# DeepVirFinder: 4-8 cores
```

### 3. Memory Optimization

```bash
# MEGAHIT memory setting (fraction of available)
--megahit_memory 0.8  # Use 80% available memory

# SPAdes memory setting (absolute value in GB)
--spades_memory 256   # Use 256GB
```

### 4. Use -resume (Optional)

```bash
# Default run (from scratch, no resume)
sbatch run_metagenome_assembly_classification_en.sh

# Enable resume manually if desired
nextflow run metagenome_assembly_classification_workflow_en.nf \
  -c metagenome_assembly_classification_en.config \
  --input samplesheet.csv \
  --outdir results \
  --virsorter2_db /path/to/virsorter2/db \
  --deepvirfinder_dir /path/to/DeepVirFinder \
  -resume

# If complete rerun needed (clear cache first):
nextflow clean -f
rm -rf work/
sbatch run_metagenome_assembly_classification_en.sh
```

---

## 📚 References

### Main Tools

1. **VirSorter2**
   - Guo, J., et al. (2021). VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome*, 9(1), 1-13.
   - DOI: [10.1186/s40168-020-00990-y](https://doi.org/10.1186/s40168-020-00990-y)
   - GitHub: https://github.com/jiarong/VirSorter2

2. **DeepVirFinder**
   - Ren, J., et al. (2020). Identifying viruses from metagenomic data using deep learning. *Quantitative Biology*, 8(1), 64-77.
   - DOI: [10.1007/s40484-019-0187-4](https://doi.org/10.1007/s40484-019-0187-4)
   - GitHub: https://github.com/jessieren/DeepVirFinder

### Assembly Tools

3. **MEGAHIT**
   - Li, D., et al. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31(10), 1674-1676.

4. **metaSPAdes**
   - Nurk, S., et al. (2017). metaSPAdes: a new versatile metagenomic assembler. *Genome research*, 27(5), 824-834.

### Workflow Management

5. **Nextflow**
   - Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature biotechnology*, 35(4), 316-319.

---

## 👥 Technical Support

### Log File Locations

| File | Location | Content |
|------|----------|---------|
| Nextflow main log | `.nextflow.log` | Overall workflow execution log |
| SLURM output | `Viral_Classification_*.out` | Standard output |
| SLURM error | `Viral_Classification_*.err` | Error messages |
| Process logs | `work/xx/xxxxx/.command.*` | Detailed logs for individual processes |

### Debugging Steps

1. **View overall process**:
   ```bash
   cat Viral_Classification_*.out | grep "ERROR\|✅\|❌"
   ```

2. **Check failed processes**:
   ```bash
   # Find work directory from log
   cd work/xx/xxxxxx/
   cat .command.err  # View errors
   cat .command.out  # View output
   bash .command.sh  # Manually reproduce issue
   ```

3. **Verify environment**:
   ```bash
   conda activate nextflow_env
   virsorter -h
   
   conda activate dvf
   python -c "import h5py, keras, theano"
   ```

### Contact Information

If you encounter issues, please provide the following information:
- Nextflow version: `nextflow -version`
- Error logs: relevant parts of `.nextflow.log`
- SLURM logs: `Viral_Classification_*.out`
- System info: memory, number of CPU cores

---

## 🔄 Version History

### v5.1.0 (2025-10-30) - Current

**New Features**:
- ✨ Added assembler comparison module (MEGAHIT vs SPAdes)
- ✨ Generate final consensus viral sequence list (four-way validation)
- ✨ Optional save clean reads and assembly contigs
- ✨ Three-level validation system
- ✨ Optional `-resume` support (script defaults to full run from scratch)

**Improvements**:
- 🔧 Optimize DeepVirFinder environment activation mechanism (absolute path + explicit PATH setting)
- 🔧 Fix merge report process dependency issues (auto-install pandas and numpy)
- 🔧 Resolve PYTHONPATH pollution issues (clean and reset environment variables)
- 🔧 Improve Keras backend configuration (force use of Theano)
- 🔧 Update h5py compatibility (downgrade to 2.10.0)
- 🔧 Improve error handling and log output
- 📊 Enhance report readability

### v5.0.0 (2025-10-24)

**New Features**:
- ✨ Integrate VirSorter2 and DeepVirFinder
- ✨ Dual assembler parallel analysis
- ✨ Automated result integration

### v4.0.0 (2025-10)

**Initial Version**:
- Basic metagenome assembly and classification pipeline

---

## 📝 Best Practices

### 1. Environment Configuration

- ✅ **Strictly follow documentation for DeepVirFinder dependencies** (strict version requirements)
- ✅ Verify all dependency versions: h5py=2.10.0, keras=2.2.4, theano=1.0.5
- ✅ Ensure Keras backend configured for theano (not TensorFlow)
- ✅ Test environment: run small dataset to ensure all tools work properly

### 2. Data Quality Control

- ✅ Use high-quality raw data (Phred score > 30)
- ✅ Remove adapters and low-quality reads
- ✅ Check for host DNA contamination
- ✅ Test workflow with subset if data is large

### 3. Parameter Selection

- ✅ Conservative thresholds: VirSorter2 max_score > 0.7, DeepVirFinder p-value < 0.01
- ✅ Permissive thresholds: VirSorter2 max_score > 0.5, DeepVirFinder p-value < 0.05
- ✅ **Recommend using final consensus sequence list** (four-way validation)
- ✅ Adjust minimum contig length based on sample type

### 4. Workflow Execution

- ✅ **Use -resume parameter** (enabled by default in script)
- ✅ Regularly check log files to monitor progress
- ✅ If errors occur, check specific process `.command.err` files
- ✅ Keep work directory until results are confirmed correct

### 5. Result Validation

- ✅ Check viral length distribution
- ✅ Compare consistency between MEGAHIT and SPAdes results
- ✅ Manually validate high-confidence sequences (BLAST comparison)
- ✅ Check consensus rate (typically should be > 50%)

### 6. Downstream Analysis

- ✅ Prioritize using `assembler_comparison/*_consensus_viral_sequences.txt`
- ✅ Perform functional annotation and classification of viruses
- ✅ Calculate relative abundance of viruses in samples
- ✅ Save all intermediate files for reanalysis

---

**Thank you for using this workflow! Feedback on issues or suggestions is welcome.** 🎉




