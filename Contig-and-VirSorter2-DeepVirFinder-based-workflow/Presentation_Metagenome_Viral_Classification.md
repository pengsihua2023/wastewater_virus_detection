% Metagenome Viral Classification Workflow
% Presentation Slides
% Based on README_EN.md

---
title: "Metagenome Viral Classification Workflow"
subtitle: "Integration of VirSorter2 and DeepVirFinder with Assembler Comparison"
author: "Workflow Documentation v5.1.0"
date: "2025-10-28"
---

# Slide 1: Title Slide

## 🦠 Metagenome Viral Classification Workflow

### Integration of VirSorter2 and DeepVirFinder

**Version 5.1.0**

- Dual-Tool Validation
- Dual-Assembler Comparison  
- High-Confidence Viral Sequence Identification

---

# Slide 2: Overview

## 📋 Workflow Overview

### Key Features

- 🔬 **Dual-Tool Validation**: VirSorter2 (ML + rules) + DeepVirFinder (Deep Learning)
- 🔄 **Dual Assemblers**: MEGAHIT (fast) + metaSPAdes (high quality)
- 🎯 **Three-Level Validation**: Single tool → Tool consensus → Assembler consensus
- ⭐ **Final Consensus**: Four-way validation (2 tools × 2 assemblers)

### Main Goal

Generate **high-confidence viral sequence list** validated by multiple independent methods

---

# Slide 3: Workflow Architecture

## 🔄 Workflow Pipeline

```
1️⃣ Quality Control (fastp)
    ↓
2️⃣ Parallel Assembly
    ├── MEGAHIT
    └── metaSPAdes
    ↓
3️⃣ Viral Identification (Parallel)
    ├── VirSorter2 (MEGAHIT & SPAdes)
    └── DeepVirFinder (MEGAHIT & SPAdes)
    ↓
4️⃣ Tool Integration (per assembler)
    ├── VirSorter2 ∩ DeepVirFinder (MEGAHIT)
    └── VirSorter2 ∩ DeepVirFinder (SPAdes)
    ↓
5️⃣ Assembler Comparison ⭐
    └── Consensus by existence
    ↓
🎯 Final Consensus Viral Sequences
```

---

# Slide 4: Tool 1 - VirSorter2

## 🔬 VirSorter2

### Characteristics

| Feature | Description |
|---------|-------------|
| **Method** | Machine Learning + Rule-based |
| **Database** | ~13GB viral reference database |
| **Min Length** | 1000 bp (default) |
| **Threshold** | max_score > 0.5 (adjustable) |
| **Advantages** | • Multiple viral types (dsDNA, ssDNA, RNA, NCLDV)<br>• Detailed confidence scores<br>• Excellent for known viruses |

### Output

- Viral scores (`*_vs2_final-viral-score.tsv`)
- Identified viral contigs (`*_vs2_final-viral-combined.fa`)

---

# Slide 5: Tool 2 - DeepVirFinder

## 🔬 DeepVirFinder

### Characteristics

| Feature | Description |
|---------|-------------|
| **Method** | Convolutional Neural Network (CNN) |
| **Model** | Pre-trained model (no database) |
| **Min Length** | 1000 bp (default) |
| **Threshold** | p-value < 0.05 (adjustable) |
| **Advantages** | • No reference database needed<br>• Suitable for novel/unknown viruses<br>• Statistical p-value significance |

### Output

- Prediction results with scores and p-values (`*_dvf_output.txt`)

---

# Slide 6: Three-Level Validation System

## 🎯 Validation Levels

### Level 1: Single Tool ⭐⭐☆☆☆
- **Location**: `virsorter2_*/`, `deepvirfinder_*/`
- **Confidence**: Low
- **Description**: Identified by one tool

### Level 2: Tool Consensus ⭐⭐⭐⭐☆
- **Location**: `merged_viral_reports_*/`
- **Confidence**: High
- **Description**: Identified by **BOTH** VirSorter2 **AND** DeepVirFinder

### Level 3: Full Consensus ⭐⭐⭐⭐⭐
- **Location**: `assembler_comparison/`
- **Confidence**: **Highest**
- **Description**: Four-way validation by existence

---

# Slide 7: Level 3 - Full Consensus

## ⭐ Level 3: Full Consensus (Recommended)

### Four-Way Validation

✅ **MEGAHIT** assembled sequences  
✅ **SPAdes** assembled sequences  
✅ **VirSorter2** identified viruses  
✅ **DeepVirFinder** identified viruses  

### Consensus Logic

- **If both assemblers detected viruses** → Output all high-confidence sequences
- **Sequence IDs differ** between assemblers (k141_XXX vs NODE_XXX)
- **Consensus = existence**, not exact sequence ID matching

### Result

Highest confidence viral sequence list for downstream analysis

---

# Slide 8: Prerequisites

## 🚀 Prerequisites

### 1. Conda Environments

**Environment 1: nextflow_env**
- Python 3.9
- Nextflow, VirSorter2, fastp

**Environment 2: dvf**  
- Python 3.6 ⚠️ (strict requirement)
- DeepVirFinder dependencies:
  - h5py=2.10.0
  - keras=2.2.4
  - theano=1.0.5

### 2. Databases

- VirSorter2 database (~13GB)
- DeepVirFinder: Built-in model

---

# Slide 9: Workflow Execution

## ⚙️ How to Run

### Step 1: Prepare Sample Sheet

```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
```

### Step 2: Configure Paths

Edit script with database paths:
- VirSorter2 database
- DeepVirFinder directory

### Step 3: Submit Job

```bash
sbatch run_metagenome_assembly_classification_en.sh
```

**Note**: `-resume` enabled by default (skips completed steps)

---

# Slide 10: Output Directory Structure

## 📊 Results Organization

```
results/
├── fastp/                          # Quality control
├── assembly_megahit/               # MEGAHIT contigs
├── assembly_spades/                # SPAdes contigs
├── virsorter2_megahit/             # VirSorter2 results
├── virsorter2_spades/
├── deepvirfinder_megahit/          # DeepVirFinder results
├── deepvirfinder_spades/
├── merged_viral_reports_megahit/   # Level 2 consensus ⭐⭐⭐⭐
├── merged_viral_reports_spades/
└── assembler_comparison/           # Level 3 consensus ⭐⭐⭐⭐⭐
    └── *_consensus_viral_sequences.txt (Final result!)
```

---

# Slide 11: Key Output Files

## 📄 Important Output Files

### Level 2: Tool Consensus (per assembler)

**`*_viral_merged_report.txt`**
- Overall statistics
- Consensus sequences (Both tools)
- Tool-specific sequences

**`*_viral_consensus.txt`**
- High-confidence sequences (Both tools agree)

### Level 3: Full Consensus ⭐

**`*_consensus_viral_sequences.txt`** 🏆
- **Most important file!**
- Final high-confidence viral sequences
- Four-way validation

---

# Slide 12: Example Results - Statistics

## 📈 Example Output Statistics

### Merge Report (Level 2)
```
VirSorter2 identified:    1,595 sequences
DeepVirFinder identified: 48,295 sequences
Consensus (Both tools):   1,595 sequences ⭐⭐⭐⭐
```

### Assembler Comparison (Level 3) ⭐
```
MEGAHIT viral sequences:  1,595
SPAdes viral sequences:   236
Total consensus:          1,831 ⭐⭐⭐⭐⭐

Consensus = Both assemblers detected viruses
All high-confidence sequences included
```

---

# Slide 13: Output File Format

## 📝 Consensus Sequence List Format

```
# High-confidence viral sequences from both MEGAHIT and SPAdes
# Both assemblers detected viruses (consensus by existence)
# Sample: sample1
# Total sequences: 1,831
#
k141_123456
NODE_1003_length_1759_cov_6.676643
k141_234567
...
```

### Key Points

- ✅ Contains all sequences when both assemblers detected viruses
- ✅ Sequence IDs differ (MEGAHIT: k141_XXX, SPAdes: NODE_XXX)
- ✅ Consensus = existence, not ID matching

---

# Slide 14: Downstream Analysis

## 🔬 Applications

### 1. Extract Viral Sequences

```bash
seqkit grep -f consensus_viral_sequences.txt \
    assembly_megahit/contigs.fa > viruses.fa
```

### 2. Functional Annotation

- Prokka (gene prediction)
- DRAM-v (viral functions)

### 3. Viral Classification

- vConTACT2 (viral clusters)

### 4. Abundance Analysis

- Map reads to viral contigs
- Calculate viral abundance

---

# Slide 15: Parameter Adjustment

## ⚙️ Customizable Parameters

### VirSorter2
- `--virsorter2_min_score`: 0.5 (default)
  - Lower = more sequences, may include false positives
  - Higher = fewer sequences, higher confidence

### DeepVirFinder
- `--deepvirfinder_pvalue`: 0.05 (default)
  - Lower = stricter (0.01)
  - Higher = more permissive (0.1, 0.2)

### Assembly
- `--megahit_min_contig_len`: 1000 bp
- `--megahit_memory`: 0.8 (fraction)

---

# Slide 16: Best Practices

## ✅ Recommendations

### 1. Threshold Selection
- **Conservative**: VirSorter2 > 0.7, DeepVirFinder p < 0.01
- **Default**: VirSorter2 > 0.5, DeepVirFinder p < 0.05
- **Use Level 3 consensus** for final analysis

### 2. Workflow Execution
- ✅ Use `-resume` (enabled by default)
- ✅ Monitor logs regularly
- ✅ Check intermediate results

### 3. Result Validation
- ✅ Compare MEGAHIT vs SPAdes consistency
- ✅ Check consensus rate (typically > 50%)
- ✅ Validate high-confidence sequences

---

# Slide 17: Common Issues & Solutions

## 🐛 Troubleshooting

### Issue 1: Consensus Sequences = 0

**Cause**: Cache issues or old logic  
**Solution**: 
- Delete `work/` directory
- Rerun workflow
- Use `generate_correct_comparison.py` script

### Issue 2: DeepVirFinder Dependencies

**Cause**: Version incompatibility  
**Solution**:
- Python 3.6 (strict)
- h5py=2.10.0, keras=2.2.4, theano=1.0.5
- Keras backend = theano

### Issue 3: Memory Issues

**Solution**: Increase SLURM memory allocation

---

# Slide 18: Performance & Optimization

## 📈 Performance Tips

### Resource Configuration
- **MEGAHIT**: 8-16 CPUs, 128GB RAM
- **SPAdes**: 16-32 CPUs, 512GB RAM
- **DeepVirFinder**: 4-8 CPUs per job

### Time Optimization
- ✅ Use `-resume` (skip completed steps)
- ✅ Parallel execution (4 tools run simultaneously)
- ✅ Caching intermediate results

### Typical Runtime
- Quality Control: ~30 min
- Assembly: 2-4 hours each (parallel)
- Viral ID: 1-3 hours each (parallel)
- Merging/Comparison: < 5 minutes

---

# Slide 19: Key Advantages

## 🎯 Why This Workflow?

### 1. Multi-Tool Validation
- Reduces false positives
- Increases confidence

### 2. Multi-Assembler Comparison
- Captures different assembly strategies
- More comprehensive viral discovery

### 3. Automated Pipeline
- No manual intervention
- Reproducible results

### 4. Three-Level Confidence
- Choose appropriate level for your analysis
- Level 3 = highest confidence

---

# Slide 20: Summary & Contact

## 🎉 Summary

### What This Workflow Provides

✅ **Quality-controlled** metagenome data  
✅ **Dual-assembler** parallel assembly  
✅ **Dual-tool** viral identification  
✅ **Integrated** results with confidence levels  
✅ **Final consensus** list for downstream analysis  

### Recommended Output

**`assembler_comparison/*_consensus_viral_sequences.txt`**

Four-way validated, high-confidence viral sequences

---

## 📚 References

- VirSorter2: Guo et al. (2021), *Microbiome*
- DeepVirFinder: Ren et al. (2020), *Quantitative Biology*
- MEGAHIT: Li et al. (2015), *Bioinformatics*
- metaSPAdes: Nurk et al. (2017), *Genome Research*

---

**Thank You!**

Questions & Feedback Welcome

