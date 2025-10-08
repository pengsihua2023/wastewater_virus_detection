# Kraken2 Database Results Integration Tool

## 📖 Overview

This tool integrates and compares Kraken2 classification results from two viral databases (RVDB and NCBI Viral RefSeq) for comprehensive viral metagenomic analysis. It provides confidence-level stratification, statistical analysis, and visualization-ready output.

**Version**: 2.0  
**Sample**: 66ce4dde  
**Author**: Cursor AI Assistant  
**Date**: October 2025

## 🎯 Key Features

### 1. Smart Classification Counting
- ✅ **Only counts direct read assignments** - Avoids parent-level duplication
- ✅ **Filters out aggregated counts** - Focuses on actual classifications
- ✅ **Accurate species enumeration** - No Kingdom/Phylum/Class redundancy

### 2. Confidence Stratification
Automatically categorizes detections using a two-dimensional matrix:

| Detection Condition | Reads Range | Confidence Level | Recommended Action |
|-------------------|-------------|------------------|-------------------|
| **Both databases** | ≥100 | **High** | Direct reporting |
| **Both databases** | 50-99 | **Medium-High** | Secondary reporting |
| **Single database** | ≥100 | **Medium** | Validation required |
| **Single database** | 50-99 | **Low-Medium** | Strong validation needed |
| **Any condition** | <50 | **Low** | Not recommended |

### 3. Comprehensive Output
Generates 6 detailed output files for analysis and reporting.

## 📁 Input Files

The script expects two Kraken2 report files:

```
66ce4dde_kraken2_RVDB.txt
66ce4dde_kraken2_NCBI.txt
```

**Kraken2 Report Format:**
```
Percentage  Total_reads  Direct_reads  Rank  TaxID  Scientific_name
10.50       1000         100           S     12345  Virus species A
```

## 🚀 Quick Start

### Prerequisites

```bash
# Python 3.6+
python3 --version

# Required package
pip install pandas
# or on cluster
pip install --user pandas
```

### Basic Usage

```bash
# Simple run (auto-detects files)
python3 integrate_66ce4dde_EN.py

# Explicit file specification
python3 integrate_66ce4dde_EN.py \
  66ce4dde_kraken2_RVDB.txt \
  66ce4dde_kraken2_NCBI.txt

# Custom thresholds
python3 integrate_66ce4dde_EN.py \
  66ce4dde_kraken2_RVDB.txt \
  66ce4dde_kraken2_NCBI.txt \
  150 75 0

# Parameters:
#   150 = High confidence threshold (default: 100)
#   75  = Medium confidence threshold (default: 50)
#   0   = Minimum direct reads (default: 0)
```

## 📊 Output Files

### 1. Complete Integration Results
**File**: `66ce4dde_integrated_full.tsv`

Contains all detected taxa with comprehensive information:
- Taxon names and ranks
- Read counts (direct and total) for both databases
- Confidence levels and descriptions
- Presence indicators (In_RVDB, In_NCBI, In_Both)

**Use case**: Detailed analysis, supplementary materials

### 2. High Confidence Species
**File**: `66ce4dde_integrated_high_confidence_species.tsv`

Species-level classifications detected in both databases with sufficient reads.

**Columns**:
- Taxon_Name
- RVDB_Reads_Direct, RVDB_Percent
- NCBI_Reads_Direct, NCBI_Percent
- Avg_Percent

**Use case**: Main results for publication

### 3. All-Level High Confidence
**File**: `66ce4dde_integrated_high_confidence_all.tsv`

High confidence detections at all taxonomic levels (Species, Genus, Family, etc.).

**Use case**: Taxonomic diversity analysis

### 4. Candidate Viruses
**File**: `66ce4dde_integrated_candidate_viruses.tsv`

Species detected in only one database but with high read counts (≥100).

**Columns**:
- Taxon_Name
- In_RVDB, In_NCBI
- Max_Direct_Reads, Max_Percent
- Description

**Use case**: Novel or rare virus discovery, validation targets

### 5. Venn Diagram Data
**File**: `66ce4dde_integrated_venn_data.txt`

Detailed breakdown of database overlap at species level:
- Statistics summary
- Intersection species list with read counts
- RVDB-unique species (top 20)
- NCBI-unique species (top 20)

**Use case**: Database comparison visualization

### 6. Summary Report
**File**: `66ce4dde_integrated_summary.txt`

Human-readable summary including:
- Database information
- Integration parameters
- Statistical overview
- Top 10 high confidence viruses
- File list

**Use case**: Quick overview, presentation slides

## 📈 Understanding the Results

### Confidence Level Interpretation

#### High Confidence (Primary findings)
```
✓ Detected in both RVDB and NCBI
✓ Direct reads ≥ 100
✓ Core virome component
→ Action: Direct reporting in main text
```

#### Medium-High Confidence (Secondary findings)
```
✓ Detected in both RVDB and NCBI
✓ Direct reads 50-99
✓ Low-abundance true virus
→ Action: Report as secondary findings
```

#### Medium Confidence (Candidate findings)
```
! Detected in only one database
! Direct reads ≥ 100
! Possibly rare or database-specific virus
→ Action: Validation required (BLAST, assembly)
```

#### Low-Medium Confidence (Borderline findings)
```
! Detected in only one database
! Direct reads 50-99
! Uncertain authenticity
→ Action: Strong validation needed
```

#### Low Confidence (Not recommended)
```
✗ Any detection condition
✗ Direct reads < 50
✗ Likely false positive or noise
→ Action: Do not report unless specifically justified
```

### Expected Results Pattern

**Clinical Sample:**
- High confidence: 5-20 species
- Medium-High confidence: 3-10 species
- Medium confidence: 5-15 species
- Low-Medium confidence: 10-25 species
- Low confidence: 20-50+ species

**Environmental Sample:**
- High confidence: 20-60 species
- Medium-High confidence: 15-40 species
- Medium confidence: 30-80 species
- Low-Medium confidence: 40-120 species
- Low confidence: 100-300+ species

## 🔬 Analysis Workflow

### Step 1: Initial Screening
```bash
# Run integration
python3 integrate_66ce4dde_EN.py

# Check summary
cat 66ce4dde_integrated_summary.txt
```

### Step 2: Examine High Confidence Results
```bash
# View top 20 high confidence viruses
head -21 66ce4dde_integrated_high_confidence_species.tsv | column -t

# Count high confidence species
tail -n +2 66ce4dde_integrated_high_confidence_species.tsv | wc -l
```

### Step 3: Investigate Candidates
```bash
# View candidate viruses
cat 66ce4dde_integrated_candidate_viruses.tsv | column -t

# Filter by reads threshold
awk -F'\t' 'NR==1 || $4>500' 66ce4dde_integrated_candidate_viruses.tsv
```

### Step 4: Database Comparison
```bash
# View Venn statistics
head -20 66ce4dde_integrated_venn_data.txt
```

## 📊 Visualization Examples

### 1. Create Venn Diagram (R)

```R
library(VennDiagram)

# Read summary statistics
venn_file <- "66ce4dde_integrated_venn_data.txt"
lines <- readLines(venn_file, n=10)

# Extract counts (adjust line numbers as needed)
rvdb_count <- as.numeric(sub(".*: ", "", lines[5]))
ncbi_count <- as.numeric(sub(".*: ", "", lines[6]))

# Create Venn diagram
venn.diagram(
  x = list(
    RVDB = 1:rvdb_count,
    NCBI_RefSeq = (rvdb_count+1):(rvdb_count+ncbi_count)
  ),
  category.names = c("RVDB", "NCBI RefSeq"),
  filename = "66ce4dde_venn.png",
  output = TRUE,
  fill = c("#440154ff", "#21908dff"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2
)
```

### 2. Top Viruses Bar Chart (R)

```R
library(ggplot2)
library(dplyr)

# Read high confidence data
data <- read.table("66ce4dde_integrated_high_confidence_species.tsv",
                   header = TRUE, sep = "\t")

# Select top 20
top20 <- head(data[order(-data$Avg_Percent),], 20)

# Create bar chart
ggplot(top20, aes(x = reorder(Taxon_Name, Avg_Percent), 
                  y = Avg_Percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 High Confidence Viruses",
       x = "Virus Species", 
       y = "Average Percentage (%)") +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("top20_viruses.pdf", width = 10, height = 8)
```

### 3. Database Correlation Plot (R)

```R
# Read high confidence data
data <- read.table("66ce4dde_integrated_high_confidence_species.tsv",
                   header = TRUE, sep = "\t")

# Correlation plot
pdf("database_correlation.pdf", width = 8, height = 8)
plot(data$RVDB_Reads_Direct, data$NCBI_Reads_Direct,
     log = "xy",
     xlab = "RVDB Direct Reads (log)",
     ylab = "NCBI Direct Reads (log)",
     main = "Database Correlation - High Confidence Species",
     pch = 19, col = rgb(0, 0, 1, 0.5))
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
legend("topleft", legend = "1:1 line", col = "red", lty = 2, lwd = 2)
dev.off()
```

## 🔍 Validation Strategies

### For Candidate Viruses (Medium and Low-Medium Confidence)

#### 1. BLAST Validation
```bash
# Extract candidate virus name
VIRUS_NAME="Candidate_virus_name"

# Search in original Kraken report
grep "$VIRUS_NAME" 66ce4dde_kraken2_RVDB.txt
grep "$VIRUS_NAME" 66ce4dde_kraken2_NCBI.txt

# Note: Extract reads and BLAST against nt database
```

#### 2. Coverage Analysis
```bash
# Check read coverage in original reports
# Higher coverage = more reliable detection
```

#### 3. Sequence Assembly
```bash
# If sufficient reads (>1000), attempt assembly
# Validate assembled contigs with BLAST
```

## 📝 Publication Reporting Template

### Methods Section

```
Viral classification was performed using Kraken2 v2.1.2 with two 
reference databases: RVDB (version X.X) and NCBI Viral RefSeq 
(downloaded YYYY-MM-DD). To integrate results and enhance reliability, 
we implemented a confidence stratification strategy. 

We implemented a two-dimensional confidence stratification strategy:

**High confidence**: Viruses detected in both databases with ≥100 direct reads  
**Medium-High confidence**: Viruses detected in both databases with 50-99 direct reads  
**Medium confidence**: Viruses detected in a single database with ≥100 direct reads  
**Low-Medium confidence**: Viruses detected in a single database with 50-99 direct reads  
**Low confidence**: Any detection with <50 direct reads (excluded from analysis)

Only classifications with direct read assignments were counted to avoid 
inflated numbers from parent-level taxonomic aggregation.
```

### Results Section

```
Viral metagenomic analysis identified X high-confidence viral species 
in sample 66ce4dde (Figure 1A, Supplementary Table S1). RVDB detected 
Y species while NCBI RefSeq detected Z species, with W species 
identified by both databases (Figure 1B). 

The most abundant virus was [Virus A], accounting for XX% of viral 
reads, followed by [Virus B] (YY%) and [Virus C] (ZZ%) (Figure 1C). 
These high-confidence detections represent the core viral community.

Additionally, N medium-high confidence viruses (detected in both databases 
but with lower read counts) and M medium-confidence candidate viruses 
(detected in a single database) were identified (Supplementary Table S2), 
potentially representing rare or recently discovered viruses requiring 
further validation.
```

## 🛠️ Troubleshooting

### Issue 1: "pandas not found"
```bash
# Solution
pip install --user pandas

# Verify
python3 -c "import pandas; print(pandas.__version__)"
```

### Issue 2: "File not found"
```bash
# Check current directory
ls -la *.txt

# Use full paths
python3 integrate_66ce4dde_EN.py \
  /full/path/to/66ce4dde_kraken2_RVDB.txt \
  /full/path/to/66ce4dde_kraken2_NCBI.txt
```

### Issue 3: No high confidence results
```bash
# Solution: Lower thresholds
python3 integrate_66ce4dde_EN.py \
  66ce4dde_kraken2_RVDB.txt \
  66ce4dde_kraken2_NCBI.txt \
  50 25 0
```

### Issue 4: Too many low confidence results
```bash
# Solution: Increase minimum direct reads
python3 integrate_66ce4dde_EN.py \
  66ce4dde_kraken2_RVDB.txt \
  66ce4dde_kraken2_NCBI.txt \
  100 50 10
```

## 📚 Understanding Database Differences

### Why Results Differ

1. **Database Size**
   - RVDB: Large, comprehensive (~5M sequences)
   - NCBI RefSeq: Curated, representative (~100K sequences)

2. **Coverage**
   - RVDB: Broad, includes uncultured viruses
   - NCBI RefSeq: Well-characterized viruses only

3. **Quality Control**
   - RVDB: Inclusive approach
   - NCBI RefSeq: Strict quality standards

### Best Practices

- ✅ **Clinical samples**: Prioritize NCBI RefSeq results
- ✅ **Environmental samples**: Leverage RVDB's broader coverage
- ✅ **Novel virus discovery**: Use RVDB for sensitivity
- ✅ **Publication**: Report intersection as high-confidence core

## 🔗 Related Tools & Resources

### Companion Scripts
- `compress_fastq_EN.sh` - Compress FASTQ files
- `find_new_data_EN.sh` - Locate new data files
- `run_taxprofiler_EN.sh` - Run TaxProfiler workflow
- `check_databases.sh` - Verify database configuration

### Documentation
- [TaxProfiler Manual](https://nf-co.re/taxprofiler)
- [Kraken2 Manual](https://github.com/DerrickWood/kraken2/wiki)
- [RVDB Database](https://rvdb-prot.pasteur.fr/)
- [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)

## 💡 Tips & Tricks

### 1. Batch Processing Multiple Samples
```bash
# Create wrapper script
for sample in sample1 sample2 sample3; do
    python3 integrate_66ce4dde_EN.py \
        ${sample}_kraken2_RVDB.txt \
        ${sample}_kraken2_NCBI.txt
    mv 66ce4dde_integrated_* ${sample}_integrated_
done
```

### 2. Extract Specific Virus Information
```bash
# Search for specific virus
VIRUS="Herpesvirus"
grep -i "$VIRUS" 66ce4dde_integrated_full.tsv
```

### 3. Compare Multiple Samples
```R
# R script to compare samples
samples <- c("sample1", "sample2", "sample3")
data_list <- list()

for (s in samples) {
    file <- paste0(s, "_integrated_high_confidence_species.tsv")
    data_list[[s]] <- read.table(file, header = TRUE, sep = "\t")
}

# Create comparison matrix
# ... (custom analysis)
```

## 📞 Support & Citation

### Getting Help
- Check documentation first
- Review troubleshooting section
- Examine output summary files

### Citation
If you use this tool in your research, please cite:

```
TaxProfiler: Costello et al. (2023) 
"Taxonomic Profiling Workflow"
DOI: 10.1101/2023.10.20.563221

Kraken2: Wood et al. (2019)
"Improved metagenomic analysis with Kraken 2"
Genome Biology 20:257

RVDB: Goodacre et al. (2018)
"A Reference Viral Database (RVDB)"
mSphere 3:e00069-18
```

## 📋 Version History

- **v2.0** (2025-10-07)
  - Only counts direct read assignments
  - Eliminates parent-level duplication
  - More accurate species enumeration
  - Enhanced output organization

- **v1.0** (2025-10-06)
  - Initial release
  - Basic integration functionality

## 📄 License

This script is provided as-is for research and educational purposes.

---

**Last Updated**: October 7, 2025  
**Script**: integrate_66ce4dde_EN.py  
**Version**: 2.0

