# Viral Metagenome Assembly and Kraken2 Classification Workflow

![Version](https://img.shields.io/badge/version-2.5.0-blue)
![Nextflow](https://img.shields.io/badge/nextflow-≥25.04-brightgreen)
![License](https://img.shields.io/badge/license-MIT-orange)

## 📋 Overview

This is an automated Nextflow-based workflow specifically designed for **viral metagenome assembly and classification**. The workflow integrates quality control (fastp), two mainstream metagenomic assemblers (MEGAHIT and metaSPAdes), and the fast taxonomic classifier Kraken2, enabling rapid assembly and taxonomic annotation of viral genomes from raw sequencing data.

### Key Features

- ✅ **Quality Control**: Fast quality checking, adapter removal, and read filtering using fastp
- ✅ **Parallel Assembly**: Simultaneous assembly using both MEGAHIT and metaSPAdes for comprehensive results
- ✅ **Fast Classification**: Rapid taxonomic annotation of assembled contigs using Kraken2
- ✅ **Comprehensive Analysis**: Automatic merging of classification results from both assembly methods with comparison report
- ✅ **Hybrid Execution Environment**: Assembly tools in Apptainer containers, fastp and Kraken2 in Conda environments
- ✅ **Resumable**: Supports Nextflow's `-resume` feature for continuing from breakpoints after failures
- ✅ **Parallel Processing**: Automatic parallel processing of multiple samples
- ✅ **Resource Optimized**: Optimized memory and CPU allocation for viral metagenomic data
- ✅ **Stable**: metaSPAdes uses `--only-assembler` mode to avoid memory issues
- ✅ **Flexible**: Optional quality control step (suitable for pre-cleaned data)

## 🔬 Workflow

```
Raw Sequencing Data (FASTQ)
         ↓
    fastp Quality Control
    (optional, enabled by default)
         ↓
   Cleaned FASTQ
         ↓
    ┌────┴────┐
    ↓         ↓
 MEGAHIT   metaSPAdes
    ↓         ↓
 contigs   contigs
    ↓         ↓
Kraken2   Kraken2
    ↓         ↓
Classification  Classification
    └─────┬──────┘
          ↓
Comprehensive Analysis Report
  (auto-merge and compare)
```

## 📦 Environment Requirements

### Required Software

- **Nextflow** ≥ 25.04.7
- **Apptainer/Singularity** (for running MEGAHIT and metaSPAdes containers)
- **Conda/Mamba** (for Kraken2 environment and managing Nextflow)

### Hardware Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU | 16 cores | 32 cores |
| Memory | 128 GB | 512 GB |
| Storage | 100 GB | 500 GB |

### Database Requirements

- **Kraken2 Database** (viral reference database)
  - Location: `/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref`
  - Size: Varies by database version

## 🚀 Quick Start

### 1. Clone or Download Workflow

```bash
cd /scratch/sp96859/Meta-genome-data-analysis/Apptainer/contig-Kraken2-based
```

### 2. Prepare Input Data

Create sample sheet `samplesheet_mag_tax.csv`:

```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

**Note**:
- First line must be header
- Paths must be absolute paths
- Compressed files (.gz) are supported

### 3. Modify Configuration (Optional)

Edit database path in `run_mag_tax_classification_en.sh`:

```bash
# Line 48: Set Kraken2 database path
KRAKEN2_DB="/path/to/your/kraken2/database"
```

### 4. Submit Job

```bash
# Submit using SLURM
sbatch run_mag_tax_classification_en.sh

# Or run Nextflow directly
nextflow run mag_tax_classification_workflow_en.nf \
    -c mag_tax_classification_en.config \
    --input samplesheet_mag_tax.csv \
    --outdir results \
    --kraken2_db /path/to/kraken2/db
```

### 5. Resume Failed Jobs

If the job fails, use `-resume` to continue from breakpoint:

```bash
nextflow run mag_tax_classification_workflow_en.nf \
    -c mag_tax_classification_en.config \
    --input samplesheet_mag_tax.csv \
    --outdir results \
    --kraken2_db /path/to/kraken2/db \
    -resume
```

## 📂 Output Results

After workflow completion, results are stored in the `results/` directory:

```
results/
├── fastp/
│   ├── sample1_fastp.html                 # fastp quality report (HTML)
│   ├── sample1_fastp.json                 # fastp quality report (JSON)
│   ├── sample2_fastp.html
│   └── sample2_fastp.json
├── kraken2_megahit/
│   ├── sample1_megahit_classification.txt  # Detailed classification of MEGAHIT contigs
│   ├── sample1_megahit_report.txt         # Classification summary of MEGAHIT contigs
│   ├── sample2_megahit_classification.txt
│   └── sample2_megahit_report.txt
├── kraken2_spades/
│   ├── sample1_spades_classification.txt   # Detailed classification of SPAdes contigs
│   ├── sample1_spades_report.txt          # Classification summary of SPAdes contigs
│   ├── sample2_spades_classification.txt
│   └── sample2_spades_report.txt
└── merged_reports/                        # ⭐ NEW: Comprehensive analysis reports
    ├── sample1_merged_report.txt          # Text format comprehensive report
    ├── sample1_merged_report.csv          # CSV format detailed data
    ├── sample2_merged_report.txt
    └── sample2_merged_report.csv
```

### Output File Descriptions

#### fastp Reports (`*_fastp.html` / `*_fastp.json`)
fastp quality control reports containing:
- **Before/After filtering statistics**: Number of reads, bases, Q20/Q30 ratios
- **Quality distribution**: Quality score distribution at each position
- **Base content**: GC content, ATCG distribution
- **Adapter detection**: Automatically detected and removed adapter sequences
- **Filtering reasons**: Statistics on low-quality, too-short, too-long filtering
- **Insert size**: Estimated insert size for paired-end reads

💡 **Tip**: Open HTML file in a browser to view visual reports

#### Classification Files (`*_classification.txt`)
Contains detailed classification information for each contig:
- C/U: Whether classified successfully
- Contig ID
- Classified taxID
- Contig length
- Classification details

#### Report Files (`*_report.txt`)
Contains classification statistical summary:
- Percentage of reads in each taxonomic unit
- Hierarchical structure of classification (from domain to species)
- Classification confidence

#### Comprehensive Analysis Reports (`*_merged_report.txt` / `*_merged_report.csv`) ⭐ Core Feature
Automatically merges classification results from MEGAHIT and SPAdes, generating easy-to-read comparison analysis reports.

**Text Report (`.txt`)** - Formatted readable report:
```
================================================================================
Kraken2 Comprehensive Analysis Report - MEGAHIT vs SPAdes Assembly
================================================================================

[Overall Statistics]
--------------------------------------------------------------------------------
SPAdes total contigs:     988,186.0
MEGAHIT total contigs:     52,329.0

SPAdes unclassified:      980,538.0 (99.23%)
MEGAHIT unclassified:      47,812.0 (91.37%)

SPAdes virus classified:      804.0 (0.08%)
MEGAHIT virus classified:     488.0 (0.93%)


[Species Level Comparison]
--------------------------------------------------------------------------------
Species Name                                               SPAdes       MEGAHIT     
--------------------------------------------------------------------------------
Shigella phage SfIV                                16           13          
Escherichia phage 500465-1                         17           9           
...

[Genus Level Comparison]
[Family Level Comparison]
[Unique Findings]
- Species found only in SPAdes (Top 10)
- Species found only in MEGAHIT (Top 10)
```

**Report Content Description**:
- **Overall Statistics**: Comparison of total contigs, unclassified ratio, virus classification ratio
- **Species Level Comparison**: Top 30 species, showing detection counts for both methods side-by-side
- **Genus Level Comparison**: Top 20 genera for higher-level classification comparison
- **Family Level Comparison**: Top 15 families to understand major viral family distribution
- **Unique Findings**: Identify species detected only by a single method, assess complementarity

**Detailed Data Table (`.csv`)**:
- Can be opened in Excel, Python, or R for further analysis
- Contains complete comparison data for all taxonomic units (TaxID, rank, reads count, percentage)
- Supports custom filtering, sorting, and visualization
- Suitable for generating publication-ready figures

💡 **Recommended Usage**:
1. **Quick Assessment**: First view the Overall Statistics section of the `.txt` file
2. **In-depth Analysis**: Review detailed comparisons at species/genus/family levels
3. **Method Evaluation**: Understand differences between assembly methods through "Unique Findings"
4. **Data Mining**: Use `.csv` file for custom analysis and visualization
5. **Result Validation**: Species detected by both methods have higher confidence

## ⚙️ Parameter Configuration

### Main Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | - | Sample sheet CSV file (required)|
| `--outdir` | `./results` | Output directory |
| `--kraken2_db` | - | Kraken2 database path (required)|
| `--skip_merge_reports` | `false` | Skip comprehensive analysis report generation |

### fastp Quality Control Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_fastp` | `false` | Skip fastp quality control (for pre-cleaned data)|
| `--fastp_qualified_quality` | 20 | Minimum quality value for qualified bases (Phred score)|
| `--fastp_unqualified_percent` | 40 | Maximum percentage of low-quality bases allowed |
| `--fastp_min_length` | 50 | Minimum read length after filtering |

**fastp Automatic Features**:
- ✅ Automatic adapter detection and removal
- ✅ Automatic paired-end read quality checking
- ✅ Generate HTML and JSON format quality reports

### MEGAHIT Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `megahit_memory` | 0.8 | Memory usage ratio |
| `megahit_min_contig_len` | 1000 | Minimum contig length |

### metaSPAdes Parameters

- Uses `--only-assembler` mode, skips error correction
- Automatically uses metaSPAdes mode (optimized for metagenomes)

## 🔧 Resource Configuration

Adjust resource allocation in `mag_tax_classification_en.config`:

```groovy
process {
    withName: 'FASTP' {
        cpus = 8
        memory = '16 GB'
        time = '4h'
    }
    
    withName: 'MEGAHIT_ASSEMBLY' {
        cpus = 16
        memory = '64 GB'
        time = '12h'
    }
    
    withName: 'SPADES_ASSEMBLY' {
        cpus = 32
        memory = '512 GB'  // Can be adjusted based on actual situation
        time = '48h'
    }
    
    withName: 'KRAKEN2_CLASSIFICATION_*' {
        cpus = 16
        memory = '48 GB'
        time = '8h'
    }
}
```

## 🐛 Troubleshooting

### Issue 1: SPAdes Out of Memory or Assertion Failure

**Symptoms**:
```
ERROR: The reads contain too many k-mers to fit into available memory
or
spades-hammer: Assertion 'data[i].count()' failed
```

**Solutions**:
1. **Already enabled by default**: Workflow uses `metaspades.py --only-assembler` to skip error correction
2. Increase memory allocation (already set to 512GB in config)
3. If still failing, check if SLURM partition supports large memory nodes

### Issue 2: Container Pull Failure

**Symptoms**:
```
Failed to pull singularity image ... manifest unknown
```

**Solutions**:
1. **Already fixed**: Workflow now uses stable container image tags
2. Check network connection to quay.io
3. Wait a moment and retry (may be temporary network issue)
4. Manually pre-pull containers (optional):
   ```bash
   apptainer pull docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1
   apptainer pull docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1
   ```

### Issue 3: Kraken2 Database Not Found

**Symptoms**:
```
Kraken2 database not found
or
Unable to find required file "hash.k2d"
```

**Solutions**:
1. Check `KRAKEN2_DB` path in `run_mag_tax_classification_en.sh`
2. Ensure database directory contains required files:
   - `hash.k2d`
   - `opts.k2d`
   - `taxo.k2d`
3. Verify database path permissions: `ls -l $KRAKEN2_DB`

### Issue 4: Sample Sheet Format Error

**Symptoms**:
```
Error parsing CSV file
or
A process input channel evaluates to null
```

**Solutions**:
1. Ensure CSV file contains header: `sample,fastq_1,fastq_2`
2. Ensure no blank lines
3. **Use absolute paths** (relative paths may cause issues)
4. Check file encoding (use UTF-8)
5. Verify FASTQ files exist: `ls -l /path/to/fastq`

### Issue 5: Nextflow Syntax Error

**Symptoms**:
```
ERROR ~ Unknown method invocation `params`
or
token recognition error at: '$'
```

**Solutions**:
- **Already fixed**: Workflow uses correct Nextflow DSL2 syntax
- If still having issues, ensure using Nextflow ≥ 25.04.7

### Issue 6: Conda Environment Creation Failed

**Symptoms**:
```
Failed to create conda environment
```

**Solutions**:
1. Check if conda cache directory has write permissions:
   ```bash
   ls -ld /scratch/sp96859/Meta-genome-data-analysis/conda_cache
   ```
2. Manually create Kraken2 environment (optional):
   ```bash
   conda create -n kraken2_env -c bioconda kraken2=2.1.3
   ```
3. Clean conda cache: `conda clean --all`

### Issue 7: Comprehensive Report Generation Failed

**Symptoms**:
```
IndentationError: unindent does not match any outer indentation level
or
ERROR ~ Error executing process > 'MERGE_KRAKEN2_REPORTS'
```

**Solutions**:
- **Already fixed** (v2.5.0): Python script indentation issues resolved
- If still having issues, ensure using latest version of workflow file
- Check for encoding issues caused by Chinese characters

**Temporary skip solution** (not recommended):
```bash
nextflow run mag_tax_classification_workflow_en.nf \
    ... \
    --skip_merge_reports
```

### Issue 8: work Directory Permission Problem

**Symptoms**:
```
Cannot create work directory
```

**Solutions**:
1. Ensure current directory has write permissions
2. Or specify work directory: `nextflow run ... -work-dir /scratch/your_work_dir`

## 🚀 Performance Optimization and Best Practices

### Resource Allocation Recommendations

Adjust resources based on different data scales:

| Data Scale | MEGAHIT Memory | SPAdes Memory | SPAdes CPU |
|------------|----------------|---------------|------------|
| Small (<10GB) | 32 GB | 256 GB | 16 cores |
| Medium (10-50GB) | 64 GB | 512 GB | 32 cores |
| Large (>50GB) | 128 GB | 1024 GB | 48 cores |

### Runtime Estimation

Based on 10GB paired-end sequencing data (~10M reads):
- **fastp quality control**: 10-30 minutes
- **MEGAHIT assembly**: 1-3 hours
- **metaSPAdes assembly**: 12-24 hours
- **Kraken2 classification**: 30 minutes-2 hours (depends on database size)

**Total**: ~15-30 hours (most time spent on metaSPAdes assembly)

### Tips for Saving Computational Resources

1. **Test with small dataset first**
   ```bash
   # Use subsample to test workflow
   seqtk sample -s100 input_R1.fastq.gz 100000 > test_R1.fastq.gz
   ```

2. **Use -resume to recover failed jobs**
   - Nextflow caches successful steps
   - Only reruns failed parts

3. **Parallel processing of multiple samples**
   - Workflow automatically processes all samples in sample sheet in parallel
   - SLURM schedules tasks based on resource availability

4. **Monitor resource usage**
   ```bash
   # Check SLURM job status
   squeue -u $USER
   
   # Check job resource usage
   sacct -j <job_id> --format=JobID,MaxRSS,MaxVMSize,CPUTime
   ```

## 📊 Result Interpretation

### Kraken2 Report File Interpretation

Example output:
```
45.21  2261    2261    U       0       unclassified
54.79  2738    145     R       1       root
54.34  2716    0       D       10239     Viruses
54.34  2716    10      D1      2731619     Duplodnaviria
...
```

Column descriptions:
1. **Percentage**: Percentage of reads assigned to this taxon
2. **Reads Count**: Number of reads classified to this taxon
3. **Direct Reads**: Number of reads directly classified to this taxon (excluding subtaxa)
4. **Rank Code**: U=unclassified, R=root, D=domain, K=kingdom, P=phylum, C=class, O=order, F=family, G=genus, S=species
5. **TaxID**: NCBI taxonomy ID
6. **Scientific Name**: Taxon name

## 📖 Technical Details

### Tools Used

| Tool | Version | Purpose | Execution Environment |
|------|---------|---------|----------------------|
| fastp | 0.23.4 | Quality control and read filtering | Conda environment |
| MEGAHIT | 1.2.9 | Fast metagenomic assembly | Apptainer container |
| metaSPAdes | 3.15.5 | High-quality metagenomic assembly | Apptainer container |
| Kraken2 | 2.1.3 | Fast taxonomic classification | Conda environment |

#### Container Images and Environments

- **fastp**: `bioconda::fastp=0.23.4` (Conda environment)
- **MEGAHIT**: `docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1`
- **SPAdes**: `docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1`
- **Kraken2**: `bioconda::kraken2=2.1.3` (Conda environment)

### Workflow Files

- `mag_tax_classification_workflow_en.nf` - Main workflow
- `mag_tax_classification_en.config` - Configuration file
- `run_mag_tax_classification_en.sh` - SLURM submission script
- `samplesheet_mag_tax.csv` - Sample input sheet

### Hybrid Execution Environment Description

This workflow adopts a **hybrid execution strategy** to improve stability and flexibility:

#### Apptainer Containers (for assembly tools)
- **Advantages**: Complete isolation, fixed versions, strong portability
- **Tools**: MEGAHIT, metaSPAdes
- **Cache location**: `/scratch/sp96859/Meta-genome-data-analysis/Apptainer/singularity`

#### Conda Environments (for quality control and classification)
- **Advantages**: Avoids large container image pull failures, faster startup speed
- **Tools**: fastp, Kraken2
- **Cache location**: `/scratch/sp96859/Meta-genome-data-analysis/conda_cache`
- **Auto-creation**: Automatically installed on first run

This hybrid strategy ensures:
- ✅ Quality control tools start quickly without pulling containers
- ✅ Assembly tools run stably in isolated environments
- ✅ Kraken2 can quickly access large database files
- ✅ Avoids container image pull timeout issues

## 🔄 Workflow Version History

### v2.5.0 (Current Version - 2025-10-14)
- ✅ **Optimized report format**: Comprehensive analysis report more beautiful and readable
  - Standardized title format (80-character separator lines)
  - Unified number format (floats + thousand separators)
  - Optimized column widths, more compact table layout
- ✅ **Added family level comparison**: Shows distribution of Top 15 viral families
- ✅ **Added unique findings**: Identifies species detected only by single assembly method
- ✅ **Expanded comparison range**:
  - Species level: Top 20 → Top 30
  - Genus level: Top 15 → Top 20
- ✅ **Fixed Python indentation error**: Resolved IndentationError in embedded Python script

### v2.4.0 (2025-10-13)
- ✅ **Added comprehensive analysis**: Automatically merges Kraken2 classification results from MEGAHIT and SPAdes
- ✅ **Intelligent comparison**: Generates text reports and CSV data tables showing differences and similarities between methods
- ✅ **Optional feature**: Supports `--skip_merge_reports` to skip comprehensive analysis
- ✅ **Automated process**: Fully integrated into workflow, no need to run scripts manually

### v2.3.0 (2025-10-13)
- ✅ **Added fastp quality control**: Automatic quality checking, adapter removal, and read filtering
- ✅ **Optional QC**: Supports `--skip_fastp` parameter to skip quality control (for pre-cleaned data)
- ✅ **HTML quality reports**: fastp automatically generates visualized quality reports
- ✅ **Optimized parameters**: Customizable quality thresholds, minimum length, and other parameters

### v2.2.0 (2025-10-13)
- ✅ **Removed CAT classification**: Simplified workflow, uses only Kraken2 for classification
- ✅ **Restored metaSPAdes**: Parallel running of MEGAHIT and metaSPAdes assembly
- ✅ **Hybrid execution environment**: Assemblers use Apptainer, Kraken2 uses Conda
- ✅ **Optimized SPAdes**: Uses `--only-assembler` to skip error correction, avoids memory issues
- ✅ **Increased resources**: SPAdes memory 512GB, global max_memory adjusted to 512GB
- ✅ **Fixed input parsing**: Correctly parses CSV format sample sheet
- ✅ **Fixed syntax errors**: Corrected params definition and variable references
- ✅ **Independent processes**: MEGAHIT and SPAdes separated into independent processes, no interference

### v2.1.0 (Deprecated)
- Simplified to use only MEGAHIT
- Temporarily removed SPAdes (memory issues)

### v2.0.0 (Deprecated)
- Attempted to use Apptainer containers
- Separated MEGAHIT and SPAdes into independent processes
- Had container image and syntax issues

### v1.0.0 (Original Version)
- Initial version
- Included MAG assembly, Kraken2 and CAT classification
- Had multiple syntax and execution issues

## 📚 References

1. **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890.

2. **MEGAHIT**: Li, D., et al. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31(10), 1674-1676.

3. **SPAdes/metaSPAdes**: Nurk, S., et al. (2017). metaSPAdes: a new versatile metagenomic assembler. *Genome Research*, 27(5), 824-834.

4. **Kraken2**: Wood, D.E., et al. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.

## 👥 Authors and Support

- **Author**: Assistant
- **Version**: 2.5.0
- **Last Updated**: 2025-10-14

### Getting Help

If you encounter issues, please check:
1. `.nextflow.log` - Nextflow log
2. `.command.log` in `work/` directory - Specific task logs
3. SLURM output files - `MAG_Tax_Classification_*.out`

## 📄 License

This workflow is open-sourced under the MIT license.

---

## ⚡ Quick Reference

### Common Commands

```bash
# 1. Submit job to SLURM
sbatch run_mag_tax_classification_en.sh

# 2. Check job status
squeue -u $USER

# 3. Cancel job
scancel <job_id>

# 4. View job output
tail -f MAG_Tax_Classification_<job_id>.out

# 5. Resume failed job
nextflow run mag_tax_classification_workflow_en.nf \
    -c mag_tax_classification_en.config \
    --input samplesheet_mag_tax.csv \
    --outdir results \
    --kraken2_db /path/to/kraken2/db \
    -resume

# 6. Check results
ls -lh results/kraken2_*/

# 7. View fastp quality report (open in browser)
firefox results/fastp/*_fastp.html

# 8. View Kraken2 classification summary
head results/kraken2_megahit/*_report.txt

# 9. View comprehensive analysis report (⭐ Highly Recommended)
cat results/merged_reports/*_merged_report.txt
# Or open CSV file in Excel
# results/merged_reports/*_merged_report.csv

# 10. Skip fastp quality control (if data already cleaned)
nextflow run mag_tax_classification_workflow_en.nf \
    -c mag_tax_classification_en.config \
    --input samplesheet_mag_tax.csv \
    --outdir results \
    --kraken2_db /path/to/kraken2/db \
    --skip_fastp

# 11. Skip comprehensive analysis report (not recommended unless only need individual results)
nextflow run mag_tax_classification_workflow_en.nf \
    -c mag_tax_classification_en.config \
    --input samplesheet_mag_tax.csv \
    --outdir results \
    --kraken2_db /path/to/kraken2/db \
    --skip_merge_reports

# 12. Clean work directory (release space after completion)
rm -rf work/
```

### File Location Quick Index

| Content | Location |
|---------|----------|
| Workflow file | `mag_tax_classification_workflow_en.nf` |
| Configuration file | `mag_tax_classification_en.config` |
| Startup script | `run_mag_tax_classification_en.sh` |
| Sample sheet | `samplesheet_mag_tax.csv` |
| Results directory | `results/` |
| Nextflow log | `.nextflow.log` |
| Work directory | `work/` |
| SLURM output | `MAG_Tax_Classification_*.out` |
| Kraken2 database | `/scratch/sp96859/.../kraken2_Viral_ref` |

### Key Parameter Quick Reference

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--input` | CSV file path | Sample sheet (required)|
| `--outdir` | Directory path | Output directory |
| `--kraken2_db` | Directory path | Kraken2 database (required)|
| `-resume` | - | Resume failed job |
| `-work-dir` | Directory path | Specify work directory |

---

**Good luck with your analysis!** 🧬🦠

If you have any questions, please refer to the troubleshooting section or check the log files.


