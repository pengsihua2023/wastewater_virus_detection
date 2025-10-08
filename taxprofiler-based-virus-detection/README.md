# TaxProfiler Metagenomic Classification Workflow

## 📖 Project Overview

This project uses the nf-core/taxprofiler workflow for metagenomic sequencing data classification to identify microbial composition in samples. The workflow is based on Nextflow and Singularity container technology, running in a SLURM cluster environment.

## 📁 File Structure

```
taxprofiler/
├── databases.csv              # Database configuration file
├── samplesheet.csv            # Sample configuration file
├── run_taxprofiler.sh         # SLURM job submission script
├── check_setup.sh             # Environment check script
├── check_results.sh           # Results check script
├── compress_fastq.sh          # FASTQ file compression script
├── find_new_data.sh           # Data file search script
├── README.md                  # Chinese documentation
├── README_EN.md               # This document (English)
└── data/
    └── reads/                 # Sequencing data directory
        ├── *_R1.fastq.gz     # Paired-end reads 1
        └── *_R2.fastq.gz     # Paired-end reads 2
```

## 🔧 Configuration Files

### 1. databases.csv

Database configuration file that defines the classification databases to use.

**Format:**
```csv
db_name,db_path,tool,db_type
kraken2_viral,/path/to/database/,kraken2,short
```

**Field Descriptions:**
- `db_name`: Database name (custom)
- `db_path`: Absolute path to database (must end with /)
- `tool`: Classification tool (kraken2, krakenuniq, metaphlan, etc.)
- `db_type`: Data type (short for short-read data)

**Current Configuration:**
- Uses Kraken2 viral database for viral classification

### 2. samplesheet.csv

Sample configuration file that defines input sequencing data.

**Format:**
```csv
sample,fastq_1,fastq_2,instrument_platform,run_accession
sample_name,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,ILLUMINA,run_001
```

**Field Descriptions:**
- `sample`: Sample name (unique identifier)
- `fastq_1`: Absolute path to Read 1 file (must be .fastq.gz or .fq.gz format)
- `fastq_2`: Absolute path to Read 2 file (for paired-end sequencing)
- `instrument_platform`: Sequencing platform (ILLUMINA, OXFORD_NANOPORE, etc.)
- `run_accession`: Run batch identifier

**Important Notes:**
- Files must be in compressed format (.fastq.gz or .fq.gz)
- Paths must be absolute paths
- Filenames cannot contain spaces

### 3. run_taxprofiler.sh

SLURM job submission script that configures computational resources and runtime parameters.

**SLURM Resource Configuration:**
- Job name: TaxProfiler
- Partition: bahl_p
- CPU cores: 32
- Memory: 256GB
- Runtime: 48 hours
- Output log: TaxProfiler_[JobID].out
- Error log: TaxProfiler_[JobID].err

**Environment Configuration:**
- Uses Miniforge3 Conda environment
- Activates nextflow_env environment
- Configures Singularity container bind paths
- Sets Singularity cache directory

**Workflow Parameters:**
- Version: 1.2.4
- Profile: singularity
- Tools: Kraken2
- QC: FastQC
- Report: MultiQC

## 🚀 Usage Guide

### Step 1: Prepare Data

#### 1.1 Ensure Sequencing Data is Compressed

If your FASTQ files are uncompressed, use the compression script:

```bash
# Add execute permission
chmod +x compress_fastq.sh

# Run compression script
./compress_fastq.sh
```

#### 1.2 Locate Data Files

If unsure of data file locations, use the search script:

```bash
# Add execute permission
chmod +x find_new_data.sh

# Run search script
./find_new_data.sh
```

### Step 2: Configure Files

#### 2.1 Configure Database File (databases.csv)

Edit `databases.csv` to ensure database paths are correct:

```bash
# Example
db_name,db_path,tool,db_type
kraken2_viral,/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref/,kraken2,short
```

#### 2.2 Configure Sample File (samplesheet.csv)

Edit `samplesheet.csv` to add your sample information:

```bash
# Example
sample,fastq_1,fastq_2,instrument_platform,run_accession
SRR30162756,/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/data/reads/SRR30162756.ss-sub4_R1.fastq.gz,/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/data/reads/SRR30162756.ss-sub4_R2.fastq.gz,ILLUMINA,run_001
```

### Step 3: Environment Check

Run the environment check script to ensure all configurations are correct:

```bash
# Add execute permission
chmod +x check_setup.sh

# Run check script
./check_setup.sh
```

**Check includes:**
- ✅ Configuration files exist
- ✅ FASTQ files exist
- ✅ Database paths are correct
- ✅ Required software is installed (Nextflow, Singularity, Conda)

### Step 4: Submit Job

```bash
# Add execute permission to run script
chmod +x run_taxprofiler.sh

# Submit SLURM job
sbatch run_taxprofiler.sh
```

### Step 5: Monitor Job

```bash
# Check job queue
squeue -u $USER

# View job output log (real-time)
tail -f TaxProfiler_[JobID].out

# View error log
tail -f TaxProfiler_[JobID].err
```

### Step 6: Check Results

After job completion, run the results check script:

```bash
# Add execute permission
chmod +x check_results.sh

# Run check script
./check_results.sh
```

## 📊 Output Files

After workflow completion, results are saved in the `results/` directory:

### Main Result Files

```
results/
├── multiqc/
│   └── multiqc_report.html          # Comprehensive analysis report (main)
├── kraken2/
│   └── kraken2_viral/
│       └── *_kraken2_viral.kraken2.kraken2.report.txt  # Kraken2 classification results
├── fastqc/
│   └── raw/
│       ├── *_raw_1_fastqc.html      # Read 1 quality report
│       └── *_raw_2_fastqc.html      # Read 2 quality report
└── pipeline_info/
    ├── execution_timeline.html      # Execution timeline
    ├── execution_report.html        # Execution report
    └── execution_trace.txt          # Execution trace
```

### Result Interpretation

#### 1. MultiQC Report (multiqc_report.html)

**View Method:**
```bash
# Open in browser
firefox results/multiqc/multiqc_report.html
```

**Contents:**
- 📊 Sequencing data quality statistics
- 🦠 Classification results visualization (by Kingdom, Phylum, Class, Order, Family, Species)
- 📈 Inter-sample comparisons
- 🔍 Unclassified sequence statistics

#### 2. Kraken2 Classification Report

**View Method:**
```bash
# View classification results
cat results/kraken2/kraken2_viral/*_kraken2_viral.kraken2.kraken2.report.txt | head -20
```

**Report Format:**
```
Percentage  Fragments  DirectlyAssigned  Rank  NCBI_TaxID  ScientificName
```

#### 3. FastQC Quality Report

**View Method:**
```bash
# Open in browser
firefox results/fastqc/raw/*_fastqc.html
```

**Contents:**
- Per base sequence quality
- Per sequence quality scores
- GC content distribution
- Sequence duplication levels
- Adapter content detection

## ⚙️ Advanced Configuration

### Modify Computational Resources

Edit SLURM parameters in `run_taxprofiler.sh`:

```bash
#SBATCH --cpus-per-task=32    # Modify CPU cores
#SBATCH --mem=256G            # Modify memory size
#SBATCH --time=48:00:00       # Modify runtime limit
```

### Add More Databases

Add additional databases in `databases.csv`:

```csv
db_name,db_path,tool,db_type
kraken2_viral,/path/to/kraken2_viral/,kraken2,short
kraken2_bacteria,/path/to/kraken2_bacteria/,kraken2,short
metaphlan4,/path/to/metaphlan4/,metaphlan4,short
```

Enable corresponding tools in `run_taxprofiler.sh`:

```bash
nextflow run nf-core/taxprofiler -r 1.2.4 \
  --run_kraken2 \
  --run_metaphlan4 \
  # Other parameters...
```

### Process Multiple Samples

Add multiple samples in `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,instrument_platform,run_accession
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,ILLUMINA,run_001
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,ILLUMINA,run_002
sample3,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz,ILLUMINA,run_003
```

### Single-End Sequencing Data

For single-end sequencing data, leave the `fastq_2` field empty in `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,instrument_platform,run_accession
sample_se,/path/to/sample_se.fastq.gz,,ILLUMINA,run_001
```

## 🛠️ Troubleshooting

### Issue 1: File Format Error

**Error Message:**
```
Error: does not match regular expression [^\S+\.f(ast)?q\.gz$]
```

**Solution:**
- Ensure FASTQ files are in compressed format (.fastq.gz or .fq.gz)
- Use `compress_fastq.sh` script to compress files

### Issue 2: File Not Found

**Error Message:**
```
Error: the file or directory 'xxx' does not exist
```

**Solution:**
1. Check if file paths are correct
2. Ensure absolute paths are used
3. Run `find_new_data.sh` to locate correct file locations

### Issue 3: Database Path Error

**Solution:**
1. Run `check_setup.sh` to check database paths
2. Ensure database paths end with a slash (/)
3. Ensure read permissions for databases

### Issue 4: Out of Memory

**Error Message:**
```
OutOfMemoryError
```

**Solution:**
- Increase SLURM job memory configuration: `#SBATCH --mem=512G`
- Reduce number of samples processed in parallel

### Issue 5: Singularity Cache Warning

**Warning Message:**
```
WARN: Singularity cache directory has not been defined
```

**Solution:**
- Already configured `NXF_SINGULARITY_CACHEDIR` in `run_taxprofiler.sh`
- Warning doesn't affect execution, but setting cache directory improves performance

### Issue 6: Job Timeout

**Solution:**
- Increase job runtime: `#SBATCH --time=96:00:00`
- Or reduce sample count and process in batches

## 📚 References

### Official Documentation

- **nf-core/taxprofiler**: https://nf-co.re/taxprofiler
- **Nextflow**: https://www.nextflow.io/docs/latest/
- **Singularity**: https://sylabs.io/docs/

### Tool Documentation

- **Kraken2**: https://github.com/DerrickWood/kraken2
- **FastQC**: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC**: https://multiqc.info/

### Publications

- **TaxProfiler Paper**: https://doi.org/10.1101/2023.10.20.563221
- **nf-core Framework**: https://doi.org/10.1038/s41587-020-0439-x

## 💡 Best Practices

### 1. Data Quality Control

- Use FastQC to check raw data quality before running
- Decide if data preprocessing is needed based on quality reports
- Low-quality data affects classification accuracy

### 2. Database Selection

- Choose appropriate databases based on research objectives
- Use viral databases for viral research
- Use bacterial databases for bacterial research
- Or use complete databases for comprehensive analysis

### 3. Resource Management

- Allocate resources reasonably based on data volume and sample count
- Test with small datasets first
- Monitor resource usage to avoid waste

### 4. Result Validation

- Check quality metrics in MultiQC report
- Compare results from multiple classification tools
- Pay attention to unclassified sequence percentages

### 5. Data Backup

- Regularly backup important result files
- Save complete configuration files
- Record analysis parameters and version information

## 📝 Version History

- **v1.0** (2025-10-06): Initial release
  - Configured Kraken2 viral database
  - Support for paired-end Illumina sequencing data
  - Added complete auxiliary scripts
  - Comprehensive documentation

## 👥 Contact

For questions or suggestions, please contact:
- User: sp96859
- Cluster: a2-20

## 📄 License

This workflow is based on nf-core/taxprofiler and follows the MIT License.

---

**Last Updated**: October 7, 2025
