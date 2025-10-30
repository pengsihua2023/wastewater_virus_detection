#!/bin/bash
#SBATCH --job-name=Viral_Classification
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Viral_Classification_%j.out
#SBATCH --error=Viral_Classification_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🦠  Metagenome Viral Classification Workflow"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0

# User's conda environment path (adjust according to your actual path)
USER_CONDA_ENV="/home/sp96859/.conda/envs/nextflow_env"

# Get conda base path
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "⚠️  Warning: conda info --base failed, trying alternative method..."
    CONDA_BASE="/usr/local/apps/eb/Miniforge3/24.11.3-0"
fi

echo "   Conda base: $CONDA_BASE"
echo "   Target env: $USER_CONDA_ENV"

# Initialize conda
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "❌ Cannot find conda.sh at $CONDA_BASE/etc/profile.d/conda.sh"
    exit 1
fi

# Check if user environment exists
if [ ! -d "$USER_CONDA_ENV" ]; then
    echo "❌ Conda environment not found: $USER_CONDA_ENV"
    echo "   Available environments:"
    conda env list
    exit 1
fi

# Activate environment using absolute path
conda activate "$USER_CONDA_ENV"

# Force update PATH - ensure user environment's bin is first
export PATH="$USER_CONDA_ENV/bin:$PATH"

# Set conda-related environment variables
export CONDA_PREFIX="$USER_CONDA_ENV"
export CONDA_DEFAULT_ENV="nextflow_env"

# Get Python version and set PYTHONPATH
PYTHON_VERSION=$("$USER_CONDA_ENV/bin/python" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "3.9")
export PYTHONPATH="$USER_CONDA_ENV/lib/python${PYTHON_VERSION}/site-packages:${PYTHONPATH:-}"

# Verify environment activation - check Python path
PYTHON_PATH=$(which python)
echo "   After PATH update:"
echo "   - Python path: $PYTHON_PATH"
echo "   - CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV"

if [[ "$PYTHON_PATH" == *"$USER_CONDA_ENV"* ]]; then
    echo "✅ Conda environment activated successfully!"
    echo "   Environment: nextflow_env"
    echo "   Python: $PYTHON_PATH"
else
    echo "❌ Failed to activate user conda environment!"
    echo "   Expected Python in: $USER_CONDA_ENV"
    echo "   Actual Python: $PYTHON_PATH"
    echo ""
    echo "Debugging PATH:"
    echo "$PATH" | tr ':' '\n' | head -5
    exit 1
fi

# Verify tools
echo "🧪 2. Verifying tools..."
echo "✅ Nextflow: $(which nextflow)"

# Check for Apptainer/Singularity (required for containers)
if command -v apptainer &> /dev/null; then
    echo "✅ Apptainer: $(which apptainer)"
elif command -v singularity &> /dev/null; then
    echo "✅ Singularity: $(which singularity)"
else
    echo "❌ Apptainer/Singularity not found (required for containers)"
    exit 1
fi

# Check VirSorter2 dependencies
echo ""
echo "🔍 Checking VirSorter2 dependencies..."
echo "   Current Python: $(which python)"
echo "   Python version: $(python --version 2>&1)"

# Check screed
if python -c "import screed" 2>/dev/null; then
    SCREED_VERSION=$(python -c "import screed; print(screed.__version__)" 2>/dev/null)
    echo "✅ screed module: installed (version: $SCREED_VERSION)"
else
    echo "❌ screed module: NOT FOUND"
    echo ""
    echo "⚠️  VirSorter2 requires the 'screed' Python module!"
    echo "   Debug info:"
    echo "   - Current conda env: $CONDA_DEFAULT_ENV"
    echo "   - Python path: $(which python)"
    echo "   - Python packages location: $(python -c 'import site; print(site.getsitepackages()[0])')"
    echo ""
    echo "   Please install screed:"
    echo "   conda activate nextflow_env"
    echo "   conda install -c bioconda screed -y"
    echo ""
    exit 1
fi

# Check VirSorter2
if command -v virsorter &> /dev/null; then
    VIRSORTER_PATH=$(which virsorter)
    echo "✅ VirSorter2: $VIRSORTER_PATH"
    # Try to get version
    VIRSORTER_VERSION=$(virsorter --version 2>&1 | grep -oP 'version \K[0-9.]+' || echo "unknown")
    echo "   Version: $VIRSORTER_VERSION"
else
    echo "⚠️  VirSorter2 command not found in PATH"
    echo "   Make sure it's installed in nextflow_env"
fi

echo ""
echo "ℹ️  Note: Workflow execution environment"
echo "   - Quality control (fastp): Conda environment"
echo "   - Assembly tools (MEGAHIT, metaSPAdes): Apptainer containers"
echo "   - Viral identification tools:"
echo "     * VirSorter2: Pre-installed in nextflow_env ✅"
echo "     * DeepVirFinder: Pre-installed in dvf environment ✅"
echo "   Container images and Conda environments will be auto-downloaded on first run"
echo ""

# Set database paths
# VirSorter2 database for viral identification (REQUIRED)
VIRSORTER2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/virsorter2/db"

# DeepVirFinder installation directory
DEEPVIRFINDER_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/DeepVirFinder"

# Verify databases and tools
echo "🗄️ 3. Verifying databases and tools..."

# VirSorter2 database (required)
if [ -d "$VIRSORTER2_DB" ]; then
    echo "✅ VirSorter2 database: $VIRSORTER2_DB"
    echo "   Database size: $(du -sh $VIRSORTER2_DB | cut -f1)"
else
    echo "❌ VirSorter2 database not found: $VIRSORTER2_DB"
    echo ""
    echo "📝 Note: VirSorter2 database setup:"
    echo "   Database: VirSorter2 reference database"
    echo "   Location: /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/virsorter2/db"
    echo ""
    echo "   To setup VirSorter2 database:"
    echo "   1. Activate conda environment with VirSorter2"
    echo "   2. Run: virsorter setup -d /path/to/db -j 4"
    echo "   3. This will download ~13GB of reference data"
    echo ""
    exit 1
fi

# DeepVirFinder installation (required)
if [ -d "$DEEPVIRFINDER_DIR" ] && [ -f "$DEEPVIRFINDER_DIR/dvf.py" ]; then
    echo "✅ DeepVirFinder: $DEEPVIRFINDER_DIR"
else
    echo "❌ DeepVirFinder not found: $DEEPVIRFINDER_DIR"
    echo "   Please ensure DeepVirFinder is installed in dvf conda environment"
    exit 1
fi

echo ""

# Verify input files
echo "📁 4. Verifying input files..."
if [ -f "samplesheet.csv" ]; then
    echo "✅ Samplesheet: samplesheet.csv"
    echo "📊 Found $(wc -l < samplesheet.csv) samples"
else
    echo "❌ Samplesheet not found: samplesheet.csv"
    exit 1
fi

# Clean previous results
echo "🧹 5. Cleaning previous results..."
if [ -d "results" ]; then
    echo "Removing previous results directory..."
    rm -rf results
fi

# Set Singularity bind paths
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases:/databases"

# Run workflow
echo "🚀 6. Running Metagenome Viral Classification workflow..."
echo "Command: nextflow run metagenome_assembly_classification_workflow_en.nf -c metagenome_assembly_classification_en.config --input samplesheet.csv --outdir results --virsorter2_db $VIRSORTER2_DB --deepvirfinder_dir $DEEPVIRFINDER_DIR"
echo ""
echo "📝 Workflow steps:"
echo "   1. fastp quality control (auto adapter removal, low-quality read filtering)"
echo "   2. MEGAHIT and metaSPAdes parallel assembly"
echo "   3. VirSorter2 viral sequence identification (machine learning + rules)"
echo "   4. DeepVirFinder viral prediction (deep learning)"
echo "   5. Tool result merging (VirSorter2 + DeepVirFinder per assembler)"
echo "   6. Assembler comparison (MEGAHIT vs SPAdes) → Final consensus viral list ⭐"
echo ""
echo "✅ Note: Using VirSorter2 from nextflow_env environment"
echo "✅ Note: Using DeepVirFinder from dvf environment"
echo "✅ Note: Full run from scratch (no resume mode)"
echo ""

nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet.csv \
    --outdir results \
    --virsorter2_db "$VIRSORTER2_DB" \
    --deepvirfinder_dir "$DEEPVIRFINDER_DIR"

# Check results
echo ""
echo "=========================================="
echo "🎯 Workflow Results"
echo "=========================================="

if [ $? -eq 0 ]; then
    echo "✅ Workflow completed successfully!"
    
    if [ -d "results" ]; then
        echo "📁 Results directory created: results/"
        echo "📊 Generated results:"
        
        # Check fastp results
        if [ -d "results/fastp" ]; then
            echo "  ✅ fastp quality reports: results/fastp/"
            FASTP_HTML=$(find results/fastp -name "*.html" | wc -l)
            echo "     - Generated $FASTP_HTML HTML quality reports"
        fi
        
        # Check clean reads
        if [ -d "results/clean_reads" ]; then
            echo "  ✅ Clean reads: results/clean_reads/"
            CLEAN_READS=$(find results/clean_reads -name "*_clean_R*.fastq.gz" | wc -l)
            echo "     - Generated $CLEAN_READS clean read files"
        fi
        
        # Check assembly results
        if [ -d "results/assembly_megahit" ]; then
            echo "  ✅ MEGAHIT assembly: results/assembly_megahit/"
            MEGAHIT_CONTIGS=$(find results/assembly_megahit -name "*_megahit_contigs.fa" | wc -l)
            echo "     - Generated $MEGAHIT_CONTIGS contig files"
        fi
        
        if [ -d "results/assembly_spades" ]; then
            echo "  ✅ SPAdes assembly: results/assembly_spades/"
            SPADES_CONTIGS=$(find results/assembly_spades -name "*_spades_contigs.fa" | wc -l)
            echo "     - Generated $SPADES_CONTIGS contig files"
        fi
        
        # Check VirSorter2 results
        if [ -d "results/virsorter2_megahit" ]; then
            echo "  ✅ VirSorter2 MEGAHIT results: results/virsorter2_megahit/"
            VS2_MEGAHIT=$(find results/virsorter2_megahit -name "*_vs2_final-viral-score.tsv" | wc -l)
            echo "     - Generated $VS2_MEGAHIT viral identification reports"
        fi
        
        if [ -d "results/virsorter2_spades" ]; then
            echo "  ✅ VirSorter2 SPAdes results: results/virsorter2_spades/"
            VS2_SPADES=$(find results/virsorter2_spades -name "*_vs2_final-viral-score.tsv" | wc -l)
            echo "     - Generated $VS2_SPADES viral identification reports"
        fi
        
        # Check DeepVirFinder results
        if [ -d "results/deepvirfinder_megahit" ]; then
            echo "  ✅ DeepVirFinder MEGAHIT results: results/deepvirfinder_megahit/"
            DVF_MEGAHIT=$(find results/deepvirfinder_megahit -name "*_dvf_output.txt" | wc -l)
            echo "     - Generated $DVF_MEGAHIT viral prediction reports"
        fi
        
        if [ -d "results/deepvirfinder_spades" ]; then
            echo "  ✅ DeepVirFinder SPAdes results: results/deepvirfinder_spades/"
            DVF_SPADES=$(find results/deepvirfinder_spades -name "*_dvf_output.txt" | wc -l)
            echo "     - Generated $DVF_SPADES viral prediction reports"
        fi
        
        # Check merged viral reports
        if [ -d "results/merged_viral_reports_megahit" ]; then
            echo "  ✅ Merged viral reports (MEGAHIT): results/merged_viral_reports_megahit/"
            MERGED_VIRAL_MEGAHIT=$(find results/merged_viral_reports_megahit -name "*_viral_merged_report.txt" | wc -l)
            echo "     - Generated $MERGED_VIRAL_MEGAHIT comprehensive viral analysis reports"
        fi
        
        if [ -d "results/merged_viral_reports_spades" ]; then
            echo "  ✅ Merged viral reports (SPAdes): results/merged_viral_reports_spades/"
            MERGED_VIRAL_SPADES=$(find results/merged_viral_reports_spades -name "*_viral_merged_report.txt" | wc -l)
            echo "     - Generated $MERGED_VIRAL_SPADES comprehensive viral analysis reports"
        fi
        
        # Check assembler comparison
        if [ -d "results/assembler_comparison" ]; then
            echo "  ✅ Assembler comparison (MEGAHIT vs SPAdes): results/assembler_comparison/"
            ASSEMBLER_COMP=$(find results/assembler_comparison -name "*_assembler_comparison.txt" | wc -l)
            echo "     - Generated $ASSEMBLER_COMP assembler comparison reports"
            CONSENSUS_SEQS=$(find results/assembler_comparison -name "*_consensus_viral_sequences.txt" | wc -l)
            echo "     - Generated $CONSENSUS_SEQS final consensus viral sequence lists"
        fi
        
        echo ""
        echo "📋 Summary of key viral identification files:"
        echo "  VirSorter2 viral scores:"
        find results/virsorter2_* -name "*_vs2_final-viral-score.tsv" 2>/dev/null | head -10
        echo ""
        echo "  DeepVirFinder predictions:"
        find results/deepvirfinder_* -name "*_dvf_output.txt" 2>/dev/null | head -10
        echo ""
        echo "  Merged viral reports:"
        find results/merged_viral_reports_* -name "*.txt" -o -name "*.csv" 2>/dev/null | head -10
        echo ""
        echo "  Assembler comparison (MEGAHIT vs SPAdes):"
        find results/assembler_comparison -name "*_assembler_comparison.txt" 2>/dev/null | head -10
        echo ""
        echo "  ⭐ Final consensus viral sequences (recommended for downstream analysis):"
        find results/assembler_comparison -name "*_consensus_viral_sequences.txt" 2>/dev/null | head -10
        echo ""
        echo "Total files: $(find results -type f | wc -l)"
        
    else
        echo "❌ Results directory not found"
    fi
    
else
    echo "❌ Workflow failed with exit code: $?"
    echo "🔍 Check the error log for details"
fi

echo ""
echo "End time: $(date)"
echo "=========================================="

