#!/bin/bash
#SBATCH --job-name=Metagenome_Assembly_Classification
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=Metagenome_Assembly_Classification_%j.out
#SBATCH --error=Metagenome_Assembly_Classification_%j.err

cd "$SLURM_SUBMIT_DIR"

echo "=========================================="
echo "🧬  fastp + Kraken2 Workflow"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load conda environment
echo "🔧 1. Setting up environment..."
module load Miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nextflow_env

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

echo ""
echo "ℹ️  Note: Workflow execution environment"
echo "   - Quality control (fastp): Conda environment"
echo "   - Assembly tools (MEGAHIT, metaSPAdes): Apptainer containers"
echo "   - Classification tool (Kraken2): Conda environment"
echo "   Container images and Conda environments will be auto-downloaded on first run"
echo ""

# Set database paths
# Only using Kraken2 database for classification
KRAKEN2_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/kraken2_Viral_ref"

# Verify databases
echo "🗄️ 3. Verifying databases..."
if [ -d "$KRAKEN2_DB" ]; then
    echo "✅ Kraken2 database: $KRAKEN2_DB"
else
    echo "❌ Kraken2 database not found: $KRAKEN2_DB"
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
echo "🚀 6. Running  fastp + Kraken2 workflow..."
echo "Command: nextflow run metagenome_assembly_classification_workflow_en.nf -c metagenome_assembly_classification_en.config --input samplesheet.csv --outdir results --kraken2_db $KRAKEN2_DB"
echo ""
echo "📝 Workflow steps:"
echo "   1. fastp quality control (auto adapter removal, low-quality read filtering)"
echo "   2. MEGAHIT and metaSPAdes parallel assembly"
echo "   3. Kraken2 taxonomic classification"
echo "   4. Comprehensive report generation"
echo ""

nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet.csv \
    --outdir results \
    --kraken2_db "$KRAKEN2_DB"

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
        
        # Check Kraken2 results
        if [ -d "results/kraken2_megahit" ]; then
            echo "  ✅ Kraken2 MEGAHIT results: results/kraken2_megahit/"
        fi
        
        if [ -d "results/kraken2_spades" ]; then
            echo "  ✅ Kraken2 SPAdes results: results/kraken2_spades/"
        fi
        
        # Check merged reports
        if [ -d "results/merged_reports" ]; then
            echo "  ✅ Comprehensive analysis reports: results/merged_reports/"
            MERGED_REPORTS=$(find results/merged_reports -name "*.txt" | wc -l)
            echo "     - Generated $MERGED_REPORTS comprehensive reports"
        fi
        
        echo ""
        echo "📋 Summary of generated files:"
        echo "  fastp reports:"
        find results/fastp -name "*.html" -o -name "*.json" 2>/dev/null | head -10
        echo ""
        echo "  Classification reports:"
        find results/kraken2_* -name "*.txt" 2>/dev/null | head -10
        echo ""
        echo "  Merged reports:"
        find results/merged_reports -name "*.txt" -o -name "*.csv" 2>/dev/null | head -10
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

