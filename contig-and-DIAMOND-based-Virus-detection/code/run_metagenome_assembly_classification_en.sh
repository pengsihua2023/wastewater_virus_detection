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
echo "🧬  fastp + Diamond Workflow"
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
echo "   - Gene prediction (Prodigal): Conda environment"
echo "   - Classification tool (Diamond): Conda environment"
echo "   Container images and Conda environments will be auto-downloaded on first run"
echo ""

# Set database paths
# Diamond database for classification (RVDB - Reference Viral DataBase)
DIAMOND_DB="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB/RVDB_prot_ref.dmnd"

# Verify databases
echo "🗄️ 3. Verifying databases..."
if [ -f "$DIAMOND_DB" ]; then
    echo "✅ Diamond database: $DIAMOND_DB"
    echo "   Database size: $(du -h $DIAMOND_DB | cut -f1)"
else
    echo "❌ Diamond database not found: $DIAMOND_DB"
    echo ""
    echo "📝 Note: RVDB database information:"
    echo "   Database: Reference Viral DataBase (RVDB)"
    echo "   Expected file: RVDB_prot_ref.dmnd"
    echo "   Location: /scratch/sp96859/Meta-genome-data-analysis/Apptainer/databases/RVDB/"
    echo ""
    echo "   To obtain RVDB database:"
    echo "   1. Download from: https://rvdb-prot.pasteur.fr/"
    echo "   2. Build Diamond index: diamond makedb --in RVDB.fasta -d RVDB_prot_ref"
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
echo "🚀 6. Running  fastp + Diamond workflow..."
echo "Command: nextflow run metagenome_assembly_classification_workflow_en.nf -c metagenome_assembly_classification_en.config --input samplesheet.csv --outdir results --diamond_db $DIAMOND_DB"
echo ""
echo "📝 Workflow steps:"
echo "   1. fastp quality control (auto adapter removal, low-quality read filtering)"
echo "   2. MEGAHIT and metaSPAdes parallel assembly"
echo "   3. Prodigal gene prediction (metagenome mode)"
echo "   4. Diamond BLASTP classification against protein database"
echo "   5. Comprehensive report generation"
echo ""

nextflow run metagenome_assembly_classification_workflow_en.nf \
    -c metagenome_assembly_classification_en.config \
    --input samplesheet.csv \
    --outdir results \
    --diamond_db "$DIAMOND_DB"

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
        
        # Check Prodigal results
        if [ -d "results/prodigal_megahit" ]; then
            echo "  ✅ Prodigal MEGAHIT gene predictions: results/prodigal_megahit/"
            MEGAHIT_GENES=$(find results/prodigal_megahit -name "*.faa" | wc -l)
            echo "     - Generated $MEGAHIT_GENES protein sequence files"
        fi
        
        if [ -d "results/prodigal_spades" ]; then
            echo "  ✅ Prodigal SPAdes gene predictions: results/prodigal_spades/"
            SPADES_GENES=$(find results/prodigal_spades -name "*.faa" | wc -l)
            echo "     - Generated $SPADES_GENES protein sequence files"
        fi
        
        # Check Diamond results
        if [ -d "results/diamond_megahit" ]; then
            echo "  ✅ Diamond MEGAHIT results: results/diamond_megahit/"
        fi
        
        if [ -d "results/diamond_spades" ]; then
            echo "  ✅ Diamond SPAdes results: results/diamond_spades/"
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
        echo "  Gene prediction results:"
        find results/prodigal_* -name "*.faa" 2>/dev/null | head -10
        echo ""
        echo "  Diamond classification results:"
        find results/diamond_* -name "*.txt" 2>/dev/null | head -10
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

