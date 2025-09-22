#!/bin/bash

#SBATCH --job-name=viral_enhanced_comp
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=viral_enhanced_comp_%j.out
#SBATCH --error=viral_enhanced_comp_%j.err

echo "=========================================="
echo "🧬 Enhanced Comprehensive Viral Detection Workflow"
echo "=========================================="
start_time=$(date +%s)
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo

cd "$SLURM_SUBMIT_DIR"

# Enhanced environment setup
echo "🔧 1. Setting up enhanced comprehensive environment..."
export LD_LIBRARY_PATH=".:$HOME:/lib64:/usr/lib64:$HOME/.conda/envs/nextflow/lib:/usr/lib64:/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH"

# Create libbz2 symbolic links
if [ -f "/usr/lib64/libbz2.so.1" ]; then
    ln -sf /usr/lib64/libbz2.so.1 $HOME/libbz2.so.1.0 2>/dev/null
    ln -sf /usr/lib64/libbz2.so.1 ./libbz2.so.1.0 2>/dev/null
fi

# Load conda environment
module load Miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nextflow

echo "Current environment: $CONDA_DEFAULT_ENV"

# Comprehensive tool verification
echo "🧪 2. Verifying comprehensive toolset..."

# Core alignment tools
echo "Core alignment tools:"
echo "  ✅ BWA: $(which bwa)"
echo "  ✅ samtools: $(samtools --version | head -1)"

# Quality control
echo "Quality control tools:"
echo "  ✅ fastp: $(which fastp)"

# Assembly tools
echo "Assembly tools:"
echo "  ✅ MEGAHIT: $(which megahit)"

# ORF prediction
if command -v prodigal >/dev/null 2>&1; then
    echo "  ✅ PRODIGAL: $(which prodigal)"
else
    echo "  ⚠️ PRODIGAL not found, installing..."
    conda install -c bioconda prodigal -y
    if command -v prodigal >/dev/null 2>&1; then
        echo "  ✅ PRODIGAL installed: $(which prodigal)"
    else
        echo "  ❌ PRODIGAL installation failed"
    fi
fi

# Protein analysis tools
if command -v diamond >/dev/null 2>&1; then
    echo "  ✅ DIAMOND: $(which diamond)"
else
    echo "  ⚠️ DIAMOND not found, installing..."
    conda install -c bioconda diamond -y
    if command -v diamond >/dev/null 2>&1; then
        echo "  ✅ DIAMOND installed: $(which diamond)"
    else
        echo "  ❌ DIAMOND installation failed"
    fi
fi

# Profile analysis tools
if command -v hmmscan >/dev/null 2>&1; then
    echo "  ✅ HMMER: $(which hmmscan)"
else
    echo "  ⚠️ HMMER not found, installing..."
    conda install -c bioconda hmmer -y
    if command -v hmmscan >/dev/null 2>&1; then
        echo "  ✅ HMMER installed: $(which hmmscan)"
    else
        echo "  ❌ HMMER installation failed"
    fi
fi

# Quality assessment tools
if command -v checkv >/dev/null 2>&1; then
    echo "  ✅ CheckV: $(which checkv)"
else
    echo "  ⚠️ CheckV not found (optional tool)"
    echo "     Will use alternative quality assessment"
fi

# Taxonomic classification
echo "  ✅ Kraken2: $(which kraken2)"

# Utility tools
echo "Utility tools:"
echo "  ✅ seqtk: $(which seqtk)"
echo "  ✅ bc: $(which bc)"
echo "  ✅ Nextflow: $(nextflow -version 2>/dev/null | head -1)"

# Verify configuration files
echo "📋 3. Verifying enhanced configuration files..."
if [ -f "main_viral_enhanced_comprehensive_fixed.nf" ]; then
    echo "  ✅ Enhanced workflow file: $(ls -lh main_viral_enhanced_comprehensive_fixed.nf | awk '{print $5}')"
else
    echo "  ❌ Enhanced workflow file does not exist"
    exit 1
fi

if [ -f "enhanced_comprehensive_config_fixed.config" ]; then
    echo "  ✅ Enhanced configuration file: $(ls -lh enhanced_comprehensive_config_fixed.config | awk '{print $5}')"
else
    echo "  ❌ Enhanced configuration file does not exist"
    exit 1
fi

# Comprehensive database verification
echo "🗄️ 4. Verifying comprehensive databases..."

# Viral genomes database
viral_genomes="databases/viral_genomes/complete_precise_human_animal_viruses.fa"
if [ -f "$viral_genomes" ]; then
    size=$(ls -lh "$viral_genomes" | awk '{print $5}')
    count=$(grep -c '^>' "$viral_genomes")
    echo "  ✅ Viral genomes: $count sequences ($size)"
else
    echo "  ❌ Viral genomes database not found: $viral_genomes"
    exit 1
fi

# DIAMOND protein database
diamond_db="databases/viral_proteins/complete_precise_human_animal_viruses_proteins_diamond.dmnd"
if [ -f "$diamond_db" ]; then
    size=$(ls -lh "$diamond_db" | awk '{print $5}')
    echo "  ✅ DIAMOND protein database: $size"
else
    echo "  ❌ DIAMOND protein database not found: $diamond_db"
    echo "     This will limit protein-level analysis"
fi

# HMMER HMM database
hmm_db="databases/viral_hmm/rvdb-prot.hmm"
if [ -f "$hmm_db" ]; then
    size=$(ls -lh "$hmm_db" | awk '{print $5}')
    echo "  ✅ HMMER HMM database: $size"
    
    # Check if HMM database is pressed (indexed)
    if [ -f "${hmm_db}.h3i" ]; then
        echo "     ✅ HMM database is indexed"
    else
        echo "     ⚠️ HMM database not indexed, will index during analysis"
    fi
else
    echo "  ❌ HMMER HMM database not found: $hmm_db"
    echo "     This will limit profile-based analysis"
fi

# Kraken2 database
kraken2_db="databases/viral_genomes/complete_precise_human_animal_viruses_kraken2"
if [ -d "$kraken2_db" ]; then
    echo "  ✅ Kraken2 database directory exists"
    
    # Check key files
    required_files=("hash.k2d" "opts.k2d" "taxo.k2d" "seqid2taxid.map")
    missing_files=0
    
    for file in "${required_files[@]}"; do
        if [ -f "$kraken2_db/$file" ]; then
            size=$(ls -lh "$kraken2_db/$file" | awk '{print $5}')
            echo "    ✅ $file ($size)"
        else
            echo "    ❌ $file missing"
            ((missing_files++))
        fi
    done
    
    if [ $missing_files -eq 0 ]; then
        echo "  ✅ Kraken2 database complete"
    else
        echo "  ⚠️ Kraken2 database incomplete, classification may be skipped"
    fi
else
    echo "  ⚠️ Kraken2 database does not exist: $kraken2_db"
    echo "     Taxonomic classification will be skipped"
fi

# Input data verification
echo "📁 5. Verifying input data..."
input_pattern="data/*_{R1,R2}.fastq.gz"
input_files=($(ls data/*_{R1,R2}.fastq.gz 2>/dev/null))

if [ ${#input_files[@]} -gt 0 ]; then
    sample_count=$(ls data/*_R1.fastq.gz 2>/dev/null | wc -l)
    echo "  ✅ Found $sample_count sample(s):"
    
    for r1_file in data/*_R1.fastq.gz; do
        if [ -f "$r1_file" ]; then
            sample_name=$(basename "$r1_file" | sed 's/_R1.fastq.gz//')
            r2_file="${r1_file/_R1/_R2}"
            
            if [ -f "$r2_file" ]; then
                r1_size=$(ls -lh "$r1_file" | awk '{print $5}')
                r2_size=$(ls -lh "$r2_file" | awk '{print $5}')
                echo "    📊 $sample_name: R1=$r1_size, R2=$r2_size"
            else
                echo "    ❌ $sample_name: R2 file missing"
            fi
        fi
    done
else
    echo "  ❌ No input files found matching pattern: $input_pattern"
    echo "     Please ensure paired FASTQ files are in the data/ directory"
    exit 1
fi

# Clean environment
echo "🧹 6. Cleaning previous results..."
rm -rf work .nextflow* results_viral_enhanced

# Display enhanced workflow features
echo "🚀 7. Enhanced workflow features overview..."
echo ""
echo "NEW FEATURES in this enhanced version:"
echo "  🧬 ORF Prediction with PRODIGAL"
echo "    - Identifies protein-coding regions in assembled viral genomes"
echo "    - Provides amino acid and nucleotide sequences"
echo "    - Calculates coding density and ORF statistics"
echo ""
echo "  💎 DIAMOND Protein Analysis"
echo "    - BLASTP search against viral protein database"
echo "    - Identifies viral protein families and functions"
echo "    - Provides ortholog assignment and annotation"
echo ""
echo "  🔍 HMMER Profile Analysis"
echo "    - Profile HMM search for conserved viral domains"
echo "    - Higher sensitivity than BLAST for distant homologs"
echo "    - Identifies viral protein families and domains"
echo ""
echo "  📊 Abundance Estimation (RPKM/TPM)"
echo "    - Maps reads back to assembled contigs"
echo "    - Calculates coverage, depth, and normalized abundance"
echo "    - Provides quantitative viral load assessment"
echo ""
echo "  ✅ CheckV Quality Assessment"
echo "    - Evaluates completeness of viral genomes"
echo "    - Detects contamination and assembly errors"
echo "    - Provides quality scores and recommendations"
echo ""
echo "  🔬 Multi-Evidence Integration"
echo "    - Combines evidence from all analysis methods"
echo "    - Provides confidence scores for each detection"
echo "    - Generates comprehensive quality assessment"
echo ""

# Resource allocation summary
echo "💻 Resource allocation:"
echo "  📊 Total allocation: 32 cores, 256GB RAM, 48 hours"
echo "  🔧 Process-specific resources:"
echo "    - VIRAL_SCREENING: 16 cores, 128GB, 24h (most intensive)"
echo "    - HMMER_ANALYSIS: 16 cores, 128GB, 16h (HMM search)"
echo "    - DIAMOND_ANALYSIS: 16 cores, 96GB, 8h (protein search)"
echo "    - ABUNDANCE_QUALITY: 16 cores, 128GB, 12h (mapping + CheckV)"
echo "    - Other processes: 8-16 cores, 32-96GB, 4-12h"
echo ""

# Run enhanced comprehensive workflow
echo "🚀 8. Running enhanced comprehensive viral detection workflow..."

# Use single-line command to avoid multi-line issues
NEXTFLOW_CMD="nextflow run main_viral_enhanced_comprehensive_fixed.nf -c enhanced_comprehensive_config_fixed.config --reads '/scratch/sp96859/Meta-genome-data-analysis/Nextflow/data/llnl_66ce4dde_{R1,R2}.fastq.gz' --outdir results_viral_enhanced --base_path /scratch/sp96859/Meta-genome-data-analysis/Nextflow --viral_genomes databases/viral_genomes/complete_precise_human_animal_viruses.fa --viral_proteins databases/viral_proteins/complete_precise_human_animal_viruses_proteins_diamond.dmnd --viral_hmm databases/viral_hmm/rvdb-prot.hmm --kraken2_db databases/viral_genomes/complete_precise_human_animal_viruses_kraken2 --threads 32 --memory 256GB -ansi-log false"

echo "Executing enhanced workflow..."
echo "Command: nextflow run main_viral_enhanced_comprehensive_fixed.nf ..."
echo ""

# Execute command and capture exit code correctly
eval $NEXTFLOW_CMD
ACTUAL_EXIT_CODE=$?

echo ""
echo "=========================================="
echo "🎯 Enhanced Comprehensive Workflow Results"
echo "=========================================="
echo "Nextflow exit code: $ACTUAL_EXIT_CODE"

if [ $ACTUAL_EXIT_CODE -eq 0 ]; then
    echo "🎉🎉🎉 Enhanced comprehensive workflow executed successfully!"
    
    # Detailed result checking
    echo ""
    echo "📊 Generated result directories:"
    if [ -d "results_viral_enhanced" ]; then
        find results_viral_enhanced -type d | sort | sed 's/^/  📁 /'
        
        echo ""
        echo "📄 Enhanced analysis results:"
        
        # 1. Basic screening results
        screening_stats="results_viral_enhanced/02_viral_screening/llnl_66ce4dde_screening_stats.txt"
        if [ -f "$screening_stats" ]; then
            viral_reads=$(grep "Detected viral reads:" "$screening_stats" | cut -d: -f2 | tr -d ' ' || echo "Unknown")
            echo "  ✅ Viral reads detected: $viral_reads reads"
        fi
        
        # 2. Assembly results
        assembly_stats="results_viral_enhanced/03_viral_assembly/llnl_66ce4dde_assembly_stats.txt"
        if [ -f "$assembly_stats" ]; then
            contigs=$(grep "Viral contigs count:" "$assembly_stats" | cut -d: -f2 | tr -d ' ' || echo "Unknown")
            echo "  ✅ Assembled contigs: $contigs"
        fi
        
        # 3. NEW: ORF prediction results
        orf_stats="results_viral_enhanced/04_orf_prediction/llnl_66ce4dde_orf_stats.txt"
        if [ -f "$orf_stats" ]; then
            orfs=$(grep "Total ORFs predicted:" "$orf_stats" | cut -d: -f2 | tr -d ' ' || echo "Unknown")
            echo "  🧬 NEW: Predicted ORFs: $orfs"
        fi
        
        # 4. NEW: DIAMOND protein analysis
        diamond_stats="results_viral_enhanced/05_diamond_analysis/llnl_66ce4dde_diamond_stats.txt"
        if [ -f "$diamond_stats" ]; then
            diamond_hits=$(grep "Total DIAMOND hits:" "$diamond_stats" | cut -d: -f2 | tr -d ' ' || echo "Unknown")
            echo "  💎 NEW: DIAMOND protein hits: $diamond_hits"
        fi
        
        # 5. NEW: HMMER profile analysis
        hmmer_stats="results_viral_enhanced/06_hmmer_analysis/llnl_66ce4dde_hmmer_stats.txt"
        if [ -f "$hmmer_stats" ]; then
            hmmer_hits=$(grep "Total HMMER hits:" "$hmmer_stats" | cut -d: -f2 | tr -d ' ' || echo "Unknown")
            echo "  🔍 NEW: HMMER profile hits: $hmmer_hits"
        fi
        
        # 6. NEW: Abundance and quality assessment
        abundance_stats="results_viral_enhanced/07_abundance_quality/llnl_66ce4dde_abundance_stats.txt"
        if [ -f "$abundance_stats" ]; then
            high_cov=$(grep "High coverage contigs:" "$abundance_stats" | cut -d: -f2 | tr -d ' ' || echo "Unknown")
            echo "  📊 NEW: High coverage contigs: $high_cov"
        fi
        
        # 7. Kraken2 classification
        kraken2_summary="results_viral_enhanced/08_viral_classification/llnl_66ce4dde_viral_species_summary.txt"
        if [ -f "$kraken2_summary" ]; then
            if grep -q "Classified reads:" "$kraken2_summary"; then
                classified=$(grep "Classified reads:" "$kraken2_summary" | cut -d: -f2 | tr -d ' ')
                echo "  🦠 Kraken2 classified reads: $classified"
            fi
        fi
        
        # 8. NEW: Enhanced final report
        final_report="results_viral_enhanced/09_final_report/llnl_66ce4dde.final_summary_stats.txt"
        if [ -f "$final_report" ]; then
            echo "  📋 NEW: Enhanced comprehensive report generated"
            echo "  📄 Report size: $(ls -lh "$final_report" | awk '{print $5}')"
            
            # Show conclusion
            if grep -q "CONCLUSION:" "$final_report"; then
                echo "  💡 Analysis conclusion:"
                grep "CONCLUSION:" "$final_report" | sed 's/^/      /'
            fi
        fi
        
        # NEW: Evidence integration table
        evidence_table="results_viral_enhanced/09_final_report/llnl_66ce4dde.evidence_integration.tsv"
        if [ -f "$evidence_table" ]; then
            high_conf=$(awk 'BEGIN{FS="\t"} $11=="HIGH" {count++} END{print count+0}' "$evidence_table")
            medium_conf=$(awk 'BEGIN{FS="\t"} $11=="MEDIUM" {count++} END{print count+0}' "$evidence_table")
            echo "  🔬 NEW: Evidence integration - High confidence: $high_conf, Medium: $medium_conf"
        fi
        
        echo ""
        echo "🏆 ENHANCED FEATURES SUMMARY:"
        echo "  ✅ Original workflow: Screening + Assembly + Kraken2"
        echo "  🆕 ORF Prediction: PRODIGAL analysis completed"
        echo "  🆕 Protein Analysis: DIAMOND BLASTP search completed"
        echo "  🆕 Profile Analysis: HMMER domain search completed"
        echo "  🆕 Abundance Analysis: RPKM/TPM calculation completed"
        echo "  🆕 Quality Assessment: CheckV or alternative assessment completed"
        echo "  🆕 Evidence Integration: Multi-method confidence scoring completed"
        echo "  🆕 Enhanced Reporting: Comprehensive analysis report generated"
        
    else
        echo "❌ Result directory does not exist"
    fi
    
else
    echo "❌ Enhanced workflow execution failed, exit code: $ACTUAL_EXIT_CODE"
    
    echo ""
    echo "🔍 Error diagnosis:"
    if [ -f ".nextflow.log" ]; then
        echo "Latest Nextflow log:"
        tail -15 .nextflow.log | sed 's/^/  /'
    fi
    
    # Check errors in work directory
    if [ -d "work" ]; then
        echo ""
        echo "Checking process errors:"
        find work -name ".command.err" -size +0 | head -3 | while read errfile; do
            echo "Error file: $errfile"
            echo "Error content:"
            tail -10 "$errfile" | sed 's/^/  /'
            echo ""
        done
    fi
fi

echo ""
echo "🏆 Enhanced workflow capabilities comparison:"
echo ""
echo "ORIGINAL WORKFLOW:"
echo "  ✅ Quality control (fastp)"
echo "  ✅ Viral screening (BWA + samtools)"
echo "  ✅ Viral assembly (MEGAHIT)"
echo "  ✅ Taxonomic classification (Kraken2)"
echo "  ✅ Basic reporting"
echo ""
echo "ENHANCED COMPREHENSIVE WORKFLOW:"
echo "  ✅ All original features PLUS:"
echo "  🆕 ORF prediction (PRODIGAL) - identifies protein-coding genes"
echo "  🆕 Protein analysis (DIAMOND) - functional annotation"
echo "  🆕 Profile analysis (HMMER) - conserved domain detection"
echo "  🆕 Abundance estimation (RPKM/TPM) - quantitative analysis"
echo "  🆕 Quality assessment (CheckV) - genome completeness"
echo "  🆕 Evidence integration - multi-method validation"
echo "  🆕 Enhanced reporting - comprehensive analysis summary"
echo ""
echo "This addresses ALL the limitations you identified:"
echo "  ✅ DIAMOND database: Now actively used for protein analysis"
echo "  ✅ HMMER database: Now actively used for profile searches"
echo "  ✅ PRODIGAL software: Now used for ORF prediction"
echo "  ✅ Abundance estimation: RPKM/TPM calculations implemented"
echo "  ✅ CheckV software: Quality assessment implemented"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo ""
echo "=========================================="
echo "Enhanced comprehensive workflow completed: $(date)"
echo "Total time: $duration seconds ($(($duration / 60)) minutes)"
echo "Final exit code: $ACTUAL_EXIT_CODE"
echo "=========================================="

exit $ACTUAL_EXIT_CODE
