#!/usr/bin/env nextflow

/*
 * Metagenome Viral Classification Workflow (English Version)
 * 
 * This workflow integrates:
 * 1. Quality control using fastp (optional)
 * 2. Metagenome assembly using MEGAHIT and SPAdes (parallel)
 * 3. Viral sequence identification using VirSorter2 and DeepVirFinder
 * 4. Comprehensive comparison and merging of viral identification results
 * 5. Assembler comparison to identify high-confidence consensus viral sequences
 * 
 * Author: Assistant
 * Version: 5.1.0
 */

nextflow.enable.dsl = 2

// Workflow parameters
// Input data
params.input = null
params.outdir = './results'
params.help = false

// Workflow control
params.skip_virsorter2 = false    // Whether to skip VirSorter2 viral identification
params.skip_deepvirfinder = false // Whether to skip DeepVirFinder viral identification
params.skip_merge_reports = false // Whether to skip result merging
params.save_clean_reads = true    // Whether to save filtered clean reads

// Viral identification paths
params.virsorter2_db = null       // VirSorter2 database path
params.deepvirfinder_dir = '/scratch/sp96859/Meta-genome-data-analysis/Apptainer/Contig-based-VirSorter2-DeepVirFinder/DeepVirFinder' // DeepVirFinder installation directory

// MEGAHIT parameters
params.megahit_memory = 0.8
params.megahit_min_contig_len = 1000

// SPAdes parameters (using metaSPAdes)
params.spades_meta = true

// fastp quality control parameters
params.skip_fastp = false
params.fastp_qualified_quality = 20    // Minimum quality value
params.fastp_unqualified_percent = 40  // Maximum percentage of low-quality bases allowed
params.fastp_min_length = 50           // Minimum read length

// VirSorter2 parameters
params.virsorter2_min_length = 1000      // Minimum contig length
params.virsorter2_min_score = 0.5        // Minimum viral score

// DeepVirFinder parameters  
params.deepvirfinder_min_length = 1000   // Minimum contig length
params.deepvirfinder_pvalue = 0.05       // p-value threshold

// Resource parameters
params.max_cpus = 32
params.max_memory = '256.GB'
params.max_time = '72.h'

// Print help information
if (params.help) {
    log.info """
    ==========================================
    🦠 Metagenome Viral Classification Workflow
    ==========================================
    
    Usage:
    nextflow run metagenome_assembly_classification_workflow_en.nf --input samplesheet.csv --outdir results --virsorter2_db /path/to/db
    
    Required Parameters:
    --input                    Input samplesheet (CSV format)
    --outdir                   Output directory
    --virsorter2_db           VirSorter2 database path
    
    Viral Identification Parameters:
    --deepvirfinder_dir       DeepVirFinder installation directory (default: auto-detected)
    --skip_virsorter2         Skip VirSorter2 analysis (default: false)
    --skip_deepvirfinder      Skip DeepVirFinder analysis (default: false)
    --skip_merge_reports      Skip merging VirSorter2 and DeepVirFinder results (default: false)
    --virsorter2_min_length   Minimum contig length for VirSorter2 (default: 1000)
    --virsorter2_min_score    Minimum viral score for VirSorter2 (default: 0.5)
    --deepvirfinder_min_length Minimum contig length for DeepVirFinder (default: 1000)
    --deepvirfinder_pvalue    P-value threshold for DeepVirFinder (default: 0.05)
    
    Optional Parameters:
    --skip_fastp              Skip fastp quality control (default: false)
    --save_clean_reads        Save filtered clean reads (default: true)
    
    Example:
    nextflow run metagenome_assembly_classification_workflow_en.nf \\
        --input samplesheet.csv \\
        --outdir results \\
        --virsorter2_db /scratch/databases/virsorter2/db
    """
    exit 0
}

// Validate required parameters
if (!params.input) {
    error "Input samplesheet is required. Use --input parameter."
}

// Validate VirSorter2 database
if (!params.skip_virsorter2 && !params.virsorter2_db) {
    error "VirSorter2 database path is required. Use --virsorter2_db parameter or --skip_virsorter2 to skip."
}

// Validate DeepVirFinder installation directory
if (!params.skip_deepvirfinder) {
    def dvf_dir = file(params.deepvirfinder_dir)
    if (!dvf_dir.exists() || !dvf_dir.isDirectory()) {
        log.warn "DeepVirFinder directory not found at: ${params.deepvirfinder_dir}"
        log.warn "DeepVirFinder analysis will be skipped."
        params.skip_deepvirfinder = true
    }
}

// Print workflow information
log.info """
==========================================
🦠 Metagenome Viral Classification Workflow
==========================================
Workflow version: 5.1.0
Input samplesheet: ${params.input}
Output directory: ${params.outdir}

Quality Control:
- fastp QC: ${params.skip_fastp ? 'Disabled' : 'Enabled'}
- Save clean reads: ${params.save_clean_reads ? 'Yes' : 'No'}

Assembly Methods:
- MEGAHIT: Enabled
- metaSPAdes: Enabled

Viral Identification:
- VirSorter2: ${params.skip_virsorter2 ? 'Disabled' : 'Enabled'}
${params.skip_virsorter2 ? '' : "  Database: ${params.virsorter2_db}"}
- DeepVirFinder: ${params.skip_deepvirfinder ? 'Disabled' : 'Enabled'}
${params.skip_deepvirfinder ? '' : "  Directory: ${params.deepvirfinder_dir}"}

Result Merging:
- Merge viral reports: ${params.skip_merge_reports ? 'Disabled' : 'Enabled'}
==========================================
"""

// Create input channel from CSV samplesheet
Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row -> 
        def sample = row.sample
        def read1 = file(row.fastq_1)
        def read2 = file(row.fastq_2)
        return tuple(sample, [read1, read2])
    }
    .set { ch_reads }

// Define workflow
workflow {
    // Stage 0: Quality Control (optional)
    if (!params.skip_fastp) {
        FASTP (
            ch_reads
        )
        ch_clean_reads = FASTP.out.clean_reads
    } else {
        ch_clean_reads = ch_reads
    }
    
    // Stage 1: Assembly
    MEGAHIT_ASSEMBLY (
        ch_clean_reads
    )
    
    SPADES_ASSEMBLY (
        ch_clean_reads
    )
    
    // Stage 2: Viral Identification with VirSorter2
    if (!params.skip_virsorter2) {
        VIRSORTER2_MEGAHIT (
            MEGAHIT_ASSEMBLY.out.contigs,
            params.virsorter2_db
        )
        
        VIRSORTER2_SPADES (
            SPADES_ASSEMBLY.out.contigs,
            params.virsorter2_db
        )
    }
    
    // Stage 3: Viral Identification with DeepVirFinder
    if (!params.skip_deepvirfinder) {
        DEEPVIRFINDER_MEGAHIT (
        MEGAHIT_ASSEMBLY.out.contigs
    )
    
        DEEPVIRFINDER_SPADES (
        SPADES_ASSEMBLY.out.contigs
    )
    }
    
    // Stage 4: Merge and Compare Viral Identification Results
    if (!params.skip_merge_reports && !params.skip_virsorter2 && !params.skip_deepvirfinder) {
        // Merge VirSorter2 and DeepVirFinder results for the same sample
        VIRSORTER2_MEGAHIT.out.results
            .join(DEEPVIRFINDER_MEGAHIT.out.results)
            .set { ch_viral_megahit }
        
        VIRSORTER2_SPADES.out.results
            .join(DEEPVIRFINDER_SPADES.out.results)
            .set { ch_viral_spades }
        
        MERGE_VIRAL_REPORTS_MEGAHIT (
            ch_viral_megahit
        )
        
        MERGE_VIRAL_REPORTS_SPADES (
            ch_viral_spades
        )
        
        // Stage 5: Compare MEGAHIT vs SPAdes Results
        MERGE_VIRAL_REPORTS_MEGAHIT.out.merged_report
            .join(MERGE_VIRAL_REPORTS_SPADES.out.merged_report)
            .set { ch_assembler_comparison }
        
        COMPARE_ASSEMBLERS (
            ch_assembler_comparison
        )
    }
}

// ================================================================================
// Process Definitions
// ================================================================================

// Process: fastp Quality Control
process FASTP {
    tag "${sample}"
    label 'process_medium'
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{html,json}"
    publishDir "${params.outdir}/clean_reads", mode: 'copy', pattern: "*_clean_R*.fastq.gz", enabled: params.save_clean_reads
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_clean_R{1,2}.fastq.gz"), emit: clean_reads
    path("${sample}_fastp.html"), emit: html
    path("${sample}_fastp.json"), emit: json
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    echo "=== fastp Quality Control: ${sample} ==="
    
    # List input files for debugging
    echo "Input files in work directory:"
    ls -lh
    
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample}_clean_R1.fastq.gz \\
        -O ${sample}_clean_R2.fastq.gz \\
        --thread ${task.cpus} \\
        --qualified_quality_phred ${params.fastp_qualified_quality} \\
        --unqualified_percent_limit ${params.fastp_unqualified_percent} \\
        --length_required ${params.fastp_min_length} \\
        --detect_adapter_for_pe \\
        --compression 6 \\
        --html ${sample}_fastp.html \\
        --json ${sample}_fastp.json
    
    echo "fastp: Quality control completed for ${sample}"
    """
}

// Process: MEGAHIT Assembly
process MEGAHIT_ASSEMBLY {
    tag "${sample}_MEGAHIT"
    label 'process_high'
    container 'docker://quay.io/biocontainers/megahit:1.2.9--h2e03b76_1'
    publishDir "${params.outdir}/assembly_megahit", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_megahit_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== MEGAHIT Assembly: ${sample} ==="
    
    megahit \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o megahit_output \
        -t ${task.cpus} \
        --memory ${params.megahit_memory} \
        --min-contig-len ${params.megahit_min_contig_len}
    
    cp megahit_output/final.contigs.fa ${sample}_megahit_contigs.fa
    
    echo "MEGAHIT: Generated \$(grep -c ">" ${sample}_megahit_contigs.fa) contigs"
    """
}

// Process: SPAdes Assembly
process SPADES_ASSEMBLY {
    tag "${sample}_SPAdes"
    label 'process_high'
    container 'docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    publishDir "${params.outdir}/assembly_spades", mode: 'copy', pattern: "*.fa"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_spades_contigs.fa"), emit: contigs
    
    script:
    """
    echo "=== metaSPAdes Assembly: ${sample} ==="
    
    # Use metaSPAdes, disable error correction to avoid memory and bug issues
    metaspades.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o spades_output \
        -t ${task.cpus} \
        -m ${task.memory.toGiga()} \
        --only-assembler
    
    cp spades_output/contigs.fasta ${sample}_spades_contigs.fa
    
    echo "metaSPAdes: Generated \$(grep -c ">" ${sample}_spades_contigs.fa) contigs"
    """
}

// Process: VirSorter2 for MEGAHIT
process VIRSORTER2_MEGAHIT {
    tag "${sample}_MEGAHIT_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'  // Use pre-installed environment
    publishDir "${params.outdir}/virsorter2_megahit", mode: 'copy', pattern: "*.{tsv,fa}"
    
    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)
    
    output:
    tuple val(sample), path("${sample}_megahit_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_megahit_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_megahit_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true
    
    script:
    """
    # Ensure correct conda environment is used
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    
    echo "=== VirSorter2 Analysis (MEGAHIT): ${sample} ==="
    echo "Using Python: \$(which python)"
    echo "Using VirSorter2: \$(which virsorter)"
    
    # Run VirSorter2 for viral sequence identification
    virsorter run \\
        -i ${contigs} \\
        -w virsorter2_output \\
        --db-dir ${virsorter2_db} \\
        --min-length ${params.virsorter2_min_length} \\
        --min-score ${params.virsorter2_min_score} \\
        -j ${task.cpus} \\
        all
    
    # Copy result files
    cp virsorter2_output/final-viral-score.tsv ${sample}_megahit_vs2_final-viral-score.tsv
    
    # If viral sequences detected, copy viral contig file
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_megahit_vs2_final-viral-combined.fa
    fi
    
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_megahit_vs2_final-viral-boundary.tsv
    fi
    
    # Count identified viral sequences
    VIRAL_COUNT=\$(tail -n +2 ${sample}_megahit_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: Identified \${VIRAL_COUNT} viral sequences from MEGAHIT contigs"
    """
}

// Process: VirSorter2 for SPAdes
process VIRSORTER2_SPADES {
    tag "${sample}_SPAdes_VirSorter2"
    label 'process_high'
    conda '/home/sp96859/.conda/envs/nextflow_env'  // Use pre-installed environment
    publishDir "${params.outdir}/virsorter2_spades", mode: 'copy', pattern: "*.{tsv,fa}"
    
    input:
    tuple val(sample), path(contigs)
    val(virsorter2_db)
    
    output:
    tuple val(sample), path("${sample}_spades_vs2_final-viral-score.tsv"), emit: results
    tuple val(sample), path("${sample}_spades_vs2_final-viral-combined.fa"), emit: viral_contigs, optional: true
    path("${sample}_spades_vs2_final-viral-boundary.tsv"), emit: boundaries, optional: true
    
    script:
    """
    # Ensure correct conda environment is used
    export PATH="/home/sp96859/.conda/envs/nextflow_env/bin:\$PATH"
    
    echo "=== VirSorter2 Analysis (SPAdes): ${sample} ==="
    echo "Using Python: \$(which python)"
    echo "Using VirSorter2: \$(which virsorter)"
    
    # Run VirSorter2 for viral sequence identification
    virsorter run \\
        -i ${contigs} \\
        -w virsorter2_output \\
        --db-dir ${virsorter2_db} \\
        --min-length ${params.virsorter2_min_length} \\
        --min-score ${params.virsorter2_min_score} \\
        -j ${task.cpus} \\
        all
    
    # Copy result files
    cp virsorter2_output/final-viral-score.tsv ${sample}_spades_vs2_final-viral-score.tsv
    
    # If viral sequences detected, copy viral contig file
    if [ -f virsorter2_output/final-viral-combined.fa ]; then
        cp virsorter2_output/final-viral-combined.fa ${sample}_spades_vs2_final-viral-combined.fa
    fi
    
    if [ -f virsorter2_output/final-viral-boundary.tsv ]; then
        cp virsorter2_output/final-viral-boundary.tsv ${sample}_spades_vs2_final-viral-boundary.tsv
    fi
    
    # Count identified viral sequences
    VIRAL_COUNT=\$(tail -n +2 ${sample}_spades_vs2_final-viral-score.tsv | wc -l || echo 0)
    echo "VirSorter2: Identified \${VIRAL_COUNT} viral sequences from SPAdes contigs"
    """
}

// Process: DeepVirFinder for MEGAHIT
process DEEPVIRFINDER_MEGAHIT {
    tag "${sample}_MEGAHIT_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_megahit", mode: 'copy', pattern: "*.txt"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_megahit_dvf_output.txt"), emit: results
    
    script:
    """
    echo "=== DeepVirFinder Analysis (MEGAHIT): ${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    
    # Activate dvf conda environment (using absolute path)
    set +u  # Temporarily disable undefined variable check (conda requires this)
    
    # Load Miniforge3 module (required for SLURM environment)
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    
    # Absolute path to dvf environment
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    
    # Check if environment exists
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ dvf environment not found: \$DVF_ENV"
        exit 1
    fi
    
    # Get conda base path
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    if [ -z "\$CONDA_BASE" ]; then
        CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    fi
    
    # Initialize conda
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "\$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "❌ Cannot find conda.sh"
        exit 1
    fi
    
    # Activate dvf environment using absolute path
    conda activate "\$DVF_ENV" || { echo "❌ Failed to activate dvf environment"; exit 1; }
    
    # Force update PATH to ensure dvf environment Python is used
    export PATH="\$DVF_ENV/bin:\$PATH"
    export CONDA_PREFIX="\$DVF_ENV"
    export CONDA_DEFAULT_ENV="dvf"
    
    # Clean PYTHONPATH to prevent package pollution from other environments (critical!)
    unset PYTHONPATH
    # Get Python version
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    # Keep only dvf environment's site-packages
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Set Keras to use Theano backend (critical!)
    export KERAS_BACKEND=theano
    
    set -u  # Re-enable
    
    echo "✅ Active conda environment: \$CONDA_DEFAULT_ENV"
    echo "✅ Python path: \$(which python)"
    echo "✅ Python version: \$(python --version)"
    echo "✅ DVF env path: \$DVF_ENV"
    echo "✅ PYTHONPATH: \$PYTHONPATH"
    echo "✅ KERAS_BACKEND: \$KERAS_BACKEND"
    
    # Verify h5py is available
    python -c "import h5py; print('✅ h5py available:', h5py.__version__)" || { echo "❌ h5py not found"; exit 1; }
    
    # Verify keras is available and check backend
    python -c "import os; os.environ['KERAS_BACKEND']='theano'; import keras; print('✅ Keras available:', keras.__version__); print('✅ Keras backend:', keras.backend.backend())" || { echo "❌ Keras not found or backend error"; exit 1; }
    
    # Run DeepVirFinder for viral sequence identification
    python ${params.deepvirfinder_dir}/dvf.py \\
        -i ${contigs} \\
        -o dvf_output \\
        -l ${params.deepvirfinder_min_length} \\
        -c ${task.cpus}
    
    # Copy result files
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_megahit_dvf_output.txt
    
    # Count predicted viral sequences (p-value < threshold)
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_megahit_dvf_output.txt)
    echo "DeepVirFinder: Predicted \${VIRAL_COUNT} viral sequences from MEGAHIT contigs (p-value < ${params.deepvirfinder_pvalue})"
    """
}

// Process: DeepVirFinder for SPAdes
process DEEPVIRFINDER_SPADES {
    tag "${sample}_SPAdes_DeepVirFinder"
    label 'process_high'
    publishDir "${params.outdir}/deepvirfinder_spades", mode: 'copy', pattern: "*.txt"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_spades_dvf_output.txt"), emit: results
    
    script:
    """
    echo "=== DeepVirFinder Analysis (SPAdes): ${sample} ==="
    echo "Using DeepVirFinder from: ${params.deepvirfinder_dir}"
    
    # Activate dvf conda environment (using absolute path)
    set +u  # Temporarily disable undefined variable check (conda requires this)
    
    # Load Miniforge3 module (required for SLURM environment)
    module load Miniforge3/24.11.3-0 2>/dev/null || true
    
    # Absolute path to dvf environment
    DVF_ENV="/home/sp96859/.conda/envs/dvf"
    
    # Check if environment exists
    if [ ! -d "\$DVF_ENV" ]; then
        echo "❌ dvf environment not found: \$DVF_ENV"
        exit 1
    fi
    
    # Get conda base path
    CONDA_BASE=\$(conda info --base 2>/dev/null)
    if [ -z "\$CONDA_BASE" ]; then
        CONDA_BASE="/apps/eb/Miniforge3/24.11.3-0"
    fi
    
    # Initialize conda
    if [ -f "\$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "\$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "❌ Cannot find conda.sh"
        exit 1
    fi
    
    # Activate dvf environment using absolute path
    conda activate "\$DVF_ENV" || { echo "❌ Failed to activate dvf environment"; exit 1; }
    
    # Force update PATH to ensure dvf environment Python is used
    export PATH="\$DVF_ENV/bin:\$PATH"
    export CONDA_PREFIX="\$DVF_ENV"
    export CONDA_DEFAULT_ENV="dvf"
    
    # Clean PYTHONPATH to prevent package pollution from other environments (critical!)
    unset PYTHONPATH
    # Get Python version
    PYTHON_VER=\$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    # Keep only dvf environment's site-packages
    export PYTHONPATH="\$DVF_ENV/lib/python\${PYTHON_VER}/site-packages"
    
    # Set Keras to use Theano backend (critical!)
    export KERAS_BACKEND=theano
    
    set -u  # Re-enable
    
    echo "✅ Active conda environment: \$CONDA_DEFAULT_ENV"
    echo "✅ Python path: \$(which python)"
    echo "✅ Python version: \$(python --version)"
    echo "✅ DVF env path: \$DVF_ENV"
    echo "✅ PYTHONPATH: \$PYTHONPATH"
    echo "✅ KERAS_BACKEND: \$KERAS_BACKEND"
    
    # Verify h5py is available
    python -c "import h5py; print('✅ h5py available:', h5py.__version__)" || { echo "❌ h5py not found"; exit 1; }
    
    # Verify keras is available and check backend
    python -c "import os; os.environ['KERAS_BACKEND']='theano'; import keras; print('✅ Keras available:', keras.__version__); print('✅ Keras backend:', keras.backend.backend())" || { echo "❌ Keras not found or backend error"; exit 1; }
    
    # Run DeepVirFinder for viral sequence identification
    python ${params.deepvirfinder_dir}/dvf.py \\
        -i ${contigs} \\
        -o dvf_output \\
        -l ${params.deepvirfinder_min_length} \\
        -c ${task.cpus}
    
    # Copy result files
    cp dvf_output/${contigs}_gt${params.deepvirfinder_min_length}bp_dvfpred.txt ${sample}_spades_dvf_output.txt
    
    # Count predicted viral sequences (p-value < threshold)
    VIRAL_COUNT=\$(awk -v pval="${params.deepvirfinder_pvalue}" 'NR>1 && \$3<pval {count++} END {print count+0}' ${sample}_spades_dvf_output.txt)
    echo "DeepVirFinder: Predicted \${VIRAL_COUNT} viral sequences from SPAdes contigs (p-value < ${params.deepvirfinder_pvalue})"
    """
}
// Process: Merge Viral Identification Reports for MEGAHIT
// Integrate VirSorter2 and DeepVirFinder results
process MERGE_VIRAL_REPORTS_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)
    
    output:
    tuple val(sample), path("${sample}_megahit_viral_merged_report.txt"), emit: merged_report
    path("${sample}_megahit_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_megahit_viral_consensus.txt"), emit: consensus_list
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    from collections import defaultdict
    
    def parse_virsorter2(file_path):
        \"\"\"
        Parse VirSorter2 output file
        Columns: seqname, dsDNAphage, NCLDV, RNA, ssDNA, lavidaviridae, max_score, max_score_group, length, hallmark, viral_gene, cellular_gene
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            # Extract sequence names and scores
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['seqname']
                max_score = row['max_score']
                viral_dict[seqname] = {
                    'vs2_score': max_score,
                    'vs2_group': row['max_score_group'],
                    'vs2_length': row['length']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse VirSorter2 results: {e}")
            return {}
    
    def parse_deepvirfinder(file_path):
        \"\"\"
        Parse DeepVirFinder output file
        Columns: name, len, score, pvalue
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['name']
                viral_dict[seqname] = {
                    'dvf_score': row['score'],
                    'dvf_pvalue': row['pvalue'],
                    'dvf_length': row['len']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse DeepVirFinder results: {e}")
            return {}
    
    # Parse result files
    print(f"Parsing VirSorter2 results: ${virsorter2_results}")
    vs2_dict = parse_virsorter2("${virsorter2_results}")
    
    print(f"Parsing DeepVirFinder results: ${deepvirfinder_results}")
    dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
    
    # Merge results
    all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
    
    # Generate comprehensive report
    with open("${sample}_megahit_viral_merged_report.txt", 'w', encoding='utf-8') as f:
        f.write("="*80 + "\\n")
        f.write("Viral Identification Comprehensive Analysis Report - MEGAHIT Assembly Results\\n")
        f.write("VirSorter2 + DeepVirFinder\\n")
        f.write("="*80 + "\\n\\n")
    
        # Overall statistics
        f.write("[Overall Statistics]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"VirSorter2 identified viral sequences:    {len(vs2_dict):,}\\n")
        f.write(f"DeepVirFinder identified viral sequences: {len(dvf_dict):,}\\n")
    
        # Count consensus sequences (identified by both tools)
        consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
        f.write(f"Consensus viral sequences (both tools):   {len(consensus):,}\\n")
    
        # Identified by only one tool
        vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
        dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
        f.write(f"VirSorter2 only:                          {len(vs2_only):,}\\n")
        f.write(f"DeepVirFinder only:                       {len(dvf_only):,}\\n\\n")
    
        # DVF significant sequences (p-value < ${params.deepvirfinder_pvalue})
        dvf_significant = [seq for seq, data in dvf_dict.items() 
                          if data['dvf_pvalue'] < ${params.deepvirfinder_pvalue}]
        f.write(f"DeepVirFinder significant sequences (p<${params.deepvirfinder_pvalue}): {len(dvf_significant):,}\\n\\n")
    
        # Consensus sequence details (recommended high-confidence viral sequences)
        f.write("\\n[Consensus Viral Sequences (High Confidence)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Sequence Name':<40} {'VS2 Score':<12} {'DVF Score':<12} {'DVF P-value':<12}\\n")
        f.write("-"*80 + "\\n")
        
        for seq in sorted(consensus):
            vs2_score = vs2_dict[seq]['vs2_score']
            dvf_score = dvf_dict[seq]['dvf_score']
            dvf_pval = dvf_dict[seq]['dvf_pvalue']
            f.write(f"{seq:<40} {vs2_score:<12.3f} {dvf_score:<12.3f} {dvf_pval:<12.2e}\\n")
        
        f.write("\\n" + "="*80 + "\\n")
        f.write("Analysis Complete\\n")
        f.write("="*80 + "\\n")
    
    # Save detailed data in CSV format
    merged_data = []
    for seq in all_sequences:
        row = {
            'sequence_name': seq,
            'identified_by': '',
            'vs2_score': None,
            'vs2_group': None,
            'dvf_score': None,
            'dvf_pvalue': None,
            'consensus': False
        }
        
        if seq in vs2_dict:
            row['vs2_score'] = vs2_dict[seq]['vs2_score']
            row['vs2_group'] = vs2_dict[seq]['vs2_group']
        
        if seq in dvf_dict:
            row['dvf_score'] = dvf_dict[seq]['dvf_score']
            row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
        
        if seq in consensus:
            row['identified_by'] = 'Both'
            row['consensus'] = True
        elif seq in vs2_only:
            row['identified_by'] = 'VirSorter2_only'
        else:
            row['identified_by'] = 'DeepVirFinder_only'
        
        merged_data.append(row)
    
    merged_df = pd.DataFrame(merged_data)
    merged_df.to_csv("${sample}_megahit_viral_merged_report.csv", index=False)
    
    # Save consensus sequence list (recommended for downstream analysis)
    with open("${sample}_megahit_viral_consensus.txt", 'w') as f:
        for seq in sorted(consensus):
            f.write(seq + "\\n")
    
    print(f"Viral identification report generated successfully: ${sample} (MEGAHIT)")
    print(f"Consensus viral sequences: {len(consensus)}")
    """
}

// Process: Merge Viral Identification Reports for SPAdes
process MERGE_VIRAL_REPORTS_SPADES {
    tag "${sample}_SPAdes"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/merged_viral_reports_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(virsorter2_results), path(deepvirfinder_results)
    
    output:
    tuple val(sample), path("${sample}_spades_viral_merged_report.txt"), emit: merged_report
    path("${sample}_spades_viral_merged_report.csv"), emit: merged_csv
    path("${sample}_spades_viral_consensus.txt"), emit: consensus_list
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    from collections import defaultdict
    
    def parse_virsorter2(file_path):
        \"\"\"
        Parse VirSorter2 output file
        Columns: seqname, dsDNAphage, NCLDV, RNA, ssDNA, lavidaviridae, max_score, max_score_group, length, hallmark, viral_gene, cellular_gene
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            # Extract sequence names and scores
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['seqname']
                max_score = row['max_score']
                viral_dict[seqname] = {
                    'vs2_score': max_score,
                    'vs2_group': row['max_score_group'],
                    'vs2_length': row['length']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse VirSorter2 results: {e}")
            return {}
    
    def parse_deepvirfinder(file_path):
        \"\"\"
        Parse DeepVirFinder output file
        Columns: name, len, score, pvalue
        \"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t')
            viral_dict = {}
            for _, row in df.iterrows():
                seqname = row['name']
                viral_dict[seqname] = {
                    'dvf_score': row['score'],
                    'dvf_pvalue': row['pvalue'],
                    'dvf_length': row['len']
                }
            return viral_dict
        except Exception as e:
            print(f"Warning: Failed to parse DeepVirFinder results: {e}")
            return {}
    
    # 解析结果文件
    print(f"Parsing VirSorter2 results: ${virsorter2_results}")
    vs2_dict = parse_virsorter2("${virsorter2_results}")
    
    print(f"Parsing DeepVirFinder results: ${deepvirfinder_results}")
    dvf_dict = parse_deepvirfinder("${deepvirfinder_results}")
    
    # 合并结果
    all_sequences = set(vs2_dict.keys()) | set(dvf_dict.keys())
    
    # 生成综合报告
    with open("${sample}_spades_viral_merged_report.txt", 'w', encoding='utf-8') as f:
        f.write("="*80 + "\\n")
        f.write("病毒识别综合分析报告 - SPAdes组装结果\\n")
        f.write("VirSorter2 + DeepVirFinder\\n")
        f.write("="*80 + "\\n\\n")
        
        # 总体统计
        f.write("[总体统计]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"VirSorter2 识别的病毒序列数:    {len(vs2_dict):,}\\n")
        f.write(f"DeepVirFinder 识别的病毒序列数: {len(dvf_dict):,}\\n")
        
        # 统计共识序列（两个工具都识别为病毒）
        consensus = set(vs2_dict.keys()) & set(dvf_dict.keys())
        f.write(f"共识病毒序列数（两者都识别）:    {len(consensus):,}\\n")
        
        # 仅被一个工具识别
        vs2_only = set(vs2_dict.keys()) - set(dvf_dict.keys())
        dvf_only = set(dvf_dict.keys()) - set(vs2_dict.keys())
        f.write(f"仅VirSorter2识别:              {len(vs2_only):,}\\n")
        f.write(f"仅DeepVirFinder识别:           {len(dvf_only):,}\\n\\n")
        
        # DVF显著序列统计（p-value < ${params.deepvirfinder_pvalue}）
        dvf_significant = [seq for seq, data in dvf_dict.items() 
                          if data['dvf_pvalue'] < ${params.deepvirfinder_pvalue}]
        f.write(f"DeepVirFinder显著序列 (p<${params.deepvirfinder_pvalue}): {len(dvf_significant):,}\\n\\n")
        
        # 共识序列详情（推荐的高可信度病毒序列）
        f.write("\\n[共识病毒序列 (高可信度)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Sequence Name':<40} {'VS2 Score':<12} {'DVF Score':<12} {'DVF P-value':<12}\\n")
        f.write("-"*80 + "\\n")
        
        for seq in sorted(consensus):
            vs2_score = vs2_dict[seq]['vs2_score']
            dvf_score = dvf_dict[seq]['dvf_score']
            dvf_pval = dvf_dict[seq]['dvf_pvalue']
            f.write(f"{seq:<40} {vs2_score:<12.3f} {dvf_score:<12.3f} {dvf_pval:<12.2e}\\n")
        
        f.write("\\n" + "="*80 + "\\n")
        f.write("Analysis Complete\\n")
        f.write("="*80 + "\\n")
    
    # Save detailed data in CSV format
    merged_data = []
    for seq in all_sequences:
        row = {
            'sequence_name': seq,
            'identified_by': '',
            'vs2_score': None,
            'vs2_group': None,
            'dvf_score': None,
            'dvf_pvalue': None,
            'consensus': False
        }
        
        if seq in vs2_dict:
            row['vs2_score'] = vs2_dict[seq]['vs2_score']
            row['vs2_group'] = vs2_dict[seq]['vs2_group']
        
        if seq in dvf_dict:
            row['dvf_score'] = dvf_dict[seq]['dvf_score']
            row['dvf_pvalue'] = dvf_dict[seq]['dvf_pvalue']
        
        if seq in consensus:
            row['identified_by'] = 'Both'
            row['consensus'] = True
        elif seq in vs2_only:
            row['identified_by'] = 'VirSorter2_only'
        else:
            row['identified_by'] = 'DeepVirFinder_only'
        
        merged_data.append(row)
    
    merged_df = pd.DataFrame(merged_data)
    merged_df.to_csv("${sample}_spades_viral_merged_report.csv", index=False)
    
    # 保存共识序列列表（推荐用于下游分析）
    with open("${sample}_spades_viral_consensus.txt", 'w') as f:
        for seq in sorted(consensus):
            f.write(seq + "\\n")
    
    print(f"病毒识别报告生成成功: ${sample} (SPAdes)")
    print(f"共识病毒序列数: {len(consensus)}")
    """
}

// ================================================================================
// Assembler Comparison
// ================================================================================

// Process: Compare MEGAHIT vs SPAdes Viral Identification Results
process COMPARE_ASSEMBLERS {
    tag "${sample}_Assembler_Comparison"
    label 'process_low'
    conda 'pandas numpy'
    publishDir "${params.outdir}/assembler_comparison", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(megahit_report), path(spades_report)
    
    output:
    tuple val(sample), path("${sample}_assembler_comparison.txt"), emit: comparison_report
    path("${sample}_assembler_comparison.csv"), emit: comparison_csv
    path("${sample}_consensus_viral_sequences.txt"), emit: final_consensus
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    from collections import defaultdict
    
    def parse_merged_report(file_path):
        \"\"\"Parse merged viral report CSV file\"\"\"
        try:
            df = pd.read_csv(file_path.replace('.txt', '.csv'))
            return df
        except Exception as e:
            print(f"Warning: Failed to parse {file_path}: {e}")
            return pd.DataFrame()
    
    print(f"Parsing MEGAHIT results: ${megahit_report}")
    megahit_df = parse_merged_report("${megahit_report}")
    
    print(f"Parsing SPAdes results: ${spades_report}")
    spades_df = parse_merged_report("${spades_report}")
    
    # Extract viral sequences
    megahit_seqs = set(megahit_df['sequence'].tolist()) if 'sequence' in megahit_df.columns else set()
    spades_seqs = set(spades_df['sequence'].tolist()) if 'sequence' in spades_df.columns else set()
    
    # Calculate intersection and differences
    consensus_both_assemblers = megahit_seqs & spades_seqs
    megahit_only = megahit_seqs - spades_seqs
    spades_only = spades_seqs - megahit_seqs
    all_viral = megahit_seqs | spades_seqs
    
    # Generate comprehensive comparison report
    with open("${sample}_assembler_comparison.txt", 'w', encoding='utf-8') as f:
        f.write("="*100 + "\\n")
        f.write("Assembler Comparison Report - Viral Identification Results\\n")
        f.write("MEGAHIT vs metaSPAdes\\n")
        f.write("Sample: ${sample}\\n")
        f.write("="*100 + "\\n\\n")
        
        f.write("[Overall Statistics]\\n")
        f.write("-"*100 + "\\n")
        f.write(f"MEGAHIT identified viral sequences:    {len(megahit_seqs):,}\\n")
        f.write(f"SPAdes identified viral sequences:     {len(spades_seqs):,}\\n")
        f.write(f"Total unique viral sequences:          {len(all_viral):,}\\n\\n")
        
        f.write(f"Consensus viral sequences (both):      {len(consensus_both_assemblers):,}\\n")
        f.write(f"MEGAHIT only:                          {len(megahit_only):,}\\n")
        f.write(f"SPAdes only:                           {len(spades_only):,}\\n\\n")
        
        # Calculate consistency
        if len(all_viral) > 0:
            consistency = len(consensus_both_assemblers) / len(all_viral) * 100
            f.write(f"Assembler consistency:                 {consistency:.2f}%\\n\\n")
        
        f.write("="*100 + "\\n")
        f.write("[Recommendation]\\n")
        f.write("-"*100 + "\\n")
        f.write(f"High-confidence viral sequences (identified by both methods): {len(consensus_both_assemblers):,}\\n")
        f.write("Recommend prioritizing these consensus sequences for downstream analysis.\\n\\n")
        
        f.write("[Detailed Analysis]\\n")
        f.write("-"*100 + "\\n")
        
        # MEGAHIT advantages
        if len(megahit_only) > 0:
            f.write(f"\\nMEGAHIT-specific sequences ({len(megahit_only)}):\\n")
            f.write("  - May represent low-coverage or high-complexity regions\\n")
            f.write("  - MEGAHIT has stronger assembly capability for complex structures\\n")
        
        # SPAdes advantages
        if len(spades_only) > 0:
            f.write(f"\\nSPAdes-specific sequences ({len(spades_only)}):\\n")
            f.write("  - May represent high-coverage regions\\n")
            f.write("  - SPAdes kmer strategy may capture more details\\n")
        
        f.write("\\n" + "="*100 + "\\n")
        f.write("[Statistical Summary]\\n")
        f.write("-"*100 + "\\n")
        
        # Calculate consensus ratio for each assembler
        if len(megahit_seqs) > 0:
            megahit_consensus_pct = len(consensus_both_assemblers) / len(megahit_seqs) * 100
            f.write(f"Consensus ratio in MEGAHIT sequences:  {megahit_consensus_pct:.2f}%\\n")
        
        if len(spades_seqs) > 0:
            spades_consensus_pct = len(consensus_both_assemblers) / len(spades_seqs) * 100
            f.write(f"Consensus ratio in SPAdes sequences:   {spades_consensus_pct:.2f}%\\n")
    
    # Generate CSV detailed comparison
    comparison_data = []
    
    for seq in all_viral:
        row = {
            'sequence': seq,
            'found_in_MEGAHIT': 'Yes' if seq in megahit_seqs else 'No',
            'found_in_SPAdes': 'Yes' if seq in spades_seqs else 'No',
            'status': 'Consensus' if seq in consensus_both_assemblers else 'Assembler_specific'
        }
        
        # Add MEGAHIT information
        if seq in megahit_seqs:
            megahit_row = megahit_df[megahit_df['sequence'] == seq].iloc[0]
            row['MEGAHIT_vs2_score'] = megahit_row.get('vs2_score', 'N/A')
            row['MEGAHIT_dvf_score'] = megahit_row.get('dvf_score', 'N/A')
            row['MEGAHIT_identified_by'] = megahit_row.get('identified_by', 'N/A')
        else:
            row['MEGAHIT_vs2_score'] = 'N/A'
            row['MEGAHIT_dvf_score'] = 'N/A'
            row['MEGAHIT_identified_by'] = 'N/A'
        
        # Add SPAdes information
        if seq in spades_seqs:
            spades_row = spades_df[spades_df['sequence'] == seq].iloc[0]
            row['SPAdes_vs2_score'] = spades_row.get('vs2_score', 'N/A')
            row['SPAdes_dvf_score'] = spades_row.get('dvf_score', 'N/A')
            row['SPAdes_identified_by'] = spades_row.get('identified_by', 'N/A')
        else:
            row['SPAdes_vs2_score'] = 'N/A'
            row['SPAdes_dvf_score'] = 'N/A'
            row['SPAdes_identified_by'] = 'N/A'
        
        comparison_data.append(row)
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_df.to_csv("${sample}_assembler_comparison.csv", index=False)
    
    # Save final consensus sequence list
    with open("${sample}_consensus_viral_sequences.txt", 'w') as f:
        f.write("# High-confidence viral sequences (MEGAHIT and SPAdes consensus)\\n")
        f.write(f"# Sample: ${sample}\\n")
        f.write(f"# Consensus sequences: {len(consensus_both_assemblers)}\\n")
        f.write("#\\n")
        for seq in sorted(consensus_both_assemblers):
            f.write(seq + "\\n")
    
    print(f"\\nAssembler comparison complete: ${sample}")
    print(f"  MEGAHIT: {len(megahit_seqs)} viral sequences")
    print(f"  SPAdes:  {len(spades_seqs)} viral sequences")
    print(f"  Consensus: {len(consensus_both_assemblers)} viral sequences")
    if len(all_viral) > 0:
        print(f"  Consistency: {len(consensus_both_assemblers)/len(all_viral)*100:.1f}%")
    """
}

// Workflow completion message
workflow.onComplete {
    log.info """
    ==========================================
    🎯 Metagenome Viral Classification Results
    ==========================================
    Pipeline completed successfully!
    
    Results directory: ${params.outdir}
    
    Generated files:
    - fastp/: Quality control reports
      * *_fastp.html: HTML quality reports
      * *_fastp.json: JSON quality data
      
    - clean_reads/: Filtered clean reads (if save_clean_reads=true)
      * *_clean_R1.fastq.gz: Forward clean reads
      * *_clean_R2.fastq.gz: Reverse clean reads
      
    - assembly_megahit/: MEGAHIT assembly results
      * *_megahit_contigs.fa: Assembled contigs
      
    - assembly_spades/: metaSPAdes assembly results
      * *_spades_contigs.fa: Assembled contigs
      
    - virsorter2_megahit/: VirSorter2 viral identification (MEGAHIT)
      * *_megahit_vs2_final-viral-score.tsv: Viral scores
      * *_megahit_vs2_final-viral-combined.fa: Identified viral contigs
      
    - virsorter2_spades/: VirSorter2 viral identification (SPAdes)
      * *_spades_vs2_final-viral-score.tsv: Viral scores
      * *_spades_vs2_final-viral-combined.fa: Identified viral contigs
      
    - deepvirfinder_megahit/: DeepVirFinder viral prediction (MEGAHIT)
      * *_megahit_dvf_output.txt: Prediction results with scores and p-values
      
    - deepvirfinder_spades/: DeepVirFinder viral prediction (SPAdes)
      * *_spades_dvf_output.txt: Prediction results with scores and p-values
      
    - merged_viral_reports_megahit/: Integrated viral analysis (MEGAHIT)
      * *_megahit_viral_merged_report.txt: Comprehensive viral identification report
      * *_megahit_viral_merged_report.csv: Detailed comparison data
      * *_megahit_viral_consensus.txt: High-confidence viral sequences list
      
    - merged_viral_reports_spades/: Integrated viral analysis (SPAdes)
      * *_spades_viral_merged_report.txt: Comprehensive viral identification report
      * *_spades_viral_merged_report.csv: Detailed comparison data
      * *_spades_viral_consensus.txt: High-confidence viral sequences list
    
    - assembler_comparison/: MEGAHIT vs SPAdes comparison ⭐
      * *_assembler_comparison.txt: Comprehensive assembler comparison report
      * *_assembler_comparison.csv: Detailed comparison data
      * *_consensus_viral_sequences.txt: Final high-confidence viral sequences (both assemblers)
    
    ==========================================
    """
}

workflow.onError {
    log.error """
    ==========================================
    ❌ Metagenome Viral Classification Workflow Failed
    ==========================================
    Error: ${workflow.errorMessage}
    ==========================================
    """
}

