#!/usr/bin/env nextflow

/*
 * Metagenome Assembly and Diamond Taxonomic Classification Workflow (English Version)
 * 
 * This workflow integrates:
 * 1. Quality control using fastp (optional)
 * 2. Metagenome assembly using MEGAHIT and SPAdes
 * 3. Gene prediction using Prodigal
 * 4. Taxonomic/functional classification using Diamond BLASTX
 * 5. Comprehensive analysis merging results from both assemblers
 * 
 * Author: Assistant
 * Version: 3.0.0
 */

nextflow.enable.dsl = 2

// Workflow parameters
// Input data
params.input = null
params.outdir = './results'
params.help = false

// MAG parameters (simplified version without these features)
// params.skip_binning = true
// params.skip_checkm = true  
// params.skip_busco = true
// params.skip_prodigal = true
// params.skip_diamond = true
// params.skip_hmmer = true

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

// Diamond classification parameters
params.diamond_db = null
params.diamond_evalue = 1e-5
params.diamond_max_target_seqs = 1
params.diamond_outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids'

// Merge analysis parameters
params.skip_merge_reports = false  // Whether to skip comprehensive report generation

// Resource parameters
params.max_cpus = 32
params.max_memory = '256.GB'
params.max_time = '72.h'

// Print help information
if (params.help) {
    log.info """
    ==========================================
    Metagenome Assembly and Diamond Taxonomic Classification Workflow
    ==========================================
    
    Usage:
    nextflow run metagenome_assembly_classification_workflow_en.nf --input samplesheet.csv --outdir results --diamond_db /path/to/db
    
    Parameters:
    --input                    Input samplesheet
    --outdir                   Output directory
    --diamond_db              Diamond database path (e.g., NCBI nr)
    --diamond_evalue          E-value threshold (default: 1e-5)
    --diamond_max_target_seqs Maximum number of target sequences (default: 1)
    
    Example:
    nextflow run metagenome_assembly_classification_workflow_en.nf \\
        --input samplesheet.csv \\
        --outdir results \\
        --diamond_db /path/to/diamond/nr.dmnd
    """
    exit 0
}

// Validate required parameters
if (!params.input) {
    error "Input samplesheet is required. Use --input parameter."
}

if (!params.diamond_db) {
    error "Diamond database path is required. Use --diamond_db parameter."
}

// Print workflow information
log.info """
==========================================
🧬 Metagenome Assembly and Diamond Taxonomic Classification Workflow
==========================================
Workflow version: 3.0.0
Input samplesheet: ${params.input}
Output directory: ${params.outdir}
Diamond database: ${params.diamond_db}
Diamond E-value: ${params.diamond_evalue}
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
    
    // Stage 2: Gene Prediction with Prodigal
    PRODIGAL_MEGAHIT (
        MEGAHIT_ASSEMBLY.out.contigs
    )
    
    PRODIGAL_SPADES (
        SPADES_ASSEMBLY.out.contigs
    )
    
    // Stage 3: Diamond Classification
    DIAMOND_CLASSIFICATION_MEGAHIT (
        PRODIGAL_MEGAHIT.out.proteins,
        params.diamond_db
    )
    
    DIAMOND_CLASSIFICATION_SPADES (
        PRODIGAL_SPADES.out.proteins,
        params.diamond_db
    )
    
    // Stage 4: Merge Reports (Comprehensive Analysis)
    if (!params.skip_merge_reports) {
        // Merge MEGAHIT and SPAdes reports for the same sample
        DIAMOND_CLASSIFICATION_MEGAHIT.out.diamond_megahit
            .join(DIAMOND_CLASSIFICATION_SPADES.out.diamond_spades)
            .set { ch_reports_to_merge }
        
        // Prepare taxonomy files
        ch_taxonomy_names = Channel.fromPath(params.taxonomy_names, checkIfExists: true)
        ch_taxonomy_nodes = Channel.fromPath(params.taxonomy_nodes, checkIfExists: true)
        
        MERGE_DIAMOND_REPORTS (
            ch_reports_to_merge,
            ch_taxonomy_names.collect(),
            ch_taxonomy_nodes.collect()
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
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("megahit_contigs.fa"), emit: contigs
    
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
    
    cp megahit_output/final.contigs.fa megahit_contigs.fa
    
    echo "MEGAHIT: Generated \$(grep -c ">" megahit_contigs.fa) contigs"
    """
}

// Process: SPAdes Assembly
process SPADES_ASSEMBLY {
    tag "${sample}_SPAdes"
    label 'process_high'
    container 'docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("spades_contigs.fa"), emit: contigs
    
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
    
    cp spades_output/contigs.fasta spades_contigs.fa
    
    echo "metaSPAdes: Generated \$(grep -c ">" spades_contigs.fa) contigs"
    """
}

// Process: Prodigal Gene Prediction for MEGAHIT
process PRODIGAL_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_medium'
    conda 'bioconda::prodigal=2.6.3'
    publishDir "${params.outdir}/prodigal_megahit", mode: 'copy', pattern: "*.faa"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_megahit_proteins.faa"), emit: proteins
    path("${sample}_megahit_genes.fna"), emit: genes
    
    script:
    """
    echo "=== Prodigal Gene Prediction (MEGAHIT): ${sample} ==="
    
    # 运行 Prodigal 进行基因预测（metagenome 模式）
    prodigal \\
        -i ${contigs} \\
        -a ${sample}_megahit_proteins.faa \\
        -d ${sample}_megahit_genes.fna \\
        -p meta \\
        -q
    
    # 统计预测的基因数量
    GENE_COUNT=\$(grep -c ">" ${sample}_megahit_proteins.faa || echo 0)
    echo "Prodigal: Predicted \${GENE_COUNT} genes from MEGAHIT contigs"
    """
}

// Process: Prodigal Gene Prediction for SPAdes
process PRODIGAL_SPADES {
    tag "${sample}_SPAdes"
    label 'process_medium'
    conda 'bioconda::prodigal=2.6.3'
    publishDir "${params.outdir}/prodigal_spades", mode: 'copy', pattern: "*.faa"
    
    input:
    tuple val(sample), path(contigs)
    
    output:
    tuple val(sample), path("${sample}_spades_proteins.faa"), emit: proteins
    path("${sample}_spades_genes.fna"), emit: genes
    
    script:
    """
    echo "=== Prodigal Gene Prediction (SPAdes): ${sample} ==="
    
    # 运行 Prodigal 进行基因预测（metagenome 模式）
    prodigal \\
        -i ${contigs} \\
        -a ${sample}_spades_proteins.faa \\
        -d ${sample}_spades_genes.fna \\
        -p meta \\
        -q
    
    # 统计预测的基因数量
    GENE_COUNT=\$(grep -c ">" ${sample}_spades_proteins.faa || echo 0)
    echo "Prodigal: Predicted \${GENE_COUNT} genes from SPAdes contigs"
    """
}

// Process: Diamond Classification for MEGAHIT
process DIAMOND_CLASSIFICATION_MEGAHIT {
    tag "${sample}_MEGAHIT"
    label 'process_high'
    conda 'bioconda::diamond=2.1.8'
    publishDir "${params.outdir}/diamond_megahit", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(proteins)
    val(diamond_db)
    
    output:
    tuple val(sample), path("${sample}_megahit_diamond.txt"), emit: diamond_megahit
    
    script:
    """
    echo "=== Diamond Classification (MEGAHIT): ${sample} ==="
    
    # 运行 Diamond BLASTX 进行分类
    diamond blastp \\
        --query ${proteins} \\
        --db ${diamond_db} \\
        --out ${sample}_megahit_diamond.txt \\
        --outfmt ${params.diamond_outfmt} \\
        --threads ${task.cpus} \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs ${params.diamond_max_target_seqs} \\
        --sensitive
    
    # 统计匹配结果
    HIT_COUNT=\$(wc -l < ${sample}_megahit_diamond.txt || echo 0)
    echo "Diamond: Found \${HIT_COUNT} hits for MEGAHIT proteins"
    """
}

// Process: Diamond Classification for SPAdes
process DIAMOND_CLASSIFICATION_SPADES {
    tag "${sample}_SPAdes"
    label 'process_high'
    conda 'bioconda::diamond=2.1.8'
    publishDir "${params.outdir}/diamond_spades", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(proteins)
    val(diamond_db)
    
    output:
    tuple val(sample), path("${sample}_spades_diamond.txt"), emit: diamond_spades
    
    script:
    """
    echo "=== Diamond Classification (SPAdes): ${sample} ==="
    
    # 运行 Diamond BLASTX 进行分类
    diamond blastp \\
        --query ${proteins} \\
        --db ${diamond_db} \\
        --out ${sample}_spades_diamond.txt \\
        --outfmt ${params.diamond_outfmt} \\
        --threads ${task.cpus} \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs ${params.diamond_max_target_seqs} \\
        --sensitive
    
    # 统计匹配结果
    HIT_COUNT=\$(wc -l < ${sample}_spades_diamond.txt || echo 0)
    echo "Diamond: Found \${HIT_COUNT} hits for SPAdes proteins"
    """
}

// Process: Merge Diamond Reports (Comprehensive Analysis)
process MERGE_DIAMOND_REPORTS {
    tag "${sample}"
    label 'process_low'
    conda 'conda-forge::pandas=2.0.3'
    publishDir "${params.outdir}/merged_reports", mode: 'copy', pattern: "*"
    
    input:
    tuple val(sample), path(megahit_report), path(spades_report)
    path(taxonomy_names)
    path(taxonomy_nodes)
    
    output:
    tuple val(sample), path("${sample}_merged_report.txt"), emit: merged_report
    path("${sample}_merged_report.csv"), emit: merged_csv
    path("${sample}_megahit_with_taxonomy.txt"), emit: megahit_enhanced
    path("${sample}_spades_with_taxonomy.txt"), emit: spades_enhanced
    
    script:
    """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    
    import pandas as pd
    import sys
    from collections import Counter, defaultdict
    
    # Taxonomy database class
    class TaxonomyDB:
        \"\"\"Parse NCBI taxonomy database\"\"\"
        
        def __init__(self, names_file, nodes_file):
            \"\"\"Initialize taxonomy database\"\"\"
            print("Loading taxonomy database...", file=sys.stderr)
            self.names = self._load_names(names_file)
            self.nodes = self._load_nodes(nodes_file)
            print(f"Loaded: {len(self.names):,} names, {len(self.nodes):,} nodes", file=sys.stderr)
        
        def _load_names(self, names_file):
            \"\"\"Load names.dmp\"\"\"
            names = {}
            with open(names_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 4:
                        taxid = parts[0]
                        name = parts[1]
                        name_class = parts[3]
                        
                        # Only keep scientific names
                        if name_class == 'scientific name':
                            names[taxid] = name
            
            return names
        
        def _load_nodes(self, nodes_file):
            \"\"\"Load nodes.dmp\"\"\"
            nodes = {}
            with open(nodes_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) >= 3:
                        taxid = parts[0]
                        parent_taxid = parts[1]
                        rank = parts[2]
                        
                        nodes[taxid] = {
                            'parent': parent_taxid,
                            'rank': rank
                        }
            
            return nodes
        
        def get_lineage(self, taxid):
            \"\"\"Get complete taxonomic lineage\"\"\"
            lineage = {
                'superkingdom': 'N/A',
                'kingdom': 'N/A',
                'phylum': 'N/A',
                'class': 'N/A',
                'order': 'N/A',
                'family': 'N/A',
                'genus': 'N/A',
                'species': 'N/A',
                'organism_name': 'N/A'
            }
            
            taxid = str(taxid)
            
            # Get current taxid name
            if taxid in self.names:
                lineage['organism_name'] = self.names[taxid]
            
            # Traverse up the taxonomy tree
            current_taxid = taxid
            visited = set()
            
            while current_taxid != '1' and current_taxid in self.nodes:
                # Prevent loops
                if current_taxid in visited:
                    break
                visited.add(current_taxid)
                
                node = self.nodes[current_taxid]
                rank = node['rank']
                
                # Only keep major ranks
                if rank in lineage and current_taxid in self.names:
                    lineage[rank] = self.names[current_taxid]
                
                # Move to parent node
                current_taxid = node['parent']
            
            return lineage
    
    def parse_diamond_output(file_path):
        \"\"\"Parse Diamond output file\"\"\"
        try:
            df = pd.read_csv(file_path, sep='\\t', header=None, 
                           names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                  'evalue', 'bitscore', 'staxids'])
            return df
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                        'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                        'evalue', 'bitscore', 'staxids'])
    
    def add_taxonomy_to_dataframe(df, taxonomy_db):
        \"\"\"Add taxonomic information to DataFrame\"\"\"
        if df.empty:
            return df
        
        print(f"Adding taxonomy info to {len(df):,} records...", file=sys.stderr)
        
        lineages = []
        for idx, row in df.iterrows():
            if (idx + 1) % 1000 == 0:
                print(f"  Processing {idx+1:,}/{len(df):,}...", file=sys.stderr)
            
            taxid = str(row['staxids']) if pd.notna(row['staxids']) else '0'
            lineage = taxonomy_db.get_lineage(taxid)
            lineages.append(lineage)
        
        # Add taxonomy columns
        for key in ['organism_name', 'superkingdom', 'kingdom', 'phylum', 
                    'class', 'order', 'family', 'genus', 'species']:
            df[key] = [lineage[key] for lineage in lineages]
        
        return df
    
    def extract_taxonomic_info(df):
        \"\"\"Extract taxonomic statistics from Diamond output\"\"\"
        if df.empty:
            return {
                'total_hits': 0,
                'unique_queries': 0,
                'unique_subjects': 0,
                'avg_identity': 0,
                'avg_length': 0,
                'taxid_counts': Counter(),
                'phylum_counts': Counter(),
                'family_counts': Counter()
            }
        
        stats = {
            'total_hits': len(df),
            'unique_queries': df['qseqid'].nunique(),
            'unique_subjects': df['sseqid'].nunique(),
            'avg_identity': df['pident'].mean(),
            'avg_length': df['length'].mean(),
            'taxid_counts': Counter(df['staxids'].dropna().astype(str))
        }
        
        # If taxonomy columns exist, count phylum and family
        if 'phylum' in df.columns:
            stats['phylum_counts'] = Counter(df['phylum'].dropna())
            stats['family_counts'] = Counter(df['family'].dropna())
        
        return stats
    
    # Initialize Taxonomy database
    taxonomy_db = TaxonomyDB("${taxonomy_names}", "${taxonomy_nodes}")
    
    # Parse Diamond output files
    print(f"\\nParsing MEGAHIT Diamond results: ${megahit_report}", file=sys.stderr)
    megahit_df = parse_diamond_output("${megahit_report}")
    
    # Add taxonomic information
    megahit_df_enhanced = add_taxonomy_to_dataframe(megahit_df.copy(), taxonomy_db)
    megahit_stats = extract_taxonomic_info(megahit_df_enhanced)
    
    # Save enhanced version (with taxonomy)
    megahit_df_enhanced.to_csv("${sample}_megahit_with_taxonomy.txt", sep='\\t', index=False)
    print(f"MEGAHIT enhanced file saved", file=sys.stderr)
    
    print(f"\\nParsing SPAdes Diamond results: ${spades_report}", file=sys.stderr)
    spades_df = parse_diamond_output("${spades_report}")
    
    # Add taxonomic information
    spades_df_enhanced = add_taxonomy_to_dataframe(spades_df.copy(), taxonomy_db)
    spades_stats = extract_taxonomic_info(spades_df_enhanced)
    
    # Save enhanced version (with taxonomy)
    spades_df_enhanced.to_csv("${sample}_spades_with_taxonomy.txt", sep='\\t', index=False)
    print(f"SPAdes enhanced file saved", file=sys.stderr)
    
    # Generate text report
    with open("${sample}_merged_report.txt", 'w', encoding='utf-8') as f:
        f.write("="*80 + "\\n")
        f.write("Diamond Comprehensive Analysis Report - MEGAHIT vs SPAdes Assembly Comparison\\n")
        f.write("="*80 + "\\n\\n")
        
        # Overall statistics
        f.write("[Overall Statistics]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"SPAdes total hits:          {spades_stats['total_hits']:,}\\n")
        f.write(f"MEGAHIT total hits:         {megahit_stats['total_hits']:,}\\n\\n")
        
        f.write(f"SPAdes unique queries:      {spades_stats['unique_queries']:,}\\n")
        f.write(f"MEGAHIT unique queries:     {megahit_stats['unique_queries']:,}\\n\\n")
        
        f.write(f"SPAdes average identity:    {spades_stats['avg_identity']:.2f}%\\n")
        f.write(f"MEGAHIT average identity:   {megahit_stats['avg_identity']:.2f}%\\n\\n")
        
        f.write(f"SPAdes average length:      {spades_stats['avg_length']:.1f} aa\\n")
        f.write(f"MEGAHIT average length:     {megahit_stats['avg_length']:.1f} aa\\n\\n")
        
        # Phylum level comparison
        f.write("\\n[Phylum Level Comparison (Top 15)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Phylum':<30} {'SPAdes Count':<15} {'MEGAHIT Count':<15} {'Total':<15}\\n")
        f.write("-"*80 + "\\n")
        
        all_phyla = set(spades_stats.get('phylum_counts', {}).keys()) | set(megahit_stats.get('phylum_counts', {}).keys())
        phylum_comparison = []
        for phylum in all_phyla:
            spades_c = spades_stats.get('phylum_counts', {}).get(phylum, 0)
            megahit_c = megahit_stats.get('phylum_counts', {}).get(phylum, 0)
            total_c = spades_c + megahit_c
            phylum_comparison.append((phylum, spades_c, megahit_c, total_c))
        
        phylum_comparison.sort(key=lambda x: x[3], reverse=True)
        for phylum, spades_c, megahit_c, total_c in phylum_comparison[:15]:
            f.write(f"{phylum:<30} {spades_c:<15} {megahit_c:<15} {total_c:<15}\\n")
        
        # Family level comparison
        f.write("\\n[Family Level Comparison (Top 15)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Family':<30} {'SPAdes Count':<15} {'MEGAHIT Count':<15} {'Total':<15}\\n")
        f.write("-"*80 + "\\n")
        
        all_families = set(spades_stats.get('family_counts', {}).keys()) | set(megahit_stats.get('family_counts', {}).keys())
        family_comparison = []
        for family in all_families:
            spades_c = spades_stats.get('family_counts', {}).get(family, 0)
            megahit_c = megahit_stats.get('family_counts', {}).get(family, 0)
            total_c = spades_c + megahit_c
            family_comparison.append((family, spades_c, megahit_c, total_c))
        
        family_comparison.sort(key=lambda x: x[3], reverse=True)
        for family, spades_c, megahit_c, total_c in family_comparison[:15]:
            f.write(f"{family:<30} {spades_c:<15} {megahit_c:<15} {total_c:<15}\\n")
        
        # Taxonomic ID comparison (Top 20)
        f.write("\\n[Taxonomic ID Comparison (Top 20)]\\n")
        f.write("-"*80 + "\\n")
        f.write(f"{'Taxonomic ID':<20} {'SPAdes Count':<15} {'MEGAHIT Count':<15}\\n")
        f.write("-"*80 + "\\n")
        
        # Merge all taxonomic IDs
        all_taxids = set(spades_stats['taxid_counts'].keys()) | set(megahit_stats['taxid_counts'].keys())
        taxid_comparison = []
        
        for taxid in all_taxids:
            spades_count = spades_stats['taxid_counts'].get(taxid, 0)
            megahit_count = megahit_stats['taxid_counts'].get(taxid, 0)
            total_count = spades_count + megahit_count
            taxid_comparison.append((taxid, spades_count, megahit_count, total_count))
        
        # Sort by total count
        taxid_comparison.sort(key=lambda x: x[3], reverse=True)
        
        for taxid, spades_c, megahit_c, total_c in taxid_comparison[:20]:
            f.write(f"{taxid:<20} {spades_c:<15} {megahit_c:<15}\\n")
        
        # Unique findings
        f.write("\\n[Unique Findings]\\n")
        f.write("-"*80 + "\\n\\n")
        
        # Taxonomic IDs found only in SPAdes
        spades_only = set(spades_stats['taxid_counts'].keys()) - set(megahit_stats['taxid_counts'].keys())
        f.write(f"Taxonomic IDs found only in SPAdes: {len(spades_only)}\\n")
        if spades_only:
            spades_only_sorted = sorted([(tid, spades_stats['taxid_counts'][tid]) 
                                        for tid in spades_only], 
                                       key=lambda x: x[1], reverse=True)
            for tid, count in spades_only_sorted[:10]:
                f.write(f"  - {tid}: {count} matches\\n")
        
        f.write("\\n")
        
        # Taxonomic IDs found only in MEGAHIT
        megahit_only = set(megahit_stats['taxid_counts'].keys()) - set(spades_stats['taxid_counts'].keys())
        f.write(f"Taxonomic IDs found only in MEGAHIT: {len(megahit_only)}\\n")
        if megahit_only:
            megahit_only_sorted = sorted([(tid, megahit_stats['taxid_counts'][tid]) 
                                         for tid in megahit_only], 
                                        key=lambda x: x[1], reverse=True)
            for tid, count in megahit_only_sorted[:10]:
                f.write(f"  - {tid}: {count} matches\\n")
        
        f.write("\\n" + "="*80 + "\\n")
        f.write("Analysis Complete\\n")
        f.write("="*80 + "\\n")
    
    # Save CSV format detailed comparison
    comparison_data = []
    for taxid, spades_c, megahit_c, total_c in taxid_comparison:
        comparison_data.append({
            'taxonomic_id': taxid,
            'spades_count': spades_c,
            'megahit_count': megahit_c,
            'total_count': total_c,
            'spades_only': spades_c > 0 and megahit_c == 0,
            'megahit_only': megahit_c > 0 and spades_c == 0,
            'both': spades_c > 0 and megahit_c > 0
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_df = comparison_df.sort_values('total_count', ascending=False)
    comparison_df.to_csv("${sample}_merged_report.csv", index=False)
    
    print(f"Reports generated successfully for ${sample}", file=sys.stderr)
    """
}

// Workflow completion message
workflow.onComplete {
    log.info """
    ==========================================
    🎯 Metagenome Assembly and Diamond Classification Results
    ==========================================
    Pipeline completed successfully!
    
    Results directory: ${params.outdir}
    
    Generated files:
    - fastp/: Quality control reports (if enabled)
      * *_fastp.html: HTML quality reports
      * *_fastp.json: JSON quality data
    - prodigal_megahit/: Gene predictions from MEGAHIT contigs
      * *_megahit_proteins.faa: Predicted protein sequences
    - prodigal_spades/: Gene predictions from SPAdes contigs
      * *_spades_proteins.faa: Predicted protein sequences
    - diamond_megahit/: Diamond classification of MEGAHIT proteins
      * *_megahit_diamond.txt: BLAST-style alignment results
    - diamond_spades/: Diamond classification of SPAdes proteins
      * *_spades_diamond.txt: BLAST-style alignment results
    - merged_reports/: Comprehensive analysis (if enabled)
      * *_merged_report.txt: Combined analysis report
      * *_merged_report.csv: Detailed comparison data
      * *_megahit_with_taxonomy.txt: MEGAHIT results with full taxonomy
      * *_spades_with_taxonomy.txt: SPAdes results with full taxonomy
    
    ==========================================
    """
}

workflow.onError {
    log.error """
    ==========================================
    ❌ Metagenome Assembly and Diamond Classification Workflow Failed
    ==========================================
    Error: ${workflow.errorMessage}
    ==========================================
    """
}

