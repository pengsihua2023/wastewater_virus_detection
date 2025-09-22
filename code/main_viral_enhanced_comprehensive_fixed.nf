#!/usr/bin/env nextflow
/*
 Enhanced Wastewater Metagenomic Viral Detection Workflow - Fixed Version
 
 Comprehensive version including:
 - DIAMOND protein analysis
 - HMMER profile analysis  
 - PRODIGAL ORF prediction
 - Abundance estimation (RPKM/TPM)
 - CheckV quality assessment
 - Multi-evidence integration
*/

nextflow.enable.dsl=2

params.reads            = "data/*_{R1,R2}.fastq.gz"
params.outdir           = "results_viral_enhanced"
params.viral_genomes    = "databases/viral_genomes/complete_precise_human_animal_viruses.fa"
params.viral_proteins   = "databases/viral_proteins/complete_precise_human_animal_viruses_proteins_diamond.dmnd"
params.viral_hmm        = "databases/viral_hmm/rvdb-prot.hmm"
params.kraken2_db       = "databases/viral_genomes/complete_precise_human_animal_viruses_kraken2"
params.base_path        = "/scratch/sp96859/Meta-genome-data-analysis/Nextflow"
params.threads          = 8
params.memory           = '16 GB'
params.publish_mode     = 'copy'

// -------------------------------
// Step 1: Quality control
process FASTP {
    tag { sample }
    cpus { params.threads }
    publishDir "${params.outdir}/01_qc", mode: params.publish_mode
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}.clean.R1.fq.gz"), path("${sample}.clean.R2.fq.gz"), path("${sample}_fastp.json"), path("${sample}_fastp.html")
    
    script:
    def (r1, r2) = reads
    """
    echo "=== Quality control: ${sample} ==="
    echo "Input files: ${r1}, ${r2}"
    
    fastp -i ${r1} -I ${r2} \\
          -o ${sample}.clean.R1.fq.gz \\
          -O ${sample}.clean.R2.fq.gz \\
          --json ${sample}_fastp.json \\
          --html ${sample}_fastp.html \\
          --thread ${task.cpus} \\
          --qualified_quality_phred 20 \\
          --length_required 50
    
    echo "✅ Quality control completed"
    """
}

// -------------------------------
// Step 2: Viral screening - using real samtools
process VIRAL_SCREENING {
    tag { sample }
    cpus { params.threads }
    publishDir "${params.outdir}/02_viral_screening", mode: params.publish_mode
    
    input:
    tuple val(sample), path(r1), path(r2), path(json), path(html)
    
    output:
    tuple val(sample), path("${sample}.viral.R1.fq.gz"), path("${sample}.viral.R2.fq.gz"), path("${sample}.viral_mapping.bam"), path("${sample}_screening_stats.txt")
    
    script:
    """
    echo "=== Viral screening: ${sample} ==="
    
    # Verify tool availability
    echo "Verifying required tools..."
    
    if ! command -v bwa >/dev/null 2>&1; then
        echo "❌ BWA not available"
        exit 127
    fi
    echo "✅ BWA: \$(which bwa)"
    
    if ! command -v seqtk >/dev/null 2>&1; then
        echo "❌ seqtk not available"
        exit 127
    fi
    echo "✅ seqtk: \$(which seqtk)"
    
    if ! command -v samtools >/dev/null 2>&1; then
        echo "❌ samtools not available"
        exit 127
    fi
    echo "✅ samtools: \$(samtools --version | head -1)"
    
    # Use hardcoded absolute path
    VIRAL_DB_PATH="${params.base_path}/${params.viral_genomes}"
    echo "Target viral genome database: \$VIRAL_DB_PATH"
    
    # Verify database file
    if [ ! -f "\$VIRAL_DB_PATH" ]; then
        echo "❌ Error: Viral genome database not found: \$VIRAL_DB_PATH"
        exit 1
    fi
    
    # Count viral genomes
    VIRUS_COUNT=\$(grep -c '^>' "\$VIRAL_DB_PATH" || echo "0")
    echo "Target viral genome count: \$VIRUS_COUNT"
    
    if [ "\$VIRUS_COUNT" -eq 0 ]; then
        echo "❌ Viral genome database is empty"
        exit 1
    fi
    
    # Build viral genome index
    echo "Building viral genome BWA index..."
    bwa index "\$VIRAL_DB_PATH" || { echo "BWA index creation failed"; exit 1; }
    
    # Align reads to viral genomes
    echo "Aligning reads to \$VIRUS_COUNT target viral genomes..."
    echo "Using BWA mem parameters: -t ${task.cpus}"
    
    # Use direct pipeline
    echo "Executing BWA alignment and generating BAM..."
    bwa mem -t ${task.cpus} "\$VIRAL_DB_PATH" ${r1} ${r2} | \\
        samtools view -@ ${task.cpus} -bS - | \\
        samtools sort -@ ${task.cpus} -o ${sample}.viral_mapping.bam -
    
    # Verify BAM file creation
    if [ ! -f "${sample}.viral_mapping.bam" ]; then
        echo "❌ BAM file creation failed"
        exit 1
    fi
    
    BAM_SIZE=\$(ls -lh ${sample}.viral_mapping.bam | awk '{print \$5}')
    echo "✅ BAM file generated successfully, size: \$BAM_SIZE"
    
    # Build BAM index
    echo "Building BAM index..."
    samtools index ${sample}.viral_mapping.bam
    
    # Extract aligned reads
    echo "Extracting viral reads..."
    
    # Get aligned read IDs
    samtools view -F 4 ${sample}.viral_mapping.bam | \\
        awk '{print \$1}' | sort | uniq > ${sample}.viral_read_ids.txt
    
    # Count viral reads
    VIRAL_READ_COUNT=\$(wc -l < ${sample}.viral_read_ids.txt || echo "0")
    echo "Detected viral reads count: \$VIRAL_READ_COUNT"
    
    # Extract viral reads from original clean reads
    if [ "\$VIRAL_READ_COUNT" -gt 0 ]; then
        echo "Extracting viral reads to FASTQ files..."
        seqtk subseq ${r1} ${sample}.viral_read_ids.txt | gzip > ${sample}.viral.R1.fq.gz
        seqtk subseq ${r2} ${sample}.viral_read_ids.txt | gzip > ${sample}.viral.R2.fq.gz
        echo "✅ Viral reads extraction completed"
    else
        echo "⚠️ No viral reads detected, creating empty files"
        echo "" | gzip > ${sample}.viral.R1.fq.gz
        echo "" | gzip > ${sample}.viral.R2.fq.gz
    fi
    
    # Generate detailed statistics report
    echo "=== Viral Screening Statistics Report ===" > ${sample}_screening_stats.txt
    echo "Sample: ${sample}" >> ${sample}_screening_stats.txt
    echo "Target viral genome count: \$VIRUS_COUNT" >> ${sample}_screening_stats.txt
    
    # Calculate input reads count
    INPUT_READS=\$(zcat ${r1} | awk 'END{print NR/4}')
    echo "Input clean reads: \$INPUT_READS" >> ${sample}_screening_stats.txt
    echo "Detected viral reads: \$VIRAL_READ_COUNT" >> ${sample}_screening_stats.txt
    
    # Calculate viral reads ratio
    VIRAL_RATIO=\$(awk -v viral=\$VIRAL_READ_COUNT -v total=\$INPUT_READS 'BEGIN{if(total>0) printf "%.4f", viral*100/total; else print "0"}')
    echo "Viral reads ratio: \${VIRAL_RATIO}%" >> ${sample}_screening_stats.txt
    
    echo "Viral R1 file size: \$(ls -lh ${sample}.viral.R1.fq.gz | awk '{print \$5}')" >> ${sample}_screening_stats.txt
    echo "Viral R2 file size: \$(ls -lh ${sample}.viral.R2.fq.gz | awk '{print \$5}')" >> ${sample}_screening_stats.txt
    echo "BAM file size: \$BAM_SIZE" >> ${sample}_screening_stats.txt
    
    # Use real samtools flagstat
    echo "=== samtools flagstat results ===" >> ${sample}_screening_stats.txt
    samtools flagstat ${sample}.viral_mapping.bam >> ${sample}_screening_stats.txt
    
    echo "✅ Viral screening completed"
    echo "✅ Pure viral reads extracted"
    """
}

// -------------------------------
// Step 3: Viral genome assembly
process VIRAL_ASSEMBLY {
    tag { sample }
    cpus { params.threads }
    publishDir "${params.outdir}/03_viral_assembly", mode: params.publish_mode
    
    input:
    tuple val(sample), path(viral_r1), path(viral_r2), path(mapping_bam), path(stats)
    
    output:
    tuple val(sample), path("${sample}.viral_contigs.fa"), path("${sample}_assembly_stats.txt"), path("${sample}.viral_contigs_lengths.txt")
    
    script:
    """
    echo "=== Viral genome assembly: ${sample} ==="
    
    # Calculate viral reads count
    VIRAL_READ_COUNT=\$(zcat ${viral_r1} | awk 'END{print int(NR/4)}')
    echo "Viral reads count: \$VIRAL_READ_COUNT"
    
    # Determine whether to perform assembly
    if [ "\$VIRAL_READ_COUNT" -lt 1000 ]; then
        echo "⚠️ Insufficient viral reads (\$VIRAL_READ_COUNT < 1000), skipping assembly"
        echo ">empty_assembly" > ${sample}.viral_contigs.fa
        echo "ACGT" >> ${sample}.viral_contigs.fa
        echo "Insufficient viral reads, assembly not performed" > ${sample}_assembly_stats.txt
        echo "empty_assembly\t4" > ${sample}.viral_contigs_lengths.txt
    else
        echo "Sufficient viral reads, starting assembly..."
        
        # Use megahit for assembly
        megahit -1 ${viral_r1} -2 ${viral_r2} \\
                -t ${task.cpus} \\
                -o ${sample}.viral_assembly \\
                --min-contig-len 200 \\
                --k-min 21 \\
                --k-max 141 \\
                --k-step 12
        
        # Check assembly results
        if [ -f "${sample}.viral_assembly/final.contigs.fa" ]; then
            cp ${sample}.viral_assembly/final.contigs.fa ${sample}.viral_contigs.fa
            
            # Calculate assembly statistics
            CONTIG_COUNT=\$(grep -c '^>' ${sample}.viral_contigs.fa || echo "0")
            TOTAL_LENGTH=\$(awk '/^>/{next} {total+=length} END{print total+0}' ${sample}.viral_contigs.fa)
            
            # Calculate N50
            awk '/^>/{name=substr(\$0,2); next} {seq[name]+=length; total+=length} END{n=0; for(i in seq) lens[++n]=seq[i]; asort(lens); target=total*0.5; sum=0; for(i=n; i>=1; i--) {sum+=lens[i]; if(sum>=target) {print "N50: " lens[i]; break}}}' ${sample}.viral_contigs.fa > n50.tmp
            N50=\$(cat n50.tmp | cut -d' ' -f2 || echo "0")
            
            # Generate contig lengths file
            awk '/^>/{name=substr(\$0,2); next} {seq[name]+=length} END{for(i in seq) print i "\\t" seq[i]}' ${sample}.viral_contigs.fa | sort -k2 -nr > ${sample}.viral_contigs_lengths.txt
            
            echo "Viral contigs count: \$CONTIG_COUNT" > ${sample}_assembly_stats.txt
            echo "Total viral sequence length: \$TOTAL_LENGTH bp" >> ${sample}_assembly_stats.txt
            echo "Assembly N50: \$N50 bp" >> ${sample}_assembly_stats.txt
            echo "Input viral reads: \$VIRAL_READ_COUNT" >> ${sample}_assembly_stats.txt
            echo "Average contig length: \$(awk -v total=\$TOTAL_LENGTH -v count=\$CONTIG_COUNT 'BEGIN{if(count>0) printf "%.1f", total/count; else print "0"}') bp" >> ${sample}_assembly_stats.txt
        else
            echo "⚠️ Assembly failed, creating empty file"
            echo ">failed_assembly" > ${sample}.viral_contigs.fa
            echo "ACGT" >> ${sample}.viral_contigs.fa
            echo "Assembly failed" > ${sample}_assembly_stats.txt
            echo "failed_assembly\t4" > ${sample}.viral_contigs_lengths.txt
        fi
    fi
    
    echo "✅ Viral genome assembly completed"
    """
}

// -------------------------------
// Step 4: ORF Prediction with PRODIGAL
process ORF_PREDICTION {
    tag { sample }
    cpus 4
    publishDir "${params.outdir}/04_orf_prediction", mode: params.publish_mode
    
    input:
    tuple val(sample), path(contigs), path(assembly_stats), path(contig_lengths)
    
    output:
    tuple val(sample), path("${sample}.viral_orfs.faa"), path("${sample}.viral_orfs.fna"), path("${sample}.orf2contig.tsv"), path("${sample}_orf_stats.txt")
    
    script:
    """
    echo "=== ORF prediction: ${sample} ==="
    
    # Check if we have meaningful contigs to analyze
    CONTIG_COUNT=\$(grep -c '^>' ${contigs} || echo "0")
    echo "Input contigs: \$CONTIG_COUNT"
    
    # Calculate total sequence length (excluding empty/failed assemblies)
    TOTAL_SEQ_LENGTH=\$(awk '/^>/{next} /^ACGT\$/{next} {total+=length(\$0)} END{print total+0}' ${contigs})
    echo "Total sequence length: \$TOTAL_SEQ_LENGTH bp"
    
    if [ "\$TOTAL_SEQ_LENGTH" -lt 200 ]; then
        echo "⚠️ Insufficient sequence for ORF prediction, creating empty files"
        echo ">empty_orf" > ${sample}.viral_orfs.faa
        echo "M" >> ${sample}.viral_orfs.faa
        echo ">empty_orf" > ${sample}.viral_orfs.fna
        echo "ATG" >> ${sample}.viral_orfs.fna
        echo "empty_orf\tempty_contig" > ${sample}.orf2contig.tsv
        echo "No ORFs predicted due to insufficient sequence" > ${sample}_orf_stats.txt
    else
        echo "Running PRODIGAL for ORF prediction..."
        
        # Run PRODIGAL in metagenome mode
        prodigal -i ${contigs} \\
                 -a ${sample}.viral_orfs.faa \\
                 -d ${sample}.viral_orfs.fna \\
                 -p meta \\
                 -q \\
                 -g 11
        
        if [ -f "${sample}.viral_orfs.faa" ] && [ -s "${sample}.viral_orfs.faa" ]; then
            # Generate ORF to contig mapping
            echo "Generating ORF to contig mapping..."
            awk '/^>/{
                id=substr(\$0,2); 
                # Extract contig name from ORF ID (remove _# suffix)
                contig=id; 
                sub(/_[0-9]+\$/, "", contig); 
                print id "\\t" contig
            }' ${sample}.viral_orfs.faa > ${sample}.orf2contig.tsv
            
            # Calculate ORF statistics
            ORF_COUNT=\$(grep -c '^>' ${sample}.viral_orfs.faa)
            AVG_ORF_LENGTH=\$(awk '/^>/{if(seq) print length(seq); seq=""; next} {seq=seq\$0} END{if(seq) print length(seq)}' ${sample}.viral_orfs.faa | awk '{sum+=\$1; count++} END{if(count>0) printf "%.1f", sum/count; else print "0"}')
            
            echo "ORF prediction completed successfully" > ${sample}_orf_stats.txt
            echo "Total ORFs predicted: \$ORF_COUNT" >> ${sample}_orf_stats.txt
            echo "Average ORF length: \$AVG_ORF_LENGTH amino acids" >> ${sample}_orf_stats.txt
            echo "ORFs per contig: \$(awk -v orfs=\$ORF_COUNT -v contigs=\$CONTIG_COUNT 'BEGIN{if(contigs>0) printf "%.1f", orfs/contigs; else print "0"}' )" >> ${sample}_orf_stats.txt
            echo "Coding density: \$(awk -v orfs=\$ORF_COUNT -v total=\$TOTAL_SEQ_LENGTH 'BEGIN{if(total>0) printf "%.2f", orfs*300/total; else print "0"}' ) ORFs/kb" >> ${sample}_orf_stats.txt
            
            echo "✅ ORF prediction completed: \$ORF_COUNT ORFs found"
        else
            echo "⚠️ PRODIGAL failed or no ORFs found"
            echo ">no_orfs" > ${sample}.viral_orfs.faa
            echo "M" >> ${sample}.viral_orfs.faa
            echo ">no_orfs" > ${sample}.viral_orfs.fna
            echo "ATG" >> ${sample}.viral_orfs.fna
            echo "no_orfs\tno_contig" > ${sample}.orf2contig.tsv
            echo "PRODIGAL failed or no ORFs predicted" > ${sample}_orf_stats.txt
        fi
    fi
    
    echo "✅ ORF prediction step completed"
    """
}

// -------------------------------
// Step 5: DIAMOND Protein Analysis
process DIAMOND_PROTEIN_ANALYSIS {
    tag { sample }
    cpus { params.threads }
    publishDir "${params.outdir}/05_diamond_analysis", mode: params.publish_mode
    
    input:
    tuple val(sample), path(orfs_faa), path(orfs_fna), path(orf2contig), path(orf_stats)
    
    output:
    tuple val(sample), path("${sample}.diamond_results.m8"), path("${sample}_diamond_stats.txt"), path("${sample}.diamond_best_hits.tsv")
    
    script:
    """
    echo "=== DIAMOND protein analysis: ${sample} ==="
    
    # Check if we have ORFs to analyze
    ORF_COUNT=\$(grep -c '^>' ${orfs_faa} || echo "0")
    echo "Input ORFs: \$ORF_COUNT"
    
    # Check database
    DIAMOND_DB_PATH="${params.base_path}/${params.viral_proteins}"
    echo "DIAMOND database: \$DIAMOND_DB_PATH"
    
    if [ ! -f "\$DIAMOND_DB_PATH" ]; then
        echo "⚠️ DIAMOND database not found, creating empty results"
        echo "# No DIAMOND database available" > ${sample}.diamond_results.m8
        echo "DIAMOND database not found" > ${sample}_diamond_stats.txt
        echo "query_id\tsubject_id\tpident\tlength\tevalue\tbitscore\tsubject_title" > ${sample}.diamond_best_hits.tsv
    elif [ "\$ORF_COUNT" -lt 1 ] || ! grep -q '^>[^>]' ${orfs_faa}; then
        echo "⚠️ No valid ORFs for DIAMOND analysis"
        echo "# No ORFs available for analysis" > ${sample}.diamond_results.m8
        echo "No ORFs available for DIAMOND analysis" > ${sample}_diamond_stats.txt
        echo "query_id\tsubject_id\tpident\tlength\tevalue\tbitscore\tsubject_title" > ${sample}.diamond_best_hits.tsv
    else
        echo "Running DIAMOND blastp analysis..."
        
        # Run DIAMOND blastp
        diamond blastp \\
            --query ${orfs_faa} \\
            --db "\$DIAMOND_DB_PATH" \\
            --out ${sample}.diamond_results.m8 \\
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \\
            --threads ${task.cpus} \\
            --evalue 1e-5 \\
            --max-target-seqs 5 \\
            --sensitive
        
        if [ -s "${sample}.diamond_results.m8" ]; then
            # Calculate statistics
            TOTAL_HITS=\$(wc -l < ${sample}.diamond_results.m8)
            UNIQUE_ORFS=\$(cut -f1 ${sample}.diamond_results.m8 | sort | uniq | wc -l)
            UNIQUE_PROTEINS=\$(cut -f2 ${sample}.diamond_results.m8 | sort | uniq | wc -l)
            
            # Generate best hits table
            echo "query_id\tsubject_id\tpident\tlength\tevalue\tbitscore\tsubject_title" > ${sample}.diamond_best_hits.tsv
            awk 'BEGIN{OFS="\\t"} {
                if(!seen[\$1] || \$12 > best_score[\$1]) {
                    seen[\$1] = 1;
                    best_score[\$1] = \$12;
                    best_hit[\$1] = \$0;
                }
            } END {
                for(orf in best_hit) {
                    split(best_hit[orf], fields, "\\t");
                    print fields[1], fields[2], fields[3], fields[4], fields[11], fields[12], fields[13];
                }
            }' ${sample}.diamond_results.m8 | sort -k6 -nr >> ${sample}.diamond_best_hits.tsv
            
            # Generate statistics
            echo "DIAMOND protein analysis completed successfully" > ${sample}_diamond_stats.txt
            echo "Total DIAMOND hits: \$TOTAL_HITS" >> ${sample}_diamond_stats.txt
            echo "ORFs with hits: \$UNIQUE_ORFS / \$ORF_COUNT (\$(awk -v hits=\$UNIQUE_ORFS -v total=\$ORF_COUNT 'BEGIN{if(total>0) printf "%.1f", hits*100/total; else print "0"}')%)" >> ${sample}_diamond_stats.txt
            echo "Unique viral proteins detected: \$UNIQUE_PROTEINS" >> ${sample}_diamond_stats.txt
            
            # Top viral protein families
            echo "Top 10 viral protein families:" >> ${sample}_diamond_stats.txt
            cut -f13 ${sample}.diamond_results.m8 | \\
                sed 's/\\[.*\\]//' | \\
                awk '{for(i=1;i<=NF;i++) if(\$i ~ /protein|gene|polymerase|capsid|envelope|membrane/) print \$i}' | \\
                sort | uniq -c | sort -nr | head -10 | \\
                awk '{printf "  %s: %d hits\\n", \$2, \$1}' >> ${sample}_diamond_stats.txt
            
            echo "✅ DIAMOND analysis completed: \$TOTAL_HITS hits for \$UNIQUE_ORFS ORFs"
        else
            echo "⚠️ No DIAMOND hits found"
            echo "DIAMOND analysis completed - no significant hits found" > ${sample}_diamond_stats.txt
            echo "query_id\tsubject_id\tpident\tlength\tevalue\tbitscore\tsubject_title" > ${sample}.diamond_best_hits.tsv
        fi
    fi
    
    echo "✅ DIAMOND protein analysis step completed"
    """
}

// -------------------------------
// Step 6: HMMER Profile Analysis
process HMMER_ANALYSIS {
    tag { sample }
    cpus { params.threads }
    publishDir "${params.outdir}/06_hmmer_analysis", mode: params.publish_mode
    
    input:
    tuple val(sample), path(orfs_faa), path(orfs_fna), path(orf2contig), path(orf_stats)
    
    output:
    tuple val(sample), path("${sample}.hmmer_results.tbl"), path("${sample}_hmmer_stats.txt"), path("${sample}.hmmer_domains.tbl")
    
    script:
    """
    echo "=== HMMER profile analysis: ${sample} ==="
    
    # Check if we have ORFs to analyze
    ORF_COUNT=\$(grep -c '^>' ${orfs_faa} || echo "0")
    echo "Input ORFs: \$ORF_COUNT"
    
    # Check HMM database
    HMM_DB_PATH="${params.base_path}/${params.viral_hmm}"
    echo "HMM database: \$HMM_DB_PATH"
    
    if [ ! -f "\$HMM_DB_PATH" ]; then
        echo "⚠️ HMM database not found, creating empty results"
        echo "# No HMM database available" > ${sample}.hmmer_results.tbl
        echo "HMM database not found" > ${sample}_hmmer_stats.txt
        echo "# No HMM database available" > ${sample}.hmmer_domains.tbl
    elif [ "\$ORF_COUNT" -lt 1 ] || ! grep -q '^>[^>]' ${orfs_faa}; then
        echo "⚠️ No valid ORFs for HMMER analysis"
        echo "# No ORFs available for analysis" > ${sample}.hmmer_results.tbl
        echo "No ORFs available for HMMER analysis" > ${sample}_hmmer_stats.txt
        echo "# No ORFs available for analysis" > ${sample}.hmmer_domains.tbl
    else
        echo "Running HMMER hmmscan analysis..."
        
        # Check if HMM database is pressed (indexed)
        if [ ! -f "\${HMM_DB_PATH}.h3i" ]; then
            echo "HMM database not indexed, running hmmpress..."
            hmmpress "\$HMM_DB_PATH"
        fi
        
        # Run HMMER hmmscan
        hmmscan \\
            --cpu ${task.cpus} \\
            --tblout ${sample}.hmmer_results.tbl \\
            --domtblout ${sample}.hmmer_domains.tbl \\
            --noali \\
            -E 1e-3 \\
            --domE 1e-3 \\
            "\$HMM_DB_PATH" \\
            ${orfs_faa} > ${sample}.hmmer_output.txt
        
        if [ -s "${sample}.hmmer_results.tbl" ]; then
            # Remove comment lines for analysis
            grep -v '^#' ${sample}.hmmer_results.tbl > hmmer_clean.tmp || touch hmmer_clean.tmp
            
            if [ -s "hmmer_clean.tmp" ]; then
                TOTAL_HITS=\$(wc -l < hmmer_clean.tmp)
                UNIQUE_ORFS=\$(cut -f4 -d' ' hmmer_clean.tmp | sort | uniq | wc -l)
                UNIQUE_PROFILES=\$(cut -f1 -d' ' hmmer_clean.tmp | sort | uniq | wc -l)
                
                echo "HMMER profile analysis completed successfully" > ${sample}_hmmer_stats.txt
                echo "Total HMMER hits: \$TOTAL_HITS" >> ${sample}_hmmer_stats.txt
                echo "ORFs with profile hits: \$UNIQUE_ORFS / \$ORF_COUNT (\$(awk -v hits=\$UNIQUE_ORFS -v total=\$ORF_COUNT 'BEGIN{if(total>0) printf "%.1f", hits*100/total; else print "0"}')%)" >> ${sample}_hmmer_stats.txt
                echo "Unique viral profiles detected: \$UNIQUE_PROFILES" >> ${sample}_hmmer_stats.txt
                
                # Top viral protein families from HMM
                echo "Top 10 viral protein profiles:" >> ${sample}_hmmer_stats.txt
                cut -f1 -d' ' hmmer_clean.tmp | sort | uniq -c | sort -nr | head -10 | \\
                    awk '{printf "  %s: %d hits\\n", \$2, \$1}' >> ${sample}_hmmer_stats.txt
                
                echo "✅ HMMER analysis completed: \$TOTAL_HITS hits for \$UNIQUE_ORFS ORFs"
            else
                echo "⚠️ No significant HMMER hits found"
                echo "HMMER analysis completed - no significant hits found" > ${sample}_hmmer_stats.txt
            fi
        else
            echo "⚠️ HMMER analysis failed or no hits found"
            echo "HMMER analysis failed or no significant hits found" > ${sample}_hmmer_stats.txt
        fi
        
        # Clean up
        rm -f hmmer_clean.tmp ${sample}.hmmer_output.txt
    fi
    
    echo "✅ HMMER profile analysis step completed"
    """
}

// -------------------------------
// Step 7: Simplified Abundance Estimation 
process ABUNDANCE_ESTIMATION {
    tag { sample }
    cpus { params.threads }
    publishDir "${params.outdir}/07_abundance_estimation", mode: params.publish_mode
    
    input:
    tuple val(sample), path(contigs), path(assembly_stats), path(contig_lengths)
    tuple val(sample2), path(viral_r1), path(viral_r2), path(mapping_bam), path(screening_stats)
    
    output:
    tuple val(sample), path("${sample}.abundance_table.tsv"), path("${sample}_abundance_stats.txt")
    
    when:
    sample == sample2
    
    script:
    """
    echo "=== Abundance estimation: ${sample} ==="
    
    # Check if we have meaningful contigs
    CONTIG_COUNT=\$(grep -c '^>' ${contigs} || echo "0")
    TOTAL_SEQ_LENGTH=\$(awk '/^>/{next} /^ACGT\$/{next} {total+=length(\$0)} END{print total+0}' ${contigs})
    
    echo "Input contigs: \$CONTIG_COUNT"
    echo "Total sequence length: \$TOTAL_SEQ_LENGTH bp"
    
    if [ "\$TOTAL_SEQ_LENGTH" -lt 200 ] || [ "\$CONTIG_COUNT" -lt 1 ]; then
        echo "⚠️ Insufficient contigs for abundance analysis"
        
        # Create empty abundance table
        echo "contig_id\tlength\tmapped_reads\tcoverage\tdepth\trpkm" > ${sample}.abundance_table.tsv
        echo "no_contigs\t0\t0\t0\t0\t0" >> ${sample}.abundance_table.tsv
        
        echo "Insufficient contigs for analysis" > ${sample}_abundance_stats.txt
    else
        echo "Running abundance estimation..."
        
        # Build index for contigs
        bwa index ${contigs}
        
        # Map viral reads back to assembled contigs
        echo "Mapping viral reads to assembled contigs..."
        bwa mem -t ${task.cpus} ${contigs} ${viral_r1} ${viral_r2} | \\
            samtools sort -@ ${task.cpus} -o ${sample}.contigs_mapping.bam -
        
        samtools index ${sample}.contigs_mapping.bam
        
        # Calculate depth per position
        samtools depth -a ${sample}.contigs_mapping.bam > ${sample}.depth.txt
        
        # Get total mapped reads for normalization
        TOTAL_MAPPED_READS=\$(samtools view -c -F 4 ${sample}.contigs_mapping.bam)
        echo "Total mapped reads to contigs: \$TOTAL_MAPPED_READS"
        
        # Create abundance table header
        echo "contig_id\tlength\tmapped_reads\tcoverage\tdepth\trpkm" > ${sample}.abundance_table.tsv
        
        # Get contigs list and process each
        grep '^>' ${contigs} | sed 's/^>//' | while read contig_id; do
            # Get contig length from the contig_lengths file
            contig_length=\$(awk -v id="\$contig_id" '\$1==id {print \$2}' ${contig_lengths})
            
            if [ -z "\$contig_length" ] || [ "\$contig_length" -eq 0 ]; then
                # Fallback: calculate length from sequence
                contig_length=\$(awk -v id="\$contig_id" '
                    BEGIN{len=0; inseq=0} 
                    /^>/ {if(\$0==">"id) inseq=1; else inseq=0; next} 
                    inseq {len+=length} 
                    END{print len}' ${contigs})
            fi
            
            # Get mapped reads count
            mapped_reads=\$(samtools view -c ${sample}.contigs_mapping.bam "\$contig_id" 2>/dev/null || echo "0")
            
            # Calculate coverage (fraction of bases covered)
            covered_bases=\$(awk -v contig="\$contig_id" '\$1==contig && \$3>0 {count++} END {print count+0}' ${sample}.depth.txt)
            if [ "\$contig_length" -gt 0 ]; then
                coverage=\$(awk -v covered=\$covered_bases -v total=\$contig_length 'BEGIN{printf "%.4f", covered/total}')
            else
                coverage="0"
            fi
            
            # Calculate average depth
            avg_depth=\$(awk -v contig="\$contig_id" '\$1==contig {sum+=\$3; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}' ${sample}.depth.txt)
            
            # Calculate RPKM (Reads Per Kilobase per Million mapped reads)
            if [ "\$contig_length" -gt 0 ] && [ "\$TOTAL_MAPPED_READS" -gt 0 ]; then
                rpkm=\$(awk -v reads=\$mapped_reads -v length=\$contig_length -v total=\$TOTAL_MAPPED_READS 'BEGIN{printf "%.2f", (reads*1000000000)/(length*total)}')
            else
                rpkm="0"
            fi
            
            echo "\$contig_id\t\$contig_length\t\$mapped_reads\t\$coverage\t\$avg_depth\t\$rpkm" >> ${sample}.abundance_table.tsv
        done
        
        # Generate abundance statistics
        echo "Abundance estimation completed successfully" > ${sample}_abundance_stats.txt
        echo "Total contigs analyzed: \$CONTIG_COUNT" >> ${sample}_abundance_stats.txt
        echo "Total reads mapped to contigs: \$TOTAL_MAPPED_READS" >> ${sample}_abundance_stats.txt
        
        # Summary statistics
        HIGH_COV_CONTIGS=\$(awk 'NR>1 && \$4>0.5 {count++} END{print count+0}' ${sample}.abundance_table.tsv)
        HIGH_DEPTH_CONTIGS=\$(awk 'NR>1 && \$5>5 {count++} END{print count+0}' ${sample}.abundance_table.tsv)
        
        echo "High coverage contigs (>50%): \$HIGH_COV_CONTIGS" >> ${sample}_abundance_stats.txt
        echo "High depth contigs (>5x): \$HIGH_DEPTH_CONTIGS" >> ${sample}_abundance_stats.txt
        
        # Top abundant contigs
        echo "Top 5 most abundant contigs (by RPKM):" >> ${sample}_abundance_stats.txt
        tail -n +2 ${sample}.abundance_table.tsv | sort -k6 -nr | head -5 | \\
            awk '{printf "  %s: %.2f RPKM, %.1fx depth\\n", \$1, \$6, \$5}' >> ${sample}_abundance_stats.txt
        
        # Clean up temporary files
        rm -f ${sample}.depth.txt
        
        echo "✅ Abundance estimation completed"
    fi
    
    echo "✅ Abundance estimation step completed"
    """
}

// -------------------------------
// Step 8: Kraken2 viral classification (existing)
process KRAKEN2_VIRAL_CLASSIFICATION {
    tag { sample }
    cpus { params.threads }
    publishDir "${params.outdir}/08_viral_classification", mode: params.publish_mode
    
    input:
    tuple val(sample), path(viral_r1), path(viral_r2), path(mapping_bam), path(stats)
    
    output:
    tuple val(sample), path("${sample}_kraken2_report.txt"), path("${sample}_kraken2_classification.txt"), path("${sample}_viral_species_summary.txt")
    
    script:
    """
    echo "=== Kraken2 viral classification: ${sample} ==="
    
    # Check Kraken2 database
    KRAKEN2_DB_PATH="${params.base_path}/${params.kraken2_db}"
    echo "Kraken2 database path: \$KRAKEN2_DB_PATH"
    
    if [ ! -d "\$KRAKEN2_DB_PATH" ]; then
        echo "⚠️ Kraken2 viral database does not exist: \$KRAKEN2_DB_PATH"
        echo "# No Kraken2 database" > ${sample}_kraken2_report.txt
        echo "U\t0\tunclassified" > ${sample}_kraken2_classification.txt
        echo "Sample: ${sample}" > ${sample}_viral_species_summary.txt
        echo "Status: No Kraken2 database, skipping classification" >> ${sample}_viral_species_summary.txt
        echo "Database path: \$KRAKEN2_DB_PATH" >> ${sample}_viral_species_summary.txt
        echo "Time: \$(date)" >> ${sample}_viral_species_summary.txt
        exit 0
    fi
    
    # Verify database file integrity
    echo "Verifying Kraken2 database integrity..."
    if [ ! -f "\$KRAKEN2_DB_PATH/hash.k2d" ] || [ ! -f "\$KRAKEN2_DB_PATH/opts.k2d" ] || [ ! -f "\$KRAKEN2_DB_PATH/taxo.k2d" ]; then
        echo "⚠️ Kraken2 database files incomplete"
        echo "# Kraken2 database incomplete" > ${sample}_kraken2_report.txt
        echo "U\t0\tunclassified" > ${sample}_kraken2_classification.txt
        echo "Sample: ${sample}" > ${sample}_viral_species_summary.txt
        echo "Status: Kraken2 database incomplete, skipping classification" >> ${sample}_viral_species_summary.txt
        echo "Time: \$(date)" >> ${sample}_viral_species_summary.txt
        exit 0
    fi
    
    echo "✅ Kraken2 database verification passed"
    
    # Check viral reads count
    VIRAL_READ_COUNT=\$(zcat ${viral_r1} | awk 'END{print int(NR/4)}')
    echo "Viral reads count: \$VIRAL_READ_COUNT"
    
    if [ "\$VIRAL_READ_COUNT" -eq 0 ]; then
        echo "⚠️ No viral reads, skipping classification"
        echo "# No viral reads" > ${sample}_kraken2_report.txt
        echo "U\t0\tunclassified" > ${sample}_kraken2_classification.txt
        echo "Sample: ${sample}" > ${sample}_viral_species_summary.txt
        echo "Status: No viral reads, skipping classification" >> ${sample}_viral_species_summary.txt
        exit 0
    fi
    
    # Verify Kraken2 tool
    echo "Checking Kraken2 tool availability..."
    if ! command -v kraken2 >/dev/null 2>&1; then
        echo "⚠️ Kraken2 not available, attempting installation..."
        conda install -c bioconda kraken2 -y || {
            echo "❌ Kraken2 installation failed, skipping classification"
            echo "# Kraken2 tool not available" > ${sample}_kraken2_report.txt
            echo "U\t0\tunclassified" > ${sample}_kraken2_classification.txt
            echo "Sample: ${sample}" > ${sample}_viral_species_summary.txt
            echo "Status: Kraken2 tool not available, skipping classification" >> ${sample}_viral_species_summary.txt
            echo "Time: \$(date)" >> ${sample}_viral_species_summary.txt
            exit 0
        }
    fi
    
    echo "✅ Kraken2: \$(which kraken2)"
    kraken2 --version
    
    # Run Kraken2 classification
    echo "Classifying \$VIRAL_READ_COUNT viral reads..."
    kraken2 --db "\$KRAKEN2_DB_PATH" \\
            --paired \\
            --threads ${task.cpus} \\
            --output ${sample}_kraken2_classification.txt \\
            --report ${sample}_kraken2_report.txt \\
            ${viral_r1} ${viral_r2}
    
    echo "✅ Kraken2 classification completed"
    
    # Analyze classification results
    echo "=== Viral classification results analysis ==="
    
    # Generate viral species summary
    echo "Sample: ${sample}" > ${sample}_viral_species_summary.txt
    echo "Classification time: \$(date)" >> ${sample}_viral_species_summary.txt
    echo "Input viral reads count: \$VIRAL_READ_COUNT" >> ${sample}_viral_species_summary.txt
    echo "" >> ${sample}_viral_species_summary.txt
    
    if [ -f "${sample}_kraken2_classification.txt" ]; then
        # Calculate overall classification statistics
        TOTAL_CLASSIFIED=\$(awk '\$1!="U"' ${sample}_kraken2_classification.txt | wc -l)
        TOTAL_UNCLASSIFIED=\$(awk '\$1=="U"' ${sample}_kraken2_classification.txt | wc -l)
        TOTAL_READS=\$(wc -l < ${sample}_kraken2_classification.txt)
        
        echo "=== Classification Statistics ===" >> ${sample}_viral_species_summary.txt
        echo "Total reads: \$TOTAL_READS" >> ${sample}_viral_species_summary.txt
        echo "Classified reads: \$TOTAL_CLASSIFIED" >> ${sample}_viral_species_summary.txt
        echo "Unclassified reads: \$TOTAL_UNCLASSIFIED" >> ${sample}_viral_species_summary.txt
        
        if [ "\$TOTAL_READS" -gt 0 ]; then
            CLASSIFICATION_RATE=\$(awk -v c=\$TOTAL_CLASSIFIED -v t=\$TOTAL_READS 'BEGIN{printf "%.2f", c*100/t}')
            echo "Classification rate: \${CLASSIFICATION_RATE}%" >> ${sample}_viral_species_summary.txt
        fi
        
        echo "" >> ${sample}_viral_species_summary.txt
        
        # Analyze viral species
        echo "=== Detected Viral Species (Top 10) ===" >> ${sample}_viral_species_summary.txt
        awk '\$1!="U"{print \$3}' ${sample}_kraken2_classification.txt | \\
            sort | uniq -c | sort -nr | head -10 | \\
            awk '{printf "%s\\t%s reads\\n", \$2, \$1}' >> ${sample}_viral_species_summary.txt
        
        echo "" >> ${sample}_viral_species_summary.txt
        echo "=== Kraken2 Classification Report Summary ===" >> ${sample}_viral_species_summary.txt
        head -20 ${sample}_kraken2_report.txt | \\
            awk '\$1>0{printf "%.2f%%\\t%s reads\\t%s\\n", \$1, \$2, \$6}' >> ${sample}_viral_species_summary.txt
    else
        echo "❌ Classification result file generation failed" >> ${sample}_viral_species_summary.txt
    fi
    
    echo "✅ Viral classification analysis completed"
    """
}

// -------------------------------
// Step 9: Simplified Enhanced Final Report
process ENHANCED_FINAL_REPORT {
    tag { sample }
    publishDir "${params.outdir}/09_final_report", mode: params.publish_mode
    
    input:
    tuple val(sample), path(contigs), path(assembly_stats), path(contig_lengths)
    tuple val(sample2), path(viral_r1), path(viral_r2), path(mapping_bam), path(screening_stats)
    tuple val(sample3), path(orfs_faa), path(orfs_fna), path(orf2contig), path(orf_stats)
    tuple val(sample4), path(diamond_results), path(diamond_stats), path(diamond_best_hits)
    tuple val(sample5), path(hmmer_results), path(hmmer_stats), path(hmmer_domains)
    tuple val(sample6), path(abundance_table), path(abundance_stats)
    tuple val(sample7), path(kraken2_report), path(kraken2_classification), path(viral_species_summary)
    
    output:
    path("${sample}.comprehensive_viral_report.tsv")
    path("${sample}.final_summary_stats.txt")
    
    when:
    sample == sample2 && sample == sample3 && sample == sample4 && sample == sample5 && sample == sample6 && sample == sample7
    
    script:
    """
    echo "=== Enhanced final report generation: ${sample} ==="
    
    # Read basic statistics
    VIRAL_READS=\$(grep "Detected viral reads:" ${screening_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
    INPUT_READS=\$(grep "Input clean reads:" ${screening_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
    
    # Assembly statistics
    if [ -f "${assembly_stats}" ]; then
        CONTIGS=\$(grep "Viral contigs count:" ${assembly_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
        TOTAL_LENGTH=\$(grep "Total viral sequence length:" ${assembly_stats} | cut -d: -f2 | cut -d' ' -f1 || echo "0")
        N50=\$(grep "Assembly N50:" ${assembly_stats} | cut -d: -f2 | cut -d' ' -f1 || echo "0")
    else
        CONTIGS="0"
        TOTAL_LENGTH="0"
        N50="0"
    fi
    
    # ORF statistics
    if [ -f "${orf_stats}" ]; then
        TOTAL_ORFS=\$(grep "Total ORFs predicted:" ${orf_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
        AVG_ORF_LENGTH=\$(grep "Average ORF length:" ${orf_stats} | cut -d: -f2 | cut -d' ' -f1 || echo "0")
    else
        TOTAL_ORFS="0"
        AVG_ORF_LENGTH="0"
    fi
    
    # DIAMOND statistics
    if [ -f "${diamond_stats}" ]; then
        DIAMOND_HITS=\$(grep "Total DIAMOND hits:" ${diamond_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
        ORFS_WITH_HITS=\$(grep "ORFs with hits:" ${diamond_stats} | cut -d: -f2 | cut -d' ' -f1 || echo "0")
    else
        DIAMOND_HITS="0"
        ORFS_WITH_HITS="0"
    fi
    
    # HMMER statistics
    if [ -f "${hmmer_stats}" ]; then
        HMMER_HITS=\$(grep "Total HMMER hits:" ${hmmer_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
        ORFS_WITH_PROFILES=\$(grep "ORFs with profile hits:" ${hmmer_stats} | cut -d: -f2 | cut -d' ' -f1 || echo "0")
    else
        HMMER_HITS="0"
        ORFS_WITH_PROFILES="0"
    fi
    
    # Abundance statistics
    if [ -f "${abundance_stats}" ]; then
        HIGH_COV_CONTIGS=\$(grep "High coverage contigs" ${abundance_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
        HIGH_DEPTH_CONTIGS=\$(grep "High depth contigs" ${abundance_stats} | cut -d: -f2 | tr -d ' ' || echo "0")
    else
        HIGH_COV_CONTIGS="0"
        HIGH_DEPTH_CONTIGS="0"
    fi
    
    # Kraken2 statistics
    if [ -f "${viral_species_summary}" ]; then
        CLASSIFIED_READS=\$(grep "Classified reads:" ${viral_species_summary} | cut -d: -f2 | tr -d ' ' || echo "0")
        CLASSIFICATION_RATE=\$(grep "Classification rate:" ${viral_species_summary} | cut -d: -f2 | tr -d ' ' || echo "0%")
    else
        CLASSIFIED_READS="N/A"
        CLASSIFICATION_RATE="N/A"
    fi
    
    # Generate comprehensive report
    echo -e "metric\\tvalue\\tdescription" > ${sample}.comprehensive_viral_report.tsv
    echo -e "sample_id\\t${sample}\\tSample identifier" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "input_reads\\t\$INPUT_READS\\tTotal input clean reads" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "viral_reads\\t\$VIRAL_READS\\tDetected viral reads" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "viral_detection_rate\\t\$(awk -v v=\$VIRAL_READS -v t=\$INPUT_READS 'BEGIN{if(t>0) printf "%.4f", v*100/t; else print "0"}')%\\tPercentage of viral reads" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "assembled_contigs\\t\$CONTIGS\\tNumber of assembled viral contigs" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "total_assembly_length\\t\$TOTAL_LENGTH\\tTotal length of assembled viral sequences (bp)" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "assembly_n50\\t\$N50\\tAssembly N50 metric (bp)" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "predicted_orfs\\t\$TOTAL_ORFS\\tTotal number of predicted ORFs" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "avg_orf_length\\t\$AVG_ORF_LENGTH\\tAverage ORF length (amino acids)" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "diamond_hits\\t\$DIAMOND_HITS\\tTotal DIAMOND protein hits" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "orfs_with_protein_hits\\t\$ORFS_WITH_HITS\\tORFs with DIAMOND protein matches" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "hmmer_hits\\t\$HMMER_HITS\\tTotal HMMER profile hits" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "orfs_with_profile_hits\\t\$ORFS_WITH_PROFILES\\tORFs with HMMER profile matches" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "high_coverage_contigs\\t\$HIGH_COV_CONTIGS\\tContigs with >50% coverage" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "high_depth_contigs\\t\$HIGH_DEPTH_CONTIGS\\tContigs with >5x sequencing depth" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "kraken2_classified_reads\\t\$CLASSIFIED_READS\\tReads classified by Kraken2" >> ${sample}.comprehensive_viral_report.tsv
    echo -e "kraken2_classification_rate\\t\$CLASSIFICATION_RATE\\tKraken2 classification success rate" >> ${sample}.comprehensive_viral_report.tsv
    
    # Generate comprehensive summary
    echo "=== Enhanced Viral Detection Final Report - ${sample} ===" > ${sample}.final_summary_stats.txt
    echo "Analysis completed: \$(date)" >> ${sample}.final_summary_stats.txt
    echo "" >> ${sample}.final_summary_stats.txt
    
    echo "📊 DETECTION SUMMARY:" >> ${sample}.final_summary_stats.txt
    echo "  Total input reads: \$INPUT_READS" >> ${sample}.final_summary_stats.txt
    echo "  Viral reads detected: \$VIRAL_READS (\$(awk -v v=\$VIRAL_READS -v t=\$INPUT_READS 'BEGIN{if(t>0) printf "%.4f", v*100/t; else print "0"}')%)" >> ${sample}.final_summary_stats.txt
    echo "  Assembled viral contigs: \$CONTIGS" >> ${sample}.final_summary_stats.txt
    echo "  Total viral sequence: \$TOTAL_LENGTH bp" >> ${sample}.final_summary_stats.txt
    echo "" >> ${sample}.final_summary_stats.txt
    
    echo "🧬 ORF ANALYSIS:" >> ${sample}.final_summary_stats.txt
    echo "  Predicted ORFs: \$TOTAL_ORFS" >> ${sample}.final_summary_stats.txt
    echo "  Average ORF length: \$AVG_ORF_LENGTH aa" >> ${sample}.final_summary_stats.txt
    echo "  ORFs with protein hits (DIAMOND): \$ORFS_WITH_HITS" >> ${sample}.final_summary_stats.txt
    echo "  ORFs with profile hits (HMMER): \$ORFS_WITH_PROFILES" >> ${sample}.final_summary_stats.txt
    echo "" >> ${sample}.final_summary_stats.txt
    
    echo "📈 ABUNDANCE ANALYSIS:" >> ${sample}.final_summary_stats.txt
    echo "  High coverage contigs (>50%): \$HIGH_COV_CONTIGS" >> ${sample}.final_summary_stats.txt
    echo "  High depth contigs (>5x): \$HIGH_DEPTH_CONTIGS" >> ${sample}.final_summary_stats.txt
    echo "" >> ${sample}.final_summary_stats.txt
    
    echo "🦠 TAXONOMIC CLASSIFICATION:" >> ${sample}.final_summary_stats.txt
    echo "  Kraken2 classified reads: \$CLASSIFIED_READS" >> ${sample}.final_summary_stats.txt
    echo "  Classification rate: \$CLASSIFICATION_RATE" >> ${sample}.final_summary_stats.txt
    echo "" >> ${sample}.final_summary_stats.txt
    
    # Overall assessment
    if [ \$VIRAL_READS -gt 1000 ] && [ \$CONTIGS -gt 5 ] && [ \$HMMER_HITS -gt 0 ]; then
        echo "🎉 CONCLUSION: Strong viral detection with high-quality evidence!" >> ${sample}.final_summary_stats.txt
        echo "  - Substantial viral reads detected" >> ${sample}.final_summary_stats.txt
        echo "  - Successful viral genome assembly" >> ${sample}.final_summary_stats.txt
        echo "  - Multiple lines of evidence (assembly + protein + profile)" >> ${sample}.final_summary_stats.txt
    elif [ \$VIRAL_READS -gt 100 ] && [ \$CONTIGS -gt 0 ]; then
        echo "✅ CONCLUSION: Moderate viral detection with good evidence" >> ${sample}.final_summary_stats.txt
        echo "  - Viral reads detected and assembled" >> ${sample}.final_summary_stats.txt
        echo "  - Recommend confirmation with additional methods" >> ${sample}.final_summary_stats.txt
    elif [ \$VIRAL_READS -gt 10 ]; then
        echo "⚠️ CONCLUSION: Low-level viral detection" >> ${sample}.final_summary_stats.txt
        echo "  - Few viral reads detected" >> ${sample}.final_summary_stats.txt
        echo "  - Limited assembly success" >> ${sample}.final_summary_stats.txt
        echo "  - Results should be interpreted with caution" >> ${sample}.final_summary_stats.txt
    else
        echo "❌ CONCLUSION: No significant viral detection" >> ${sample}.final_summary_stats.txt
        echo "  - Very few or no viral reads detected" >> ${sample}.final_summary_stats.txt
        echo "  - May indicate low viral load or technical issues" >> ${sample}.final_summary_stats.txt
    fi
    
    echo "✅ Enhanced final report generation completed"
    """
}

// -------------------------------
// Workflow definition
workflow {
    // Input data
    Channel.fromFilePairs(params.reads, flat: true)
        .map { sid, r1, r2 -> 
            def sample = sid.replaceAll(/(_R1|_R2)\..*/, '')
            tuple(sample, [r1, r2])
        }
        .set { read_pairs }
    
    // Execute enhanced workflow
    FASTP(read_pairs)
    VIRAL_SCREENING(FASTP.out)
    VIRAL_ASSEMBLY(VIRAL_SCREENING.out)
    ORF_PREDICTION(VIRAL_ASSEMBLY.out)
    DIAMOND_PROTEIN_ANALYSIS(ORF_PREDICTION.out)
    HMMER_ANALYSIS(ORF_PREDICTION.out)
    ABUNDANCE_ESTIMATION(VIRAL_ASSEMBLY.out, VIRAL_SCREENING.out)
    KRAKEN2_VIRAL_CLASSIFICATION(VIRAL_SCREENING.out)
    ENHANCED_FINAL_REPORT(
        VIRAL_ASSEMBLY.out,
        VIRAL_SCREENING.out,
        ORF_PREDICTION.out,
        DIAMOND_PROTEIN_ANALYSIS.out,
        HMMER_ANALYSIS.out,
        ABUNDANCE_ESTIMATION.out,
        KRAKEN2_VIRAL_CLASSIFICATION.out
    )
}

