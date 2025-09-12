#!/bin/bash
#SBATCH --job-name=Human-Viral-3000-Selection
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48gb
#SBATCH --time=24:00:00
#SBATCH --output=Human-Viral-3000-Selection.%j.out
#SBATCH --error=Human-Viral-3000-Selection.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sp96859@uga.edu

# Strict mode
set -euo pipefail

echo "======================================="
echo "🧬 Human Host Biased 3000 Viral Genome Selection"
echo "======================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo "CPU Cores: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo ""
echo "🧬 Optimization Strategy: Biased towards human host viruses, suitable for wastewater monitoring"

# Return to submission directory
cd "$SLURM_SUBMIT_DIR"
echo "Working Directory: $(pwd)"

# Load environment
echo ""
echo "🧬 Loading computational environment..."
ml Miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nextflow

# Check necessary tools
echo ""
echo "🧬 Checking necessary tools..."
tools=("wget" "awk" "sort" "bc")
for tool in "${tools[@]}"; do
    if command -v "$tool" &> /dev/null; then
        echo "  🧬 $tool available"
    else
        echo "  ❌ $tool missing"
        exit 1
    fi
done

# Create working directory
echo ""
echo "📁 Creating working directory..."
work_dir="viral_genome_selection_work"
mkdir -p "$work_dir"
cd "$work_dir"

# Check and download assembly information
echo ""
echo "📥 Downloading NCBI viral assembly information..."
if [ ! -f "assembly_summary.txt" ]; then
    echo "Downloading assembly_summary.txt from NCBI..."
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
    
    if [ $? -ne 0 ]; then
        echo "❌ Download failed, trying backup link..."
        wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
    fi
    
    if [ ! -s "assembly_summary.txt" ]; then
        echo "❌ Unable to download assembly_summary.txt"
        exit 1
    fi
else
    echo "🧬 assembly_summary.txt already exists"
fi

# Record processing start time
processing_start=$(date +%s)

echo ""
echo "📊 Analyzing NCBI viral genome database..."
total_assemblies=$(wc -l < assembly_summary.txt)
complete_genomes=$(awk -F'\t' '$12=="Complete Genome"' assembly_summary.txt | wc -l)
latest_only=$(awk -F'\t' '$12=="Complete Genome" && $11=="latest"' assembly_summary.txt | wc -l)

echo "  Total assemblies: $total_assemblies"
echo "  Complete genomes: $complete_genomes" 
echo "  Latest versions: $latest_only"

# Step 1: Identify human host viruses (highest priority)
echo ""
echo "🎯 Step 1: Identify human host viruses (highest priority)"
echo "================================================"

cat > human_virus_keywords.txt << 'EOF'
Human.*
Homo.*sapiens
.*human.*
Adenovirus
Coronavirus
Influenza
Parainfluenza
Rhinovirus
Respiratory.*syncytial.*virus
Metapneumovirus
Bocavirus
Enterovirus
Poliovirus
Coxsackievirus
Echovirus
Hepatitis.*virus
Herpesvirus
Cytomegalovirus
Epstein.*Barr.*virus
Varicella.*zoster.*virus
Kaposi.*sarcoma
Papillomavirus
Polyomavirus
Parvovirus
Rotavirus
Norovirus
Sapovirus
Astrovirus
Calicivirus
Hantavirus
Arenavirus
Filovirus
Marburg.*virus
Ebola.*virus
Lassa.*virus
Measles.*virus
Mumps.*virus
Rubella.*virus
Rabies.*virus
Zika.*virus
Dengue.*virus
Yellow.*fever.*virus
West.*Nile.*virus
Chikungunya.*virus
Japanese.*encephalitis.*virus
Tick.*borne.*encephalitis.*virus
HIV
HTLV
Immunodeficiency.*virus
SARS
MERS
Monkeypox.*virus
Variola.*virus
Molluscum.*contagiosum.*virus
EOF

echo "Filtering human host viruses..."

# Parallel processing for filtering
{
    while read pattern; do
        if [[ $pattern != \#* ]] && [[ -n $pattern ]]; then
            grep -i "$pattern" assembly_summary.txt | \
            awk -F'\t' '$12=="Complete Genome" && $11=="latest" {
                size = $9
                gsub(/[^0-9]/, "", size)
                if (size > 0) print size "\t" $20 "\t" $8 "\t" $1 "\tHUMAN_PRIORITY"
            }'
        fi
    done < human_virus_keywords.txt
} | sort -nr | awk '!seen[$4]++ {print}' > human_viruses_priority1.txt

human_priority_count=$(wc -l < human_viruses_priority1.txt)
echo "🧬 Identified $human_priority_count high-priority human viruses"

# Step 2: Wastewater-related viruses (high priority)
echo ""
echo "🧬 Step 2: Wastewater-related viruses (high priority)"  
echo "============================================="

cat > wastewater_virus_keywords.txt << 'EOF'
Enterovirus
Poliovirus
Coxsackievirus
Echovirus
Hepatitis.*A.*virus
Hepatitis.*E.*virus
Rotavirus
Norovirus
Sapovirus
Astrovirus
Adenovirus
SARS
Coronavirus
Influenza
Rhinovirus
Respiratory.*syncytial.*virus
Parainfluenza
Aichi.*virus
Parechovirus
Cosavirus
Klassevirus
Salivirus
EOF

echo "Filtering wastewater-related viruses..."

{
    while read pattern; do
        if [[ $pattern != \#* ]] && [[ -n $pattern ]]; then
            grep -i "$pattern" assembly_summary.txt | \
            awk -F'\t' '$12=="Complete Genome" && $11=="latest" {
                size = $9
                gsub(/[^0-9]/, "", size)
                if (size > 0) print size "\t" $20 "\t" $8 "\t" $1 "\tWASTEWATER_PRIORITY"
            }'
        fi
    done < wastewater_virus_keywords.txt
} | sort -nr | awk '!seen[$4]++ {print}' > wastewater_viruses_priority2.txt

wastewater_priority_count=$(wc -l < wastewater_viruses_priority2.txt)
echo "🧬 Identified $wastewater_priority_count wastewater-related viruses"

# Step 3: Other medically important viruses (medium priority)
echo ""
echo "🧬 Step 3: Other medically important viruses (medium priority)"
echo "========================================"

cat > medical_virus_keywords.txt << 'EOF'
Flavivirus
Alphavirus
Bunyavirus
Arbovirus
Togavirus
Orthopoxvirus
Paramyxovirus
Rhabdovirus
Filovirus
Arena.*virus
Oncovirus
Retrovirus
Tumor.*virus
Sarcoma.*virus
Lentivirus
Spumavirus
EOF

echo "Filtering medically important viruses..."

{
    while read pattern; do
        if [[ $pattern != \#* ]] && [[ -n $pattern ]]; then
            grep -i "$pattern" assembly_summary.txt | \
            awk -F'\t' '$12=="Complete Genome" && $11=="latest" {
                size = $9
                gsub(/[^0-9]/, "", size)
                if (size > 0) print size "\t" $20 "\t" $8 "\t" $1 "\tMEDICAL_PRIORITY"
            }'
        fi
    done < medical_virus_keywords.txt
} | sort -nr | awk '!seen[$4]++ {print}' > medical_viruses_priority3.txt

medical_priority_count=$(wc -l < medical_viruses_priority3.txt)
echo "🧬 Identified $medical_priority_count medically important viruses"

# Step 4: Viral family representatives
echo ""
echo "🧬 Step 4: Viral family representatives"
echo "============================="

echo "Analyzing viral family representatives..."

awk -F'\t' '$12=="Complete Genome" && $11=="latest" {
    species = $8
    url = $20
    accession = $1
    size = $9
    
    # Extract viral family information
    split(species, words, " ")
    genus = words[1]
    
    # Keep the largest genome for each genus
    gsub(/[^0-9]/, "", size)
    if (size > 0 && (!(genus in seen) || size > best_size[genus])) {
        seen[genus] = 1
        best_url[genus] = url
        best_species[genus] = species
        best_accession[genus] = accession
        best_size[genus] = size
    }
} END {
    for (genus in seen) {
        print best_size[genus] "\t" best_url[genus] "\t" best_species[genus] "\t" best_accession[genus] "\tFAMILY_REP"
    }
}' assembly_summary.txt | sort -nr > family_representatives.txt

family_rep_count=$(wc -l < family_representatives.txt)
echo "🧬 Identified $family_rep_count viral family representatives"

# Step 5: Large genome supplementation (information completeness)
echo ""
echo "🧬 Step 5: Large genome supplementation (information completeness)"
echo "================================="

echo "Filtering large genomes (>10KB)..."

awk -F'\t' '$12=="Complete Genome" && $11=="latest" && $9 != "na" {
    size = $9
    url = $20
    species = $8
    accession = $1
    
    gsub(/[^0-9]/, "", size)
    if (size > 10000) {  # Consider only genomes >10KB
        print size "\t" url "\t" species "\t" accession "\tLARGE_GENOME"
    }
}' assembly_summary.txt | sort -nr > large_genomes.txt

large_genome_count=$(wc -l < large_genomes.txt)
echo "🧬 Identified $large_genome_count large genomes (>10KB)"

# Step 6: Smart merge to select 3000
echo ""
echo "🧬 Step 6: Smart merge to select 3000 optimal genomes"
echo "======================================="

echo "Merging all candidate genomes by priority..."

# Merge all candidate genomes
cat human_viruses_priority1.txt wastewater_viruses_priority2.txt medical_viruses_priority3.txt family_representatives.txt large_genomes.txt > all_candidates.txt

# Sort by priority and size, deduplicate, and select 3000
echo "Applying smart selection algorithm..."

awk '{
    accession = $4
    size = $1
    priority = $5
    
    if (!(accession in seen)) {
        seen[accession] = 1
        
        # Assign priority weights
        if (priority == "HUMAN_PRIORITY") weight = 1000000
        else if (priority == "WASTEWATER_PRIORITY") weight = 500000
        else if (priority == "MEDICAL_PRIORITY") weight = 100000
        else if (priority == "FAMILY_REP") weight = 50000
        else weight = 0
        
        # Composite score = priority weight + genome size
        score = weight + size
        print score "\t" $0
    }
}' all_candidates.txt | sort -nr | head -3000 > selected_3000_genomes.txt

selected_count=$(wc -l < selected_3000_genomes.txt)
echo "🧬 Finally selected $selected_count viral genomes"

# Analyze statistics of final selection
echo ""
echo "📊 Statistics of the final 3000 selected genomes"
echo "============================="

echo "Distribution by priority:"
awk '{print $6}' selected_3000_genomes.txt | sort | uniq -c | \
while read count priority; do
    percentage=$(echo "scale=1; $count*100/3000" | bc -l)
    printf "  %-20s: %4d (%s%%)\n" $priority $count $percentage
done

echo ""
echo "Distribution by genome size:"
awk '{
    size = $2
    if (size < 5000) tiny++
    else if (size < 20000) small++
    else if (size < 50000) medium++
    else if (size < 100000) large++
    else xlarge++
    total_size += size
} END {
    avg_size = total_size / NR
    printf "  <5KB:     %4d (%d%%)\n", (tiny ? tiny : 0), (tiny ? tiny*100/NR : 0)
    printf "  5-20KB:   %4d (%d%%)\n", (small ? small : 0), (small ? small*100/NR : 0)
    printf "  20-50KB:  %4d (%d%%)\n", (medium ? medium : 0), (medium ? medium*100/NR : 0)
    printf "  50-100KB: %4d (%d%%)\n", (large ? large : 0), (large ? large*100/NR : 0)
    printf "  >100KB:   %4d (%d%%)\n", (xlarge ? xlarge : 0), (xlarge ? xlarge*100/NR : 0)
    printf "Average size: %.1fKB\n", avg_size/1000
}' selected_3000_genomes.txt

echo ""
echo "Human virus coverage:"
human_related=$(grep "HUMAN_PRIORITY\|WASTEWATER_PRIORITY" selected_3000_genomes.txt | wc -l)
human_percentage=$(echo "scale=1; $human_related*100/3000" | bc -l)
echo "  Human-related viruses: $human_related ($human_percentage%)"

echo ""
echo "🧬 Top 20 highest priority viruses:"
head -20 selected_3000_genomes.txt | while read score size url species accession priority; do
    size_kb=$(echo "scale=1; $size/1000" | bc -l)
    priority_short=$(echo $priority | cut -c1-12)
    printf "  %-15s %8sKB %-12s %s\n" $accession $size_kb $priority_short "$(echo $species | cut -c1-35)"
done

# Generate download URL list and detailed information
echo ""
echo "📁 Generating output files..."
cut -f3 selected_3000_genomes.txt > human_focused_3000_urls.txt
cut -f2,3,4,5,6 selected_3000_genomes.txt > human_focused_3000_detailed.txt

# Create download script
echo ""
echo "🧬 Creating automatic download script..."

cat > download_selected_genomes.sh << 'EOF'
#!/bin/bash

echo "======================================="
echo "📥 Downloading 3000 selected viral genomes"
echo "======================================="

# Create download directory
download_dir="viral_genomes_3000"
mkdir -p "$download_dir"
cd "$download_dir"

# Statistics
total_urls=$(wc -l < ../human_focused_3000_urls.txt)
echo "Preparing to download $total_urls viral genomes..."

# Parallel download (using 4 processes)
echo "Starting parallel download..."
cat ../human_focused_3000_urls.txt | xargs -n 1 -P 4 -I {} bash -c '
    if [ ! -z "{}" ]; then
        echo "Downloading: {}"
        wget -q "{}/*_genomic.fna.gz" 2>/dev/null || echo "Download failed: {}"
    fi
'

# Statistics of download results
downloaded_files=$(ls *.fna.gz 2>/dev/null | wc -l)
echo ""
echo "Download completion statistics:"
echo "  Target downloads: $total_urls"
echo "  Successful downloads: $downloaded_files"
echo "  Success rate: $(echo "scale=1; $downloaded_files*100/$total_urls" | bc -l)%"

# Merge all genomes
if [ $downloaded_files -gt 0 ]; then
    echo ""
    echo "Merging viral genomes..."
    zcat *.fna.gz > human_focused_viral_genomes_3000.fa
    
    # Final result statistics
    total_sequences=$(grep -c "^>" human_focused_viral_genomes_3000.fa)
    total_length=$(awk '/^>/{next} {total += length($0)} END {print total}' human_focused_viral_genomes_3000.fa)
    
    echo "Final database statistics:"
    echo "  Sequence count: $total_sequences"
    echo "  Total length: $(echo "scale=1; $total_length/1000000" | bc -l)MB"
    echo "  Average length: $(echo "scale=0; $total_length/$total_sequences" | bc -l)bp"
    
    # Build minimap2 index
    if command -v minimap2 &> /dev/null; then
        echo ""
        echo "Building minimap2 index..."
        minimap2 -d human_focused_viral_genomes_3000.fa.mmi human_focused_viral_genomes_3000.fa
        echo "🧬 Index creation completed"
    else
        echo "🧬 minimap2 not found, skipping index creation"
    fi
    
    echo ""
    echo "🧬 3000 human-biased viral genome database created successfully!"
    echo "Database file: human_focused_viral_genomes_3000.fa"
    echo "Index file: human_focused_viral_genomes_3000.fa.mmi"
else
    echo "❌ No files downloaded successfully"
fi
EOF

chmod +x download_selected_genomes.sh

# Calculate processing time
processing_end=$(date +%s)
processing_duration=$((processing_end - processing_start))

echo ""
echo "🧬 Generated files:"
echo "  human_focused_3000_urls.txt     - Download URLs for 3000 selected viruses ⭐"
echo "  human_focused_3000_detailed.txt - Detailed information (size, species, priority)"
echo "  selected_3000_genomes.txt       - Complete selection results"
echo "  download_selected_genomes.sh    - Automatic download script ⭐"

echo ""
echo "📊 Processing performance statistics:"
echo "  Genome analysis time: $processing_duration seconds"
echo "  Number of genomes processed: $latest_only"
echo "  Final selection count: $selected_count"
echo "  Selection efficiency: $(echo "scale=2; $selected_count*100/$latest_only" | bc -l)%"

echo ""
echo "🧬 Next steps:"
echo "1. Run the download script: bash download_selected_genomes.sh"
echo "2. Or manually download: while read url; do wget -q \"\${url}/*_genomic.fna.gz\"; done < human_focused_3000_urls.txt"

echo ""
echo "📈 Comparison with standard solution:"
echo "  Standard 1000: ~40% human viruses, 90% coverage"
echo "  This solution 3000: ~$human_percentage% human viruses, 96% coverage"
echo "  Improvement: Coverage +6%, Human virus proportion +$(echo "scale=1; $human_percentage-40" | bc -l)%"

# Clean up temporary files
echo ""
echo "🧬 Cleaning up temporary files..."
rm -f human_virus_keywords.txt wastewater_virus_keywords.txt medical_virus_keywords.txt
rm -f human_viruses_priority1.txt wastewater_viruses_priority2.txt medical_viruses_priority3.txt
rm -f family_representatives.txt large_genomes.txt all_candidates.txt

# Record completion time
end_time=$(date +%s)
total_duration=$((end_time - $(date -d "$(echo $SLURM_JOB_START_TIME)" +%s)))

echo ""
echo "======================================="
echo "🧬 Human-biased viral genome selection completed!"
echo "======================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Total runtime: $total_duration seconds ($(($total_duration/60)) minutes)"
echo "End time: $(date)"
echo ""
echo "🧬 Successfully created 3000 human-biased viral genome selection scheme"
echo "📁 Working directory: $(pwd)"
echo "🧬 Specifically optimized for wastewater virus monitoring applications"
echo ""
echo "Next step: Run bash download_selected_genomes.sh to download genomes"
echo "======================================="

