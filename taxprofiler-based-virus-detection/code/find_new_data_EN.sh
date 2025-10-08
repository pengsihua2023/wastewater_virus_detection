#!/bin/bash
# Script to Find New Data Files

echo "=== Find New FASTQ Files ==="
echo "Check time: $(date)"
echo

DATA_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/data"

echo "🔍 Searching for FASTQ files..."
echo "Search directory: $DATA_DIR"
echo

# Find all possible FASTQ files
echo "📁 FASTQ files found:"
find "$DATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | sort

echo
echo "📁 Grouped by directory:"
find "$DATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | while read file; do
    dir=$(dirname "$file")
    echo "  Directory: $dir"
    echo "    File: $(basename "$file")"
done

echo
echo "🔍 Searching for paired R1 and R2 files:"
# Find paired R1 and R2 files
find "$DATA_DIR" -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" 2>/dev/null | while read r1_file; do
    # Get corresponding R2 file
    r2_file=$(echo "$r1_file" | sed 's/_R1_/_R2_/')
    if [ -f "$r2_file" ]; then
        echo "  ✅ Paired files found:"
        echo "    R1: $r1_file"
        echo "    R2: $r2_file"
        
        # Extract sample name
        sample_name=$(basename "$r1_file" | sed 's/_R1_.*//')
        echo "    Suggested sample name: $sample_name"
        echo
    fi
done

echo "=== Search Complete ==="
echo
echo "💡 Tip: If new files are found, please update samplesheet.csv"

