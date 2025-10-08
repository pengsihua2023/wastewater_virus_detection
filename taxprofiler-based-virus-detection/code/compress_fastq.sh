#!/bin/bash
# FASTQ File Compression Script

echo "=== Compress FASTQ Files ==="
echo "Check time: $(date)"
echo

DATA_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/data/reads"

echo "🔍 Searching for FASTQ files to compress..."
echo "Search directory: $DATA_DIR"
echo

# Find uncompressed FASTQ files
find "$DATA_DIR" -name "*.fastq" -not -name "*.gz" 2>/dev/null | while read fastq_file; do
    compressed_file="${fastq_file}.gz"
    
    if [ -f "$fastq_file" ] && [ ! -f "$compressed_file" ]; then
        echo "📦 Compressing file: $(basename "$fastq_file")"
        echo "   Input: $fastq_file"
        echo "   Output: $compressed_file"
        
        # Compress file using gzip
        gzip -c "$fastq_file" > "$compressed_file"
        
        if [ $? -eq 0 ]; then
            echo "   ✅ Compression successful"
            
            # Display file size comparison
            original_size=$(du -h "$fastq_file" | cut -f1)
            compressed_size=$(du -h "$compressed_file" | cut -f1)
            echo "   Original size: $original_size"
            echo "   Compressed size: $compressed_size"
            
            # Verify compressed file
            if gzip -t "$compressed_file" 2>/dev/null; then
                echo "   ✅ Compressed file verification passed"
            else
                echo "   ❌ Compressed file verification failed"
            fi
        else
            echo "   ❌ Compression failed"
        fi
        echo
    fi
done

echo "🔍 Searching for compressed files..."
find "$DATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | sort

echo
echo "=== Compression Complete ==="

