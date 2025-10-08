#!/bin/bash
# 压缩FASTQ文件脚本

echo "=== 压缩FASTQ文件 ==="
echo "检查时间: $(date)"
echo

DATA_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/data/reads"

echo "🔍 查找需要压缩的FASTQ文件..."
echo "搜索目录: $DATA_DIR"
echo

# 查找未压缩的FASTQ文件
find "$DATA_DIR" -name "*.fastq" -not -name "*.gz" 2>/dev/null | while read fastq_file; do
    compressed_file="${fastq_file}.gz"
    
    if [ -f "$fastq_file" ] && [ ! -f "$compressed_file" ]; then
        echo "📦 压缩文件: $(basename "$fastq_file")"
        echo "   输入: $fastq_file"
        echo "   输出: $compressed_file"
        
        # 使用gzip压缩文件
        gzip -c "$fastq_file" > "$compressed_file"
        
        if [ $? -eq 0 ]; then
            echo "   ✅ 压缩成功"
            
            # 显示文件大小对比
            original_size=$(du -h "$fastq_file" | cut -f1)
            compressed_size=$(du -h "$compressed_file" | cut -f1)
            echo "   原始大小: $original_size"
            echo "   压缩后: $compressed_size"
            
            # 验证压缩文件
            if gzip -t "$compressed_file" 2>/dev/null; then
                echo "   ✅ 压缩文件验证通过"
            else
                echo "   ❌ 压缩文件验证失败"
            fi
        else
            echo "   ❌ 压缩失败"
        fi
        echo
    fi
done

echo "🔍 查找压缩后的文件..."
find "$DATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | sort

echo
echo "=== 压缩完成 ==="
