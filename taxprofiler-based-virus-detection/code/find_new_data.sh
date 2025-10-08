#!/bin/bash
# 查找新数据文件的脚本

echo "=== 查找新的FASTQ文件 ==="
echo "检查时间: $(date)"
echo

DATA_DIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/data"

echo "🔍 搜索FASTQ文件..."
echo "搜索目录: $DATA_DIR"
echo

# 查找所有可能的FASTQ文件
echo "📁 找到的FASTQ文件:"
find "$DATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | sort

echo
echo "📁 按目录分组:"
find "$DATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | while read file; do
    dir=$(dirname "$file")
    echo "  目录: $dir"
    echo "    文件: $(basename "$file")"
done

echo
echo "🔍 查找成对的R1和R2文件:"
# 查找成对的R1和R2文件
find "$DATA_DIR" -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" 2>/dev/null | while read r1_file; do
    # 获取对应的R2文件
    r2_file=$(echo "$r1_file" | sed 's/_R1_/_R2_/')
    if [ -f "$r2_file" ]; then
        echo "  ✅ 找到配对文件:"
        echo "    R1: $r1_file"
        echo "    R2: $r2_file"
        
        # 提取样本名
        sample_name=$(basename "$r1_file" | sed 's/_R1_.*//')
        echo "    建议样本名: $sample_name"
        echo
    fi
done

echo "=== 搜索完成 ==="
echo
echo "💡 提示: 如果找到了新的文件，请更新samplesheet.csv"
